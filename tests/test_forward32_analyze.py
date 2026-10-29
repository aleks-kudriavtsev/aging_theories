"""Synthetic orchestration tests; no public participant data or network."""
from contextlib import redirect_stdout
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from research.mortality.forward32 import analyze as a
from research.mortality.forward32 import cohort


def fixture():
    rng = np.random.default_rng(320)
    early_n, late_n = 175, 120  # first early stratum has3PSUs; all others2
    n = early_n + late_n
    frame = pd.DataFrame(rng.normal(size=(n, 11)), columns=a.CORE)
    frame['SSCYST'] = np.exp2(rng.normal(size=n)); frame['SSB2M'] = np.exp2(rng.normal(size=n))
    frame['partition'] = ['train'] * early_n + ['late'] * late_n
    frame['cycle'] = [1999] * 90 + [2001] * (early_n - 90) + [2003] * late_n
    frame['cluster'] = ([f'early_synthetic_psu_{i // 5}' for i in range(early_n)]
                        + [f'late_synthetic_psu_{i // 5}' for i in range(late_n)])
    frame['stratum'] = ([f'early_{0 if i // 5 < 3 else 1 + (i // 5 - 3) // 2}' for i in range(early_n)]
                        + [f'late_{i // 10}' for i in range(late_n)])
    w = rng.lognormal(0, .4, n)
    frame['WTSSCB4Y'] = np.r_[w[:early_n], np.full(late_n, np.nan)]
    frame['WTSSCB2Y'] = w * 1.3
    frame['weight'] = np.r_[frame.WTSSCB4Y[:early_n], frame.WTSSCB2Y[early_n:]]
    frame['weight6y_original'] = np.where(frame.partition.eq('train'), 2 * frame.weight / 3, frame.weight / 3)
    frame['dead'] = np.arange(n) % 4 == 1
    frame['t_exam'] = np.where(frame.dead, 60., 120.)
    frame['not_for_export'] = 'SYNTHETIC_PERSON_RECORD_MUST_STAY_LOCAL'
    common = np.ones(n, bool); common[::5] = False
    clusters = frame.loc[frame.partition.eq('train'), 'cluster'].to_numpy()
    folds = a.model.numerical.psu_partition(clusters, a.model.INNER_SALT, 3)
    mapping = {key: int(folds[np.flatnonzero(clusters == key)[0]]) for key in set(clusters)}
    audit = dict(gate=dict(temporal_evaluation_allowed=True, blockers=[]),
                 eligible_derived_early_inner_fold_mapping=mapping,
                 periods={'train': {'native_weight': 'WTSSCB4Y'}, 'late': {'native_weight': 'WTSSCB2Y'}})
    manifest = dict(schema='forward32-source-manifest/1', intake={'synthetic_only': True},
        source=dict(files=[dict(file='synthetic.xpt', url='https://example.invalid/source',
            resolved_url='https://example.invalid/release', bytes=10, sha256='not-real-data')]))
    return frame, common, audit, manifest


def prepared():
    frame, common, audit, _ = fixture()
    return (a._prepare_period(frame, common, 'train'), a._prepare_period(frame, common, 'late'), audit)


def fast_metrics(y, p, w):
    if np.sum(w) <= 0:
        raise ValueError('synthetic zero-weight evaluation')
    observed, predicted = float(np.average(y, weights=w)), float(np.average(p, weights=w))
    return dict(auc=.6 if np.unique(y[w > 0]).size == 2 else None,
        brier=float(np.average((y - p)**2, weights=w)), observed=observed, predicted=predicted,
        bias_pp=100 * (predicted - observed), observed_expected=observed / predicted,
        calibration_in_large=None, calibration_intercept=None, calibration_slope=None,
        calibration_status='not_fitted_in_synthetic_orchestration_test')


def fast_arm(core, proteins, y, weights, clusters, late_core, late_proteins, late_clusters, arm):
    if np.sum(weights) <= 0:
        raise ValueError('synthetic zero-weight training')
    p = np.full(len(late_core), np.clip(np.average(y, weights=weights) + .01 * a.ARMS.index(arm), .01, .99))
    return p, dict(selected_penalty=10., synthetic_only=True)


def point_predictions(train, late):
    predictions = {arm: fast_arm(train['core'], train['protein'], train['y'], train['weight'],
        train['clusters'], late['core'], late['protein'], late['clusters'], arm)[0] for arm in a.ARMS}
    return predictions, {arm: {'status': 'computed'} for arm in a.ARMS}


def small_config(replicates=2):
    return dict(replicates_per_mode=replicates, late_seed=20261103, train_seed=20261104)


class NativeSplitTests(unittest.TestCase):
    def test_native_weights_are_used_without_requiring_late_W4(self):
        frame, common, audit, manifest = fixture()
        a._validate_schema(frame, common, audit, manifest)
        train = a._prepare_period(frame, common, 'train')
        late = a._prepare_period(frame, common, 'late')
        self.assertTrue(late['eligible'].WTSSCB4Y.isna().all())
        np.testing.assert_array_equal(train['eligible'].weight, train['eligible'].WTSSCB4Y)
        np.testing.assert_array_equal(late['eligible'].weight, late['eligible'].WTSSCB2Y)
        self.assertGreater(len(train['eligible']), len(train['y']))
        self.assertGreater(len(late['eligible']), len(late['y']))
        a._check_inner_mapping(train, late, audit)

    def test_six_year_weights_wrong_split_and_changed_mapping_are_rejected(self):
        frame, common, audit, manifest = fixture()
        wrong = frame.copy(); wrong['weight'] = wrong.weight6y_original
        with self.assertRaisesRegex(ValueError, 'Native-period'):
            a._prepare_period(wrong, common, 'train')
        wrong = frame.copy(); wrong.loc[0, 'cycle'] = 2003
        with self.assertRaisesRegex(ValueError, 'calendar split'):
            a._validate_schema(wrong, common, audit, manifest)
        train, late, audit = prepared()
        audit['eligible_derived_early_inner_fold_mapping'] = {}
        with self.assertRaisesRegex(ValueError, 'eligible-derived'):
            a._check_inner_mapping(train, late, audit)

    def test_common_inventory_loss_and_source_schema_errors_are_not_fit_failures(self):
        frame, common, audit, manifest = fixture()
        common[:5] = False
        with self.assertRaisesRegex(ValueError, 'eligible PSU inventory'):
            a._prepare_period(frame, common, 'train')
        with self.assertRaisesRegex(ValueError, 'source schema'):
            a._validate_schema(frame.drop(columns='SSCYST'), common, audit, manifest)
        with self.assertRaisesRegex(ValueError, 'boolean mask'):
            a._validate_schema(frame, common.astype(int), audit, manifest)


class PairedBootstrapTests(unittest.TestCase):
    def draw(self, train_seed=20261104, late_seed=20261103, replicates=3):
        train, late, _ = prepared(); predictions, status = point_predictions(train, late)
        with patch.object(a.model, 'fit_training_arm', side_effect=fast_arm), \
             patch.object(a, 'metrics', side_effect=fast_metrics), redirect_stdout(io.StringIO()):
            draws, selection = a.bootstrap_modes(train, late, predictions, status,
                replicates=replicates, late_seed=late_seed, train_seed=train_seed)
        return draws, selection

    def test_same_late_weights_both_modes_and_independent_early_stream(self):
        first, _ = self.draw(); other_train, _ = self.draw(train_seed=20261105)
        other_late, _ = self.draw(late_seed=20261106)
        self.assertEqual(first['fixed_training'], other_train['fixed_training'])
        for rep in range(3):
            fixed, joint = first['fixed_training'][rep], first['joint_refit'][rep]
            for key in fixed:
                if key.startswith('late_'): self.assertEqual(fixed[key], joint[key])
            for key in joint:
                if key.startswith('late_'):
                    self.assertEqual(joint[key], other_train['joint_refit'][rep][key])
                if key.startswith('early_'):
                    self.assertEqual(joint[key], other_late['joint_refit'][rep][key])
            for arm in a.ARMS:
                self.assertEqual(fixed[arm + '|observed'], joint[arm + '|observed'])
        self.assertTrue(any(first['joint_refit'][rep]['early_weight_sum']
                            != other_train['joint_refit'][rep]['early_weight_sum'] for rep in range(3)))

    def test_full_eligible_frame_sampled_and_zero_rows_kept_in_fit(self):
        train, late, _ = prepared(); predictions, statuses = point_predictions(train, late)
        actual = a.design_bootstrap_weights; frames = []; fit_sizes = []
        def capture(frame, rng):
            frames.append((len(frame), frame.partition.unique().item(), frame.weight.sum()))
            return actual(frame, rng)
        def capture_fit(*args):
            fit_sizes.append((len(args[2]), int((args[3] == 0).sum())))
            return fast_arm(*args)
        with patch.object(a, 'design_bootstrap_weights', side_effect=capture), \
             patch.object(a.model, 'fit_training_arm', side_effect=capture_fit), \
             patch.object(a, 'metrics', side_effect=fast_metrics):
            a.bootstrap_modes(train, late, predictions, statuses, replicates=2)
        self.assertEqual([row[:2] for row in frames], [(120, 'late'), (175, 'train')] * 2)
        self.assertEqual(len(fit_sizes), 8)
        self.assertTrue(all(n == len(train['y']) and zeros > 0 for n, zeros in fit_sizes))
        self.assertTrue(all(value[2] > 0 for value in frames))

    def test_exact_500_attempts_per_mode_no_redraw_and_partial_failures(self):
        train, late, _ = prepared(); predictions, statuses = point_predictions(train, late)
        actual_bootstrap = a.design_bootstrap_weights
        def cys_failure(*args):
            if args[-1] == 'routine_cys': raise ArithmeticError('synthetic arm failure')
            return fast_arm(*args)
        with patch.object(a.model, 'fit_training_arm', side_effect=cys_failure) as fitting, \
             patch.object(a, 'metrics', side_effect=fast_metrics), \
             patch.object(a, 'design_bootstrap_weights', side_effect=actual_bootstrap) as resampling, \
             redirect_stdout(io.StringIO()):
            draws, selection = a.bootstrap_modes(train, late, predictions, statuses, replicates=500)
        self.assertEqual(fitting.call_count, 4 * 500)
        self.assertEqual(resampling.call_count, 2 * 500)
        for mode in a.MODES:
            self.assertEqual(len(draws[mode]), 500)
            self.assertEqual([row['replicate'] for row in draws[mode]], list(range(500)))
        self.assertTrue(all(row['primary|status'] == 'computed' for row in draws['joint_refit']))
        self.assertTrue(all(row['secondary|status'] == 'undefined' for row in draws['joint_refit']))
        self.assertEqual(selection['routine_cys']['10.0'], 0)
        self.assertEqual(selection['routine']['10.0'], 500)

    def test_mode_intervals_are_separate_and_metric_local(self):
        draws, _ = self.draw()
        draws['joint_refit'][0]['primary|delta_auc'] = None
        for mode in a.MODES:
            summary = a._mode_summary(draws[mode], 3, mode)
            self.assertEqual(len(summary['intervals']), 40)
            self.assertEqual(summary['intervals']['primary|delta_brier']['valid'], 3)
            self.assertEqual(summary['intervals']['routine|calibration_slope']['undefined'], 3)
            self.assertEqual(summary['is_primary_uncertainty'], mode == 'joint_refit')
            self.assertEqual(summary['intervals']['primary|delta_brier']['early_training_refitted'], mode == 'joint_refit')
        self.assertEqual(a._mode_summary(draws['joint_refit'], 3, 'joint_refit')['intervals']['primary|delta_auc']['undefined'], 1)


class ExportTests(unittest.TestCase):
    def run_fixture(self, directory, supplied=None, failure=None):
        supplied = fixture() if supplied is None else supplied
        fit = a.model.fit_training_arm
        def fitting(*args):
            if failure is not None and args[-1] == failure:
                raise ArithmeticError('synthetic point/replicate arm failure')
            return fit(*args)
        with patch.object(cohort, 'load', return_value=supplied), \
             patch.object(a, '_configuration', return_value=small_config()), \
             patch.object(a.model, 'fit_training_arm', side_effect=fitting), \
             patch.object(a, 'metrics', side_effect=fast_metrics), redirect_stdout(io.StringIO()) as stdout:
            result = a.run('/synthetic-no-source-directory', directory)
        return result, stdout.getvalue()

    def test_export_exact_bytes_hashes_and_two_separate_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'synthetic-output'
            result, log = self.run_fixture(out)
            self.assertEqual(json.loads((out / 'RESULTS.json').read_bytes()), result)
            self.assertEqual(json.loads((out / 'SOURCE_MANIFEST.json').read_bytes()), fixture()[-1])
            for filename, digest in result['outputs_sha256'].items():
                self.assertEqual(hashlib.sha256((out / filename).read_bytes()).hexdigest(), digest)
            fixed = pd.read_csv(out / 'fixed_training_bootstrap.csv')
            joint = pd.read_csv(out / 'joint_refit_bootstrap.csv')
            paired = pd.read_csv(out / 'paired_bootstrap.csv')
            self.assertEqual(len(fixed), 2); self.assertEqual(len(joint), 2); self.assertEqual(len(paired), 4)
            self.assertEqual(set(paired['mode']), set(a.MODES))
            for key in fixed.columns:
                if key.startswith('late_'): np.testing.assert_array_equal(fixed[key], joint[key])
            self.assertEqual(result['primary_uncertainty_mode'], 'joint_refit')
            self.assertFalse(result['uncertainty_modes_pooled'])
            self.assertIn('FORWARD32_PAIRED_BOOTSTRAP_CSV_BEGIN', log)
            self.assertIn('https://example.invalid/release', log)
            self.assertNotIn('SYNTHETIC_PERSON_RECORD_MUST_STAY_LOCAL', log)
            artifact = json.loads((out / 'MODEL_ARTIFACTS.json').read_bytes())
            self.assertNotIn('point', artifact)
            self.assertTrue(artifact['inner_scores_are_training_diagnostics_not_validation'])

    def test_point_secondary_arm_failure_preserves_primary_and_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'synthetic-output'
            result, _ = self.run_fixture(out, failure='routine_cys')
            self.assertEqual(result['contrasts']['primary|status'], 'computed')
            self.assertEqual(result['contrasts']['secondary|status'], 'undefined')
            self.assertEqual(result['uncertainty']['joint_refit']['contrast_failure_counts']['primary'], 0)
            self.assertEqual(result['uncertainty']['joint_refit']['contrast_failure_counts']['secondary'], 2)
            artifact = json.loads((out / 'MODEL_ARTIFACTS.json').read_bytes())
            self.assertEqual(set(artifact['models']), set(a.ARMS) - {'routine_cys'})
            self.assertEqual(artifact['model_statuses']['routine_cys']['status'], 'undefined')

    def test_late_outcomes_cannot_change_models_between_passing_cases(self):
        supplied = fixture(); changed = copy.deepcopy(supplied)
        is_late = changed[0].partition.eq('late')
        changed[0].loc[is_late, 'dead'] = ~changed[0].loc[is_late, 'dead']
        changed[0].loc[is_late, 't_exam'] = np.where(changed[0].loc[is_late, 'dead'], 60., 120.)
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / 'first', Path(directory) / 'second'
            result_a, _ = self.run_fixture(first, supplied)
            result_b, _ = self.run_fixture(second, changed)
            self.assertEqual((first / 'MODEL_ARTIFACTS.json').read_bytes(), (second / 'MODEL_ARTIFACTS.json').read_bytes())
            self.assertNotEqual(result_a['point']['routine']['brier'], result_b['point']['routine']['brier'])

    def test_blocked_gate_exports_diagnostic_only_without_fits_or_draws(self):
        supplied = fixture()
        supplied[2]['gate'] = dict(temporal_evaluation_allowed=False,
                                  blockers=['late:event_support_below_predeclared_minimum'])
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(a, '_fit_arms') as fitting, patch.object(a, 'bootstrap_modes') as resampling:
            out = Path(directory) / 'blocked'
            result, log = self.run_fixture(out, supplied)
            self.assertEqual(result['status'], 'blocked_feasibility_no_model_fits')
            self.assertEqual(result['model_fits_attempted'], 0)
            self.assertEqual(set(p.name for p in out.iterdir()), {'RESULTS.json', 'SOURCE_MANIFEST.json'})
            self.assertNotIn('point', result); self.assertNotIn('uncertainty', result)
            self.assertIn('FORWARD32_MANIFEST_BEGIN', log)
            fitting.assert_not_called(); resampling.assert_not_called()

    def test_invalid_source_schema_raises_instead_of_becoming_blocked_result(self):
        frame, common, audit, manifest = fixture()
        audit['gate']['temporal_evaluation_allowed'] = False
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 'source schema'):
                self.run_fixture(Path(directory) / 'invalid', (frame.drop(columns='SSCYST'), common, audit, manifest))

    def test_empty_eligible_late_part_is_feasibility_diagnostic_without_fits(self):
        frame, common, audit, manifest = fixture()
        keep = frame.partition.eq('train').to_numpy()
        audit['gate'] = dict(temporal_evaluation_allowed=False, blockers=['late:empty_primary_domain'])
        supplied = frame.loc[keep].copy(), common[keep], audit, manifest
        with tempfile.TemporaryDirectory() as directory, patch.object(a, '_fit_arms') as fitting:
            result, _ = self.run_fixture(Path(directory) / 'blocked-empty', supplied)
            self.assertEqual(result['status'], 'blocked_feasibility_no_model_fits')
            fitting.assert_not_called()


if __name__ == '__main__': unittest.main()
