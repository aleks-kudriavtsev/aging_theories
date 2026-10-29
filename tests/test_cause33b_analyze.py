"""Synthetic orchestration/mocked-fit tests; no empirical data or retrieval."""
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

from research.mortality.cause33b import analyze as a
from research.mortality.cause33b import cohort


def fixture():
    rng = np.random.default_rng(330)
    n = 175  #35PSUs: one3-PSU stratum and16two-PSU strata
    d = pd.DataFrame(rng.normal(size=(n, 11)), columns=a.CORE)
    d['SSCYST'] = np.exp2(rng.normal(size=n)); d['SSB2M'] = np.exp2(rng.normal(size=n))
    d['weight'] = rng.lognormal(0, .4, n)
    d['cluster'] = [f'synthetic_psu_{i // 5}' for i in range(n)]
    d['stratum'] = [f'stratum_{0 if i // 5 < 3 else 1 + (i // 5 - 3) // 2}' for i in range(n)]
    code = np.array([2., 1., 3., 10., np.nan, np.nan])[np.arange(n) % 6]
    d['cause'] = code; d['dead'] = np.isfinite(code).astype(float)
    d['t_exam'] = np.where(d.dead.eq(1), 96., 120.)
    d.loc[np.arange(n) % 11 == 0, 't_exam'] = 120.  #some late deaths
    d.loc[1, 't_exam'] = 0.
    d['person_secret'] = 'SYNTHETIC_PERSON_VALUE_NOT_FOR_EXPORT'
    common = np.ones(n, bool); common[::5] = False
    targets = cohort.cause_targets(d, common)
    chosen = d.loc[common]
    arrays = (chosen[list(a.CORE)].to_numpy(float),
        np.log2(chosen[list(a.PROTEINS)].to_numpy(float)), targets,
        chosen.weight.to_numpy(float), chosen.cluster.to_numpy(),
        a.assign_outer(d.cluster.to_numpy())[common], d, common)
    gate = dict(gate=dict(source_endpoint_qualified=True, later_comparisons_allowed=True),
        proceed_codes=list(a.CAUSES), causes={c: {'proceed_to_later_protocol': True} for c in a.CAUSES})
    audit = dict(scientific_gate_equality=dict(exact_structure_and_counts=True,
        expected_proceed_codes=list(a.CAUSES)), reproduced33A_gate=gate,
        comparative_models_fitted=False)
    manifest = dict(schema='synthetic-manifest/1', intake={'synthetic_only': True},
        source={'files': [{'file': 'synthetic.xpt', 'resolved_url': 'https://example.invalid/source', 'sha256': 'a' * 64}]})
    return d, common, audit, manifest, arrays


def fast_metrics(y, p, w):
    if w.sum() <= 0: raise ValueError('synthetic zero-weight evaluation')
    observed, predicted = float(np.average(y, weights=w)), float(np.average(p, weights=w))
    return dict(auc=.6 if np.unique(y[w > 0]).size == 2 else None,
        brier=float(np.average((y - p)**2, weights=w)), observed=observed, predicted=predicted,
        bias_pp=100 * (predicted - observed), observed_expected=observed / predicted,
        calibration_in_large=None, calibration_intercept=None, calibration_slope=None,
        calibration_status='not_fitted_in_synthetic_orchestration_test')


def mock_nested(core, proteins, y, w, clusters, folds, arm):
    """Mock model outputs, preserving real pure preprocessing for invariants."""
    p = np.full(len(y), np.clip(np.average(y, weights=w) + .01 * a.ARMS.index(arm), .01, .99))
    rows = []
    for fold in range(5):
        train, score = folds != fold, folds == fold
        prep = a.model.fit_preprocessing(core[train], proteins[train], w[train], arm)
        rows.append(dict(model=arm, outer_fold=fold, **{'lambda': 10.},
            training_mean=prep['mean'], training_scale=prep['scale'], constant_features=prep['constant'],
            feature_names=list(a.model.CORE_NAMES) + [a.model.PROTEIN_NAMES[j] for j in a.model.MODEL_INDICES[arm]],
            core_feature_indices=list(range(11)), protein_feature_indices=list(a.model.MODEL_INDICES[arm]),
            training_n=int(train.sum()), scoring_n=int(score.sum()),
            training_positive_weight_n=int((train & (w > 0)).sum()), training_weight_sum=float(w[train].sum()),
            training_psu_n=len(set(clusters[train])), scoring_psu_n=len(set(clusters[score])), PSU_overlap_n=0))
    return p, rows


class ArraysTests(unittest.TestCase):
    def analyze(self, replicates=3, fit=mock_nested, score=fast_metrics):
        args = fixture()[-1]
        with patch.object(a.model, 'nested_crossfit', side_effect=fit), \
             patch.object(a, 'metrics', side_effect=score), redirect_stdout(io.StringIO()):
            return a.analyze_arrays(*args, replicates=replicates)

    def test_targets_include_competing_early_late_deaths_and_horizon_boundary(self):
        d, common, _, _, arrays = fixture()
        y = arrays[2]; s = d.loc[common]
        for cause in a.CAUSES:
            np.testing.assert_array_equal(y[cause], (s.dead.eq(1) & s.t_exam.le(96) & s.cause.eq(int(cause))).to_numpy())
            self.assertFalse(y[cause][(s.cause.eq(3) & s.t_exam.le(96)).to_numpy()].any())
            self.assertFalse(y[cause][s.t_exam.gt(96).to_numpy()].any())
        self.assertEqual(int(y['001'][0]), 1)  #month0 heart death retained
        self.assertTrue(any((s.dead == 1) & (s.cause == 2) & (s.t_exam == 96)))

    def test_shared_draw_weights_and_folds_across_all_eight_refits(self):
        args = fixture()[-1]; calls = []; frames = []
        actual = a.design_bootstrap_weights
        def draw(frame, rng):
            frames.append(len(frame)); return actual(frame, rng)
        def fitting(*values):
            calls.append((values[2].copy(), values[3].copy(), values[4].copy(), values[5].copy()))
            return mock_nested(*values)
        with patch.object(a.model, 'nested_crossfit', side_effect=fitting), \
             patch.object(a, 'metrics', side_effect=fast_metrics), \
             patch.object(a, 'design_bootstrap_weights', side_effect=draw):
            result, parameters, bins, wide, paired = a.analyze_arrays(*args, replicates=3)
        self.assertEqual(frames, [len(args[6])] * 3)
        self.assertEqual(len(calls), 8 * 4)
        for group in range(4):
            rows = calls[group * 8:(group + 1) * 8]
            for row in rows:
                np.testing.assert_array_equal(row[1], rows[0][1])
                np.testing.assert_array_equal(row[2], args[4]); np.testing.assert_array_equal(row[3], args[5])
            if group:
                self.assertTrue((rows[0][1] == 0).any())
                self.assertEqual(wide[group - 1]['shared_common_weight_sha256'], a._weight_sha(rows[0][1]))
            np.testing.assert_array_equal(rows[0][0], args[2]['002'])
            np.testing.assert_array_equal(rows[4][0], args[2]['001'])
        self.assertEqual(len(parameters), 40); self.assertEqual(len(paired), 12)
        self.assertTrue(result['shared_preprocessing_verified_for_available_matched_fits'])

    def test_exact_500_attempts_and_2000_paired_rows_no_redraw(self):
        args = fixture()[-1]
        def failed_arm(*values):
            if values[-1] == 'routine_cys' and np.array_equal(values[2], args[2]['002']):
                raise ArithmeticError('synthetic secondary failure')
            return mock_nested(*values)
        with patch.object(a.model, 'nested_crossfit', side_effect=failed_arm) as fitting, \
             patch.object(a, 'metrics', side_effect=fast_metrics), \
             patch.object(a, 'design_bootstrap_weights', wraps=a.design_bootstrap_weights) as drawing, \
             redirect_stdout(io.StringIO()):
            result, parameters, bins, wide, paired = a.analyze_arrays(*args, replicates=500)
        self.assertEqual(drawing.call_count, 500); self.assertEqual(fitting.call_count, 8 * 501)
        self.assertEqual(len(wide), 500); self.assertEqual(len(paired), 2000)
        self.assertEqual([row['replicate'] for row in wide], list(range(500)))
        self.assertTrue(all(row['002|primary|status'] == 'computed' for row in wide))
        self.assertTrue(all(row['002|secondary|status'] == 'undefined' for row in wide))
        self.assertTrue(all(row['001|secondary|status'] == 'computed' for row in wide))
        self.assertEqual(result['bootstrap_penalty_selection_counts']['002']['routine']['10.0'], 2500)
        self.assertEqual(result['bootstrap_penalty_selection_counts']['002']['routine_cys']['10.0'], 0)
        self.assertEqual(len(parameters), 35)  #point fits only, neverbootstrap vectors

    def test_metric_local_auc_and_calibration_undefined_keep_brier(self):
        def undefined_auc(y, p, w):
            result = fast_metrics(y, p, w); result['auc'] = None; return result
        result, _, _, wide, paired = self.analyze(score=undefined_auc)
        for cause in a.CAUSES:
            self.assertEqual(len(result['intervals'][cause]), 40)
            self.assertEqual(result['intervals'][cause]['primary|delta_brier']['valid'], 3)
            self.assertEqual(result['intervals'][cause]['primary|delta_auc']['undefined'], 3)
            self.assertEqual(result['intervals'][cause]['routine|calibration_slope']['undefined'], 3)
            self.assertEqual(result['bootstrap_failures'][cause]['model_failure_counts']['routine'], 0)
            self.assertEqual(result['bootstrap_failures'][cause]['calibration_undefined_reasons']['routine'],
                             {'not_fitted_in_synthetic_orchestration_test': 3})
        self.assertTrue(all(row['status'] == 'computed' and row['auc_status'] == 'undefined' for row in paired))

    def test_preprocessing_cannot_depend_on_cause(self):
        args = fixture()[-1]
        def contaminated(*values):
            p, rows = mock_nested(*values)
            if np.array_equal(values[2], args[2]['001']): rows[0]['training_mean'][0] += 1.
            return p, rows
        with self.assertRaisesRegex(ArithmeticError, 'cause-dependent'):
            self.analyze(fit=contaminated)

    def test_structural_errors_stop_before_any_fit(self):
        good = fixture()[-1]
        bad = []
        args = list(good); args[2] = {'002': good[2]['002']}; bad.append(args)
        args = list(good); args[2] = {c: y.copy() for c, y in good[2].items()}
        args[2]['001'][args[2]['002'] == 1] = 1; bad.append(args)
        args = list(good); args[3] = args[3] * 2; bad.append(args)
        args = list(good); args[5] = args[5].copy(); args[5][0] = (args[5][0] + 1) % 5; bad.append(args)
        args = list(good); args[0] = args[0][:, :10]; bad.append(args)
        for args in bad:
            with patch.object(a.model, 'nested_crossfit') as fitting:
                with self.assertRaises(ValueError): a.analyze_arrays(*args, replicates=1)
                fitting.assert_not_called()

    def test_seed_replay_aggregate_only_and_unchanged_inputs(self):
        args = fixture()[-1]
        arrays = [args[i].copy() for i in (0, 1, 3, 4, 5, 7)]
        registry = copy.deepcopy(a.model.MODEL_INDICES)
        first = self.analyze(); second = self.analyze()
        self.assertEqual(a._bytes(first), a._bytes(second))
        self.assertEqual(registry, a.model.MODEL_INDICES)
        for index, saved in zip((0, 1, 3, 4, 5, 7), arrays): np.testing.assert_array_equal(args[index], saved)
        text = a._bytes(first).decode()
        self.assertNotIn('SYNTHETIC_PERSON_VALUE_NOT_FOR_EXPORT', text)
        self.assertNotIn('synthetic_psu_', text)
        self.assertFalse(first[0]['between_cause_effect_comparisons_performed'])
        self.assertFalse(first[0]['independent_binary_probabilities_summed'])


class ExportTests(unittest.TestCase):
    def run_fixture(self, out, supplied=None):
        d, common, audit, manifest, _ = fixture() if supplied is None else supplied
        targets = cohort.cause_targets(d, common)
        contract = json.loads((a.ROOT / 'CONTRACT.json').read_text())
        contract.update(common_n=int(common.sum()), eligible_n=len(d),
            cause_events8={c: int(y.sum()) for c, y in targets.items()},
            controls8={c: int(common.sum() - y.sum()) for c, y in targets.items()})
        with patch.object(cohort, 'load', return_value=(d, common, audit, manifest)), \
             patch.object(a, '_contract', return_value=contract), patch.object(a, 'REPLICATES', 2), \
             patch.object(a.model, 'nested_crossfit', side_effect=mock_nested), \
             patch.object(a, 'metrics', side_effect=fast_metrics), redirect_stdout(io.StringIO()) as stdout:
            result = a.run('/synthetic-no-source-files', out)
        return result, stdout.getvalue()

    def test_six_files_exact_marked_bytes_hashes_and_pair_shape(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'result'
            result, stdout = self.run_fixture(out)
            names = {'RESULTS_SUMMARY.json', 'FOLD_PARAMETERS.json', 'CALIBRATION_BINS.json',
                     'SOURCE_MANIFEST.json', 'bootstrap_aggregates.csv', 'paired_bootstrap.csv'}
            self.assertEqual({p.name for p in out.iterdir()}, names)
            for name in names:
                content = (out / name).read_bytes()
                begin, end = 'CA33B_FILE_BEGIN ' + name + '\n', 'CA33B_FILE_END ' + name + '\n'
                self.assertEqual(stdout.split(begin)[1].split(end)[0].encode(), content)
                if name != 'RESULTS_SUMMARY.json': self.assertEqual(hashlib.sha256(content).hexdigest(), result['outputs_sha256'][name])
            self.assertEqual(json.loads((out / 'SOURCE_MANIFEST.json').read_bytes()), fixture()[3])
            self.assertEqual(len(pd.read_csv(out / 'bootstrap_aggregates.csv')), 2)
            pairs = pd.read_csv(out / 'paired_bootstrap.csv', dtype={'cause': str})
            self.assertEqual(len(pairs), 8); self.assertEqual(set(pairs.cause), set(a.CAUSES))
            self.assertNotIn('SYNTHETIC_PERSON_VALUE_NOT_FOR_EXPORT', stdout)
            self.assertTrue(result['comparative_models_fitted'])
            self.assertFalse(result['source_audit']['comparative_models_fitted'])
            self.assertEqual(len(json.loads((out / 'FOLD_PARAMETERS.json').read_bytes())), 40)

    def test_gate_failure_or_changed_counts_stop_before_fit_and_export(self):
        for mutation in ('gate', 'counts'):
            supplied = list(fixture())
            if mutation == 'gate': supplied[2]['scientific_gate_equality']['exact_structure_and_counts'] = False
            else: supplied[2]['reproduced33A_gate']['proceed_codes'] = ['002']
            with tempfile.TemporaryDirectory() as directory, patch.object(a, 'analyze_arrays') as computation:
                out = Path(directory) / 'blocked'
                with self.assertRaisesRegex(ValueError, 'reproduced33A'):
                    self.run_fixture(out, tuple(supplied))
                computation.assert_not_called(); self.assertFalse(out.exists())

    def test_contract_and_existing_output_guard(self):
        actual = a._contract()
        self.assertEqual(actual['bootstrap_seed'], 20261105)
        self.assertEqual(actual['bootstrap_replicates'], 500)
        with tempfile.TemporaryDirectory() as directory, patch.object(cohort, 'load') as loading:
            with self.assertRaisesRegex(ValueError, 'New analysis output'): a.run('unused', directory)
            loading.assert_not_called()


if __name__ == '__main__': unittest.main()
