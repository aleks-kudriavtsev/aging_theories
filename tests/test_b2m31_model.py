"""Synthetic-only Cys/B2M selection, nested preprocessing and PSU tests."""
import copy
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from research.mortality.binary28 import model as legacy
from research.mortality.routine29 import model as prior
from research.mortality.b2m31 import model as m
from research.mortality.b2m31 import analyze as analysis


def synthetic(seed=31, n=300):
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(n, 11)); proteins = rng.normal(size=(n, 2))
    y = (rng.random(n) < .3).astype(int); w = rng.lognormal(0, .4, n)
    clusters = np.array([f'psu_{i // 5}' for i in range(n)])
    folds = m.psu_partition(clusters, 'b2m31:outer', 5)
    return core, proteins, y, w, clusters, folds


class PreprocessingTests(unittest.TestCase):
    def test_exact_registered_protein_columns_and_layout(self):
        core, proteins, y, w, *_ = synthetic(n=100)
        expected = {'routine': (), 'routine_b2m': (1,), 'routine_cys': (0,),
                    'routine_cys_b2m': (0, 1)}
        self.assertEqual(m.MODEL_INDICES, expected)
        for name, columns in expected.items():
            with self.subTest(name=name):
                spec = m.fit_preprocessing(core, proteins, w, name)
                z = m.transform_features(core, proteins, spec)
                self.assertEqual(z.shape, (100, 11 + len(columns)))
                self.assertEqual(spec['protein_indices'], list(columns))
                raw = np.column_stack([core, proteins[:, columns]])
                np.testing.assert_allclose(spec['mean'], np.average(raw, axis=0, weights=w))

    def test_unused_protein_does_not_change_preprocessing(self):
        core, proteins, _, w, *_ = synthetic(n=100)
        for name, unused in [('routine', (0, 1)), ('routine_b2m', (0,)), ('routine_cys', (1,))]:
            other = proteins.copy(); other[:, unused] += 1e6
            self.assertEqual(m.fit_preprocessing(core, proteins, w, name),
                             m.fit_preprocessing(core, other, w, name))

    def test_constant_columns_and_zero_weight_support(self):
        core = np.zeros((4, 11)); proteins = np.ones((4, 2))
        core[:, 8] = [1., 2., 2., 1e6]
        w = np.array([1., 2., 1., 0.])
        fit = m.fit_preprocessing(core, proteins, w, 'routine_cys_b2m')
        self.assertAlmostEqual(fit['mean'][8], 1.75)
        self.assertFalse(fit['constant'][8]); self.assertTrue(fit['constant'][0])
        self.assertEqual(fit['scale'][0], 1.)
        self.assertEqual(fit['mean'][11:], [1., 1.])

    def test_invalid_shapes_nonfinite_weights_and_registration_fail(self):
        core, proteins, _, w, *_ = synthetic(n=100)
        for c, p, weight, name in [(core[:, :6], proteins, w, 'routine'),
                                  (core, proteins[:, :1], w, 'routine'),
                                  (core, proteins, w, 'routine_all3'),
                                  (core, proteins, np.zeros(100), 'routine')]:
            with self.assertRaises(ValueError): m.fit_preprocessing(c, p, weight, name)
        core[0, 0] = np.nan
        with self.assertRaises(ValueError): m.fit_preprocessing(core, proteins, w, 'routine')
        core[0, 0] = 0.; proteins[0, 0] = np.inf
        with self.assertRaises(ValueError): m.fit_preprocessing(core, proteins, w, 'routine_b2m')

    def test_tampered_preprocessing_indices_rejected(self):
        core, proteins, _, w, *_ = synthetic(n=100)
        spec = m.fit_preprocessing(core, proteins, w, 'routine_b2m')
        spec['protein_indices'] = [0]
        with self.assertRaisesRegex(ValueError, 'registered'):
            m.transform_features(core, proteins, spec)


class NestedTests(unittest.TestCase):
    def test_unused_protein_does_not_change_nested_prediction(self):
        core, proteins, y, w, clusters, folds = synthetic(n=200)
        for name, unused in [('routine', (0, 1)), ('routine_b2m', (0,)), ('routine_cys', (1,))]:
            other = proteins.copy(); other[:, unused] += 1e4
            p, rows = m.nested_crossfit(core, proteins, y, w, clusters, folds, name)
            q, other_rows = m.nested_crossfit(core, other, y, w, clusters, folds, name)
            np.testing.assert_array_equal(p, q); self.assertEqual(rows, other_rows)

    def test_heldout_outcome_cannot_change_own_prediction_or_selection(self):
        args = synthetic(); changed = list(args); changed[2] = args[2].copy()
        changed[2][args[-1] == 0] = 1 - changed[2][args[-1] == 0]
        p, rows = m.nested_crossfit(*args, model_name='routine_cys_b2m')
        q, other = m.nested_crossfit(*changed, model_name='routine_cys_b2m')
        np.testing.assert_array_equal(p[args[-1] == 0], q[args[-1] == 0])
        self.assertEqual(rows[0], other[0] | {'scoring_events': rows[0]['scoring_events']})

    def test_all_preprocessing_calls_exclude_outer_heldout(self):
        core, proteins, y, w, clusters, folds = synthetic()
        core[:, 0] = np.arange(len(core)); calls = []; actual = m.fit_preprocessing
        def capture(c, p, weight, name):
            calls.append(set(c[:, 0].astype(int)))
            return actual(c, p, weight, name)
        with patch.object(m, 'fit_preprocessing', side_effect=capture):
            m.nested_crossfit(core, proteins, y, w, clusters, folds, 'routine_cys_b2m')
        heldout = set(np.flatnonzero(folds == 0))
        self.assertEqual(len(calls), 20)
        for call in calls[:4]: self.assertFalse(call & heldout)

    def test_shared_outer_support_and_feature_names(self):
        args = synthetic(n=200); supports = []
        for name, indices in m.MODEL_INDICES.items():
            p, rows = m.nested_crossfit(*args, model_name=name)
            self.assertTrue(np.isfinite(p).all()); self.assertTrue(((p >= 0) & (p <= 1)).all())
            supports.append([(r['outer_fold'], r['training_n'], r['scoring_n'], r['training_events'],
                              r['scoring_events']) for r in rows])
            for row in rows:
                self.assertEqual(row['protein_feature_indices'], list(indices))
                self.assertEqual(row['feature_names'], list(m.CORE_NAMES) + [m.PROTEIN_NAMES[i] for i in indices])
                self.assertEqual(len(row['core_coefficients_standardized']), 11)
                self.assertEqual(len(row['protein_coefficients_standardized']), len(indices))
                self.assertEqual(row['PSU_overlap_n'], 0); self.assertLessEqual(row['gradient_max'], 1e-4)
        self.assertTrue(all(s == supports[0] for s in supports))

    def test_psu_crossing_missing_folds_and_wrong_outcomes_rejected(self):
        args = list(synthetic()); args[-1][1] = (args[-1][1] + 1) % 5
        with self.assertRaisesRegex(ValueError, 'PSU crosses'): m.nested_crossfit(*args)
        args = list(synthetic()); args[-1][args[-1] == 4] = 3
        with self.assertRaises(ValueError): m.nested_crossfit(*args)
        args = list(synthetic()); args[2][0] = 2
        with self.assertRaises(ValueError): m.nested_crossfit(*args)

    def test_zero_weight_psus_keep_original_folds_and_raw_count_normalization(self):
        args = list(synthetic()); args[3][:25] = 0
        _, rows = m.nested_crossfit(*args, model_name='routine_b2m')
        for row in rows:
            self.assertLess(row['training_positive_weight_n'], row['training_n'])
            self.assertAlmostEqual(row['normalized_weight_sum'], row['training_n'])

    def test_global_weight_scale_does_not_change_prediction(self):
        args = list(synthetic(n=200)); p, rows = m.nested_crossfit(*args, model_name='routine_b2m')
        args[3] = 17 * args[3]; q, other = m.nested_crossfit(*args, model_name='routine_b2m')
        np.testing.assert_allclose(p, q, atol=1e-12, rtol=1e-11)
        self.assertEqual([r['lambda'] for r in rows], [r['lambda'] for r in other])

    def test_constant_features_tie_to_stronger_ridge(self):
        args = list(synthetic()); args[0][:] = 0; args[1][:] = 0
        _, rows = m.nested_crossfit(*args, model_name='routine_cys_b2m')
        self.assertTrue(all(r['lambda'] == 100. for r in rows))

    def test_legacy_registries_and_inputs_not_mutated(self):
        registries = copy.deepcopy((legacy.MODEL_INDICES, prior.MODEL_INDICES))
        args = synthetic(n=200); originals = [a.copy() for a in args]
        for name in m.MODEL_INDICES: m.nested_crossfit(*args, model_name=name)
        self.assertEqual((legacy.MODEL_INDICES, prior.MODEL_INDICES), registries)
        for actual, original in zip(args, originals): np.testing.assert_array_equal(actual, original)

    def test_routine_numerics_match_prior_core_with_identical_inner_salt(self):
        core, proteins, y, w, clusters, folds = synthetic(n=200)
        original = m.psu_partition
        def previous_salt(c, salt, folds=3): return original(c, salt.replace('b2m31:', 'routine29:'), folds)
        with patch.object(m, 'psu_partition', side_effect=previous_salt):
            p, _ = m.nested_crossfit(core, proteins, y, w, clusters, folds, 'routine')
        q, _ = prior.nested_crossfit(core, np.column_stack([proteins, np.zeros(len(y))]),
                                    y, w, clusters, folds, 'routine')
        np.testing.assert_allclose(p, q, atol=1e-14, rtol=0)


class AggregatePipelineTests(unittest.TestCase):
    def test_outer_partition_is_outcome_blind_and_whole_psu(self):
        clusters = synthetic()[4]
        first = analysis.assign_outer(clusters)
        order = np.arange(len(clusters))[::-1]
        second = analysis.assign_outer(clusters[order])
        np.testing.assert_array_equal(first[order], second)
        self.assertEqual(set(first), set(range(5)))
        for unit in set(clusters): self.assertEqual(len(set(first[clusters == unit])), 1)
        with self.assertRaises(ValueError): analysis.assign_outer(clusters[:20])

    def test_rescaled_weights_share_multiplier_inside_each_psu(self):
        _, _, _, weights, clusters, _ = synthetic(n=200)
        strata = np.array([f'stratum_{i // 10}' for i in range(200)])
        frame = pd.DataFrame(dict(weight=weights, cluster=clusters, stratum=strata))
        first = analysis.design_bootstrap_weights(frame, np.random.default_rng(31))
        second = analysis.design_bootstrap_weights(frame, np.random.default_rng(31))
        np.testing.assert_array_equal(first, second)
        for stratum in set(strata):
            values = first[strata == stratum] / weights[strata == stratum]
            np.testing.assert_array_equal(np.unique(values), [0., 2.])
        for unit in set(clusters):
            self.assertEqual(len(set(first[clusters == unit] / weights[clusters == unit])), 1)

    def test_three_psu_stratum_uses_m_minus_one_draws_and_rescaling(self):
        frame = pd.DataFrame(dict(weight=np.ones(12), cluster=np.repeat(['a', 'b', 'c'], 4),
                                  stratum=['three_psu_stratum'] * 12))
        w = analysis.design_bootstrap_weights(frame, np.random.default_rng(31))
        values = [w[frame.cluster.eq(unit)].iloc[0] if isinstance(w, pd.Series)
                  else w[frame.cluster.eq(unit).to_numpy()][0] for unit in ('a', 'b', 'c')]
        self.assertTrue(set(values) <= {0., 1.5, 3.})
        self.assertEqual(sum(values), 3.)

    def test_invalid_or_singleton_design_is_not_silently_resampled(self):
        _, _, _, weights, clusters, _ = synthetic(n=100)
        frame = pd.DataFrame(dict(weight=weights, cluster=clusters,
                                  stratum=np.array([f'stratum_{i // 10}' for i in range(100)])))
        for mutation in ('singleton', 'cross_stratum', 'negative_weight', 'missing'):
            broken = frame.copy()
            if mutation == 'singleton': broken.loc[:4, 'stratum'] = 'singleton'
            if mutation == 'cross_stratum': broken.loc[0, 'stratum'] = 'other'
            if mutation == 'negative_weight': broken.loc[0, 'weight'] = -1.
            if mutation == 'missing': broken.loc[0, 'stratum'] = None
            with self.subTest(mutation=mutation):
                with self.assertRaises(ValueError):
                    analysis.design_bootstrap_weights(broken, np.random.default_rng(31))

    def test_model_failure_is_local_to_its_prespecified_contrast(self):
        args = synthetic(n=100)
        for failed_arm, failed_pair, remaining_pair in [('routine_cys', 'secondary', 'primary'),
                                                       ('routine_b2m', 'primary', 'secondary')]:
            def fake_fit(*values):
                if values[-1] == failed_arm: raise ArithmeticError('synthetic fit failure')
                return np.linspace(.15, .45, len(values[2])), []
            with patch.object(analysis, 'nested_crossfit', side_effect=fake_fit):
                _, points, statuses, _ = analysis.fit_models(*args, 'bootstrap', tolerate_fit_failure=True)
            flat = analysis.flatten_scores(points, statuses)
            self.assertEqual(flat[remaining_pair + '|status'], 'computed')
            self.assertIsNotNone(flat[remaining_pair + '|delta_brier'])
            self.assertEqual(flat[failed_pair + '|status'], 'undefined')
            self.assertIsNone(flat[failed_pair + '|delta_brier'])
            self.assertEqual(flat[failed_arm + '|status'], 'undefined')

    def test_metric_local_intervals_retain_failures_and_all_500_attempts(self):
        points = {arm: dict.fromkeys(analysis.METRICS, .2) | {'calibration_status': 'computed'}
                  for arm in analysis.ARMS}
        statuses = {arm: {'status': 'computed'} for arm in analysis.ARMS}
        first = analysis.flatten_scores(points, statuses)
        second = dict(first)
        second['routine_cys|calibration_slope'] = None
        second['primary|delta_brier'] = .1
        result = analysis.intervals([first, second], 2)
        self.assertEqual(result['routine_cys|calibration_slope']['undefined'], 1)
        self.assertTrue(result['routine_cys|calibration_slope']['defined_draws_only'])
        self.assertEqual(result['primary|delta_brier']['valid'], 2)
        self.assertAlmostEqual(result['primary|delta_brier']['lower'], .0025)
        self.assertAlmostEqual(result['primary|delta_brier']['upper'], .0975)
        with self.assertRaises(ValueError): analysis.intervals([first], 2)

    def test_synthetic_aggregate_export_binds_full_source_manifest_bytes(self):
        from research.mortality.b2m31.features import CORE
        core, proteins, y, weights, clusters, folds = synthetic(n=200)
        frame = pd.DataFrame(core, columns=CORE)
        frame['SSCYST'] = np.exp2(proteins[:, 0]); frame['SSB2M'] = np.exp2(proteins[:, 1])
        frame['dead'] = y; frame['t_exam'] = np.where(y == 1, 60., 120.)
        frame['weight'] = weights; frame['cluster'] = clusters
        frame['stratum'] = [f'stratum_{i // 10}' for i in range(len(y))]
        mask = np.ones(len(frame), dtype=bool)
        manifest = dict(retrieved_utc='synthetic-no-network', files=[dict(file='synthetic.xpt',
            url='https://example.invalid/source', resolved_url='https://example.invalid/release',
            bytes=10, sha256='synthetic-source-fingerprint')])
        module = types.ModuleType('research.mortality.b2m31.cohort')
        module.load = lambda _: (frame, mask, mask, {'synthetic_only': True}, manifest)
        def aggregate_fits(c, p, outcome, w, units, outer, domain, tolerate_fit_failure=False):
            predictions = {arm: np.linspace(.15, .45, len(outcome)) for arm in analysis.ARMS}
            points = {arm: analysis.metrics(outcome, values, w) for arm, values in predictions.items()}
            statuses = {arm: {'status': 'computed'} for arm in analysis.ARMS}
            return predictions, points, statuses, [dict(model=arm, **{'lambda': 10.}) for arm in analysis.ARMS]
        with tempfile.TemporaryDirectory() as directory, \
             patch.dict(sys.modules, {module.__name__: module}), \
             patch.object(analysis, 'fit_models', side_effect=aggregate_fits), \
             patch.object(analysis, 'REPLICATES', 2), redirect_stdout(io.StringIO()) as log:
            out = Path(directory) / 'synthetic-aggregate-output'
            result = analysis.run(Path(directory) / 'no-source-files', out)
            source_bytes = (out / 'SOURCE_MANIFEST.json').read_bytes()
            self.assertEqual(json.loads(source_bytes), manifest)
            self.assertEqual(result['source_manifest_export_sha256'], hashlib.sha256(source_bytes).hexdigest())
            self.assertEqual(json.loads((out / 'RESULTS.json').read_bytes()), result)
            for name in ('primary', 'secondary'):
                payload = (out / (name + '_bootstrap.csv')).read_bytes()
                self.assertEqual(result[name + '_bootstrap_csv_sha256'], hashlib.sha256(payload).hexdigest())
                self.assertEqual(len(pd.read_csv(io.BytesIO(payload))), 2)
            wide = (out / 'bootstrap_aggregates.csv').read_bytes()
            self.assertEqual(result['bootstrap_aggregates_csv_sha256'], hashlib.sha256(wide).hexdigest())
            self.assertEqual(len(pd.read_csv(io.BytesIO(wide))), 2)
            self.assertIn('B2M31_MANIFEST_BEGIN', log.getvalue())
            self.assertIn('https://example.invalid/release', log.getvalue())
            self.assertNotIn('psu_', log.getvalue().replace('psu_n', 'design_unit_count'))


if __name__ == '__main__': unittest.main()
