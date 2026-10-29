"""Synthetic final-training and source-identity tests; no cohort records."""
import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.mortality.binary28 import model as legacy
from research.mortality.routine29 import model
from research.mortality.freeze30 import freeze


def synthetic(seed=30, n=180):
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(n, 11))
    proteins = rng.normal(size=(n, 3))
    y = (rng.random(n) < .3).astype(int)
    w = rng.lognormal(0, .4, n)
    clusters = np.array([f'synthetic_psu_{i // 5}' for i in range(n)])
    return core, proteins, y, w, clusters


class FinalFitTests(unittest.TestCase):
    def test_two_arms_refit_with_fresh_full_training_scaling(self):
        core, proteins, y, w, clusters = synthetic()
        actual = model.fit_preprocessing
        calls = []
        def capture(c, p, weights, arm):
            calls.append((arm, len(c)))
            return actual(c, p, weights, arm)
        with patch.object(model, 'fit_preprocessing', side_effect=capture):
            result = freeze.fit_final(core, proteins, y, w, clusters)
        self.assertEqual(set(result['models']), {'routine', 'routine_all3'})
        self.assertEqual(len(calls), 8)
        for arm in freeze.ARMS:
            row = result['models'][arm]
            self.assertEqual(row['preprocessing'], actual(core, proteins, w, arm))
            self.assertEqual([n for name, n in calls if name == arm][-1], len(y))
            self.assertAlmostEqual(row['fit']['normalized_weight_sum'], len(y))
            self.assertLessEqual(row['fit']['gradient_max'], 1e-4)
            self.assertEqual(row['fit']['penalty'], row['selected_penalty'])
            self.assertEqual(len(row['inner_fit_diagnostics']), 9)
        self.assertEqual(result['selection_scope'], 'training_only_not_validation')

    def test_inner_preprocessing_and_support_follow_whole_psu_partition(self):
        core, proteins, y, w, clusters = synthetic()
        core[:, 0] = np.arange(len(core))
        actual = model.fit_preprocessing
        calls = []
        def capture(c, p, weights, arm):
            calls.append((arm, set(c[:, 0].astype(int))))
            return actual(c, p, weights, arm)
        with patch.object(model, 'fit_preprocessing', side_effect=capture):
            result = freeze.fit_final(core, proteins, y, w, clusters)
        folds = model.psu_partition(clusters, freeze.INNER_SALT, folds=3)
        for arm in freeze.ARMS:
            rows = result['models'][arm]['inner_folds_aggregate']
            observed = [indices for name, indices in calls if name == arm]
            self.assertEqual(len(observed), 4)
            for fold, row in enumerate(rows):
                training = folds != fold
                scoring = ~training
                self.assertEqual(observed[fold], set(np.flatnonzero(training)))
                self.assertEqual(row['training_n'], int(training.sum()))
                self.assertEqual(row['scoring_n'], int(scoring.sum()))
                self.assertEqual(row['training_events'], int(y[training].sum()))
                self.assertEqual(row['scoring_events'], int(y[scoring].sum()))
                self.assertEqual(row['training_psu_n'], len(set(clusters[training])))
                self.assertEqual(row['scoring_psu_n'], len(set(clusters[scoring])))
                self.assertEqual(row['PSU_overlap_n'], 0)
                self.assertFalse(set(clusters[training]) & set(clusters[scoring]))
        self.assertEqual(result['models']['routine']['inner_folds_aggregate'],
                         result['models']['routine_all3']['inner_folds_aggregate'])

    def test_penalty_scores_are_pooled_weighted_inner_predictions(self):
        core, proteins, y, w, clusters = synthetic(n=120)
        fitted = freeze.fit_final(core, proteins, y, w, clusters)
        folds = model.psu_partition(clusters, freeze.INNER_SALT, folds=3)
        for arm in freeze.ARMS:
            expected = {}
            for penalty in model.PENALTIES:
                residual_sum = 0.
                for fold in range(3):
                    train = folds != fold
                    heldout = ~train
                    prep = model.fit_preprocessing(core[train], proteins[train], w[train], arm)
                    fit = model.fit_ridge(model.transform_features(core[train], proteins[train], prep),
                                          y[train], w[train], penalty)
                    p = model.predict_ridge(model.transform_features(core[heldout], proteins[heldout], prep), fit)
                    residual_sum += float(np.dot(w[heldout], (y[heldout] - p)**2))
                expected[str(penalty)] = residual_sum / w.sum()
            row = fitted['models'][arm]
            for penalty, score in expected.items():
                self.assertAlmostEqual(row['inner_brier_by_penalty'][penalty], score, places=14)
            minimum = min(expected.values())
            self.assertEqual(row['selected_penalty'],
                max(float(p) for p, v in expected.items() if v <= minimum + 1e-12))

    def test_constant_features_choose_stronger_tied_ridge(self):
        core, proteins, y, w, clusters = synthetic()
        core[:] = 3.; proteins[:] = -1.
        result = freeze.fit_final(core, proteins, y, w, clusters)
        expected = np.average(y, weights=w)
        for row in result['models'].values():
            self.assertEqual(row['selected_penalty'], 100.)
            np.testing.assert_array_equal(row['fit']['coefficients'], np.zeros(len(row['feature_names'])))
            self.assertTrue(all(row['preprocessing']['constant']))
            self.assertAlmostEqual(1 / (1 + np.exp(-row['fit']['intercept'])), expected, places=10)

    def test_global_weight_scaling_preserves_parameters_and_selection(self):
        core, proteins, y, w, clusters = synthetic()
        first = freeze.fit_final(core, proteins, y, w, clusters)
        second = freeze.fit_final(core, proteins, y, w * 13., clusters)
        for arm in freeze.ARMS:
            a, b = first['models'][arm], second['models'][arm]
            self.assertEqual(a['selected_penalty'], b['selected_penalty'])
            for key in ('mean', 'scale'):
                np.testing.assert_allclose(a['preprocessing'][key], b['preprocessing'][key], atol=1e-14, rtol=1e-13)
            np.testing.assert_allclose(a['fit']['coefficients'], b['fit']['coefficients'], atol=1e-10, rtol=1e-9)
            self.assertAlmostEqual(a['fit']['intercept'], b['fit']['intercept'], places=10)

    def test_zero_weight_rows_do_not_determine_scaling_but_keep_partition(self):
        core, proteins, y, w, clusters = synthetic()
        w[:10] = 0
        core[:10, 8] = 1e5
        result = freeze.fit_final(core, proteins, y, w, clusters)
        self.assertEqual(result['training']['positive_weight_n'], len(y) - 10)
        self.assertEqual(result['training']['psu_n'], len(set(clusters)))
        for row in result['models'].values():
            self.assertAlmostEqual(row['preprocessing']['mean'][8], np.average(core[:, 8], weights=w))
            self.assertAlmostEqual(row['fit']['normalized_weight_sum'], len(y))
            self.assertTrue(any(r['training_positive_weight_n'] < r['training_n']
                                for r in row['inner_folds_aggregate']))

    def test_inputs_and_legacy_registries_unchanged_and_bytes_deterministic(self):
        args = synthetic(n=120)
        originals = [a.copy() for a in args]
        registry = copy.deepcopy(legacy.MODEL_INDICES)
        first = freeze.fit_final(*args)
        second = freeze.fit_final(*args)
        self.assertEqual(freeze._canonical_bytes(first), freeze._canonical_bytes(second))
        self.assertEqual(legacy.MODEL_INDICES, registry)
        for actual, original in zip(args, originals):
            np.testing.assert_array_equal(actual, original)

    def test_no_participant_vectors_or_training_performance_exported(self):
        result = freeze.fit_final(*synthetic())
        self.assertEqual(set(result), {'models', 'training', 'selection_scope'})
        encoded = freeze._canonical_bytes(result).decode()
        for forbidden in ('synthetic_psu_', 'predictions', 'outcomes', 'training_auc',
                          'training_brier', 'calibration', 'records'):
            self.assertNotIn(forbidden, encoded)
        for row in result['models'].values():
            self.assertEqual(set(row), {'preprocessing', 'fit', 'feature_names',
                'inner_brier_by_penalty', 'selected_penalty', 'inner_folds_aggregate',
                'inner_fit_diagnostics'})

    def test_invalid_shapes_values_and_cluster_support_fail(self):
        good = synthetic(n=90)
        bad = []
        for index, value in ((0, good[0][:, :10]), (1, good[1][:, :2]),
                             (2, np.ones(90)), (3, np.zeros(90)),
                             (4, good[4][:-1]), (4, np.array(['only'] * 90)),
                             (4, np.array([None] * 90))):
            args = list(good); args[index] = value; bad.append(args)
        for index, position, value in ((0, (0, 1), np.nan), (1, (0, 0), np.inf),
                                       (2, 0, 2), (3, 0, -1)):
            args = [a.copy() for a in good]; args[index][position] = value; bad.append(args)
        for args in bad:
            with self.subTest(shape=[a.shape for a in args]):
                with self.assertRaises((ValueError, ArithmeticError)):
                    freeze.fit_final(*args)


class SourceIdentityTests(unittest.TestCase):
    def fixtures(self):
        names = ['index.htm'] + [f'source_{i}' for i in range(13)]
        files = [dict(file=name, sha256=f'old_{name}', url=f'https://example.invalid/{name}',
                      resolved_url=f'https://example.invalid/{name}') for name in names]
        prior = {'scientific_python_sha256': {'old/model.py': 'old_code'},
                 'source_audit': {'routine_source_fingerprints':
                                  {r['file']: r['sha256'] for r in files}}}
        return prior, {'files': files}

    def invoke(self, old, manifest, code_digest='old_code'):
        def metadata(path):
            return old if Path(path).name == 'RESULTS_SUMMARY.json' else self.fixtures()[1]
        sentinel = object()
        with patch.object(freeze.adapter, 'check_study'), \
             patch.object(freeze.adapter, 'strict_json', side_effect=metadata), \
             patch.object(freeze.adapter, 'digest', return_value=code_digest), \
             patch.object(freeze, 'load', return_value=(sentinel, 'common', 'strict', {}, manifest)) as loader:
            result = freeze._check_training_sources(Path('/not-a-real-source'))
        return result, loader

    def test_only_dynamic_index_change_allowed_and_recorded(self):
        old, manifest = self.fixtures()
        manifest['files'][0]['sha256'] = 'new_index'
        result, loader = self.invoke(old, manifest)
        policy = result[-1]
        self.assertEqual(policy['pinned_non_index_files'], 13)
        self.assertEqual(policy['prior_index_sha256'], 'old_index.htm')
        self.assertEqual(policy['actual_index_sha256'], 'new_index')
        self.assertTrue(policy['index_changed'])
        self.assertTrue(policy['every_non_index_source_sha_unchanged'])
        loader.assert_called_once()

    def test_historical_code_change_blocks_before_loading(self):
        old, manifest = self.fixtures()
        with patch.object(freeze.adapter, 'check_study'), \
             patch.object(freeze.adapter, 'strict_json', return_value=old), \
             patch.object(freeze.adapter, 'digest', return_value='changed'), \
             patch.object(freeze, 'load') as loader:
            with self.assertRaisesRegex(ValueError, 'Historical scientific code changed'):
                freeze._check_training_sources('/not-a-real-source')
            loader.assert_not_called()

    def test_non_index_bytes_file_set_and_release_url_changes_block(self):
        for mutation in ('bytes', 'missing', 'extra', 'release'):
            old, manifest = self.fixtures()
            if mutation == 'bytes': manifest['files'][1]['sha256'] = 'new_data'
            if mutation == 'missing': manifest['files'].pop()
            if mutation == 'extra': manifest['files'].append(dict(file='unexpected', sha256='new'))
            if mutation == 'release': manifest['files'][0]['resolved_url'] += '?release=2'
            with self.subTest(mutation=mutation):
                with self.assertRaisesRegex(ValueError, 'Historical'):
                    self.invoke(old, manifest)

    def test_existing_output_is_rejected_before_source_or_training(self):
        with tempfile.TemporaryDirectory() as existing, \
             patch.object(freeze, '_check_training_sources') as source, \
             patch.object(freeze, 'fit_final') as fitting:
            with self.assertRaisesRegex(ValueError, 'New freeze output directory required'):
                freeze.run('/not-a-real-source', existing)
            source.assert_not_called(); fitting.assert_not_called()


if __name__ == '__main__':
    unittest.main()
