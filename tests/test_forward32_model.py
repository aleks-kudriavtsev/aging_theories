"""Synthetic-only early-training/late-scoring leakage and numerical tests."""
import copy
import inspect
import json
import unittest
from unittest.mock import patch

import numpy as np

from research.mortality.b2m31 import model as numerical
from research.mortality.forward32 import model as m


def synthetic(seed=32, n_train=180, n_late=80):
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(n_train, 11)); proteins = rng.normal(size=(n_train, 2))
    y = (rng.random(n_train) < .3).astype(int); w = rng.lognormal(0, .4, n_train)
    clusters = np.array([f'early_synthetic_psu_{i // 5}' for i in range(n_train)])
    late_core = rng.normal(size=(n_late, 11)); late_proteins = rng.normal(size=(n_late, 2))
    late_clusters = np.array([f'late_synthetic_psu_{i // 5}' for i in range(n_late)])
    return core, proteins, y, w, clusters, late_core, late_proteins, late_clusters


class ForwardFitTests(unittest.TestCase):
    def test_api_accepts_no_late_outcomes_or_late_weights(self):
        expected = ['core_train', 'protein_train', 'y_train', 'w_train', 'clusters_train',
                    'core_late', 'protein_late', 'clusters_late']
        self.assertEqual(list(inspect.signature(m.fit_forward).parameters), expected)
        self.assertEqual(list(inspect.signature(m.fit_training_arm).parameters), expected + ['arm'])
        for name in ('y_late', 'w_late', 'late_outcomes', 'late_weights'):
            with self.assertRaises(TypeError): m.fit_forward(*synthetic(), **{name: np.zeros(80)})

    def test_four_arms_fresh_full_early_transform_and_weight_metadata(self):
        core, proteins, y, w, clusters, *late = synthetic()
        predictions, artifact = m.fit_forward(core, proteins, y, w, clusters, *late)
        self.assertEqual(set(predictions), set(m.ARMS)); self.assertEqual(set(artifact['models']), set(m.ARMS))
        self.assertEqual(artifact['training']['weight_name'], 'WTSSCB4Y')
        self.assertEqual(artifact['late_domain']['evaluation_weight_name'], 'WTSSCB2Y')
        self.assertEqual(artifact['training']['cycles'], [1999, 2001])
        self.assertEqual(artifact['late_domain']['cycle'], 2003)
        self.assertEqual(artifact['late_domain']['training_psu_overlap_n'], 0)
        for arm, row in artifact['models'].items():
            self.assertEqual(row['preprocessing'], numerical.fit_preprocessing(core, proteins, w, arm))
            self.assertEqual(row['fit']['penalty'], row['selected_penalty'])
            self.assertAlmostEqual(row['fit']['normalized_weight_sum'], len(y))
            self.assertLessEqual(row['fit']['gradient_max'], 1e-4)
            self.assertEqual(len(row['inner_fit_diagnostics']), 9)
            self.assertEqual(row['feature_names'], list(numerical.CORE_NAMES)
                             + [numerical.PROTEIN_NAMES[j] for j in numerical.MODEL_INDICES[arm]])
            self.assertEqual(predictions[arm].shape, (len(late[0]),))
            self.assertTrue(((predictions[arm] >= 0) & (predictions[arm] <= 1)).all())

    def test_late_measurements_cannot_change_training_artifacts(self):
        args = list(synthetic())
        first, a = m.fit_forward(*args)
        args[5] = args[5] + 1000.; args[6] = args[6] - 1000.
        second, b = m.fit_forward(*args)
        self.assertEqual(a, b)
        self.assertTrue(any(not np.array_equal(first[arm], second[arm]) for arm in m.ARMS))

    def test_inner_preprocessing_excludes_inner_heldout_and_all_late_rows(self):
        args = list(synthetic())
        args[0][:, 0] = np.arange(len(args[0])); args[5][:, 0] = 1e6 + np.arange(len(args[5]))
        calls = []; actual = numerical.fit_preprocessing
        def capture(core, proteins, w, arm):
            calls.append((arm, set(core[:, 0].astype(int))))
            return actual(core, proteins, w, arm)
        with patch.object(numerical, 'fit_preprocessing', side_effect=capture):
            _, artifact = m.fit_forward(*args)
        folds = numerical.psu_partition(args[4], m.INNER_SALT, 3)
        self.assertEqual(len(calls), 16)
        for arm in m.ARMS:
            observed = [rows for name, rows in calls if name == arm]
            for fold in range(3):
                self.assertEqual(observed[fold], set(np.flatnonzero(folds != fold)))
            self.assertEqual(observed[-1], set(range(len(args[0]))))
        for row in artifact['inner_folds_aggregate']:
            train, score = folds != row['fold'], folds == row['fold']
            self.assertEqual(row['training_n'], int(train.sum()))
            self.assertEqual(row['scoring_n'], int(score.sum()))
            self.assertEqual(row['training_events'], int(args[2][train].sum()))
            self.assertEqual(row['scoring_events'], int(args[2][score].sum()))
            self.assertEqual(row['PSU_overlap_n'], 0)
            self.assertFalse(set(args[4][train]) & set(args[4][score]))

    def test_penalty_scores_are_pooled_weighted_inner_brier(self):
        core, proteins, y, w, clusters, *late = synthetic(n_train=120)
        _, artifact = m.fit_forward(core, proteins, y, w, clusters, *late)
        folds = numerical.psu_partition(clusters, m.INNER_SALT, 3)
        for arm in m.ARMS:
            scores = {}
            for penalty in numerical.PENALTIES:
                total = 0.
                for fold in range(3):
                    train, score = folds != fold, folds == fold
                    prep = numerical.fit_preprocessing(core[train], proteins[train], w[train], arm)
                    fit = numerical.fit_ridge(numerical.transform_features(core[train], proteins[train], prep),
                                              y[train], w[train], penalty)
                    p = numerical.predict_ridge(numerical.transform_features(core[score], proteins[score], prep), fit)
                    total += float(np.dot(w[score], (y[score] - p)**2))
                scores[penalty] = total / w.sum()
            actual = artifact['models'][arm]
            for penalty, score in scores.items():
                self.assertAlmostEqual(actual['inner_brier_by_penalty'][str(penalty)], score, places=14)
            best = min(scores.values())
            self.assertEqual(actual['selected_penalty'], max(p for p, v in scores.items() if v <= best + 1e-12))

    def test_psu_overlap_is_rejected_even_for_zero_weight_training_psu(self):
        args = list(synthetic())
        args[3][:5] = 0.; args[-1][0] = args[4][0]
        with patch.object(numerical, 'fit_preprocessing') as fitting:
            with self.assertRaisesRegex(ValueError, 'share PSU'): m.fit_forward(*args)
            fitting.assert_not_called()
            with self.assertRaisesRegex(ValueError, 'share PSU'):
                m.fit_training_arm(*args, 'routine')

    def test_single_arm_bootstrap_helper_matches_baseline_each_arm(self):
        args = synthetic()
        predictions, artifact = m.fit_forward(*args)
        for arm in m.ARMS:
            p, row = m.fit_training_arm(*args, arm)
            np.testing.assert_array_equal(p, predictions[arm])
            self.assertEqual(row, artifact['models'][arm])
        with self.assertRaisesRegex(ValueError, 'Unregistered'):
            m.fit_training_arm(*args, 'best_posthoc_arm')

    def test_constant_features_tie_to_larger_ridge_and_weighted_prevalence(self):
        args = list(synthetic()); args[0][:] = 0.; args[1][:] = 1.
        predictions, artifact = m.fit_forward(*args)
        prevalence = np.average(args[2], weights=args[3])
        for arm in m.ARMS:
            row = artifact['models'][arm]
            self.assertEqual(row['selected_penalty'], 100.)
            np.testing.assert_array_equal(row['fit']['coefficients'], np.zeros(len(row['feature_names'])))
            np.testing.assert_allclose(predictions[arm], prevalence, atol=1e-10, rtol=0)

    def test_weight_rescaling_preserves_selection_and_predictions(self):
        args = list(synthetic()); p, a = m.fit_forward(*args)
        args[3] = 13. * args[3]; q, b = m.fit_forward(*args)
        for arm in m.ARMS:
            np.testing.assert_allclose(p[arm], q[arm], atol=1e-12, rtol=1e-11)
            self.assertEqual(a['models'][arm]['selected_penalty'], b['models'][arm]['selected_penalty'])

    def test_zero_weight_values_do_not_determine_full_training_scaling(self):
        args = list(synthetic()); args[3][:10] = 0.; args[0][:10, 8] = 1e6
        _, artifact = m.fit_forward(*args)
        self.assertEqual(artifact['training']['positive_weight_n'], len(args[2]) - 10)
        for row in artifact['models'].values():
            self.assertAlmostEqual(row['preprocessing']['mean'][8], np.average(args[0][:, 8], weights=args[3]))
            self.assertAlmostEqual(row['fit']['normalized_weight_sum'], len(args[2]))

    def test_unused_protein_ignores_late_values_and_single_arm_training_values(self):
        args = synthetic()
        for arm, unused in [('routine', (0, 1)), ('routine_b2m', (0,)), ('routine_cys', (1,))]:
            p, row = m.fit_training_arm(*args, arm)
            changed = [a.copy() for a in args]
            changed[1][:, unused] += 1e3; changed[6][:, unused] -= 1e3
            q, other = m.fit_training_arm(*changed, arm)
            np.testing.assert_array_equal(p, q); self.assertEqual(row, other)

    def test_inputs_old_registries_and_deterministic_artifact_unchanged(self):
        args = synthetic(); originals = [a.copy() for a in args]
        registry = copy.deepcopy(numerical.MODEL_INDICES)
        p, a = m.fit_forward(*args); q, b = m.fit_forward(*args)
        self.assertEqual(json.dumps(a, allow_nan=False), json.dumps(b, allow_nan=False))
        self.assertEqual(numerical.MODEL_INDICES, registry)
        for actual, original in zip(args, originals): np.testing.assert_array_equal(actual, original)
        for arm in m.ARMS: np.testing.assert_array_equal(p[arm], q[arm])

    def test_invalid_feature_rows_outcomes_weights_and_clusters_fail(self):
        good = synthetic(); invalid = []
        for index, value in [(0, good[0][:, :10]), (1, good[1][:, :1]),
                             (2, np.ones(len(good[2]))), (3, np.zeros(len(good[3]))),
                             (4, good[4][:-1]), (4, np.array(['only'] * len(good[4]))),
                             (5, good[5][:, :10]), (6, good[6][:, :1]),
                             (7, np.array([None] * len(good[7])))]:
            args = list(good); args[index] = value; invalid.append(args)
        for index, position, value in [(0, (0, 0), np.nan), (1, (0, 1), np.inf),
                                       (2, 0, 2), (3, 0, -1.), (5, (0, 0), np.nan)]:
            args = [a.copy() for a in good]; args[index][position] = value; invalid.append(args)
        for args in invalid:
            with self.subTest(shapes=[a.shape for a in args]):
                with self.assertRaises((ValueError, ArithmeticError)): m.fit_forward(*args)


class FixedPredictTests(unittest.TestCase):
    def test_predict_only_replays_without_any_fit_or_preprocessing_learning(self):
        args = synthetic(); expected, artifact = m.fit_forward(*args)
        artifact = json.loads(json.dumps(artifact, allow_nan=False))
        with patch.object(numerical, 'fit_ridge', side_effect=AssertionError('fit forbidden')), \
             patch.object(numerical, 'fit_preprocessing', side_effect=AssertionError('learning forbidden')):
            actual = m.predict_forward(args[5], args[6], artifact)
        for arm in m.ARMS: np.testing.assert_array_equal(actual[arm], expected[arm])

    def test_incomplete_or_mismatched_artifact_rejected(self):
        args = synthetic(); _, original = m.fit_forward(*args)
        for mutation in ('schema', 'missing_arm', 'model_name', 'protein_indices', 'penalty'):
            artifact = copy.deepcopy(original)
            if mutation == 'schema': artifact['schema'] = 'other/1'
            if mutation == 'missing_arm': del artifact['models']['routine']
            if mutation == 'model_name': artifact['models']['routine']['preprocessing']['model_name'] = 'routine_b2m'
            if mutation == 'protein_indices': artifact['models']['routine_b2m']['preprocessing']['protein_indices'] = [0]
            if mutation == 'penalty': artifact['models']['routine']['fit']['penalty'] = -1.
            with self.subTest(mutation=mutation):
                with self.assertRaises(ValueError): m.predict_forward(args[5], args[6], artifact)

    def test_artifact_contains_only_parameters_and_aggregate_diagnostics(self):
        _, artifact = m.fit_forward(*synthetic())
        encoded = json.dumps(artifact, allow_nan=False)
        for forbidden in ('early_synthetic_psu_', 'late_synthetic_psu_', 'predictions',
                          'training_auc', 'training_brier', 'late_auc', 'late_brier', 'records'):
            if forbidden == 'records':
                self.assertNotIn('"records"', encoded)
            else:
                self.assertNotIn(forbidden, encoded)
        self.assertFalse(artifact['late_domain']['outcomes_used_for_fitting'])
        self.assertFalse(artifact['late_domain']['weights_used_for_fitting'])
        self.assertTrue(artifact['no_outer_cross_validation'])


if __name__ == '__main__': unittest.main()
