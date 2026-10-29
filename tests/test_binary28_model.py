"""Analytic and synthetic-only tests for new binary28 numerical methods."""
import unittest
from unittest.mock import patch
import numpy as np
from research.mortality.binary28 import model as m


def synthetic(seed=2828, n=300):
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(n, 6)); protein = rng.normal(size=(n, 3))
    outcome = (rng.random(n) < .3).astype(int)
    weights = rng.lognormal(0, .5, n)
    clusters = np.array([f'psu{i//5}' for i in range(n)])
    outer = m.psu_partition(clusters, 'binary28:outer', 5)
    return core, protein, outcome, weights, clusters, outer


class LogisticTests(unittest.TestCase):
    def test_gradient_and_hessian_agree_with_independent_central_differences(self):
        rng = np.random.default_rng(28)
        x = np.column_stack([np.ones(25), rng.normal(size=(25, 3))])
        y = np.arange(25) % 2; w = rng.uniform(.5, 2, 25)
        beta = np.array([.2, -.5, .1, .3]); penalty = 3.
        loss, gradient, hessian = m.logistic_objective(beta, x, y, w, penalty)
        epsilon = 1e-5; numeric_gradient = []; numeric_hessian = []
        for j in range(len(beta)):
            shift = np.zeros(len(beta)); shift[j] = epsilon
            plus = m.logistic_objective(beta+shift, x, y, w, penalty)
            minus = m.logistic_objective(beta-shift, x, y, w, penalty)
            numeric_gradient.append((plus[0]-minus[0])/(2*epsilon))
            numeric_hessian.append((plus[1]-minus[1])/(2*epsilon))
        np.testing.assert_allclose(gradient, numeric_gradient, atol=1e-8, rtol=1e-8)
        np.testing.assert_allclose(hessian, np.array(numeric_hessian).T, atol=1e-8, rtol=1e-8)
        self.assertTrue(np.linalg.eigvalsh(hessian).min() > 0)
        self.assertTrue(np.isfinite(loss))

    def test_intercept_is_unpenalized_and_recovers_weighted_prevalence(self):
        y = np.array([0., 0., 1., 1., 0.]); w = np.array([1., 3., 4., 2., 1.])
        x = np.zeros((5, 4))
        for penalty in m.PENALTIES:
            fit = m.fit_ridge(x, y, w, penalty)
            np.testing.assert_allclose(m.predict_ridge(x, fit), np.average(y, weights=w), atol=1e-12)
            np.testing.assert_array_equal(fit['coefficients'], np.zeros(4))
            self.assertLessEqual(fit['gradient_max'], 1e-6)

    def test_global_weight_scale_cannot_change_fit(self):
        core, _, y, w, _, _ = synthetic(n=100)
        a = m.fit_ridge(core, y, w, 10.)
        b = m.fit_ridge(core, y, 1000*w, 10.)
        np.testing.assert_allclose(a['coefficients'], b['coefficients'], atol=1e-12)
        self.assertAlmostEqual(a['normalized_weight_sum'], len(y))

    def test_zero_weight_rows_remain_in_declared_normalization(self):
        core, _, y, w, _, _ = synthetic(n=100); w[:20] = 0.
        fit = m.fit_ridge(core, y, w, 10.)
        self.assertAlmostEqual(fit['normalized_weight_sum'], 100.)

    def test_single_class_or_invalid_weights_do_not_create_events(self):
        x = np.zeros((5, 2))
        for y, w in [(np.zeros(5), np.ones(5)), (np.array([0, 0, 0, 1, 1]), np.array([1, 1, 1, 0, 0])),
                     (np.array([0, 0, 0, 1, 1]), np.array([1, 1, 1, -1, 1]))]:
            with self.assertRaises(ValueError): m.fit_ridge(x, y, w, 1.)

    def test_missing_features_are_not_imputed(self):
        x = np.zeros((5, 2)); x[1, 1] = np.nan
        with self.assertRaises(ValueError): m.fit_ridge(x, [0, 0, 0, 1, 1], np.ones(5), 1.)

    def test_nonfinite_fitted_intercept_is_not_a_valid_probability_model(self):
        with self.assertRaises(ValueError):
            m.predict_ridge(np.zeros((2, 1)), {'intercept': np.inf, 'coefficients': [0.]})

    def test_nonfinite_objective_is_rejected_explicitly(self):
        with np.errstate(over='ignore', invalid='ignore'):
            with self.assertRaises(ArithmeticError):
                m.logistic_objective(np.array([0., 1e308]), np.array([[1., 1e308]]), np.array([1.]), np.array([1.]), 1.)


class PreprocessingTests(unittest.TestCase):
    def data(self):
        core = np.zeros((4, 6)); core[:, 0] = [0., 1., 2., 1000.]
        proteins = np.column_stack([np.array([1., 2., 2., 1000.])]*3)
        return core, proteins, np.array([1., 2., 1., 0.])

    def test_weighted_ecdf_is_right_continuous_and_ignores_zero_weight(self):
        core, proteins, weights = self.data()
        fit = m.fit_preprocessing(core, proteins, weights, 'rank_all3')
        support = fit['protein_transforms'][0]
        np.testing.assert_array_equal(support['values'], [1., 2.])
        np.testing.assert_array_equal(support['cdf'], [.25, 1.])
        new_core = np.zeros((5, 6)); new_proteins = np.tile(np.array([.5, 1., 1.5, 2., 5000.])[:, None], (1, 3))
        unscaled = m._unscaled(new_core, new_proteins, fit)
        np.testing.assert_array_equal(unscaled[:, 6], [0., .25, .25, 1., 1.])

    def test_rank_representation_is_invariant_to_strict_monotonic_recalibration(self):
        core, proteins, w, *_ = synthetic(n=100)
        fit = m.fit_preprocessing(core, proteins, w, 'rank_all3')
        recalibrated = np.exp(.4*proteins)+10
        other = m.fit_preprocessing(core, recalibrated, w, 'rank_all3')
        np.testing.assert_allclose(m.transform_features(core, proteins, fit),
                                   m.transform_features(core, recalibrated, other), atol=0, rtol=0)

    def test_clip_limits_use_weighted_midpoint_cdf_and_no_scoring_values(self):
        core = np.zeros((100, 6)); proteins = np.tile(np.arange(100.)[:, None], (1, 3)); w = np.ones(100)
        fit = m.fit_preprocessing(core, proteins, w, 'clip_all3')
        spec = fit['protein_transforms'][0]
        self.assertAlmostEqual(spec['lower'], .5); self.assertAlmostEqual(spec['upper'], 98.5)
        test = np.tile(np.array([-1e6, 1e6])[:, None], (1, 3))
        values = m._unscaled(np.zeros((2, 6)), test, fit)[:, 6]
        np.testing.assert_allclose(values, [.5, 98.5])
        self.assertEqual(spec, fit['protein_transforms'][0])

    def test_mean_scale_are_training_weighted_and_constants_flagged(self):
        core, proteins, w = self.data()
        fit = m.fit_preprocessing(core, proteins, w, 'core')
        self.assertEqual(fit['mean'][0], 1.)
        self.assertAlmostEqual(fit['scale'][0], np.sqrt(.5))
        self.assertEqual(fit['constant'], [False, True, True, True, True, True])
        self.assertEqual(fit['scale'][1:], [1.]*5)
        np.testing.assert_array_equal(m.transform_features(core, proteins, fit)[:, 1:], np.zeros((4, 5)))

    def test_preprocessing_does_not_change_input_matrices(self):
        core, proteins, w = self.data(); old_core = core.copy(); old_proteins = proteins.copy()
        for name in m.MODEL_INDICES:
            spec = m.fit_preprocessing(core, proteins, w, name)
            m.transform_features(core, proteins, spec)
        np.testing.assert_array_equal(core, old_core); np.testing.assert_array_equal(proteins, old_proteins)


class NestedCVTests(unittest.TestCase):
    def test_cluster_partition_is_outcome_blind_and_row_order_invariant(self):
        _, _, _, _, clusters, _ = synthetic()
        a = m.psu_partition(clusters, 'binary28:outer', 5)
        np.testing.assert_array_equal(m.psu_partition(clusters[::-1], 'binary28:outer', 5), a[::-1])
        for cluster in set(clusters): self.assertEqual(len(set(a[clusters == cluster])), 1)

    def test_psu_cannot_cross_outer_boundary(self):
        core, proteins, y, w, clusters, folds = synthetic()
        folds[1] = (folds[1]+1) % 5
        with self.assertRaisesRegex(ValueError, 'PSU crosses'):
            m.nested_crossfit(core, proteins, y, w, clusters, folds)

    def test_heldout_outcomes_cannot_affect_own_predictions_or_lambda(self):
        core, proteins, y, w, clusters, folds = synthetic()
        p, fits = m.nested_crossfit(core, proteins, y, w, clusters, folds, 'all3')
        changed = y.copy(); changed[folds == 0] = 1-changed[folds == 0]
        q, other = m.nested_crossfit(core, proteins, changed, w, clusters, folds, 'all3')
        np.testing.assert_allclose(p[folds == 0], q[folds == 0], atol=0, rtol=0)
        for key in ('selected_penalty', 'training_mean', 'training_scale', 'inner_brier_by_penalty'):
            self.assertEqual(fits[0][key], other[0][key])
        self.assertEqual(len(fits), 5)
        self.assertTrue(all(r['PSU_overlap_n'] == 0 for r in fits))

    def test_heldout_features_never_enter_first_outer_training_preprocessing(self):
        core, proteins, y, w, clusters, folds = synthetic()
        calls = []; actual = m.fit_preprocessing
        # Distinct integer row-key feature enables inspecting selected training rows.
        core[:, 0] = np.arange(len(core))
        def capture(c, p, weights, name):
            calls.append(set(c[:, 0].astype(int)))
            return actual(c, p, weights, name)
        with patch.object(m, 'fit_preprocessing', side_effect=capture):
            m.nested_crossfit(core, proteins, y, w, clusters, folds, 'all3')
        heldout = set(np.flatnonzero(folds == 0))
        # Three inner preprocessing fits + one outer fit, with cached transforms.
        self.assertEqual(len(calls), 20)
        for call in calls[:4]: self.assertFalse(call & heldout)

    def test_constant_added_proteins_cannot_artificially_change_risks(self):
        core, protein, y, w, clusters, folds = synthetic()
        protein[:] = 2.
        p, _ = m.nested_crossfit(core, protein, y, w, clusters, folds, 'core')
        q, _ = m.nested_crossfit(core, protein, y, w, clusters, folds, 'all3')
        np.testing.assert_allclose(p, q, atol=1e-12, rtol=0)

    def test_tied_inner_brier_prefers_larger_penalty(self):
        core, protein, y, w, clusters, folds = synthetic()
        core[:] = 0.; protein[:] = 0.
        _, fits = m.nested_crossfit(core, protein, y, w, clusters, folds, 'core')
        self.assertTrue(all(row['selected_penalty'] == max(m.PENALTIES) for row in fits))

    def test_all_prespecified_arms_export_valid_probabilities_and_scalar_diagnostics(self):
        args = synthetic(n=200)
        for name in m.MODEL_INDICES:
            with self.subTest(model=name):
                p, fits = m.nested_crossfit(*args, model_name=name)
                self.assertTrue(np.isfinite(p).all()); self.assertTrue(((p >= 0) & (p <= 1)).all())
                self.assertTrue(all(row['gradient_max'] <= 1e-4 for row in fits))
                self.assertTrue(all(row['training_positive_weight_events'] > 0 for row in fits))
                self.assertTrue(all(row['feature_names'][:6] == list(m.CORE_NAMES) for row in fits))


if __name__ == '__main__': unittest.main()
