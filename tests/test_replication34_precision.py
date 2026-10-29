"""Pure arithmetic/toy fixtures, never production scenarios or model results."""
import copy
import json
import math
import unittest

import numpy as np
from scipy.stats import binom

from research.mortality.replication34 import precision as p


def scenario(**changes):
    value = dict(validation_mode='algorithm_replication', mu_d=0., sd_d=.05,
        alpha=.05, half_width=.002, n_evaluable=1000, deff_d=1.5,
        deff_basis='total_relative_to_iid_unweighted_n', q_joint=.75,
        q_known_conditional=1., censoring_fraction=0.,
        raw_evaluable_event_probability=.03, event_count_model='independent_binomial',
        deff_event=1., technical_event_floor=50, id='SYNTHETIC_TECHNICAL_FIXTURE')
    value.update(changes)
    return value


def atoms():
    return [dict(q=.5, r=.1, p_new=.2, p_ref=.1),
            dict(q=.5, r=.8, p_new=.7, p_ref=.6)]


class PrecisionTests(unittest.TestCase):
    def test_hand_calculated_conditional_precision_and_ceiling(self):
        result = p.plan(scenario())
        self.assertAlmostEqual(result['conditional_mean_variance'], .00000375)
        self.assertAlmostEqual(result['conditional_SE'], .0019364916731037084)
        self.assertAlmostEqual(result['normal_quantile'], 1.959963984540054)
        self.assertEqual(result['conditional_n_for_half_width'], 3602)
        self.assertEqual(result['candidate_recruitment_expectation'], 4803)
        self.assertEqual(result['event_count']['expected_raw_events'], 30.)

    def test_precision_uses_paired_variance_not_mean_effect_or_event_rate(self):
        reference = p.plan(scenario())
        for mean in (-.01, 0., .01):
            for incidence in (.001, .03, .7):
                result = p.plan(scenario(mu_d=mean, raw_evaluable_event_probability=incidence))
                self.assertEqual(result['conditional_n_for_half_width'], reference['conditional_n_for_half_width'])
                self.assertFalse(result['mean_effect_used_to_choose_n'])
        tight = p.plan(scenario(half_width=.001))
        doubled = p.plan(scenario(deff_d=3.))
        self.assertTrue(4 * 3602 - 3 <= tight['conditional_n_for_half_width'] <= 4 * 3602)
        self.assertTrue(2 * 3602 - 1 <= doubled['conditional_n_for_half_width'] <= 2 * 3602)

    def test_anchored_SE_recovers_individual_loss_variance(self):
        value = scenario()
        del value['sd_d']
        value.update(se_reference=.001, n_reference=2000, deff_reference=2.)
        result = p.plan(value)
        self.assertAlmostEqual(result['paired_loss_variance'], .001)
        value['sd_d'] = .05
        with self.assertRaisesRegex(ValueError, 'not both'): p.plan(value)
        del value['sd_d']; del value['deff_reference']
        with self.assertRaisesRegex(ValueError, 'both reference'): p.plan(value)

    def test_moment_bounds_and_zero_variance_are_explicit(self):
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            p.plan(scenario(mu_d=.9, sd_d=.8))
        zero = p.plan(scenario(mu_d=1., sd_d=0.))
        self.assertTrue(zero['degenerate_zero_loss_variance'])
        self.assertEqual(zero['conditional_SE'], 0.)
        self.assertEqual(zero['conditional_n_for_half_width'], 1)
        self.assertFalse(zero['clinical_adequacy_or_power_established'])

    def test_no_algorithm_CIF_AUC_or_calibration_adequacy_gate(self):
        for mode in ('algorithm_replication', 'fixed_scorer'):
            result = p.plan(scenario(validation_mode=mode))
            self.assertEqual(result['algorithm_precision_status'], 'adequacy_unqualified')
            self.assertFalse(result['full_CIF_validation_eligible'])
            self.assertFalse(result['algorithm_training_uncertainty_included'])
            self.assertFalse(result['calibration_or_AUC_precision_established'])
            self.assertFalse(result['clinical_adequacy_or_power_established'])
            self.assertEqual(result['frozen_source_artifact_required_before_external_outcomes'], mode == 'fixed_scorer')
        for censoring in (.000001, .1, 1.):
            result = p.plan(scenario(censoring_fraction=censoring))
            self.assertTrue(result['IPCW_method_extension_required'])
            self.assertFalse(result['bounded_uncensored_loss_formula_applies_to_full_CIF'])
            self.assertFalse(result['full_CIF_validation_eligible'])

    def test_endpoint_and_near_endpoint_variance_cannot_be_clipped_into_feasibility(self):
        for mean in (-1., 1., float(np.nextafter(1., 0.))):
            with self.subTest(mean=mean), self.assertRaisesRegex(ValueError, 'incompatible'):
                p.plan(scenario(mu_d=mean, sd_d=1e-7))
            value = scenario(mu_d=mean)
            del value['sd_d']
            value.update(se_reference=1e-7, n_reference=1000, deff_reference=1.)
            with self.assertRaisesRegex(ValueError, 'incompatible'):
                p.plan(value)
        with self.assertRaisesRegex(ValueError, 'underflows'):
            p.plan(scenario(sd_d=1e-300))
        feasible = p.plan(scenario(mu_d=.6, sd_d=.8))
        self.assertAlmostEqual(feasible['paired_loss_variance'], .64)

    def test_recruitment_is_joint_conditional_expectation_not_assurance(self):
        result = p.plan(scenario(q_joint=.5, q_known_conditional=.8))
        self.assertEqual(result['joint_evaluable_yield_probability'], .4)
        self.assertEqual(result['candidate_recruitment_expectation'], math.ceil(3602 / .4))
        self.assertTrue(result['recruitment_is_expectation_not_assurance'])
        self.assertEqual(p.plan(scenario(q_joint=1., q_known_conditional=1.))['candidate_recruitment_expectation'], 3602)
        for field in ('q_joint', 'q_known_conditional'):
            unavailable = p.plan(scenario(**{field: 0.}))
            self.assertIsNone(unavailable['candidate_recruitment_expectation'])
            self.assertFalse(unavailable['recruitment_target_reachable_under_assumed_availability'])

    def test_binomial_quantiles_use_raw_N_and_separate_event_DEFF(self):
        a = p.plan(scenario(deff_d=1.))['event_count']
        b = p.plan(scenario(deff_d=3.))['event_count']
        self.assertEqual(a, b)
        self.assertEqual(a['n_raw_trials'], 1000)
        self.assertEqual(a['predictive_range'], [int(binom.ppf(.025, 1000, .03)), int(binom.ppf(.975, 1000, .03))])
        self.assertAlmostEqual(a['probability_below_technical_floor'], binom.cdf(49, 1000, .03))
        self.assertFalse(a['technical_floor_establishes_power_or_adequacy'])
        with self.assertRaisesRegex(ValueError, 'independence'): p.plan(scenario(deff_event=2.))
        dependent = p.plan(scenario(event_count_model='dependent_unspecified', deff_event=2.))['event_count']
        self.assertEqual(dependent['expected_raw_events'], 30.)
        self.assertAlmostEqual(dependent['variance_proxy'], 58.2)
        self.assertIsNone(dependent['predictive_range'])
        self.assertIsNone(dependent['probability_below_technical_floor'])

    def test_binomial_probability_boundaries_and_CDF_inverse(self):
        for incidence, expected in ((0., [0, 0]), (1., [1000, 1000])):
            result = p.plan(scenario(raw_evaluable_event_probability=incidence))['event_count']
            self.assertEqual(result['predictive_range'], expected)
        result = p.plan(scenario(technical_event_floor=0))['event_count']
        self.assertEqual(result['probability_below_technical_floor'], 0.)
        result = p.plan(scenario())['event_count']
        for k, probability in zip(result['predictive_range'], (.025, .975)):
            self.assertGreaterEqual(binom.cdf(k, 1000, .03), probability - 1e-12)
            self.assertLess(binom.cdf(k - 1, 1000, .03), probability + 1e-12)

    def test_schema_nonfinite_bool_counts_and_double_counting_rejected(self):
        mutations = [('mu_d', True), ('sd_d', float('nan')), ('alpha', 0.),
            ('half_width', -1.), ('n_evaluable', 1000.), ('n_evaluable', True),
            ('deff_d', 0.), ('deff_basis', 'multiply_Kish_and_total_DEFF'),
            ('q_joint', 1.1), ('q_known_conditional', -.1), ('censoring_fraction', float('inf')),
            ('technical_event_floor', True), ('event_count_model', 'binomial_using_n_over_DEFF')]
        for key, value in mutations:
            with self.subTest(key=key), self.assertRaises(ValueError): p.plan(scenario(**{key: value}))
        with self.assertRaises(ValueError): p.plan(scenario(kish_effective_n=100.))
        value = scenario(); del value['mu_d']
        with self.assertRaises(ValueError): p.plan(value)
        with self.assertRaises(ValueError): p.plan(scenario(half_width=1e-300))


class AtomsTests(unittest.TestCase):
    def test_hand_calculated_moments_and_total_variance(self):
        result = p.finite_prediction_atoms_moments(atoms())
        self.assertAlmostEqual(result['raw_event_probability'], .45)
        self.assertAlmostEqual(result['mean_d'], -.01)
        self.assertAlmostEqual(result['second_moment_d'], .0055)
        self.assertAlmostEqual(result['variance_d'], .0054)
        self.assertAlmostEqual(result['mean_conditional_Y_variance'], .005)
        self.assertAlmostEqual(result['variance_of_conditional_means'], .0004)
        self.assertLess(result['variance_decomposition_error'], 1e-15)
        self.assertTrue(result['predictions_fixed_before_Y'])

    def test_sign_symmetry_identical_rules_and_prediction_before_Y(self):
        original = atoms(); switched = [dict(row, p_new=row['p_ref'], p_ref=row['p_new']) for row in original]
        a, b = p.finite_prediction_atoms_moments(original), p.finite_prediction_atoms_moments(switched)
        self.assertAlmostEqual(a['mean_d'], -b['mean_d'])
        self.assertAlmostEqual(a['variance_d'], b['variance_d'])
        identical = [dict(row, p_new=row['p_ref']) for row in original]
        zero = p.finite_prediction_atoms_moments(identical)
        self.assertEqual(zero['mean_d'], 0.); self.assertEqual(zero['variance_d'], 0.)

    def test_rare_incidence_does_not_bound_overall_loss_variance(self):
        values = [dict(q=.5, r=.001, p_new=1., p_ref=0.),
                  dict(q=.5, r=.001, p_new=0., p_ref=1.)]
        result = p.finite_prediction_atoms_moments(values)
        self.assertAlmostEqual(result['raw_event_probability'], .001)
        self.assertAlmostEqual(result['mean_d'], 0.)
        self.assertAlmostEqual(result['variance_d'], 1.)
        self.assertGreater(result['variance_d'], 4 * .001 * .999)

    def test_invalid_atom_fields_weights_or_predictions_fail(self):
        bad = [[], [dict(q=.5, r=.1, p_new=.2, p_ref=.1)],
               [dict(q=1., r=True, p_new=.2, p_ref=.1)],
               [dict(q=1., r=.1, p_new=1.01, p_ref=.1)],
               [dict(q=1., r=.1, p_new=.2, p_ref=.1, y=1)]]
        for values in bad:
            with self.assertRaises(ValueError): p.finite_prediction_atoms_moments(values)


class SimulationTests(unittest.TestCase):
    def test_toy_cluster_moments_agree_with_predeclared_MC_tolerance(self):
        #Small independent UNIT fixture, not any committed20k production case.
        result = p.simulate_case(atoms(), 25, 4, .2, 2000, 341)
        self.assertAlmostEqual(result['deff'], 1.6)
        self.assertAlmostEqual(result['exact_variance_mean_d'], 1.6 * .0054 / 100)
        self.assertLess(abs(result['MC_mean_d'] - result['exact_mean_d']), 6 * result['MCSE_mean_d'])
        self.assertLess(abs(result['MC_variance_mean_d'] - result['exact_variance_mean_d']),
                        6 * result['approximate_fourth_moment_MCSE_variance'])
        self.assertLess(abs(result['MC_mean_event_count'] - result['exact_expected_event_count']),
                        6 * result['MCSE_mean_event_count'])
        self.assertFalse(result['refitted_algorithm_precision_validated'])
        self.assertFalse(result['external_cohort_or_clinical_validation_performed'])

    def test_coverage_is_reported_not_forced_to95_and_has_MC_uncertainty(self):
        values = [dict(q=1., r=.001, p_new=.2, p_ref=.1)]
        result = p.simulate_case(values, 1, 1, 0., 2000, 342)
        coverage = result['known_toy_variance_normal_CI_coverage']
        self.assertEqual(coverage['replicates'], 2000)
        self.assertTrue(0 <= coverage['proportion'] <= 1)
        self.assertNotEqual(coverage['proportion'], .95)
        self.assertTrue(0 <= coverage['MC_95pct_Wilson_interval'][0] <= coverage['MC_95pct_Wilson_interval'][1] <= 1)
        self.assertIsNone(result['estimated_cluster_mean_SE_normal_CI_coverage'])

    def test_degenerate_rule_and_seed_replay_return_only_aggregates(self):
        values = [dict(q=1., r=.2, p_new=.1, p_ref=.1)]
        first = p.simulate_case(values, 3, 2, 1., 100, 343)
        second = p.simulate_case(values, 3, 2, 1., 100, 343)
        self.assertEqual(first, second)
        self.assertEqual(first['MC_mean_d'], 0.); self.assertEqual(first['MC_variance_mean_d'], 0.)
        coverage = first['known_toy_variance_normal_CI_coverage']
        self.assertEqual(coverage['proportion'], 1.)
        self.assertLess(coverage['MC_95pct_Wilson_interval'][0], 1.)
        self.assertFalse(first['participant_simulation_arrays_exported'])
        json.dumps(first, allow_nan=False)

    def test_tiny_positive_variance_has_no_absolute_coverage_padding(self):
        #Independent unit law, not a committed production scenario. The known
        #normal interval covers exactly K=1,2,3 among4Bernoulli(.5) outcomes.
        values = [dict(q=1., r=.5, p_new=.5000000000001, p_ref=.5)]
        moments = p.finite_prediction_atoms_moments(values)
        delta = values[0]['p_new'] - values[0]['p_ref']
        d0 = delta * (values[0]['p_new'] + values[0]['p_ref'])
        d1 = delta * (values[0]['p_new'] + values[0]['p_ref'] - 2)
        half_width = 1.959963984540054 * math.sqrt(moments['variance_d'] / 4)
        exact_coverage = sum(math.comb(4, k) for k in range(5)
            if abs(((4-k)*d0+k*d1)/4 - moments['mean_d']) <= half_width) / 16
        self.assertEqual(exact_coverage, .875)
        result = p.simulate_case(values, 4, 1, 0., 2000, 345)
        mcse_at_truth = math.sqrt(.875 * .125 / 2000)
        for field in ('known_toy_variance_normal_CI_coverage',
                      'estimated_cluster_mean_SE_normal_CI_coverage'):
            self.assertLess(abs(result[field]['proportion'] - exact_coverage), 6 * mcse_at_truth)
        self.assertTrue(result['positive_variance_CI_coverage_uses_no_added_tolerance'])

    def test_nonzero_constant_loss_uses_only_relative_degenerate_roundoff(self):
        values = [dict(q=1., r=0., p_new=.2, p_ref=.1)]
        result = p.simulate_case(values, 7, 3, .4, 20, 346)
        self.assertEqual(result['exact_variance_mean_d'], 0.)
        self.assertEqual(result['known_toy_variance_normal_CI_coverage']['proportion'], 1.)
        self.assertIn('assert_constant_losses', result['exact_zero_variance_handling'])

    def test_simulation_input_guard_and_no_mutation(self):
        values = atoms(); original = copy.deepcopy(values)
        p.simulate_case(values, 2, 2, .1, 3, 344)
        self.assertEqual(values, original)
        cases = [(0, 2, .1, 3, 1), (2, True, .1, 3, 1), (2, 2, 1.1, 3, 1),
                 (2, 2, .1, 1, 1), (2, 2, .1, 3, -1)]
        for args in cases:
            with self.assertRaises(ValueError): p.simulate_case(values, *args)


if __name__ == '__main__': unittest.main()
