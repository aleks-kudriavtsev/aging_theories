"""Independent numerical and censoring checks for the frozen temporal assessment."""
import copy
import json
import tempfile
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from research.mortality.temporal26 import analyze as a


class AssayBridgeTests(unittest.TestCase):
    def source(self):
        return pd.DataFrame({'LBXSAL': [4., 3.], 'LBXSAPSI': [100., 50.],
                             'LBXSCR': [1., .7], 'LBXSGTSI': [30., 20.],
                             'LBXSTB': [.8, .5], 'untouched': [7., 8.]})

    def test_backward_bridge_uses_published_direction(self):
        result = a.bridge(self.source())
        self.assertAlmostEqual(result.LBXSAL.iloc[0], 4.18728)
        self.assertAlmostEqual(result.LBXSCR.iloc[0], .98155)
        self.assertAlmostEqual(result.LBXSGTSI.iloc[0], 26.489)
        self.assertAlmostEqual(result.LBXSAPSI.iloc[0], 10 ** 1.95906)
        np.testing.assert_array_equal(result.LBXSTB, [.8, .5])
        np.testing.assert_array_equal(result.untouched, [7., 8.])

    def test_bridge_does_not_mutate_released_measurements(self):
        source = self.source(); original = source.copy(deep=True)
        a.bridge(source)
        pd.testing.assert_frame_equal(source, original)

    def test_nonpositive_source_and_corrected_values_remain_missing(self):
        source = self.source()
        source.loc[0, ['LBXSAL', 'LBXSAPSI', 'LBXSGTSI']] = [0., -1., 0.]
        source.loc[0, 'LBXSCR'] = .02  # Positive source, negative corrected result.
        result = a.bridge(source)
        self.assertTrue(result.loc[0, ['LBXSAL', 'LBXSAPSI', 'LBXSCR', 'LBXSGTSI']].isna().all())


class FrozenHazardTests(unittest.TestCase):
    def model(self, logarithmic=False):
        model = {'preprocessing': {'features': ['x'], 'lower': [-2.], 'upper': [2.],
                                 'mean': [0.], 'sd': [1.]},
                 'interval_ends_years': [5., 8.],
                 'causes': {'1': {'coefficients': [np.log(2.)]},
                            '10': {'coefficients': [-np.log(2.)]}}}
        for key, rates in [('1', [.01, .04]), ('10', [.02, .08])]:
            model['causes'][key]['baseline_log_hazards' if logarithmic else 'baseline_rates'] = (
                np.log(rates).tolist() if logarithmic else rates)
        return model

    def test_total_risk_matches_analytic_first_and_second_intervals(self):
        frame = pd.DataFrame({'x': [-1., 0., 1.]})
        for horizon in (1., 6.):
            exposure = min(horizon, 5.) + 4 * max(horizon - 5., 0.)
            total = exposure * (.01 * 2 ** frame.x.to_numpy() + .02 * 2 ** -frame.x.to_numpy())
            expected = 1 - np.exp(-total)
            np.testing.assert_allclose(a.total_risk(frame, self.model(), horizon), expected, atol=1e-15, rtol=1e-13)

    def test_rate_and_log_hazard_encodings_are_numerically_identical(self):
        frame = pd.DataFrame({'x': [-1., 0., 1.]})
        np.testing.assert_allclose(a.total_risk(frame, self.model(), 1),
                                   a.total_risk(frame, self.model(True), 1), atol=1e-15)

    def test_cause_dictionary_order_cannot_change_total_risk(self):
        frame = pd.DataFrame({'x': [-1., 0., 1.]})
        model = self.model(); reverse = copy.deepcopy(model)
        reverse['causes'] = dict(reversed(list(model['causes'].items())))
        np.testing.assert_allclose(a.total_risk(frame, model, 1), a.total_risk(frame, reverse, 1), atol=1e-15)

    def test_frozen_winsorization_is_retained(self):
        frame = pd.DataFrame({'x': [-20., 20.]})
        bounded = pd.DataFrame({'x': [-2., 2.]})
        np.testing.assert_allclose(a.total_risk(frame, self.model(), 1), a.total_risk(bounded, self.model(), 1), atol=0)

    def test_missing_predictor_cannot_be_imputed_silently(self):
        with self.assertRaises(ValueError):
            a.total_risk(pd.DataFrame({'x': [np.nan]}), self.model(), 1)


class IPCWAssessmentTests(unittest.TestCase):
    def data(self):
        # Independent weighted KM: .25 + .75*(1/5) = .40 by one year.
        return pd.DataFrame({'dead': [1., 0., 1., 0.], 'cause': [1., np.nan, 10., np.nan],
                             'time_years': [.25, .5, 1., 1.]})

    def predictions(self, n=4):
        return {'expanded13_bridge': np.full(n, .15), 'compact4_bridge': np.full(n, .2),
                'expanded13_raw': np.full(n, .13), 'compact4_raw': np.full(n, .18)}

    def test_ipcw_observed_incidence_matches_manual_weighted_km(self):
        result = a.assess(self.data(), self.predictions(), np.array([2., 1., 1., 4.]))
        for metrics in result['metrics'].values():
            self.assertAlmostEqual(metrics['observed'], .4)
            self.assertEqual(metrics['early_censored_n'], 1)
            self.assertEqual(metrics['events'], 2)

    def test_censored_record_is_not_scored_as_survivor(self):
        predictions = self.predictions(); altered = {k: v.copy() for k, v in predictions.items()}
        for v in altered.values(): v[1] = .99
        weights = np.array([2., 1., 1., 4.])
        first = a.assess(self.data(), predictions, weights)
        second = a.assess(self.data(), altered, weights)
        for key in predictions:
            self.assertAlmostEqual(first['metrics'][key]['brier'], second['metrics'][key]['brier'])
            self.assertNotEqual(first['metrics'][key]['predicted'], second['metrics'][key]['predicted'])

    def test_primary_contrast_uses_same_domain_and_fixed_bridge_labels(self):
        result = a.assess(self.data(), self.predictions(), np.array([2., 1., 1., 4.]))
        # E[(Y-.15)^2]-E[(Y-.2)^2] = .1*.4 - .0175 = .0225.
        self.assertAlmostEqual(result['delta_brier'], .0225)
        self.assertEqual(result['metrics']['expanded13_bridge']['n'], 4)
        self.assertEqual(result['metrics']['compact4_bridge']['n'], 4)

    def test_unequal_prediction_lengths_cannot_change_comparison_population(self):
        predictions = self.predictions(); predictions['compact4_bridge'] = np.array([.2, .2])
        with self.assertRaises(ValueError):
            a.assess(self.data(), predictions, np.ones(4))

    def test_unobserved_horizon_fails_instead_of_inventing_event_labels(self):
        data = pd.DataFrame({'dead': [0., 0.], 'cause': [np.nan, np.nan], 'time_years': [.3, .5]})
        with self.assertRaises(ValueError):
            a.assess(data, self.predictions(2), np.ones(2))

    def test_no_events_has_undefined_auc_and_defined_brier(self):
        data = pd.DataFrame({'dead': [0., 0.], 'cause': [np.nan, np.nan], 'time_years': [1., 2.]})
        result = a.assess(data, self.predictions(2), np.ones(2))
        self.assertIsNone(result['delta_auc'])
        self.assertIsNone(result['metrics']['expanded13_bridge']['auc'])
        self.assertAlmostEqual(result['metrics']['expanded13_bridge']['brier'], .15 ** 2)


class BootstrapSummaryTests(unittest.TestCase):
    def test_failures_and_auc_undefined_are_counted_separately(self):
        draws = [{'delta_brier': .01, 'delta_auc': .02},
                 {'delta_brier': .03, 'delta_auc': None},
                 {'status': 'undefined'}]
        result = a.summarize(draws, 3)
        self.assertEqual(result['delta_brier']['valid'], 2)
        self.assertEqual(result['delta_brier']['undefined'], 1)
        self.assertEqual(result['delta_auc']['valid'], 1)
        self.assertEqual(result['delta_auc']['undefined'], 2)
        self.assertTrue(result['delta_brier']['defined_draws_only'])
        self.assertAlmostEqual(result['delta_brier']['lower'], .0105)
        self.assertAlmostEqual(result['delta_brier']['upper'], .0295)

    def test_all_failures_have_no_numerical_interval(self):
        result = a.summarize([{'status': 'undefined'}] * 3, 3)
        for key in ('delta_brier', 'delta_auc'):
            self.assertEqual(result[key]['valid'], 0)
            self.assertEqual(result[key]['undefined'], 3)
            self.assertIsNone(result[key]['lower'])
            self.assertIsNone(result[key]['upper'])


class SourceGateTests(unittest.TestCase):
    def manifest(self):
        return {'analysis_lock_commit': '5fbb68c',
                'lock_sha256': a.base.sha(a.ROOT / 'LOCK.json'), 'files': []}

    def test_missing_required_sources_prevent_any_result_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp); out = source / 'result'
            (source / 'manifest.json').write_text(json.dumps(self.manifest()))
            with self.assertRaisesRegex(ValueError, 'Missing required source'):
                a.run(source, out)
            self.assertFalse(out.exists())

    def test_manifest_cannot_rebind_sources_to_another_protocol(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp); manifest = self.manifest()
            manifest['lock_sha256'] = '0' * 64
            (source / 'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'locked protocol'):
                a.verify_sources(source)

    def test_recorded_hash_is_verified_before_any_data_parsing(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp); manifest = self.manifest()
            name, url = next(iter(a.source_urls().items()))
            (source / name).write_bytes(b'invalid source')
            manifest['files'] = [dict(file=name, url=url, status='retrieved',
                                      bytes=14, sha256='0' * 64)]
            (source / 'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'Source changed'):
                a.verify_sources(source)

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp); marker = output / 'old.txt'; marker.write_text('original')
            with self.assertRaisesRegex(ValueError, 'New output'):
                a.run(output / 'missing', output)
            self.assertEqual(marker.read_text(), 'original')


if __name__ == '__main__':
    unittest.main()
