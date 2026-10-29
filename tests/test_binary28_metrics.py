"""Analytic weighted scoring/calibration tests; no empirical participant data."""
import unittest
import numpy as np
from scipy.special import expit, logit
from unittest.mock import patch
from research.mortality.binary28.metrics import metrics, calibration_bins


class ScoreTests(unittest.TestCase):
    def test_weighted_brier_and_tied_auc_match_hand_calculation(self):
        y=np.array([0,1,0,1]);p=np.array([.2,.2,.8,.8]);w=np.array([1.,2.,3.,4.])
        r=metrics(y,p,w)
        self.assertAlmostEqual(r['brier'],(1*.04+2*.64+3*.64+4*.04)/10)
        # positive weight6, negative4: concordance4*1; tied .5*(2*1+4*3)
        self.assertAlmostEqual(r['auc'],11/24)
        self.assertAlmostEqual(r['observed'],.6)
        self.assertAlmostEqual(r['predicted'],.62)
        self.assertAlmostEqual(r['bias_pp'],2.)
        self.assertAlmostEqual(r['observed_expected'],.6/.62)

    def test_zero_weight_records_do_not_create_class_support(self):
        r=metrics([0,0,1],[.2,.3,.9],[1.,2.,0.])
        self.assertIsNone(r['auc']);self.assertEqual(r['n_positive_weight'],2)
        self.assertEqual(r['events'],1);self.assertEqual(r['events_positive_weight'],0)
        self.assertEqual(r['calibration_status'],'undefined_single_class')
        self.assertIsNone(r['calibration_slope'])
        self.assertIsNone(r['calibration_in_large'])

    def test_zero_prediction_has_no_infinite_observed_expected(self):
        r=metrics([0,1],[0.,0.],[1.,1.])
        self.assertEqual(r['brier'],.5);self.assertIsNone(r['observed_expected'])
        self.assertTrue(np.isfinite(r['calibration_in_large']))

    def test_inputs_are_unchanged_and_scaled_weights_same_scores(self):
        y=np.array([0.,1.,0.,1.]);p=np.array([.1,.7,.6,.5]);w=np.array([1.,2.,3.,4.])
        old=[a.copy() for a in (y,p,w)];r=metrics(y,p,w);q=metrics(y,p,1e7*w)
        for a,b in zip((y,p,w),old):np.testing.assert_array_equal(a,b)
        for key in ('auc','brier','observed','predicted','bias_pp'):
            self.assertAlmostEqual(r[key],q[key])

    def test_invalid_inputs_fail(self):
        for y,p,w in [([],[],[]),([0,1],[.2],[1,1]),([0,2],[.2,.4],[1,1]),
                      ([0,1],[np.nan,.4],[1,1]),([0,1],[-.1,.4],[1,1]),
                      ([0,1],[.2,1.4],[1,1]),([0,1],[.2,.4],[-1,1]),
                      ([0,1],[.2,.4],[0,0])]:
            with self.subTest(y=y,p=p,w=w):
                with self.assertRaises(ValueError):metrics(y,p,w)


class CalibrationTests(unittest.TestCase):
    def test_constant_prediction_identifies_offset_not_full_slope(self):
        r=metrics([0,0,0,1],[.4]*4,[1.]*4)
        self.assertAlmostEqual(r['calibration_in_large'],logit(.25)-logit(.4),places=7)
        self.assertIsNone(r['calibration_slope']);self.assertIsNone(r['calibration_intercept'])
        self.assertEqual(r['calibration_status'],'flat_prediction_slope_unidentified')

    def test_grouped_known_logistic_calibration_and_score_equations(self):
        p=np.repeat([.2,.5,.8],10)
        y=np.r_[np.r_[np.ones(1),np.zeros(9)],np.r_[np.ones(5),np.zeros(5)],np.r_[np.ones(9),np.zeros(1)]]
        r=metrics(y,p,np.ones(30))
        self.assertEqual(r['calibration_status'],'computed')
        self.assertAlmostEqual(r['calibration_in_large'],0.,places=7)
        self.assertAlmostEqual(r['calibration_intercept'],0.,places=7)
        self.assertAlmostEqual(r['calibration_slope'],logit(.9)/logit(.8),places=6)
        X=np.column_stack([np.ones(30),logit(p)])
        b=np.array([r['calibration_intercept'],r['calibration_slope']])
        self.assertLess(np.max(abs(X.T@(expit(X@b)-y))),1e-5)

    def test_separation_is_not_reported_as_finite_calibration(self):
        r=metrics([0,0,1,1],[.1,.3,.7,.9],[1.]*4)
        self.assertIsNone(r['calibration_slope']);self.assertIsNone(r['calibration_intercept'])
        self.assertEqual(r['calibration_status'],'separated_calibration_slope_unidentified')
        self.assertTrue(np.isfinite(r['calibration_in_large']))

    def test_nonfinite_optimizer_result_does_not_leak_nan_diagnostic(self):
        class Bad:
            x=np.array([np.nan])
        with patch('research.mortality.binary28.metrics.minimize',return_value=Bad()):
            r=metrics([0,1,0,1],[.2,.3,.6,.8],[1.]*4)
        self.assertIsNone(r['calibration_in_large'])
        self.assertEqual(r['calibration_status'],'undefined_offset_fit_nonconverged')
        self.assertTrue(np.isfinite(r['brier']))


class BinTests(unittest.TestCase):
    def test_weighted_quintiles_match_simple_hand_groups(self):
        y=np.array([0,1,0,1,0,1,0,1,0,1])
        p=np.linspace(.05,.95,10);w=np.ones(10)
        rows=calibration_bins(y,p,w)
        self.assertEqual(len(rows),5)
        self.assertEqual([r['n'] for r in rows],[2]*5)
        self.assertEqual(sum(r['events'] for r in rows),5)
        for j,row in enumerate(rows):
            self.assertAlmostEqual(row['observed'],.5)
            self.assertAlmostEqual(row['predicted'],float(p[2*j:2*j+2].mean()))

    def test_tied_risks_are_never_split_to_force_five_bins(self):
        rows=calibration_bins([0,1,0,1],[.4]*4,[1,2,3,4])
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['n'],4)
        self.assertAlmostEqual(rows[0]['predicted'],.4)
        self.assertAlmostEqual(rows[0]['observed'],.6)

    def test_nonuniform_weights_determine_boundaries_and_zero_mass_explicit(self):
        rows=calibration_bins([0,1,1],[.1,.9,1.],[9.,1.,0.])
        self.assertEqual([r['n'] for r in rows],[1,2])
        self.assertEqual(sum(r['n_positive_weight'] for r in rows),2)
        self.assertAlmostEqual(rows[0]['weight_sum'],9.)
        self.assertAlmostEqual(rows[1]['predicted'],.9)
        self.assertTrue(all(not any('id'==k or 'predictions'==k for k in r) for r in rows))

if __name__=='__main__':unittest.main()
