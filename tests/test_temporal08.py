"""Synthetic arithmetic tests plus optional actual-data integrity checks."""
import json,math,os,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from research.mortality.validate_temporal08 import (fit_hazard_factor,cumulative_hazard,weighted_auc,
 audit_mortality_reader,verify_new_sources,design_bootstrap_weights,engineer_new,MODEL_SHA,sha)
from research.mortality.replay_temporal08 import replay,recalibrate_probability,SPEC_PATH,SPEC_SHA256
from research.mortality.workbench.runtime import example,BUNDLE_PATH


class Mathematics(unittest.TestCase):
    def test_factor_hand_calculation(self):
        self.assertEqual(fit_hazard_factor([1,0],[.2,.3],[1,1]),2.)
    def test_factor_score_equation(self):
        y=np.array([1,0,1]);h=np.array([.2,.6,.3]);w=np.array([2,1,3]);r=fit_hazard_factor(y,h,w)
        self.assertAlmostEqual(float(w@(y-r*h)),0.)
    def test_factor_weight_scaling(self):
        self.assertEqual(fit_hazard_factor([1,0],[.2,.3],[2,2]),fit_hazard_factor([1,0],[.2,.3],[1,1]))
    def test_factor_invalid_inputs(self):
        for y,h,w in [([],[],[]),([0,0],[.1,.2],[1,1]),([1],[0],[1]),([1],[-1],[1]),([1],[.1],[-1]),([2],[.2],[1])]:
            with self.subTest(y=y),self.assertRaises(ValueError):fit_hazard_factor(y,h,w)
    def test_piecewise_hazard(self):
        model={'interval_ends_years':[1,5,10],'baseline_log_hazards':np.log([.1,.2,.3]).tolist(),'coefficients':[0.]}
        self.assertTrue(np.allclose(cumulative_hazard(np.zeros((4,1)),model,[0,1,5,10]),[0,.1,.9,2.4]))
    def test_hazard_no_extrapolation(self):
        with self.assertRaises(ValueError):cumulative_hazard(np.zeros((1,1)),{},[11])
    def test_auc_hand_tie(self):
        self.assertEqual(weighted_auc([0,1],[.5,.5],[1,1]),.5)
    def test_auc_independent_implementation(self):
        rng=np.random.default_rng(6)
        for _ in range(25):
            y=rng.integers(0,2,50);p=rng.integers(0,10,50)/10;w=rng.random(50)
            self.assertAlmostEqual(weighted_auc(y,p,w),roc_auc_score(y,p,sample_weight=w),places=13)
    def test_resample_design(self):
        d=pd.DataFrame({'weight':[1.,1.,1.,1.],'stratum':['a','a','b','b'],'SDMVPSU':[1,2,1,2]})
        w=design_bootstrap_weights(d,np.random.default_rng(5));self.assertEqual(w.sum(),4);self.assertEqual(set(w),{0,2})
    def test_singleton_refused(self):
        d=pd.DataFrame({'weight':[1.],'stratum':['a'],'SDMVPSU':[1]})
        with self.assertRaises(ValueError):design_bootstrap_weights(d,np.random.default_rng(5))
    def test_recalibration_identity(self):
        self.assertAlmostEqual(recalibrate_probability(.2,1),.2)
    def test_recalibration_squared_survival(self):
        self.assertAlmostEqual(recalibrate_probability(.2,2),.36)
    def test_recalibration_endpoints_and_monotonicity(self):
        self.assertEqual(recalibrate_probability(0,.7),0);self.assertEqual(recalibrate_probability(1,.7),1)
        p=[recalibrate_probability(v,.7) for v in [0,.1,.5,.9,1]];self.assertEqual(p,sorted(p))
    def test_invalid_probability(self):
        for p,f in [(-.1,1),(1.1,1),(.1,0),(.1,-1),(.1,math.nan),(True,1)]:
            with self.subTest(p=p,f=f),self.assertRaises(ValueError):recalibrate_probability(p,f)


class ReplayTests(unittest.TestCase):
    def test_frozen_hash_unchanged(self):self.assertEqual(sha(BUNDLE_PATH),MODEL_SHA)
    def test_spec_hash(self):self.assertEqual(sha(SPEC_PATH),SPEC_SHA256)
    def test_acknowledgement(self):
        with self.assertRaises(ValueError):replay(example())
    def test_model_restriction(self):
        with self.assertRaises(ValueError):replay(example('M4_cystatinC'),acknowledge_temporal_research=True)
    def test_full_profile_replay(self):
        r=replay(example(),acknowledge_temporal_research=True)
        self.assertEqual(r['status'],'calculated_research_only');self.assertFalse(r['clinical_use_ready'])
        self.assertEqual(r['recalibration']['evaluation_n'],5287)
        self.assertTrue(all(a['probability']<b['probability'] for a,b in zip(r['probabilities'],r['frozen_probabilities'])))
    def test_missing_value_blocked(self):
        p=example();p['measurements'].pop('crp');r=replay(p,acknowledge_temporal_research=True)
        self.assertEqual(r['status'],'blocked');self.assertIsNone(r['probabilities'])
    def test_country_and_horizon_remain_blocked(self):
        for k,v in [('country','RU'),('horizons_years',[20])]:
            p=example();p[k]=v;self.assertEqual(replay(p,acknowledge_temporal_research=True)['status'],'blocked')
    def test_spec_tamper_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json';p.write_text('{}')
            with self.assertRaises(ValueError):replay(example(),acknowledge_temporal_research=True,spec_path=p)


class ActualDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=os.environ.get('NHANES_NEW_SOURCE_DIR')
        if not raw:raise unittest.SkipTest('Requires separately retrieved official NHANES2005-2008 data')
        cls.source=Path(raw)
    def test_new_sources(self):self.assertEqual(verify_new_sources(self.source)['required_data_files'],30)
    def test_independent_readers(self):
        for year in [2005,2007]:self.assertTrue(audit_mortality_reader(self.source/f'{year}_mortality.dat')['exact_match'])

if __name__=='__main__':unittest.main()
