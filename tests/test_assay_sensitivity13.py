"""Perturbation contracts and optional repeated-cohort checks, not new validation."""
import json,os,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from research.mortality.assay_sensitivity13 import (scenarios,perturb,weighted_quantile,metric_row,
    SENTINELS,MODEL_SHA,SECONDARY_SHA,MODEL_PATH,sha)
from research.mortality.method_bridge13 import deming
from research.mortality.transport_no_crp12 import load_evaluation,predict,labels

class ContractTests(unittest.TestCase):
    def test_grid(self):
        s=scenarios();self.assertEqual(len(s),45);self.assertEqual(len({r['id'] for r in s}),45)
        self.assertTrue(set(SENTINELS)<={r['id'] for r in s})
    def test_invalid_changes(self):
        with self.assertRaises(ValueError):perturb(pd.DataFrame(),{'unknown':{'factor':1,'offset':0}})
    def test_invalid_magnitude(self):
        for v in [True,0,-1,float('inf')]:
            with self.subTest(v=v),self.assertRaises(ValueError):perturb(pd.DataFrame(),{'rdw':{'factor':v,'offset':0}})
    def test_missing_transform_part(self):
        with self.assertRaises(ValueError):perturb(pd.DataFrame(),{'rdw':{'factor':1}})
    def test_uacr_offset_forbidden(self):
        with self.assertRaises(ValueError):perturb(pd.DataFrame(),{'uacr':{'factor':1,'offset':2}})
    def test_weighted_quantile(self):self.assertEqual(weighted_quantile([1,10],[99,1],.95),1)
    def test_zero_weight_quantile(self):self.assertEqual(weighted_quantile([0,1,2],[0,1,1],0),1)
    def test_negative_weight(self):
        with self.assertRaises(ValueError):weighted_quantile([1,2],[-1,2],.5)
    def test_frozen_models(self):
        self.assertEqual(sha(MODEL_PATH),MODEL_SHA)
        self.assertEqual(sha(MODEL_PATH.parent/'model_sensitivity12.json'),SECONDARY_SHA)
    def test_deming_independent_profile_loss(self):
        x=np.array([1.,2.,4.,5.,8.,11.,13.,17.]);y=np.array([3.,4.,8.,9.,15.,18.,23.,31.]);lam=4.
        r=deming(x.tolist(),y.tolist(),variance_ratio=lam)
        # Profile out latent true reference values; independent numerical minimization.
        def loss(b):return np.sum((y-b[0]-b[1]*x)**2/(lam+b[1]**2))
        fit=minimize(loss,[0.,1.],method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-12})
        self.assertTrue(fit.success);self.assertAlmostEqual(r['slope'],fit.x[1],places=6)
        self.assertAlmostEqual(r['intercept'],fit.x[0],places=6)

class ActualProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        src=os.environ.get('ASSAY13_SOURCE_DIR')
        if not src:raise unittest.SkipTest('Needs original, previously examined NHANES input files')
        cls.d,cls.audit=load_evaluation(Path(src),2013,MODEL_SHA)
        cls.primary=json.loads(MODEL_PATH.read_text())['models']['routine_no_crp']
        cls.secondary=json.loads((MODEL_PATH.parent/'model_sensitivity12.json').read_text())['model']
    def test_reused_cohort_size(self):self.assertEqual(len(self.d),3013)
    def test_baseline_identity(self):
        p=predict(self.d,self.primary,4);q=predict(perturb(self.d,{}),self.primary,4)
        np.testing.assert_array_equal(p,q)
    def test_uacr_recomputed(self):
        q=perturb(self.d,{'uacr':{'factor':1.1,'offset':0}})
        np.testing.assert_allclose(q.log_uacr-self.d.log_uacr,np.log2(1.1),rtol=1e-12)
    def test_creatinine_recomputes_egfr(self):
        q=perturb(self.d,{'creatinine':{'factor':1.1,'offset':0}})
        self.assertTrue((q.egfr_high<=self.d.egfr_high).all());self.assertTrue((q.egfr_low<=self.d.egfr_low).all())
    def test_albumin_offset_units(self):
        q=perturb(self.d,{'albumin':{'factor':1,'offset':-2}})
        np.testing.assert_allclose(q.albumin-self.d.albumin,-2)
    def test_do_not_change_covariates_or_events(self):
        q=perturb(self.d,{'rdw':{'factor':1,'offset':.5}})
        for c in ('SEQN','RIDAGEYR','RIAGENDR','dead','cause','time_years','weight','age','sbp','log_uacr'):
            pd.testing.assert_series_equal(q[c],self.d[c])
    def test_no_cbc_negative_control(self):
        q=perturb(self.d,{'rdw':{'factor':1.1,'offset':0},'wbc':{'factor':.9,'offset':0}})
        np.testing.assert_array_equal(predict(q,self.secondary,4),predict(self.d,self.secondary,4))
    def test_probability_mass_all_scenarios(self):
        d=self.d.iloc[:20]
        for s in scenarios():
            p=predict(perturb(d,s['changes']),self.primary,4)
            self.assertLess(np.max(np.abs(p.sum(axis=1)-1)),1e-12)
    def test_no_silent_row_dropping(self):
        for s in scenarios():self.assertEqual(len(perturb(self.d,s['changes'])),len(self.d))
    def test_baseline_metrics_reproduce(self):
        p=predict(self.d,self.primary,4);r=metric_row(labels(self.d,4),p,p,self.d.weight)
        self.assertAlmostEqual(r['six_state_brier'],.06377063237897933,places=12)
        self.assertEqual(r['mean_shift_pp'],0)
    def test_rdw_fixed_scenario_regression(self):
        p=predict(self.d,self.primary,4);q=predict(perturb(self.d,{'rdw':{'factor':1,'offset':.5}}),self.primary,4)
        r=metric_row(labels(self.d,4),q,p,self.d.weight)
        self.assertAlmostEqual(r['mean_shift_pp'],.484579842,places=6)
    def test_negative_result_is_not_dropped(self):
        q=perturb(self.d,{'rdw':{'factor':1,'offset':-.5}})
        self.assertLess(np.average(1-predict(q,self.primary,4)[:,0],weights=self.d.weight),
                        np.average(1-predict(self.d,self.primary,4)[:,0],weights=self.d.weight))

if __name__=='__main__':unittest.main()
