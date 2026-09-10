"""Synthetic mathematical tests; actual cohort replay is an optional separate group."""
from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.integrate import quad
from research.mortality.competing_risks09 import (
    classify,exposure_design,preprocess_fit,transform,fit_one_cause,
    integrate_hazards,multiclass_brier,labels_at,predict,load_test,GROUPS,sha)
from research.mortality.replay_competing09 import (
    example,run,score_engineered,load_bundle,BUNDLE_PATH,BUNDLE_SHA256)


class CauseTests(unittest.TestCase):
    def test_partition_and_unknown_death(self):
        np.testing.assert_array_equal(classify([0,1,1,1,1],[np.nan,1,2,4,np.nan]),[-1,0,1,2,2])
    def test_all_ten_source_codes(self):
        np.testing.assert_array_equal(classify(np.ones(10),np.arange(1,11)),[0,1]+[2]*8)
    def test_infection_and_external_remain_events(self):
        np.testing.assert_array_equal(classify([1,1],[4,8]),[2,2])
    def test_invalid_cause_and_vital_status(self):
        for status,cause in [([0],[1]),([1],[0]),([1],[11]),([1],[np.inf]),([.5],[1]),([1,0],[1])]:
            with self.subTest(status=status,cause=cause),self.assertRaises(ValueError):classify(status,cause)
    def test_exposure_endpoints(self):
        exposure,events=exposure_design([1,5,8,10],[0,1,2,0])
        np.testing.assert_allclose(exposure.sum(axis=1),[1,5,8,8])
        self.assertEqual(events[0,0,0],1);self.assertEqual(events[1,1,1],1);self.assertEqual(events[2,2,2],1)
        self.assertEqual(events[3].sum(),0)
    def test_competing_deaths_contribute_person_time(self):
        exposure,events=exposure_design([.5,3,7],[0,1,2])
        np.testing.assert_allclose(exposure.sum(axis=1),[.5,3,7]);self.assertEqual(events.sum(),3)
    def test_survivor_has_no_event(self):
        exposure,events=exposure_design([9],[-1]);self.assertEqual(exposure.sum(),8);self.assertEqual(events.sum(),0)
    def test_invalid_exposure(self):
        for time,group in [([0],[0]),([-1],[0]),([1],[3]),([np.inf],[0])]:
            with self.subTest(time=time),self.assertRaises(ValueError):exposure_design(time,group)
    def test_competing_event_is_control_for_cause_cif(self):
        frame=pd.DataFrame({'dead':[1,1,0,1],'cause':[1,2,np.nan,8],'time_years':[2,3,9,7]})
        np.testing.assert_array_equal(labels_at(frame,5),[1,2,0,0])
    def test_early_censoring_not_ignored(self):
        frame=pd.DataFrame({'dead':[0],'cause':[np.nan],'time_years':[3]})
        with self.assertRaises(ValueError):labels_at(frame,5)


class IntegrationTests(unittest.TestCase):
    def test_constant_hazards_hand_calculation(self):
        rates=np.tile([.1,.2,.3],(1,3,1));actual=integrate_hazards(rates,5)[0]
        mass=1-math.exp(-3)
        np.testing.assert_allclose(actual,[math.exp(-3),mass/6,mass/3,mass/2],atol=1e-15)
    def test_zero_hazards(self):
        np.testing.assert_array_equal(integrate_hazards(np.zeros((2,3,3)),8),[[1,0,0,0]]*2)
    def test_zero_time(self):
        np.testing.assert_array_equal(integrate_hazards(np.ones((1,3,3)),0),[[1,0,0,0]])
    def test_single_cause_reduces_to_survival_formula(self):
        rates=np.array([[[.01],[.02],[.03]]]);p=integrate_hazards(rates,8)[0]
        self.assertAlmostEqual(p[1],-math.expm1(-(.01+4*.02+3*.03)))
    def test_piecewise_independent_quadrature(self):
        rates=np.array([[[.02,.03,.01],[.04,.01,.05],[.03,.06,.02]]]);p=integrate_hazards(rates,7)[0]
        def survival(t):
            spans=[min(t,1),max(0,min(t,5)-1),max(0,t-5)]
            return math.exp(-sum(a*sum(r) for a,r in zip(spans,rates[0])))
        for k in range(3):
            expected=sum(quad(lambda t:survival(t)*rates[0,j,k],a,b)[0] for j,(a,b) in enumerate([(0,1),(1,5),(5,7)]))
            self.assertAlmostEqual(p[k+1],expected,places=14)
    def test_random_mass_and_monotonicity(self):
        rng=np.random.default_rng(902);rates=rng.uniform(0,.7,(30,3,3));last=None
        for h in [0,.2,1,4,5,7,8]:
            p=integrate_hazards(rates,h);np.testing.assert_allclose(p.sum(axis=1),1,atol=1e-14)
            if last is not None:self.assertTrue(np.all(p[:,1:]>=last[:,1:]));self.assertTrue(np.all(p[:,0]<=last[:,0]))
            last=p
    def test_no_extrapolation(self):
        for h in [9,-1,True,np.inf,np.nan]:
            with self.subTest(h=h),self.assertRaises(ValueError):integrate_hazards(np.ones((1,3,3)),h)
    def test_no_negative_or_nonfinite_hazard(self):
        for value in [-.1,np.nan,np.inf]:
            with self.subTest(value=value),self.assertRaises(ValueError):integrate_hazards(np.full((1,3,3),value),5)
    def test_malformed_intervals(self):
        with self.assertRaises(ValueError):integrate_hazards(np.ones((1,3,3)),5,[1,1,8])
    def test_multiclass_brier_perfect_and_worst(self):
        self.assertEqual(multiclass_brier([0,1],[[1,0],[0,1]],[1,1]),0)
        self.assertEqual(multiclass_brier([0,1],[[0,1],[1,0]],[1,1]),2)
    def test_multiclass_brier_uniform(self):
        self.assertEqual(multiclass_brier([0],[[.25]*4],[1]),.75)
    def test_bad_probability_partition(self):
        for p in [[[.3,.3]],[[np.nan,.5]],[[-.1,1.1]]]:
            with self.subTest(p=p),self.assertRaises(ValueError):multiclass_brier([0],p,[1])


class FittingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(901);cls.z=rng.normal(size=(1800,2));cls.w=rng.uniform(.5,2,1800)
        t=rng.exponential(1/(.15*np.exp(.2*cls.z[:,0]-.25*cls.z[:,1])))
        group=np.where(t<=8,0,-1);cls.exposure,ev=exposure_design(np.minimum(t,10),group);cls.events=ev[:,:,0]
    def test_profiled_likelihood_against_full_poisson_GLM(self):
        fitted=fit_one_cause(self.z,self.exposure,self.events,self.w,0)
        people,interval=np.nonzero(self.exposure>0)
        x=np.column_stack([np.eye(3)[interval],self.z[people]])
        weights=self.w/self.w.mean()
        independent=sm.GLM(self.events[people,interval],x,family=sm.families.Poisson(),
            offset=np.log(self.exposure[people,interval]),freq_weights=weights[people]).fit()
        np.testing.assert_allclose(np.r_[fitted['baseline_log_hazards'],fitted['coefficients']],independent.params,atol=2e-6)
    def test_ridge_joint_score_equations(self):
        fitted=fit_one_cause(self.z,self.exposure,self.events,self.w,10)
        w=self.w/self.w.mean();beta=np.array(fitted['coefficients'])
        hazard=np.exp((self.z@beta)[:,None]+fitted['baseline_log_hazards'])
        residual=w[:,None]*(self.events-self.exposure*hazard)
        np.testing.assert_allclose(residual.sum(axis=0),0,atol=1e-9)
        np.testing.assert_allclose(self.z.T@residual.sum(axis=1)-10*beta,0,atol=1e-4)
    def test_weight_scaling_invariance(self):
        a=fit_one_cause(self.z,self.exposure,self.events,self.w,10)
        b=fit_one_cause(self.z,self.exposure,self.events,self.w*7,10)
        np.testing.assert_allclose(a['coefficients'],b['coefficients'],atol=1e-6)
    def test_zero_event_cell_blocks_no_pseudocount(self):
        events=self.events.copy();events[:,0]=0
        with self.assertRaises(ValueError):fit_one_cause(self.z,self.exposure,events,self.w,10)
    def test_bad_likelihood_input(self):
        for penalty in [-1,True,np.nan]:
            with self.subTest(penalty=penalty),self.assertRaises(ValueError):fit_one_cause(self.z,self.exposure,self.events,self.w,penalty)
    def test_preprocessing_is_train_only(self):
        frame=pd.DataFrame({'a':[0,1,2,3],'b':[7]*4});pp=preprocess_fit(frame,['a','b']);before=deepcopy(pp)
        transform(pd.DataFrame({'a':[1000],'b':[-99]}),pp)
        self.assertEqual(pp,before);self.assertEqual(pp['sd'][1],1.)
    def test_missing_never_imputed(self):
        frame=pd.DataFrame({'a':[1.,2.]});pp=preprocess_fit(frame,['a'])
        with self.assertRaises(ValueError):transform(pd.DataFrame({'a':[np.nan]}),pp)
        with self.assertRaises(ValueError):preprocess_fit(pd.DataFrame({'a':[np.nan]}),['a'])


class ReplayTests(unittest.TestCase):
    def test_fixed_bundle_and_source_endpoint(self):
        self.assertEqual(sha(BUNDLE_PATH),BUNDLE_SHA256)
        self.assertFalse(load_bundle()['noninfectious_risk_identified'])
    def test_both_panels(self):
        for panel in ['clinical','routine']:
            r=run(example(panel));self.assertEqual(r['status'],'calculated_research_only')
            self.assertFalse(r['clinical_use_ready']);self.assertIsNone(r['individual_uncertainty_interval'])
            self.assertIsNone(r['noninfectious_natural_probability']);self.assertEqual(len(r['probabilities']),3)
    def test_no_silent_routing(self):
        for key,value in [('country','RU'),('country','DE'),('age_years',80),('sex',None),
            ('panel','M1_routine'),('endpoint','noninfectious_natural_mortality'),('schema_version','0.6.0'),
            ('horizons_years',[10]),('horizons_years',[20]),('horizons_years',[5,5]),('horizons_years',[True]),
            ('model','M1_routine'),('patient_name','SYNTHETIC'),('research_only',False)]:
            with self.subTest(key=key,value=value):
                p=example();p[key]=value;r=run(p);self.assertEqual(r['status'],'blocked');self.assertIsNone(r['probabilities'])
    def test_complete_profile_required_even_comparator(self):
        for panel in ['clinical','routine']:
            p=example(panel);p['measurements'].pop('uacr');self.assertEqual(run(p)['status'],'blocked')
    def test_assay_and_units_preserved(self):
        p=example();p['measurements']['creatinine']['assay']='unknown';self.assertEqual(run(p)['status'],'blocked')
    def test_input_is_not_mutated(self):
        p=example();original=deepcopy(p);run(p);self.assertEqual(p,original)
    def test_invalid_structure(self):
        for payload in [[],None,{},'x']:
            with self.subTest(payload=payload):self.assertEqual(run(payload)['status'],'blocked')
    def test_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'model.json';path.write_text('{}')
            with self.assertRaises(ValueError):run(example(),bundle_path=path)
    def test_stdlib_replay_matches_numpy(self):
        bundle=load_bundle()
        from research.mortality.workbench.runtime import normalize,engineer
        for panel in ['clinical','routine']:
            p=example(panel);adapter=deepcopy(p);adapter.pop('panel');adapter.update(schema_version='0.6.0',endpoint='all_cause_death',model='M1_routine',horizons_years=[5])
            values,_,errors,_=normalize(adapter);self.assertEqual(errors,[])
            features=engineer(p['age_years'],p['sex'],p['clinical'],values)
            rows,_=score_engineered(features,bundle['models'][panel],[1,5,8])
            for h,row in zip([1,5,8],rows):
                vector=[row['survival']]+[row['cause_specific_cif'][g] for g in GROUPS]
                np.testing.assert_allclose(vector,predict(pd.DataFrame([features]),bundle['models'][panel],h)[0],atol=1e-14)
    def test_cli_stdlib_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'demo.json';cmd=[sys.executable,'-S','-m','research.mortality.replay_competing09','--demo','--out',str(path)]
            first=subprocess.run(cmd,capture_output=True,text=True);self.assertEqual(first.returncode,0,first.stderr)
            self.assertFalse(json.loads(path.read_text())['clinical_use_ready'])
            second=subprocess.run(cmd,capture_output=True,text=True);self.assertNotEqual(second.returncode,0)


class ActualDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=os.environ.get('NHANES_TEST09_SOURCE_DIR')
        if not path:raise unittest.SkipTest('Needs official 2009-2010 inputs, not synthetic substitutes')
        cls.bundle=load_bundle();cls.data,cls.flow=load_test(Path(path),cls.bundle['models']['routine']['preprocessing']['features'])
    def test_frozen_sample_and_followup(self):
        self.assertEqual(len(self.data),3182);self.assertEqual((labels_at(self.data,5)>0).sum(),150)
        self.assertEqual((labels_at(self.data,8)>0).sum(),302);self.assertGreaterEqual(self.flow['minimum_survivor_followup'],8)
    def test_reader_exact_match(self):
        self.assertTrue(self.flow['reader_audit']['exact_match']);self.assertEqual(self.flow['reader_audit']['cell_comparisons'],84296)
    def test_all_predictions_partition_and_monotonicity(self):
        for model in self.bundle['models'].values():
            a=predict(self.data,model,5);b=predict(self.data,model,8)
            np.testing.assert_allclose(a.sum(axis=1),1,atol=1e-14);self.assertTrue(np.all(b[:,1:]>=a[:,1:]))

if __name__=='__main__':unittest.main()
