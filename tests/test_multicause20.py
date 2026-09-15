"""Synthetic invariants and explicit actual-source integration for multicause20."""
import json,math,os,tempfile,unittest
from pathlib import Path
from copy import deepcopy
import numpy as np
import pandas as pd
from research.mortality.multicause20.cohort import classify,load,TARGETS,IDS,CORE,HEPATIC,CLINICAL,read_fwf
from research.mortality.multicause20.model import fit,predict
from research.mortality.multicause20.scoring import outcomes,metrics,aalen_johansen,auc,resample
from research.mortality.multicause20 import replay

class Taxonomy(unittest.TestCase):
 def test_twenty_targets(self):self.assertEqual(len(TARGETS),20);self.assertEqual(len(set(IDS)),20)
 def test_ihd_not_all_heart(self):self.assertEqual(classify('4109'),'ihd');self.assertEqual(classify('4280'),'heart_failure')
 def test_individual_cancers(self):
  for code,target in [('1629','lung_cancer'),('1539','colorectal_cancer'),('1579','pancreatic_cancer'),('185','prostate_cancer'),('1519','gastric_cancer')]:self.assertEqual(classify(code),target)
 def test_infections_not_counted_as_targets(self):self.assertEqual(classify('486'),'infection');self.assertEqual(classify('0389'),'infection');self.assertNotIn('infection',IDS)
 def test_external_implicit_E(self):self.assertEqual(classify('8889'),'external');self.assertNotIn('external',IDS)
 def test_unknown_not_zero(self):self.assertEqual(classify('0000'),'unknown_cause')
 def test_no_ICD10_substitution(self):
  with self.assertRaises(ValueError):classify('I219')
 def test_no_reinterpretation_broad_renal(self):self.assertEqual(classify('5849'),'renal');self.assertEqual(classify('585'),'renal')
 def test_decimal_forms(self):self.assertEqual(classify('331.0'),'dementia')
 def test_all_possible_prefixes_unambiguous(self):
  for c in range(1,1000):classify(str(c).zfill(3))

class Scoring(unittest.TestCase):
 def test_death_other_is_not_censor(self):
  y,w,_=outcomes([1,3,4],[1,1,0],[2,1,0],[1,1,1],3,2);np.testing.assert_equal(y,[2,1,0]);np.testing.assert_equal(w,[1,1,1])
 def test_missing_death_cause_refused(self):
  with self.assertRaises(ValueError):outcomes([1,4],[1,0],[0,0],[1,1],3,2)
 def test_early_censor_zero_weight(self):
  y,w,i=outcomes([1,2,3],[1,0,1],[1,0,2],[1,1,1],3,2);np.testing.assert_equal(w,[1,0,2]);self.assertEqual(i['early_censored_n'],1)
 def test_ties_AJ_reference(self):
  t=[1,1,2,4];d=[1,0,1,0];c=[1,0,2,0];w=np.array([2.,3.,4.,5.])
  y,ew,_=outcomes(t,d,c,w,3,2);np.testing.assert_allclose(np.bincount(y,weights=ew,minlength=3)/w.sum(),aalen_johansen(t,d,c,w,3,2),atol=1e-14)
 def test_random_multicause_AJ(self):
  rng=np.random.default_rng(20)
  for _ in range(50):
   t=rng.integers(1,10,80);d=rng.integers(0,2,80);c=np.where(d,rng.integers(1,8,80),0);w=rng.uniform(.1,5,80)
   y,ew,_=outcomes(t,d,c,w,4,7);np.testing.assert_allclose(np.bincount(y,weights=ew,minlength=8)/w.sum(),aalen_johansen(t,d,c,w,4,7),atol=1e-13)
 def test_no_support(self):
  with self.assertRaises(ValueError):outcomes([1],[0],[0],[1],3,1)
 def test_auc_ties(self):self.assertEqual(auc([0,1],[.5,.5],[1,1]),.5)
 def test_auc_undefined(self):self.assertIsNone(auc([0,0],[.3,.4],[1,1]))
 def test_brier_hand(self):
  p=np.array([[.8,.2],[.4,.6]]);r=metrics(np.array([0,1]),p,np.ones(2),np.ones(2));self.assertAlmostEqual(r[0]['brier'],.20)
 def test_singleton_resampling_refused(self):
  with self.assertRaises(ValueError):resample(pd.DataFrame({'stratum':[1],'psu':[1],'weight':[1]}),np.random.default_rng(1))

class Runtime(unittest.TestCase):
 def test_primary_demo(self):self.assertEqual(replay.run(replay.example())['status'],'calculated_historical_research_only')
 def test_hepatic_demo(self):self.assertEqual(len(replay.run(replay.example('hepatic'))['selected_noninfectious_families']),17)
 def test_sums(self):
  for r in replay.run(replay.example())['probabilities']:self.assertAlmostEqual(r['survival']+sum(r['causes'].values()),1,places=14)
 def test_cause_count_not_inflated(self):self.assertEqual(len(replay.run(replay.example())['selected_noninfectious_families']),17)
 def test_noninfectious_bounds_not_CI(self):
  for r in replay.run(replay.example())['probabilities']:
   b=r['unresolved_noninfectious_probability_bounds'];self.assertLessEqual(b['lower'],b['upper']);self.assertLessEqual(b['upper'],r['all_cause']);self.assertIn('not_confidence',b['kind'])
 def test_no_imputation(self):
  p=replay.example();p['measurements'].pop('albumin')
  with self.assertRaises(ValueError):replay.run(p)
 def test_no_unmeasured_TNF(self):
  p=replay.example();p['measurements']['TNF']={'value':2}
  with self.assertRaises(ValueError):replay.run(p)
 def test_country_block(self):
  p=replay.example();p['country']='RU'
  with self.assertRaises(ValueError):replay.run(p)
 def test_method_block(self):
  p=replay.example();p['measurements']['albumin']['assay_context']='new_aptamer_assay'
  with self.assertRaises(ValueError):replay.run(p)
 def test_unsupported_horizon(self):
  for h in [[30],[True],[15,15]]:
   p=replay.example();p['horizons_years']=h
   with self.assertRaises(ValueError):replay.run(p)
 def test_input_nonfinite(self):
  for v in [None,'<3',True,float('inf'),10**400]:
   p=replay.example();p['measurements']['albumin']['value']=v
   with self.assertRaises(ValueError):replay.run(p)
 def test_female_prostate_structural_zero(self):
  p=replay.example();p['clinical']['sex']='female'
  for r in replay.run(p)['probabilities']:self.assertEqual(r['causes']['prostate_cancer'],0)
 def test_mutation(self):
  p=replay.example();old=deepcopy(p);replay.run(p);self.assertEqual(p,old)
 def test_duplicate_JSON(self):
  with self.assertRaises(ValueError):replay.strict_json('{"x":1,"x":2}')
 def test_no_clinical_release(self):self.assertFalse(replay.run(replay.example())['clinical_use_ready'])

class ActualSources(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  source=os.environ.get('NHEFS20_SOURCE')
  if not source:raise unittest.SkipTest('Requires official source files')
  cls.d,cls.audit=load(source)
 def test_actual_size(self):self.assertEqual(self.audit['vital_n'],14407);self.assertEqual(len(self.d),10900)
 def test_primary_complete(self):self.assertEqual(int(np.isfinite(self.d[CORE]).all(axis=1).sum()),9395)
 def test_cause_crosscheck(self):self.assertEqual(self.audit['certificate_mismatches'],0);self.assertEqual(self.audit['certificate_causes_compared'],4497)
 def test_nonoverlap_renal_albumin(self):self.assertEqual(int(self.d[['albumin','creatinine']].notna().all(axis=1).sum()),0)
 def test_full_replay_agreement(self):
  from research.mortality.multicause20.run import arr
  compared=0;maximum=0
  for study in ['primary','hepatic']:
   doc,_=replay.load(study);d=self.d.copy()
   if study=='hepatic':d=d[d.subsample.eq(2)]
   d=d[np.isfinite(d[doc['models']['laboratory']['preprocessing']['features']]).all(axis=1)&d.stratum.isin(doc['heldout_strata'])]
   for panel,model in doc['models'].items():
    native={h:predict(d,model,h) for h in [10,15]}
    for idx,(_,row) in enumerate(d.iterrows()):
     rows,_=replay.score(row.to_dict(),'male' if row.sex==1 else 'female',model,[10,15])
     for r in rows:
      ref=native[r['horizon_years']][idx];val=np.array([r['survival'],*r['causes'].values()]);maximum=max(maximum,float(np.max(abs(val-ref))));compared+=len(val)
  self.assertEqual(compared,375648);self.assertLess(maximum,1e-12)
  print('NUMERICAL_REPLAY_AUDIT',compared,maximum)

if __name__=='__main__':unittest.main()
