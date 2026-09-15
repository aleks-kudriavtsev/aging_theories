"""Artificial fixtures test logic; source-dependent tests are explicitly separate."""
from copy import deepcopy
import json,os,unittest
from pathlib import Path
import numpy as np
from research.mortality.joint23.robustness import partition,crossfit,contrast
from research.mortality.joint23.joint_readiness import audit,REQUIRED,PROTEINS
from research.mortality.nfl22.analyze import (cohort,load_model,baseline_hazard,
    assign_folds,fit_update,predict_update)


def fixture():
    v={k:{'released':True,'variable':'test_'+k,'source_release':'test_release',
       'baseline_window':'test_same_visit','source_locator':'synthetic_fixture',
       'uniprot':PROTEINS.get(k),'matrix':'serum','assay':'synthetic_test_assay',
       'scale':'pg/mL','quantification_policy':'test_quantified','measured_form':'NT-proBNP'} for k in REQUIRED}
    return {'variables':v,'access_approved':True,'raw_records_received':True,'permission_locator':'synthetic_not_real_permission',
       'linkage':{'same_person_verified':True,'same_baseline_verified':True,'outcome_after_baseline_verified':True,
          'mortality_followup_verified':True,'weights_and_sampling_documented':True,'assay_transport_adjudicated':True,
          'joint_complete_participants':100,'deaths_after_baseline':10}}

class ContractTests(unittest.TestCase):
    def test_complete_only_design_review(self):
        r=audit(fixture());self.assertEqual(r['status'],'ready_for_design_review');self.assertFalse(r['clinical_use_ready']);self.assertFalse(r['coefficient_training_authorized_by_this_audit'])
    def test_TNFR1_never_TNF(self):
        c=fixture();c['variables']['TNF']['uniprot']='P19438';self.assertEqual(audit(c)['status'],'blocked_before_training')
    def test_NPPB_not_always_NTproBNP(self):
        c=fixture();c['variables']['NT_PROBNP']['measured_form']='NPPB';self.assertTrue(audit(c)['identity_errors'])
    def test_no_permission_no_ingestion(self):
        c=fixture();c['permission_locator']=None;self.assertFalse(audit(c)['gates']['access_approved'])
    def test_metadata_is_not_raw_records(self):
        c=fixture();c['raw_records_received']=False;self.assertEqual(audit(c)['status'],'blocked_before_training')
    def test_no_same_time_assertion_by_year(self):
        c=fixture();c['linkage']['same_baseline_verified']=False;self.assertEqual(audit(c)['status'],'blocked_before_training')
    def test_joint_count_not_public_sample_size(self):
        c=fixture();c['linkage']['joint_complete_participants']=None;self.assertIsNone(audit(c)['joint_participant_count'])
    def test_missing_UACR_not_imputed(self):
        c=fixture();c['variables'].pop('uacr');self.assertIn('uacr',audit(c)['missing_components'])
    def test_clinical_missing(self):
        c=fixture();c['variables'].pop('smoking');self.assertEqual(audit(c)['status'],'blocked_before_training')
    def test_relative_units_no_conversion(self):
        c=fixture();c['variables']['GDF15']['scale']='NPX';self.assertTrue(audit(c)['method_notes']);self.assertFalse(audit(c)['existing_NHANES_model_applicable'])
    def test_no_events_not_training(self):
        c=fixture();c['linkage']['deaths_after_baseline']=0;self.assertEqual(audit(c)['status'],'blocked_before_training')
    def test_released_not_assumed(self):
        c=fixture();c['variables']['GFAP']['released']=None;self.assertTrue(audit(c)['unresolved_metadata'])
    def test_bad_type(self):
        with self.assertRaises(ValueError):audit([])
        with self.assertRaises(ValueError):audit({'variables':[]})
    def test_booleans_not_counts(self):
        c=fixture();c['linkage']['joint_complete_participants']=True;self.assertIsNone(audit(c)['joint_participant_count'])
    def test_no_mutation(self):
        c=fixture();old=deepcopy(c);audit(c);self.assertEqual(c,old)

class RefitTests(unittest.TestCase):
    def setUp(self):
        n=100;self.f=np.arange(n)%5;self.y=(np.arange(n)%7==0).astype(int)
        self.h=np.linspace(.03,.3,n);self.x=np.linspace(1,5,n);self.w=np.ones(n)
    def test_partition_keeps_cluster(self):
        cl=np.repeat(np.arange(30),3);f=partition(cl,'test')
        self.assertTrue(all(len(set(f[cl==k]))==1 for k in set(cl)))
    def test_partition_deterministic(self):
        a=list(range(30));self.assertTrue(np.array_equal(partition(a,'x'),partition(a,'x')))
    def test_partition_reorder_invariance(self):
        a=list(range(30));self.assertTrue(np.array_equal(partition(a,'x'),partition(a[::-1],'x')[::-1]))
    def test_too_few_clusters(self):
        with self.assertRaises(ValueError):partition([1,2,3],'test')
    def test_original_update_exact(self):
        p,_=crossfit(self.y,self.h,self.h,self.w,self.x,self.f)
        tr=self.f!=0;te=~tr
        fit=fit_update(self.y[tr],self.h[tr],self.w[tr],self.x[tr],True)
        np.testing.assert_array_equal(p['nfl'][te],predict_update(self.h[te],self.x[te],fit))
    def test_test_outcomes_do_not_train_own_predictions(self):
        p,_=crossfit(self.y,self.h,self.h,self.w,self.x,self.f)
        y=self.y.copy();y[self.f==0]=1-y[self.f==0]
        q,_=crossfit(y,self.h,self.h,self.w,self.x,self.f)
        np.testing.assert_array_equal(p['nfl'][self.f==0],q['nfl'][self.f==0])
    def test_zero_weights_supported(self):
        w=self.w.copy();w[::3]=0;p,_=crossfit(self.y,self.h,self.h,w,self.x,self.f)
        self.assertTrue(np.isfinite(p['nfl']).all())
    def test_no_events_refused(self):
        with self.assertRaises(ValueError):crossfit(self.y*0,self.h,self.h,self.w,self.x,self.f)
    def test_constant_nfl_same_calibration(self):
        p,_=crossfit(self.y,self.h,self.h,self.w,self.x*0,self.f)
        np.testing.assert_allclose(p['nfl'],p['calibration'],atol=1e-12)
    def test_contrast_direction(self):
        p={'nfl':self.y*.8+.1,'calibration':np.full(100,.5)}
        self.assertLess(contrast(self.y,p,self.w)['delta_brier'],0)
    def test_fold_bounds(self):
        with self.assertRaises(ValueError):crossfit(self.y,self.h,self.h,self.w,self.x,self.f*2)
    def test_shape_check(self):
        with self.assertRaises(ValueError):crossfit(self.y,self.h[:-1],self.h,self.w,self.x,self.f)

class SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=os.getenv('JOINT23_CDC_SOURCE')
        if not p:raise unittest.SkipTest('Actual CDC inputs required')
        cls.d,cls.a=cohort(Path(p));cls.m=load_model()
    def test_domain_counts(self):
        mask=np.isfinite(self.d[self.m['preprocessing']['features']]).all(axis=1)
        self.assertEqual(int(mask.sum()),1270);self.assertEqual(int(((self.d.dead==1)&(self.d.time_years<=4)&mask).sum()),37)
    def test_source_integrity(self):self.assertEqual(self.a['files_verified'],29)
    def test_reproduce_actual_point(self):
        mask=np.isfinite(self.d[self.m['preprocessing']['features']]).all(axis=1).to_numpy();d=self.d.loc[mask]
        y=((d.dead==1)&(d.time_years<=4)).to_numpy(int);w=d.weight.to_numpy();x=np.log2(d.SSSNFL.to_numpy())
        p,_=crossfit(y,baseline_hazard(d,self.m,np.minimum(d.time_years,4)),baseline_hazard(d,self.m,4),w,x,assign_folds(self.d.cluster)[mask])
        r=contrast(y,p,w);self.assertAlmostEqual(r['delta_brier'],.000324702278134,places=12)
if __name__=='__main__':unittest.main()
