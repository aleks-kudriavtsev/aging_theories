"""Synthetic mathematical tests and explicit optional real-source audit."""
from copy import deepcopy
from pathlib import Path
import json,math,os,tempfile,unittest
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from research.mortality.nfl22 import analyze as a, joint_audit as j


class UpdateTests(unittest.TestCase):
    def setUp(self):
        self.y=np.array([1,0,0,1,0,0,1,0],float)
        self.h=np.array([.2,.3,.1,.4,.05,.2,.3,.25])
        self.w=np.array([1,2,3,1,4,2,1,2],float)
        self.x=np.array([2,3,4,4,2,3,6,5],float)
    def test_intercept_analytic(self):
        f=a.fit_update(self.y,self.h,self.w,self.x,False)
        self.assertAlmostEqual(math.exp(f['intercept']),float(self.w@self.y/(self.w@self.h)))
    def test_ridge_matches_independent_profile_likelihood(self):
        fit=a.fit_update(self.y,self.h,self.w,self.x,True)
        w=self.w/self.w.mean();x=self.x-np.average(self.x,weights=w);D=w@self.y
        def f(b):
            intercept=np.log(D/(w@(self.h*np.exp(b*x))))
            lp=intercept+b*x
            return float(w@(self.h*np.exp(lp)-self.y*lp)+.5*b*b)
        ref=minimize_scalar(f,bounds=(-2,2),method='bounded',options={'xatol':1e-12})
        self.assertAlmostEqual(fit['log2_nfl_coefficient'],ref.x,places=6)
    def test_global_weight_scale_invariant(self):
        p=a.fit_update(self.y,self.h,self.w,self.x,True);q=a.fit_update(self.y,self.h,self.w*100,self.x,True)
        self.assertAlmostEqual(p['log2_nfl_coefficient'],q['log2_nfl_coefficient'],places=8)
    def test_training_log_shift_changes_center_not_prediction(self):
        p=a.fit_update(self.y,self.h,self.w,self.x,True);q=a.fit_update(self.y,self.h,self.w,self.x+4,True)
        np.testing.assert_allclose(a.predict_update(self.h,self.x,p),a.predict_update(self.h,self.x+4,q),atol=1e-9)
    def test_no_events_rejected(self):
        with self.assertRaises(ValueError):a.fit_update(self.y*0,self.h,self.w,self.x,True)
    def test_bad_inputs(self):
        for key in ['h','w','x','y']:
            data=dict(y=self.y.copy(),h=self.h.copy(),w=self.w.copy(),x=self.x.copy());data[key][0]=np.nan
            with self.subTest(key=key),self.assertRaises(ValueError):a.fit_update(data['y'],data['h'],data['w'],data['x'],True)
    def test_negative_exposure_and_weight_rejected(self):
        with self.assertRaises(ValueError):a.fit_update(self.y,-self.h,self.w,self.x,True)
        with self.assertRaises(ValueError):a.fit_update(self.y,self.h,-self.w,self.x,True)
    def test_constant_nfl_returns_zero_slope(self):
        f=a.fit_update(self.y,self.h,self.w,np.ones(8)*4,True);self.assertAlmostEqual(f['log2_nfl_coefficient'],0)
    def test_hazard_probability_relation(self):
        f={'intercept':0,'log2_nfl_coefficient':0,'log2_nfl_center':0}
        np.testing.assert_allclose(a.predict_update([0,.2,1],[0,0,0],f),1-np.exp(-np.array([0,.2,1])))
    def test_fold_cluster_integrity(self):
        clusters=[f'{i}:1' for i in range(10)]*3;folds=a.assign_folds(clusters)
        self.assertEqual(len(set(folds)),5)
        for c in set(clusters):self.assertEqual(len(set(folds[np.array(clusters)==c])),1)
    def test_fold_row_order_invariant(self):
        clusters=[f'{i}:1' for i in range(10)]
        np.testing.assert_array_equal(a.assign_folds(clusters),a.assign_folds(clusters[::-1])[::-1])
    def test_fold_missing_clusters(self):
        with self.assertRaises(ValueError):a.assign_folds(['1:1','1:2'])
    def test_frozen_model_hash(self):self.assertEqual(a.sha(a.MODEL_PATH),a.MODEL_SHA)
    def test_tampered_model_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json';p.write_text('{}')
            with self.assertRaises(ValueError):a.load_model(p)
    def test_baseline_constant_hand_comparison(self):
        model={'preprocessing':{'features':['x'],'lower':[0.],'upper':[1.],'mean':[0.],'sd':[1.]},
               'causes':{'1':{'coefficients':[0.],'baseline_rates':[.1,.2]},'2':{'coefficients':[0.],'baseline_rates':[.2,.4]}}}
        got=a.baseline_hazard(pd.DataFrame({'x':[0.,0.,0.]}),model,[1,4,8])
        np.testing.assert_allclose(got,[.3,1.2,3.3])
    def test_baseline_time_restriction(self):
        m=a.load_model();d=pd.DataFrame({k:[m['preprocessing']['mean'][i]] for i,k in enumerate(m['preprocessing']['features'])})
        with self.assertRaises(ValueError):a.baseline_hazard(d,m,9)
    def test_no_result_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):a.run(Path(d)/'source',Path(d))


def fixture(oid='OID1',uniprot='P14136',sample='SSNA-001-PR1',npx='3',lod='1'):
    return {'SampleID':sample,'OlinkID':oid,'UniProt':uniprot,'Assay':'GFAP','Panel':'test',
      'Panel_Lot_Nr':'lot','PlateID':'plate','QC_Warning':'PASS','Assay_Warning':'PASS',
      'NPX':npx,'LOD':lod,'Normalization':'Plate control'}


class JointTests(unittest.TestCase):
    def test_subject_visit(self):self.assertEqual(j.parse_sample('SSNA-001B-PR1'),('001B','PR1'))
    def test_controls_not_people(self):self.assertIsNone(j.parse_sample('CONTROL_SAMPLE'))
    def test_four_visits_not_four_people(self):
        rows=[fixture(sample='SSNA-001-'+v) for v in ['PR1','PR2','PT1','PT2']]
        s,_,_=j.audit_rows(rows);self.assertEqual(s['biological_subjects'],1);self.assertEqual(s['biological_timepoints'],4)
    def test_negative_NPX_not_missing(self):
        self.assertTrue(j.quantitative_qc(fixture(npx='-.4',lod='-.8'))['usable_above_LOD'])
    def test_below_LOD_not_quantitative_zero(self):
        q=j.quantitative_qc(fixture(npx='-.9',lod='-.8'));self.assertTrue(q['below_LOD']);self.assertFalse(q['usable_above_LOD'])
    def test_equal_LOD_retained(self):self.assertTrue(j.quantitative_qc(fixture(npx='1',lod='1'))['usable_above_LOD'])
    def test_QC_and_assay_flags(self):
        for key in ['QC_Warning','Assay_Warning']:
            row=fixture();row[key]='WARN';self.assertFalse(j.quantitative_qc(row)['usable_above_LOD'])
    def test_missing_or_nonfinite_not_usable(self):
        for x in ['NaN','Inf','<3','']:
            self.assertFalse(j.quantitative_qc(fixture(npx=x))['usable_above_LOD'])
    def test_duplicate_not_averaged(self):
        with self.assertRaises(ValueError):j.audit_rows([fixture(),fixture()])
    def test_OlinkID_identity_conflict(self):
        row=fixture(sample='SSNA-002-PR1');row['Panel_Lot_Nr']='different'
        with self.assertRaises(ValueError):j.audit_rows([fixture(),row])
    def test_distinct_assays_preserved(self):
        s,t,_=j.audit_rows([fixture(),fixture(oid='OID2')]);self.assertEqual(len(t),2)
    def test_TNFR_not_TNF(self):
        s,t,_=j.audit_rows([fixture(uniprot='P19438')]);self.assertFalse(s['requested_four_measured']);self.assertEqual(t,[])
    def test_NPPB_not_NTproBNP(self):
        s,t,_=j.audit_rows([fixture(uniprot='P16860')]);self.assertFalse(s['NT_PROBNP_proteoform_confirmed'])
        self.assertEqual(t[0]['target'],'NPPB_unresolved_proteoform')
    def test_no_mortality_invented(self):
        s,_,_=j.audit_rows([fixture()]);self.assertFalse(s['mortality_training_ready']);self.assertFalse(s['concentration_mass_units'])
    def test_missing_header_refused(self):
        with self.assertRaises(ValueError):j.audit_rows([{'NPX':1}])
    def test_no_new_output_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):j.run('absent',d)


class RealSources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=os.environ.get('NFL22_SOURCE');joint=os.environ.get('JOINT22_ARCHIVE')
        if not raw or not joint:raise unittest.SkipTest('Requires separately retrieved public source data')
        cls.raw=Path(raw);cls.joint=Path(joint);cls.d,cls.audit=a.cohort(cls.raw)
    def test_real_files_and_counts(self):
        self.assertEqual(self.audit['files_verified'],29);self.assertEqual(len(self.d),1355)
    def test_real_complete_and_events(self):
        m=a.load_model();d=self.d[np.isfinite(self.d[m['preprocessing']['features']]).all(axis=1)]
        self.assertEqual(len(d),1270);self.assertEqual(int((d.dead.eq(1)&d.time_years.le(4)).sum()),37)
    def test_real_reader_and_followup(self):
        self.assertEqual(self.audit['reader_audit']['cell_comparisons'],81400)
        self.assertGreaterEqual(self.audit['minimum_survivor_followup'],4)
    def test_real_Gfap_QC_and_joint(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'qc';s=j.run(self.joint,out)
            self.assertEqual(s['biological_subjects'],10);self.assertEqual(s['max_joint_four_PASS_above_LOD_timepoints'],1)
            table=pd.read_csv(out/'target_assays.csv');r=table[table.target=='GFAP'].iloc[0]
            self.assertEqual(int(r.below_LOD),39);self.assertEqual(int(r.PASS_and_at_or_above_LOD),1)
    def test_baseline_vector_vs_individual_reference(self):
        m=a.load_model();d=self.d[np.isfinite(self.d[m['preprocessing']['features']]).all(axis=1)].head(64)
        pp=m['preprocessing'];vector=a.baseline_hazard(d,m,4)
        manual=[]
        for _,row in d.iterrows():
            z=[(min(max(float(row[name]),pp['lower'][i]),pp['upper'][i])-pp['mean'][i])/pp['sd'][i]
               for i,name in enumerate(pp['features'])]
            manual.append(math.fsum(math.exp(math.fsum(b*x for b,x in zip(c['coefficients'],z)))*4*c['baseline_rates'][0]
                                    for c in m['causes'].values()))
        np.testing.assert_allclose(vector,manual,atol=1e-13,rtol=1e-13)

if __name__=='__main__':unittest.main()
