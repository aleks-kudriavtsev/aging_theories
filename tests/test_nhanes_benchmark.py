import hashlib, json, math, os, tempfile, unittest
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from research.mortality.nhanes_benchmark import (
    INTERVAL_ENDS, PANELS, Preprocessor, read_mortality, weighted_quantile,
    interval_design, fit_piecewise, predict, evaluate, calibration,
    make_cohort, research_predict, verify_sources)

class ArithmeticTests(unittest.TestCase):
    def test_constant_hazard_closed_form(self):
        m={'baseline_log_hazards':[math.log(.02)]*3,'coefficients':[0.]}
        for h in [0,1,3,5,10]:
            np.testing.assert_allclose(predict([[4.]],m,h),[1-math.exp(-.02*h)],atol=1e-15)
    def test_piecewise_partial_interval(self):
        m={'baseline_log_hazards':np.log([.01,.02,.03]),'coefficients':[math.log(2)]}
        np.testing.assert_allclose(predict([[1.]],m,7),[1-math.exp(-2*(.01+4*.02+2*.03))])
    def test_horizon_boundaries(self):
        m={'baseline_log_hazards':[-4]*3,'coefficients':[0]}
        for h in [20,-1,True,np.nan,np.inf]:
            with self.subTest(h=h),self.assertRaises(ValueError):predict([[0]],m,h)
    def test_invalid_model_predictors(self):
        m={'baseline_log_hazards':[-4]*3,'coefficients':[0]}
        for z in [[np.nan],[[1,2]],[[np.inf]]]:
            with self.subTest(z=z),self.assertRaises(ValueError):predict(z,m,5)
    def test_interval_exposure_and_one_event(self):
        t=np.array([.5,1,5,7,10,15]); event=np.ones(6)
        x,exp,d,ids=interval_design(np.zeros((6,1)),t,event)
        np.testing.assert_allclose(np.bincount(ids,weights=exp),np.minimum(t,10))
        np.testing.assert_allclose(np.bincount(ids,weights=d),[1,1,1,1,1,0])
        self.assertTrue(np.all(exp>0))
    def test_censor_not_event(self):
        _,_,d,_=interval_design(np.zeros((2,1)),[2,9],[0,0])
        self.assertEqual(d.sum(),0)
    def test_invalid_event_and_dimensions(self):
        for t,event in [([1],[2]),([1],[.5]),([1,2],[1]),([0],[1]),([math.nan],[1])]:
            with self.subTest(t=t,event=event),self.assertRaises(ValueError):
                interval_design(np.zeros((len(t),1)),t,event)
    def test_brier_and_auc_manual(self):
        r=evaluate([0,1],[.25,.75],[1,1])
        self.assertEqual(r['auc'],1);self.assertEqual(r['brier'],.0625)
    def test_weights_change_prevalence(self):
        r=evaluate([0,1],[.1,.8],[3,1]);self.assertEqual(r['weighted_observed'],.25)
    def test_calibration_detects_overprediction(self):
        p=np.repeat([.1,.3,.6,.8],100)
        y=np.concatenate([np.r_[np.ones(n),np.zeros(100-n)]for n in [5,15,30,40]])
        r=calibration(y,p,np.ones(400));self.assertLess(r['calibration_in_the_large_logit'],0)

class PreprocessingTests(unittest.TestCase):
    def test_train_only_statistics(self):
        tr=pd.DataFrame({'age':[0.,1.,2.,np.nan]}); pp=Preprocessor().fit(tr,['age'],np.ones(4))
        before=json.dumps(pp.stats,sort_keys=True)
        pp.transform(pd.DataFrame({'age':[1000.,-999.,np.nan]}))
        self.assertEqual(before,json.dumps(pp.stats,sort_keys=True))
    def test_input_not_mutated(self):
        d=pd.DataFrame({'age':[0.,1.,np.inf,np.nan]});before=d.copy(deep=True)
        Preprocessor().fit(d,['age'],np.ones(4));pd.testing.assert_frame_equal(d,before)
    def test_missing_indicator(self):
        d=pd.DataFrame({'age':[0.,1.,np.nan]});pp=Preprocessor().fit(d,['age'],[1,1,1])
        self.assertEqual(pp.columns,['age','age__missing'])
        np.testing.assert_equal(pp.transform(d)[:,-1],[0,0,1])
    def test_all_missing_rejected(self):
        with self.assertRaises(ValueError):Preprocessor().fit(pd.DataFrame({'age':[np.nan]}),['age'],[1])
    def test_quantile_weight_scale_invariance(self):
        self.assertEqual(weighted_quantile([1,2,9],[1,2,3],.5),weighted_quantile([1,2,9],[10,20,30],.5))
    def test_outcome_fields_absent(self):
        forbidden={'dead','cause','t_exam','time_years','death_diabetes','death_hypertension','elig','SEQN'}
        for fields in PANELS.values():self.assertFalse(set(fields)&forbidden)

class FitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(99);cls.z=rng.normal(size=(650,2))
        true_t=rng.exponential(1/(.03*np.exp(.35*cls.z[:,0]-.25*cls.z[:,1])))
        cls.t=np.minimum(true_t,12);cls.e=(true_t<=12).astype(int)
        cls.w=rng.uniform(.5,2,size=650)
        cls.model=fit_piecewise(cls.z,cls.t,cls.e,cls.w,penalty=0)
    def test_independent_statsmodels_poisson(self):
        x,exp,d,ids=interval_design(self.z,self.t,self.e)
        w=self.w/self.w.mean()
        fit=sm.GLM(d,x,family=sm.families.Poisson(),offset=np.log(exp),freq_weights=w[ids]).fit()
        actual=np.r_[self.model['baseline_log_hazards'],self.model['coefficients']]
        np.testing.assert_allclose(actual,fit.params,rtol=1e-4,atol=2e-5)
    def test_weight_scale_invariance(self):
        m=fit_piecewise(self.z,self.t,self.e,self.w*10,penalty=0)
        np.testing.assert_allclose(predict(self.z,m,10),predict(self.z,self.model,10),atol=1e-6)
    def test_json_roundtrip(self):
        replay=json.loads(json.dumps(self.model))
        np.testing.assert_array_equal(predict(self.z,self.model,5),predict(self.z,replay,5))
    def test_invalid_weights(self):
        for w in [np.zeros(650),np.ones(3),np.full(650,np.nan)]:
            with self.subTest(w=w.shape),self.assertRaises(ValueError):fit_piecewise(self.z,self.t,self.e,w)
    def test_invalid_penalty(self):
        for p in [-1,True,np.nan]:
            with self.subTest(p=p),self.assertRaises(ValueError):fit_piecewise(self.z,self.t,self.e,self.w,penalty=p)
    def test_risk_monotonic(self):
        p=np.array([predict(self.z,self.model,t)for t in [0,1,5,10]])
        self.assertTrue(np.all(np.diff(p,axis=0)>=0));self.assertTrue(np.all((p>=0)&(p<=1)))

class ReaderTests(unittest.TestCase):
    def test_fixed_width_and_missingness(self):
        line=list(' '*48)
        for lo,hi,text in [(0,6,'123'),(14,15,'1'),(15,16,'1'),(16,19,'001'),(19,20,'0'),(20,21,'0'),(42,45,'101'),(45,48,'99')]:
            line[lo:hi]=text.rjust(hi-lo)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'lmf.dat';p.write_text(''.join(line)+'\n')
            r=read_mortality(p).iloc[0]
            self.assertEqual(r.SEQN,123);self.assertEqual(r.dead,1);self.assertEqual(r.t_exam,99);self.assertEqual(r.t_int,101)
    def test_bad_record_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'lmf.dat';p.write_text('1234\n')
            with self.assertRaises(ValueError):read_mortality(p)
    def test_integrity_failure_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'x').write_bytes(b'bad')
            (p/'manifest.json').write_text(json.dumps({'files':[{'file':'x','bytes':3,'sha256':'0'*64}]}))
            with self.assertRaises(ValueError):verify_sources(p)

class SourceIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=os.environ.get('NHANES_SOURCE_DIR');results=os.environ.get('NHANES_RESULT_DIR')
        if not path or not results:raise unittest.SkipTest('Set NHANES_SOURCE_DIR and NHANES_RESULT_DIR')
        cls.source=Path(path);cls.results=Path(results)
        cls.data,cls.flow=make_cohort(cls.source)
        cls.models=json.loads((cls.results/'fitted_models.json').read_text())
        cls.test=cls.data[cls.data.cycle.eq(2003)].copy()
    def test_all_source_hashes(self):self.assertEqual(verify_sources(self.source)['checked_files'],89)
    def test_locked_split_and_events(self):
        tr=self.data[self.data.cycle<2003];te=self.test
        self.assertEqual(len(tr),4552);self.assertEqual(len(te),2384);self.assertFalse(set(tr.SEQN)&set(te.SEQN))
        self.assertEqual(((tr.dead==1)&(tr.time_years<=10)).sum(),756)
        self.assertEqual(((te.dead==1)&(te.time_years<=10)).sum(),395)
    def test_weighted_prevalence(self):
        te=self.test;y=(te.dead==1)&(te.time_years<=10)
        self.assertAlmostEqual(np.average(y,weights=te.weight),.11506015719576836)
    def test_replayed_metrics(self):
        m=self.models['M4_cystatinC'];te=self.test
        pred=research_predict(te,m,10,country='US',acknowledge_research_only=True)
        y=(te.dead==1)&(te.time_years<=10);actual=evaluate(y,pred,te.weight)
        self.assertAlmostEqual(actual['auc'],.8715652965094911)
        self.assertAlmostEqual(actual['weighted_predicted'],.132666645277662)
    def test_country_and_research_gate(self):
        m=self.models['M4_cystatinC']
        for kwargs in [dict(country='RU',acknowledge_research_only=True),dict(country='US')]:
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):research_predict(self.test,m,10,**kwargs)
    def test_age_and_sex_gate(self):
        m=self.models['M4_cystatinC']
        for col,value in [('RIDAGEYR',80),('male',np.nan),('age',999)]:
            d=self.test.head(1).copy();d[col]=value
            with self.subTest(col=col),self.assertRaises(ValueError):research_predict(d,m,10,country='US',acknowledge_research_only=True)
    def test_units_and_creatinine_calibration(self):
        d=self.data;old=d[d.cycle.eq(1999)]
        np.testing.assert_allclose(old.creatinine_standardized,1.013*old.LBXSCR+.147,equal_nan=True)
        np.testing.assert_allclose(d.uacr,100*d.URXUMA/d.URXUCR,equal_nan=True)
        np.testing.assert_allclose(d.albumin,d.LBXSAL*10,equal_nan=True)
    def test_no_twenty_year_prediction(self):
        with self.assertRaises(ValueError):research_predict(self.test,self.models['M1_routine'],20,country='US',acknowledge_research_only=True)
    def test_same_assay_sample_and_followup(self):
        d=self.data;self.assertTrue(d[['SSBNP','SSTNT','SSCYST']].gt(0).all().all())
        self.assertGreater(d.loc[d.dead.eq(0),'time_years'].min(),10)
    def test_no_individual_output(self):
        for path in self.results.glob('*.csv'):
            self.assertNotIn('SEQN',pd.read_csv(path,nrows=1).columns)
        for m in self.models.values():self.assertIs(m['clinical_use'],False)

if __name__=='__main__':unittest.main()
