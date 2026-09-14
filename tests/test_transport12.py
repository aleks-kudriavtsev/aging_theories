"""Mathematical and input-control tests; synthetic unless labelled frozen artifact."""
import json,math,tempfile,unittest
from copy import deepcopy
from pathlib import Path
import numpy as np
import pandas as pd
from research.mortality.transport_no_crp12 import (classify,exposure,fit_cause,labels,ROUTINE,GROUPS,
    predict,engineer_no_crp,sha)
from research.mortality.competing_risks09 import integrate_hazards,multiclass_brier,transform
from research.mortality.validate_temporal08 import engineer_new
from research.mortality.replay_transport12 import (load_bundle,MODEL_PATH,MODEL_SHA,example,run,score_engineered)

class EndpointTests(unittest.TestCase):
    def test_partition(self):
        self.assertEqual(classify([0,1,1,1,1,1,1],[np.nan,1,2,3,5,4,np.nan]).tolist(),[-1,0,1,2,3,4,4])
    def test_unknown_is_death(self):self.assertEqual(classify([1],[np.nan])[0],4)
    def test_unsupported_code(self):
        for code in (0,11,1.5,np.inf):
            with self.subTest(code=code),self.assertRaises(ValueError):classify([1],[code])
    def test_survivor_cause(self):
        with self.assertRaises(ValueError):classify([0],[1])
    def test_bad_status_shape(self):
        for dead,cause in [([2],[1]),([1,0],[1]),([[1]],[[1]])]:
            with self.assertRaises(ValueError):classify(dead,cause)
    def test_event_boundary(self):
        dt,e=exposure([5.,8.,9.],[2,3,4]);self.assertEqual(e.sum(),2)
        self.assertEqual(e[0,0,2],1);self.assertEqual(e[1,1,3],1);self.assertTrue(np.allclose(dt,[[5,0],[5,3],[5,3]]))
    def test_nonpositive_time(self):
        for t in (0,-1,np.nan):
            with self.assertRaises(ValueError):exposure([t],[0])
    def test_early_censoring_not_alive(self):
        d=pd.DataFrame({'dead':[0],'cause':[np.nan],'time_years':[4.]})
        with self.assertRaises(ValueError):labels(d,5)
    def test_later_death_alive_at_horizon(self):
        d=pd.DataFrame({'dead':[1],'cause':[3],'time_years':[7.]});self.assertEqual(labels(d,5)[0],0)
    def test_five_cause_mass(self):
        rates=np.ones((3,2,5))*.01;p=integrate_hazards(rates,4,[5,8])
        self.assertTrue(np.allclose(p.sum(axis=1),1));self.assertAlmostEqual(p[0,0],math.exp(-.2))
    def test_zero_rates(self):self.assertEqual(integrate_hazards(np.zeros((1,2,5)),5,[5,8])[0,0],1)
    def test_monotone_cif(self):
        r=np.random.default_rng(5).uniform(0,.05,(10,2,5))
        self.assertTrue((integrate_hazards(r,5,[5,8])[:,1:]>=integrate_hazards(r,4,[5,8])[:,1:]).all())
    def test_brier_hand(self):self.assertAlmostEqual(multiclass_brier(np.array([0]),[[.5,.1,.1,.1,.1,.1]],[1]),.3)
    def test_missing_CRp_feature(self):self.assertNotIn('log_crp',ROUTINE);self.assertEqual(len(ROUTINE),22)

class LikelihoodTests(unittest.TestCase):
    def test_unpenalized_intercepts(self):
        z=np.zeros((4,1));dt=np.array([[1,0],[5,1],[5,3],[5,3]],float);ev=np.array([[1,0],[0,1],[0,0],[0,0]])
        m=fit_cause(z,dt,ev,[1,1,1,1],10.)
        self.assertTrue(np.allclose(np.exp(m['baseline_log_hazards']),[1/16,1/7]))
    def test_no_pseudoevents(self):
        with self.assertRaises(ValueError):fit_cause(np.zeros((2,1)),np.ones((2,2)),np.zeros((2,2)),np.ones(2),10.)
    def test_invalid_weights(self):
        with self.assertRaises(ValueError):fit_cause(np.zeros((2,1)),np.ones((2,2)),np.eye(2),[-1,1],10.)
    def test_normalized_weight_scale(self):
        z=np.arange(4.)[:,None];dt=np.array([[1,0],[5,1],[5,3],[5,3]],float);ev=np.array([[1,0],[0,1],[0,0],[0,0]])
        a=fit_cause(z,dt,ev,[1,2,3,4],10.);b=fit_cause(z,dt,ev,[10,20,30,40],10.)
        self.assertTrue(np.allclose(a['coefficients'],b['coefficients']))
    def test_independent_poisson_fit(self):
        import statsmodels.api as sm
        rng=np.random.default_rng(44);n=150;z=rng.normal(size=(n,2));dt=np.column_stack([np.ones(n)*5,np.ones(n)*3])
        ev=np.zeros((n,2));ev[:25,0]=1;ev[25:45,1]=1
        m=fit_cause(z,dt,ev,np.ones(n),0.)
        X=np.column_stack([np.tile(np.eye(2),(n,1)),np.repeat(z,2,axis=0)])
        ref=sm.GLM(ev.reshape(-1),X,family=sm.families.Poisson(),offset=np.log(dt.reshape(-1))).fit()
        self.assertTrue(np.allclose(m['coefficients'],ref.params[2:],atol=1e-7))
        self.assertTrue(np.allclose(m['baseline_log_hazards'],ref.params[:2],atol=1e-7))

class RuntimeTests(unittest.TestCase):
    def test_hash(self):self.assertEqual(sha(MODEL_PATH),MODEL_SHA)
    def test_example_no_crp(self):self.assertNotIn('crp',example()['measurements'])
    def test_valid_both_panels(self):
        for panel in ('clinical','routine_no_crp'):self.assertEqual(run(example(panel))['status'],'calculated_research_only')
    def test_no_crp_placeholder(self):
        r=run(example());self.assertFalse(r['CRP_used']);self.assertFalse(r['imputation_applied'])
        self.assertNotIn('log_crp',[x['feature'] for x in r['feature_trace']])
    def test_crp_rejected(self):
        p=example();p['measurements']['crp']={'value':1};self.assertEqual(run(p)['status'],'blocked')
    def test_each_missing_measurement_blocks(self):
        for field in example()['measurements']:
            p=example();p['measurements'].pop(field)
            with self.subTest(field=field):self.assertEqual(run(p)['status'],'blocked')
    def test_country_endpoint_horizon(self):
        for field,value in [('country','RU'),('country','DE'),('endpoint','IHD'),('horizons_years',[20]),('horizons_years',[8]),('horizons_years',[True])]:
            p=example();p[field]=value
            with self.subTest(field=field,value=value):self.assertEqual(run(p)['status'],'blocked')
    def test_invalid_values(self):
        for v in (None,True,'<1',0,math.nan,math.inf,10**500):
            p=example();p['measurements']['albumin']['value']=v
            with self.subTest(value=str(v)[:12]):self.assertEqual(run(p)['status'],'blocked')
    def test_unknown_identifier(self):
        p=example();p['name']='synthetic';self.assertEqual(run(p)['status'],'blocked')
    def test_clinical_missing_blocks(self):
        p=example();p['clinical'].pop('cancer_history');self.assertEqual(run(p)['status'],'blocked')
    def test_wrong_unit_matrix_assay(self):
        for field,value in [('unit','ng/L'),('matrix','urine'),('assay','unknown')]:
            p=example();p['measurements']['creatinine'][field]=value;self.assertEqual(run(p)['status'],'blocked')
    def test_no_mutation(self):
        p=example();old=deepcopy(p);run(p);self.assertEqual(old,p)
    def test_equivalent_units(self):
        p=example();p['measurements']['albumin'].update(value=4.3,unit='g/dL')
        self.assertEqual(run(p)['probabilities'],run(example())['probabilities'])
    def test_checksum_failure(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'bad.json';f.write_text('{}')
            with self.assertRaises(ValueError):run(example(),bundle_path=f)
    def test_no_clinical_upgrade(self):
        r=run(example());self.assertFalse(r['clinical_use_ready']);self.assertIsNone(r['noninfectious_natural_probability']);self.assertIsNone(r['individual_uncertainty_interval'])
    def test_independent_vector_replay(self):
        b=load_bundle();rng=np.random.default_rng(32)
        for panel,m in b['models'].items():
            pp=m['preprocessing'];x=rng.uniform(pp['lower'],pp['upper'],(50,len(pp['features'])));df=pd.DataFrame(x,columns=pp['features'])
            reference={h:predict(df,m,h) for h in (1,4,5)}
            for i,rec in enumerate(df.to_dict('records')):
                out,_=score_engineered(rec,m,[1,4,5])
                for h,row in zip((1,4,5),out):
                    actual=[row['survival']]+[row['cause_specific_cif'][g] for g in GROUPS]
                    self.assertTrue(np.allclose(actual,reference[h][i],atol=1e-13,rtol=0))
    def test_changed_dimensions(self):
        m=deepcopy(load_bundle()['models']['clinical']);m['causes'][GROUPS[0]]['coefficients'].pop()
        pp=m['preprocessing'];f=dict(zip(pp['features'],pp['mean']))
        with self.assertRaises(ValueError):score_engineered(f,m,[5])

if __name__=='__main__':unittest.main()
