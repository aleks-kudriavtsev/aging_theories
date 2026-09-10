"""Synthetic tests; no participant data and no clinical validation claims."""
from copy import deepcopy
from hashlib import sha256
import json,math,subprocess,sys,tempfile,unittest
from pathlib import Path
from research.mortality.replay_compact14 import (run,example,score_engineered,load_bundle,MODEL_PATH,MODEL_SHA)
from research.mortality.workbench.runtime import strict_json,engineer

class CompactRuntimeTests(unittest.TestCase):
    def test_checksum(self):self.assertEqual(sha256(MODEL_PATH.read_bytes()).hexdigest(),MODEL_SHA)
    def test_actual_frozen_scope(self):
        b=load_bundle();self.assertEqual(b['development_n'],15280);self.assertFalse(b['clinical_use_ready'])
    def test_no_extra_assays_in_example(self):
        self.assertEqual(set(example()['measurements']),{'bmi','sbp','hba1c','creatinine','uacr','albumin'})
    def test_example_runs(self):self.assertEqual(run(example())['status'],'calculated_research_only')
    def test_not_a_person(self):self.assertEqual(run(example())['declaration'],'synthetic_not_a_person')
    def test_clinical_readiness_false(self):self.assertFalse(run(example())['clinical_use_ready'])
    def test_no_false_four_physical_assays(self):self.assertEqual(run(example())['minimum_underlying_lab_determinations'],5)
    def test_three_causes_only(self):
        for row in run(example())['probabilities']:self.assertEqual(set(row['cause_specific_cif']),{'heart_diseases','malignant_neoplasms','other_or_unknown'})
    def test_mass_conservation(self):
        for row in run(example())['probabilities']:self.assertAlmostEqual(row['survival']+sum(row['cause_specific_cif'].values()),1,places=14)
    def test_monotonicity(self):
        p=run(example())['probabilities'];self.assertLess(p[0]['all_cause_death_probability'],p[1]['all_cause_death_probability'])
    def test_missing_laboratory_inputs(self):
        for key in ['hba1c','creatinine','uacr','albumin']:
            p=example();p['measurements'].pop(key)
            with self.subTest(key=key):self.assertEqual(run(p)['status'],'blocked')
    def test_no_unused_values_accepted(self):
        for key in ['crp','rdw','wbc','total_cholesterol','hdl','TNF_alpha']:
            p=example();p['measurements'][key]={'value':1}
            with self.subTest(key=key):self.assertEqual(run(p)['status'],'blocked')
    def test_no_missing_history(self):
        p=example();p['clinical'].pop('smoking');self.assertEqual(run(p)['status'],'blocked')
    def test_horizons(self):
        for h in [[5],[8],[20],[True],[1,1],[]]:
            p=example();p['horizons_years']=h
            with self.subTest(h=h):self.assertEqual(run(p)['status'],'blocked')
    def test_country_restricted(self):
        for c in ['RU','DE',None]:
            p=example();p['country']=c
            with self.subTest(c=c):self.assertEqual(run(p)['status'],'blocked')
    def test_age_sex_limits(self):
        for k,v in [('age_years',39),('age_years',80),('sex','unknown'),('sex',None)]:
            p=example();p[k]=v;self.assertEqual(run(p)['status'],'blocked')
    def test_endpoint_not_renamed(self):
        p=example();p['endpoint']='noninfectious';self.assertEqual(run(p)['status'],'blocked')
    def test_ack_required(self):
        p=example();p['acknowledge_limitations']=False;self.assertEqual(run(p)['status'],'blocked')
    def test_tampered_model(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json';p.write_text('{}')
            with self.assertRaises(ValueError):run(example(),bundle_path=p)
    def test_units(self):
        a=example();b=deepcopy(a);b['measurements']['albumin'].update(value=4.2,unit='g/dL')
        a['measurements']['albumin'].update(value=42.,unit='g/L')
        self.assertEqual(run(a)['probabilities'],run(b)['probabilities'])
    def test_wrong_units_matrix_method(self):
        for key,field,value in [('uacr','matrix','serum'),('albumin','unit','mg/L'),('creatinine','assay','unknown')]:
            p=example();p['measurements'][key][field]=value;self.assertEqual(run(p)['status'],'blocked')
    def test_invalid_numbers(self):
        for v in [None,True,math.nan,math.inf,'<3',0,10**400]:
            p=example();p['measurements']['uacr']['value']=v
            with self.subTest(v=str(v)[:10]):self.assertEqual(run(p)['status'],'blocked')
    def test_input_not_changed(self):
        p=example();old=deepcopy(p);run(p);self.assertEqual(p,old)
    def test_no_false_individual_CI(self):self.assertIsNone(run(example())['individual_uncertainty_interval'])
    def test_clipping_visible(self):
        p=example();p['measurements']['uacr']['value']=10000
        self.assertTrue(any(r['clipped'] for r in run(p)['feature_trace']))
    def test_missing_engineered_feature(self):
        with self.assertRaises(ValueError):score_engineered({},load_bundle()['model'],[3])
    def test_analytic_constant_hazard(self):
        m={'preprocessing':{'features':['x'],'lower':[0],'upper':[1],'mean':[0],'sd':[1]},'interval_ends_years':[5.,8.],
           'causes':{k:{'coefficients':[0.],'baseline_log_hazards':[math.log(.01),math.log(.02)]} for k in load_bundle()['groups']}}
        r,_=score_engineered({'x':0},m,[3]);self.assertAlmostEqual(r[0]['all_cause_death_probability'],1-math.exp(-.15))
        self.assertAlmostEqual(r[0]['cause_specific_cif']['other_or_unknown'],(1-math.exp(-.15))*3/5)
    def test_cli_without_dependencies(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'demo.json';r=subprocess.run([sys.executable,'-S','-m','research.mortality.replay_compact14','--demo','--out',str(p)],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr);self.assertEqual(json.loads(p.read_text())['status'],'calculated_research_only')
    def test_no_output_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'present';p.write_text('unchanged')
            r=subprocess.run([sys.executable,'-S','-m','research.mortality.replay_compact14','--demo','--out',str(p)],capture_output=True)
            self.assertNotEqual(r.returncode,0);self.assertEqual(p.read_text(),'unchanged')
    def test_structural_errors(self):
        for p in [None,[],1,'input']:
            self.assertEqual(run(p)['status'],'blocked')

if __name__=='__main__':unittest.main()
