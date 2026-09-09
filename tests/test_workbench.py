"""Offline application tests. Synthetic examples are not clinical observations."""
from copy import deepcopy
from hashlib import sha256
import http.client
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import threading
import unittest
from research.mortality.workbench.runtime import (BUNDLE_PATH,BUNDLE_SHA256,FIELDS,MODELS,
    DISCLAIMER,MAX_BODY,example,normalize,run,load_bundle,strict_json,engineer,score_features)
from research.mortality.workbench.report import render,write_package
from research.mortality.workbench.server import create_server
from research.mortality.workbench.__main__ import batch,doctor
ROOT=Path(__file__).resolve().parents[1]

class InputTests(unittest.TestCase):
    def block(self,mutate,field=None):
        p=example('M4_cystatinC');mutate(p);r=run(p)
        self.assertEqual(r['status'],'blocked');self.assertIsNone(r['probabilities'])
        self.assertIs(r['clinical_use_ready'],False)
        if field:self.assertTrue(any(e['field']==field for e in r['errors']),r['errors'])
    def test_example_is_synthetic(self):
        self.assertEqual(run(example())['declaration'],'synthetic_not_a_person')
    def test_all_six_models(self):
        for model in MODELS:
            with self.subTest(model=model):self.assertEqual(run(example(model))['status'],'calculated_research_only')
    def test_no_clinical_use_flag(self):
        self.block(lambda p:p.update(research_only=False),'research_only')
    def test_ack_required(self):
        self.block(lambda p:p.update(acknowledge_limitations=1),'acknowledge_limitations')
    def test_country_restriction(self):
        for country in ['RU','DE','USA',None]:
            with self.subTest(country=country):self.block(lambda p:p.update(country=country),'country')
    def test_no_noninfectious_endpoint_rename(self):
        self.block(lambda p:p.update(endpoint='noninfectious_natural_mortality'),'endpoint')
    def test_twenty_year_and_intermediate_horizons(self):
        for h in [[20],[2],[1,1],[],[True],[1,float('nan')],None,[[1]]]:
            with self.subTest(h=h):self.block(lambda p:p.update(horizons_years=h),'horizons_years')
    def test_age_limits_and_types(self):
        for age in [39,80,True,'55',None,math.inf,[],{}]:
            with self.subTest(age=age):self.block(lambda p:p.update(age_years=age),'age_years')
    def test_age_boundary_allowed(self):
        for a in [40,79]:
            p=example();p['age_years']=a;self.assertEqual(run(p)['status'],'calculated_research_only')
    def test_sex_not_imputed(self):
        for sex in [None,'unknown',1,[]]:
            with self.subTest(sex=sex):self.block(lambda p:p.update(sex=sex),'sex')
    def test_model_not_inferred(self):
        for m in [None,[],{},'custom']:
            with self.subTest(m=m):self.block(lambda p:p.update(model=m),'model')
    def test_no_missing_measurements(self):
        for key in example()['measurements']:
            with self.subTest(key=key):self.block(lambda p:p['measurements'].pop(key),'measurements.'+key)
    def test_no_LOD_number_substitution(self):
        for x in [None,'<3',0,True,float('nan'),float('inf')]:
            with self.subTest(x=x):self.block(lambda p:p['measurements']['hs_ctnt'].update(value=x))
    def test_no_silent_clinical_imputation(self):
        for value in [None,0,1,'false']:
            with self.subTest(value=value):self.block(lambda p:p['clinical'].update(cvd_history=value),'clinical.cvd_history')
    def test_smoking_code_required(self):
        self.block(lambda p:p['clinical'].update(smoking='sometimes'),'clinical.smoking')
    def test_rdw_SD_not_CV(self):
        self.block(lambda p:p['measurements']['rdw'].update(unit='fL'),'measurements.rdw.unit')
    def test_troponin_I_not_T(self):
        self.block(lambda p:p['measurements'].update(hs_ctni=p['measurements'].pop('hs_ctnt')))
    def test_BNP_not_NTproBNP(self):
        self.block(lambda p:p['measurements'].update(bnp=p['measurements'].pop('ntprobnp')))
    def test_hba1c_IFCC_requires_explicit_conversion_outside_app(self):
        self.block(lambda p:p['measurements']['hba1c'].update(unit='mmol/mol'),'measurements.hba1c.unit')
    def test_matrix_required(self):
        self.block(lambda p:p['measurements']['uacr'].update(matrix='serum'),'measurements.uacr.matrix')
    def test_idms_required(self):
        self.block(lambda p:p['measurements']['creatinine'].update(assay='not_reported'),'measurements.creatinine.assay')
    def test_protein_assay_not_interchangeable(self):
        self.block(lambda p:p['measurements']['cystatin_c'].update(assay='new_TRFIA'),'measurements.cystatin_c.assay')
    def test_all_sources_structurally_checked(self):
        for key in ['clinical','measurements']:
            with self.subTest(key=key):self.block(lambda p:p.update({key:[]}))
    def test_identifiers_not_allowed(self):
        self.block(lambda p:p.update(patient_name='SYNTHETIC'),'patient_name')
    def test_postmortem_predictor_not_allowed(self):
        self.block(lambda p:p.update(dead=True),'dead')
    def test_hdl_total_consistency(self):
        self.block(lambda p:p['measurements']['hdl'].update(value=250),'measurements.hdl')
    def test_no_unsupported_dataset(self):
        self.block(lambda p:p.update(dataset_context='modern_population_clinical'),'dataset_context')
    def test_not_an_object(self):
        for x in [None,[],1,'x']:
            self.assertEqual(run(x)['status'],'blocked')
    def test_input_not_mutated(self):
        p=example();saved=deepcopy(p);run(p);self.assertEqual(p,saved)
    def test_ignored_values_explicit(self):
        p=example();p['measurements']['ntprobnp']={'value':'unused'};r=run(p)
        self.assertEqual(r['status'],'calculated_research_only')
        self.assertTrue(any(w['code']=='measurements_not_used_by_selected_model' for w in r['warnings']))
    def test_M0_ignores_extraneous_partial_clinical(self):
        p=example('M0_age_sex');p['clinical']={'smoking':'current'}
        self.assertEqual(run(p)['status'],'calculated_research_only')

class NumericTests(unittest.TestCase):
    def test_bundle_hash(self):self.assertEqual(sha256(BUNDLE_PATH.read_bytes()).hexdigest(),BUNDLE_SHA256)
    def test_tampered_bundle_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'model.json';p.write_bytes(BUNDLE_PATH.read_bytes()+b' ')
            with self.assertRaises(ValueError):run(example(),bundle_path=p)
    def test_risks_finite_monotonic(self):
        for m in MODELS:
            r=run(example(m));p=[x['probability'] for x in r['probabilities']]
            self.assertEqual(p,sorted(p));self.assertTrue(all(0<=x<=1 for x in p))
    def test_explicit_unit_equivalence(self):
        orig=example('M4_cystatinC');expected=run(orig)['probabilities']
        for k,spec in FIELDS.items():
            for unit,mult in spec['conversions'].items():
                with self.subTest(k=k,unit=unit):
                    p=deepcopy(orig);p['measurements'][k]['unit']=unit;p['measurements'][k]['value']/=mult
                    actual=run(p)['probabilities']
                    for a,b in zip(actual,expected):self.assertAlmostEqual(a['probability'],b['probability'],places=13)
    def test_creatinine_no_double_1999_calibration(self):
        from research.mortality.engine import egfr_ckd_epi_2021
        r=run(example());self.assertEqual(r['egfr_ckd_epi_2021'],egfr_ckd_epi_2021(55,'female',.9))
    def test_no_imputation_feature_indicators(self):
        r=run(example());self.assertFalse(r['imputation_applied']);self.assertNotIn('median',str(r['feature_trace']))
    def test_winsorization_visible(self):
        p=example();p['measurements']['crp']['value']=900
        r=run(p);self.assertIn('log_crp',r['clipped_features'])
        row=next(x for x in r['feature_trace'] if x['feature']=='log_crp')
        self.assertNotEqual(row['raw_value'],row['used_value'])
    def test_no_false_individual_CI_or_cause_fractions(self):
        r=run(example());self.assertIsNone(r['individual_uncertainty_interval']);self.assertIsNone(r['cause_specific_probabilities'])
    def test_strict_json_duplicate(self):
        with self.assertRaises(ValueError):strict_json('{"age_years":40,"age_years":55}')
    def test_strict_json_nonfinite(self):
        for s in ['{"a":NaN}','{"a":Infinity}','{"a":-Infinity}']:
            with self.assertRaises(ValueError):strict_json(s)
    def test_json_too_big(self):
        with self.assertRaises(ValueError):strict_json(b' '*(MAX_BODY+1))
    def test_json_invalid_encoding(self):
        with self.assertRaises(ValueError):strict_json(b'\xff')
    def test_doctor_does_not_claim_clinical_readiness(self):
        r=doctor();self.assertTrue(r['runtime_ready']);self.assertFalse(r['clinical_use_ready']);self.assertFalse(r['cause_specific_models_ready'])

class DeliveryTests(unittest.TestCase):
    def test_standalone_report_escapes_strings(self):
        r=run(example());r['warnings'].append({'code':'<script>alert(1)</script>'});html=render(r)
        self.assertNotIn('<script>alert',html);self.assertIn('&lt;script&gt;',html)
    def test_report_has_no_remote_assets(self):
        h=render(run(example()));self.assertNotIn('<script',h);self.assertNotIn('src="http',h)
    def test_no_overwrite_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'r';write_package(run(example()),p)
            before=(p/'report.json').read_bytes()
            with self.assertRaises(ValueError):write_package(run(example('M4_cystatinC')),p)
            self.assertEqual(before,(p/'report.json').read_bytes())
    def test_blocked_report_no_probability(self):
        p=example();p['country']='RU';r=run(p);h=render(r)
        self.assertIn('Расчёт заблокирован',h);self.assertNotIn('class="card"',h)
    def test_batch_retains_invalid_rows(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.jsonl';bad=example();bad['country']='RU'
            p.write_text(json.dumps(example())+'\n'+json.dumps(bad)+'\n{invalid}\n')
            r=batch(p,Path(d)/'results');self.assertEqual(r['rows'],3);self.assertEqual(r['blocked'],2);self.assertEqual(r['calculated'],1)
    def test_batch_refuses_empty(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'empty';p.write_text('')
            with self.assertRaises(ValueError):batch(p,Path(d)/'out')
    def test_no_site_packages_required(self):
        r=subprocess.run([sys.executable,'-S','-m','research.mortality.workbench','doctor'],cwd=ROOT,capture_output=True,text=True,timeout=10)
        self.assertEqual(r.returncode,0,r.stderr);self.assertTrue(json.loads(r.stdout)['runtime_ready'])
    def test_cli_demo_then_predict(self):
        with tempfile.TemporaryDirectory() as d:
            r=subprocess.run([sys.executable,'-S','-m','research.mortality.workbench','demo','--out',d+'/one'],cwd=ROOT,capture_output=True,text=True,timeout=10)
            self.assertEqual(r.returncode,0,r.stderr)
            r=subprocess.run([sys.executable,'-S','-m','research.mortality.workbench','predict',d+'/one/synthetic_input.json','--out',d+'/two'],cwd=ROOT,capture_output=True,text=True,timeout=10)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual((Path(d)/'one/report.json').read_bytes(),(Path(d)/'two/report.json').read_bytes())
    def test_cli_blocks_wrong_country_with_exit_code(self):
        with tempfile.TemporaryDirectory() as d:
            p=example();p['country']='DE';src=Path(d)/'input.json';src.write_text(json.dumps(p))
            r=subprocess.run([sys.executable,'-S','-m','research.mortality.workbench','predict',str(src),'--out',d+'/out'],cwd=ROOT,capture_output=True,text=True,timeout=10)
            self.assertEqual(r.returncode,2);self.assertIsNone(json.loads((Path(d)/'out/report.json').read_text())['probabilities'])

class ServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server=create_server(0);cls.thread=threading.Thread(target=cls.server.serve_forever,daemon=True);cls.thread.start()
        cls.port=cls.server.server_port
        conn=http.client.HTTPConnection('127.0.0.1',cls.port,timeout=10);conn.request('GET','/')
        res=conn.getresponse();html=res.read().decode();conn.close()
        cls.token=re.search(r'name="workbench-token" content="([^"]+)"',html)[1]
    @classmethod
    def tearDownClass(cls):cls.server.shutdown();cls.server.server_close();cls.thread.join(5)
    def req(self,path='/',method='GET',body=None,headers=None):
        c=http.client.HTTPConnection('127.0.0.1',self.port,timeout=10)
        c.request(method,path,body=body,headers=headers or {});r=c.getresponse();out=(r.status,dict(r.getheaders()),r.read());c.close();return out
    def headers(self):return {'Content-Type':'application/json','X-Workbench-Token':self.token}
    def test_loopback_binding(self):self.assertEqual(self.server.server_address[0],'127.0.0.1')
    def test_home_and_static(self):
        for p in ['/','/app.js','/style.css','/api/schema']:
            with self.subTest(p=p):self.assertEqual(self.req(p)[0],200)
    def test_no_filesystem_routes(self):
        for p in ['/../../etc/passwd','/research/mortality/nhanes04/model_bundle.json','/.git/config','//example.org']:
            with self.subTest(p=p):self.assertEqual(self.req(p)[0],404)
    def test_bad_host(self):self.assertEqual(self.req(headers={'Host':'attacker.test'})[0],403)
    def test_token_required(self):self.assertEqual(self.req('/api/predict','POST',json.dumps(example()),{'Content-Type':'application/json'})[0],403)
    def test_cross_origin_refused(self):
        h=self.headers();h['Origin']='https://attacker.test'
        self.assertEqual(self.req('/api/predict','POST',json.dumps(example()),h)[0],403)
    def test_content_type_required(self):
        h=self.headers();h['Content-Type']='text/plain'
        self.assertEqual(self.req('/api/predict','POST','{}',h)[0],415)
    def test_json_duplicates_rejected(self):self.assertEqual(self.req('/api/predict','POST','{"a":1,"a":2}',self.headers())[0],400)
    def test_negative_length(self):
        h=self.headers();h['Content-Length']='-1'
        self.assertEqual(self.req('/api/predict','POST','',h)[0],413)
    def test_oversized_body_before_read(self):
        h=self.headers();h['Content-Length']=str(MAX_BODY+1)
        self.assertEqual(self.req('/api/predict','POST','{}',h)[0],413)
    def test_security_headers_and_no_cors(self):
        status,h,b=self.req();self.assertEqual(h['Cache-Control'],'no-store');self.assertEqual(h['X-Frame-Options'],'DENY');self.assertNotIn('Access-Control-Allow-Origin',h)
    def test_post_result_matches_cli_runtime(self):
        p=example();status,h,b=self.req('/api/predict','POST',json.dumps(p),self.headers())
        self.assertEqual(status,200);self.assertEqual(json.loads(b)['result'],run(p))
    def test_API_blocked_semantics(self):
        p=example();p['country']='DE';s,h,b=self.req('/api/predict','POST',json.dumps(p),self.headers())
        self.assertEqual(s,422);self.assertIsNone(json.loads(b)['result']['probabilities'])
    def test_no_port_range_bug(self):
        for value in [True,-1,65536,'8765']:
            with self.assertRaises(ValueError):create_server(value)

if __name__=='__main__':unittest.main()
