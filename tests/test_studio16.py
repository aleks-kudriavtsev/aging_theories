"""Synthetic integration tests, not new clinical validation."""
import hashlib
import http.client
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
from copy import deepcopy
from research.mortality.studio16 import registry, engine, report, server
from research.mortality.studio16.__main__ import batch
from research.mortality import replay_compact14 as compact, replay_transport12 as transport


class RegistryTests(unittest.TestCase):
    def test_five_views_three_artifacts(self):
        self.assertEqual(len(registry.catalog()['models']),5);self.assertEqual(len(registry.ARTIFACTS),3)
    def test_hard_fail_tampered_artifact(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/registry.ARTIFACTS['compact'][0];p.parent.mkdir();p.write_text('{}')
            with patch.object(registry,'ROOT',Path(d)),self.assertRaises(ValueError):registry.load_artifact('compact')
    def test_missing_not_ready(self):
        with patch.object(registry,'ROOT',Path('/nonexistent')):self.assertEqual(server.doctor()['status'],'not_ready')
    def test_no_clinical_upgrade(self):
        self.assertTrue(all(not e['clinical_use_ready'] for e in registry.catalog()['models']))
    def test_no_automatic_selection(self):self.assertFalse(registry.catalog()['default_model_selected'])
    def test_three_year_metrics_not_attached_to_one_year(self):
        self.assertNotIn('all_cause_auc',registry.evidence('compact4',[1])[0])
        self.assertIn('all_cause_auc',registry.evidence('compact4',[3])[0])
    def test_detailed_cohort_matches_horizon(self):
        a=registry.evidence('full8_detailed',[4,5]);self.assertEqual([r['n'] for r in a],[3013,2666])
    def test_coarse_never_validates_five_heads(self):
        self.assertEqual(registry.evidence('full8',[3])[0]['evaluated_death_groups'],3)
    def test_minimum_physical_determinations(self):self.assertEqual(registry.REGISTRY['compact4']['determinations'],5)
    def test_clinical_same_domain(self):self.assertEqual(registry.REGISTRY['clinical3']['required'],registry.FULL)
    def test_copy_not_global_mutation(self):
        a=registry.catalog();a['models'][0]['required'].clear();self.assertTrue(registry.REGISTRY['compact4']['required'])
    def test_unknown_model(self):
        with self.assertRaises(ValueError):registry.model_spec('latest_best')


class RuntimeTests(unittest.TestCase):
    def test_all_models_run(self):
        for mid in registry.REGISTRY:
            with self.subTest(mid=mid):self.assertEqual(engine.run(engine.example(mid))['status'],'calculated_research_only')
    def test_legacy_compact_identical(self):
        p=engine.example('compact4');expected=compact.run(compact.example())
        self.assertEqual(engine.run(p)['probabilities'],expected['probabilities'])
    def test_legacy_transport_identical(self):
        self.assertEqual(engine.run(engine.example('full8_detailed'))['probabilities'],transport.run(transport.example())['probabilities'])
    def test_many_synthetic_profiles_exact(self):
        rng=random.Random(20260916)
        for _ in range(60):
            p=engine.example('compact4');p['age_years']=rng.randint(40,79);p['sex']=rng.choice(['male','female'])
            for k,lo,hi in [('bmi',18,38),('sbp',95,185),('hba1c',4.5,10),('creatinine',.5,2.5),('uacr',1,300),('albumin',30,50)]:p['measurements'][k]['value']=rng.uniform(lo,hi)
            q=deepcopy(p);q.pop('model_id');q.update(schema_version='0.14.0',panel='compact4')
            self.assertEqual(engine.run(p)['probabilities'],compact.run(q)['probabilities'])
    def test_mass_and_monotonicity(self):
        for mid in registry.REGISTRY:
            r=engine.run(engine.example(mid));vals=[]
            for row in r['probabilities']:
                self.assertAlmostEqual(row['survival']+sum(row['cause_specific_cif'].values()),1)
                vals.append(row['all_cause_death_probability'])
            self.assertEqual(vals,sorted(vals))
    def test_missing_each_marker(self):
        for mid,e in registry.REGISTRY.items():
            for field in e['required']:
                p=engine.example(mid);p['measurements'].pop(field)
                with self.subTest(mid=mid,field=field):self.assertEqual(engine.run(p)['status'],'blocked')
    def test_extra_analyte_refused(self):
        p=engine.example('compact4');p['measurements']['tnf_alpha']={'value':2}
        self.assertEqual(engine.run(p)['status'],'blocked')
    def test_no_unrelated_inputs(self):
        p=engine.example('compact4');p['patient_name']='private'
        self.assertEqual(engine.run(p)['status'],'blocked')
    def test_biological_age_not_inferred(self):self.assertIsNone(engine.run(engine.example('compact4'))['biological_age'])
    def test_country_and_purpose(self):
        for key,value in [('country','RU'),('country','DE'),('dataset_context','current_patient'),('research_only',False),('acknowledge_limitations',False)]:
            p=engine.example('compact4');p[key]=value;self.assertEqual(engine.run(p)['status'],'blocked')
    def test_wrong_endpoint(self):
        p=engine.example('compact4');p['endpoint']='five_competing_death_groups';self.assertEqual(engine.run(p)['status'],'blocked')
    def test_wrong_horizon(self):
        for v in [[5],[20],[True],[1,1],[],[math.nan],None]:
            p=engine.example('compact4');p['horizons_years']=v;self.assertEqual(engine.run(p)['status'],'blocked')
    def test_invalid_numbers(self):
        for v in [None,True,math.nan,math.inf,'<3',-1,10**400]:
            p=engine.example('compact4');p['measurements']['creatinine']['value']=v
            self.assertEqual(engine.run(p)['status'],'blocked')
    def test_units_converted(self):
        p=engine.example('compact4');q=deepcopy(p);p['measurements']['albumin']['value']=40
        q['measurements']['albumin'].update(value=4,unit='g/dL')
        self.assertEqual(engine.run(p)['probabilities'],engine.run(q)['probabilities'])
    def test_urine_not_plasma(self):
        p=engine.example('compact4');p['measurements']['uacr']['matrix']='plasma';self.assertEqual(engine.run(p)['status'],'blocked')
    def test_serum_not_plasma(self):
        p=engine.example('compact4');p['measurements']['albumin']['matrix']='plasma';self.assertEqual(engine.run(p)['status'],'blocked')
    def test_creatinine_standardization(self):
        p=engine.example('compact4');p['measurements']['creatinine']['assay']='unknown';self.assertEqual(engine.run(p)['status'],'blocked')
    def test_HDL_exceeds_total(self):
        p=engine.example('full8');p['measurements']['hdl']['value']=250;self.assertEqual(engine.run(p)['status'],'blocked')
    def test_explicit_history(self):
        p=engine.example('compact4');p['clinical']['diabetes_history']=None;self.assertEqual(engine.run(p)['status'],'blocked')
    def test_no_mutation(self):
        p=engine.example('six6');copy=deepcopy(p);engine.run(p);self.assertEqual(p,copy)
    def test_request_hash_changes(self):
        p=engine.example('compact4');a=engine.run(p)['request_sha256'];p['age_years']=60
        self.assertNotEqual(a,engine.run(p)['request_sha256'])
    def test_clipping_exposed(self):
        p=engine.example('compact4');p['measurements']['uacr']['value']=10000;self.assertIn('log_uacr',engine.run(p)['clipped_features'])
    def test_bad_root(self):
        for p in [None,[],1,'abc',{}]:self.assertEqual(engine.run(p)['status'],'blocked')
    def test_never_fit_or_impute(self):
        r=engine.run(engine.example('compact4'));self.assertFalse(r['imputation_applied']);self.assertFalse(r['refitting_on_request'])
    def test_no_target_overclaim(self):
        r=engine.run(engine.example('full8_detailed'));self.assertIsNone(r['noninfectious_natural_probability']);self.assertIsNone(r['individual_uncertainty_interval'])


class ComparisonTests(unittest.TestCase):
    def test_same_profile_comparison(self):
        r=engine.compare(engine.example('full8'),['clinical3','compact4','six6','full8'])
        self.assertEqual(r['status'],'compared_research_only');self.assertEqual(len(r['results']),4)
    def test_no_average_no_winner(self):
        r=engine.compare(engine.example('full8'),['compact4','full8']);self.assertIsNone(r['automatic_winner']);self.assertIsNone(r['averaged_risk'])
    def test_projected_features_explicit(self):
        r=engine.compare(engine.example('full8'),['compact4','full8'])
        self.assertEqual(r['results'][0]['explicitly_unused_input_measurements'],['hdl','rdw','total_cholesterol','wbc'])
    def test_single_equals_comparison(self):
        p=engine.example('full8');r=engine.compare(p,['compact4','full8'])
        self.assertEqual(r['results'][1]['probabilities'],engine.run(p)['probabilities'])
    def test_no_missing_full_domain(self):
        self.assertEqual(engine.compare(engine.example('compact4'),['compact4','full8'])['status'],'blocked')
    def test_different_endpoint_refused(self):
        self.assertEqual(engine.compare(engine.example('full8'),['compact4','full8_detailed'])['status'],'blocked')
    def test_duplicate_and_unknown_models(self):
        for models in [['compact4'],['compact4','compact4'],['foo','compact4'],None,'all']:
            self.assertEqual(engine.compare(engine.example('full8'),models)['status'],'blocked')


class OutputTests(unittest.TestCase):
    def test_escape_report(self):
        s=report.render(engine.blocked([{'field':'<script>alert(1)</script>','code':'bad'}]))
        self.assertNotIn('<script>alert(1)</script>',s);self.assertIn('&lt;script&gt;',s)
    def test_save(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'out';report.save(engine.run(engine.example('compact4')),p)
            self.assertTrue((p/'report.html').is_file());self.assertTrue((p/'result.json').is_file())
    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):report.save(engine.blocked([]),Path(d))
    def test_batch_invalid_kept(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'in.jsonl';p.write_text(json.dumps(engine.example('compact4'))+'\n{bad\n')
            r=batch(p,Path(d)/'out.json');self.assertEqual(r['calculated'],1);self.assertEqual(r['blocked'],1)
    def test_cli_no_dependencies(self):
        with tempfile.TemporaryDirectory() as d:
            r=subprocess.run([sys.executable,'-S','-m','research.mortality.studio16','demo','--model','compact4','--out',str(Path(d)/'demo')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)


class HTTPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.srv=server.create_server(0);cls.thread=threading.Thread(target=cls.srv.serve_forever,daemon=True);cls.thread.start();cls.port=cls.srv.server_port
        cls.origin=f'http://127.0.0.1:{cls.port}'
        c=http.client.HTTPConnection('127.0.0.1',cls.port);c.request('GET','/');r=c.getresponse();html=r.read().decode();c.close()
        cls.token=re.search('name="workbench-token" content="([^"]+)"',html).group(1)
    @classmethod
    def tearDownClass(cls):cls.srv.shutdown();cls.srv.server_close();cls.thread.join()
    def request(self,path,method='GET',payload=None,headers=None,raw=None):
        c=http.client.HTTPConnection('127.0.0.1',self.port,timeout=10)
        h={'Content-Type':'application/json','Origin':self.origin,'X-Workbench-Token':self.token};h.update(headers or {})
        body=raw if raw is not None else json.dumps(payload) if payload is not None else None
        c.request(method,path,body=body,headers=h);r=c.getresponse();data=r.read();result=(r.status,dict(r.getheaders()),data);c.close();return result
    def test_actual_api(self):
        code,_,data=self.request('/api/predict','POST',engine.example('compact4'));self.assertEqual(code,200);self.assertIn('probabilities',json.loads(data)['result'])
    def test_actual_compare(self):
        code,_,data=self.request('/api/compare','POST',{'input':engine.example('full8'),'models':['compact4','full8']});self.assertEqual(code,200)
    def test_unsupported_country422(self):
        p=engine.example('compact4');p['country']='RU';self.assertEqual(self.request('/api/predict','POST',p)[0],422)
    def test_cross_origin(self):self.assertEqual(self.request('/api/predict','POST',{},headers={'Origin':'https://example.com'})[0],403)
    def test_bad_host(self):self.assertEqual(self.request('/api/health',headers={'Host':'localhost'})[0],403)
    def test_no_token(self):self.assertEqual(self.request('/api/predict','POST',{},headers={'X-Workbench-Token':''})[0],403)
    def test_invalid_json400(self):self.assertEqual(self.request('/api/predict','POST',raw='{bad')[0],400)
    def test_wrong_content_type(self):self.assertEqual(self.request('/api/predict','POST',{},headers={'Content-Type':'text/plain'})[0],415)
    def test_unknown_route(self):self.assertEqual(self.request('/etc/passwd')[0],404)
    def test_no_cache(self):
        code,h,data=self.request('/api/catalog');self.assertEqual(code,200);self.assertEqual(h['Cache-Control'],'no-store')
    def test_no_external_CSP(self):self.assertIn("connect-src 'self'",self.request('/')[1]['Content-Security-Policy'])
    def test_health(self):self.assertEqual(json.loads(self.request('/api/health')[2])['registered_views'],5)
    def test_unknown_example(self):self.assertEqual(self.request('/api/example?model_id=unknown')[0],400)
    def test_duplicate_json_key(self):self.assertEqual(self.request('/api/predict','POST',raw='{"model_id":"compact4","model_id":"full8"}')[0],400)
    def test_nonfinite_json(self):self.assertEqual(self.request('/api/predict','POST',raw='{"age_years":NaN}')[0],400)
    def test_no_external_binding(self):self.assertEqual(self.srv.server_address[0],'127.0.0.1')

if __name__=='__main__':unittest.main()
