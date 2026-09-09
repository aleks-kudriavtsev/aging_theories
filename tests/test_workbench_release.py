"""Release regression tests; fixtures are synthetic, not clinical evidence."""
from copy import deepcopy
import http.client
import json
from pathlib import Path
import re
import tempfile
import threading
import unittest
from unittest.mock import patch
from research.mortality.workbench import __version__
from research.mortality.workbench.runtime import example, run, strict_json, finite, schema, BUNDLE_SHA256
from research.mortality.workbench.server import create_server, readiness, check_assets
from research.mortality.workbench.__main__ import doctor, batch


class NumericHardeningTests(unittest.TestCase):
    def test_huge_positive_integer(self):
        self.assertFalse(finite(10**400))
    def test_huge_negative_integer(self):
        self.assertFalse(finite(-(10**400)))
    def test_direct_huge_age_blocks_without_exception(self):
        p=example();p['age_years']=10**400
        self.assertEqual(run(p)['status'],'blocked')
    def test_direct_huge_measurement_blocks_without_exception(self):
        p=example();p['measurements']['sbp']['value']=10**400
        self.assertEqual(run(p)['status'],'blocked')
    def test_direct_huge_horizon_blocks_without_exception(self):
        p=example();p['horizons_years']=[10**400]
        self.assertEqual(run(p)['status'],'blocked')
    def test_parser_rejects_large_integer(self):
        with self.assertRaises(ValueError):strict_json('{"n":'+'1'+'0'*400+'}')
    def test_parser_rejects_negative_large_integer(self):
        with self.assertRaises(ValueError):strict_json('{"n":-'+'1'+'0'*400+'}')
    def test_parser_rejects_exponent_overflow(self):
        with self.assertRaises(ValueError):strict_json('{"n":1e999}')
    def test_parser_preserves_finite_numbers(self):
        self.assertEqual(strict_json('{"n":1e-10,"i":2026}'),{'n':1e-10,'i':2026})
    def test_batch_preserves_invalid_numeric_record(self):
        with tempfile.TemporaryDirectory() as d:
            p=example();p['age_years']=10**400
            source=Path(d)/'batch.jsonl';source.write_text(json.dumps(p)+'\n'+json.dumps(example())+'\n')
            result=batch(source,Path(d)/'out')
            self.assertEqual((result['rows'],result['blocked'],result['calculated']),(2,1,1))
    def test_input_contract_not_silently_versioned(self):
        self.assertEqual(example()['schema_version'],'0.6.0')
        self.assertEqual(schema()['input_schema_version'],'0.6.0')
        self.assertEqual(run(example())['application_version'],__version__)
    def test_frozen_model_unchanged(self):
        self.assertEqual(BUNDLE_SHA256,'7b287f3f4778b2127f2226f79d2e503e1e39087306a172e5f524a561964b07d0')


class ReadyTests(unittest.TestCase):
    def test_ready_remains_research_only(self):
        self.assertEqual(readiness()['status'],'ready_research_only')
        self.assertFalse(readiness()['clinical_use_ready'])
    def test_integrity_failure_is_unavailable(self):
        with patch('research.mortality.workbench.server.load_bundle',side_effect=ValueError('internal path')):
            result=readiness()
        self.assertEqual(result['status'],'not_ready')
        self.assertNotIn('internal path',json.dumps(result))
    def test_missing_ui_file_blocks_readiness(self):
        with tempfile.TemporaryDirectory() as d,patch('research.mortality.workbench.server.STATIC',Path(d)):
            self.assertEqual(readiness()['status'],'not_ready')
            with self.assertRaises(ValueError):check_assets()
    def test_doctor_does_not_hide_missing_ui(self):
        with tempfile.TemporaryDirectory() as d,patch('research.mortality.workbench.server.STATIC',Path(d)):
            with self.assertRaises(ValueError):doctor()
    def test_complete_static_files_reported(self):
        self.assertTrue(doctor()['static_assets_present'])


class HTTPReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server=create_server(0);cls.thread=threading.Thread(target=cls.server.serve_forever,daemon=True);cls.thread.start()
        cls.port=cls.server.server_port
        c=http.client.HTTPConnection('127.0.0.1',cls.port,timeout=10);c.request('GET','/')
        h=c.getresponse().read().decode();c.close();cls.token=re.search(r'name="workbench-token" content="([^"]+)"',h)[1]
    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown();cls.server.server_close();cls.thread.join(5)
    def request(self,path,body=None):
        c=http.client.HTTPConnection('127.0.0.1',self.port,timeout=10)
        c.request('POST' if body is not None else 'GET',path,body=body,
          headers={'Content-Type':'application/json','X-Workbench-Token':self.token})
        r=c.getresponse();result=(r.status,json.loads(r.read()));c.close();return result
    def test_health_endpoint(self):
        status,body=self.request('/api/health');self.assertEqual(status,200);self.assertFalse(body['clinical_use_ready'])
    def test_unavailable_health_status(self):
        with patch('research.mortality.workbench.server.load_bundle',side_effect=ValueError('broken')):
            status,body=self.request('/api/health')
        self.assertEqual(status,503);self.assertEqual(body['status'],'not_ready')
    def test_unavailable_model_is_not_reported_as_bad_patient_data(self):
        with patch('research.mortality.workbench.server.load_bundle',side_effect=ValueError('broken')):
            status,body=self.request('/api/predict',json.dumps(example()))
        self.assertEqual(status,503);self.assertEqual(body['error'],'application_not_ready')
    def test_huge_json_does_not_crash_http_handler(self):
        p=example();p['age_years']=10**400
        status,body=self.request('/api/predict',json.dumps(p));self.assertEqual(status,400)
        status,body=self.request('/api/predict',json.dumps(example()));self.assertEqual(status,200)
        self.assertEqual(body['result']['status'],'calculated_research_only')


if __name__=='__main__':unittest.main()
