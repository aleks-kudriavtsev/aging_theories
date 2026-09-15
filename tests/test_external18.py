from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from research.mortality.external18 import runner as r
from research.mortality.external18 import demo
from research.mortality.external18.precision import components, analyse
from research.mortality.studio16.engine import run, example

class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name)/'demo';self.c=demo.generate(self.p)
    def check(self):return r.preflight(self.c,self.p/'protocol.txt',self.p/'dictionary.json')
    def test_synthetic_ready_not_clinical(self):
        v=self.check();self.assertTrue(v['ready']);self.assertFalse(v['clinical_use_ready']);self.assertFalse(v['participant_data_opened'])
    def test_real_no_scope_blocks(self):
        self.c['data_origin']='authorized_cohort';v=self.check();self.assertFalse(v['ready']);self.assertIn('approval_document_missing_or_changed',v['errors'])
    def test_protocol_mutation(self):
        (self.p/'protocol.txt').write_text('changed');self.assertIn('protocol_integrity',self.check()['errors'])
    def test_dictionary_mutation(self):
        (self.p/'dictionary.json').write_text('{}');self.assertIn('dictionary_integrity',self.check()['errors'])
    def test_wrong_country(self):
        self.c['country']='US';self.assertFalse(self.check()['ready'])
    def test_no_endpoint_substitution(self):
        self.c['endpoint']='cardiovascular_events';self.assertFalse(self.check()['ready'])
    def test_delayed_entry_not_silently_ignored(self):
        self.c['delayed_entry']=True;self.assertFalse(self.check()['ready'])
    def test_model_hash(self):
        self.c['model_sha256']='0'*64;self.assertFalse(self.check()['ready'])
    def test_real_attestation_required_individually(self):
        self.c.update(data_origin='authorized_cohort',bootstrap_replicates=1000)
        ap=self.p/'permission.txt';ap.write_text('SYNTHETIC TEST PERMISSION NOT A REAL APPROVAL')
        self.c['approval_sha256']=r.digest(ap)
        self.c.update(owner_review_reference='TEST_ONLY',mapping_review_reference='TEST_ONLY',execution_environment_reference='TEST_ONLY')
        self.c['attestations']=dict.fromkeys(r.ATTESTATIONS,True)
        self.assertTrue(r.preflight(self.c,self.p/'protocol.txt',self.p/'dictionary.json',ap)['ready'])
        for name in r.ATTESTATIONS:
            self.c['attestations'][name]=False
            self.assertFalse(r.preflight(self.c,self.p/'protocol.txt',self.p/'dictionary.json',ap)['ready'])
            self.c['attestations'][name]=True
    def test_late_horizon_not_allowed(self):
        self.c['horizon_years']=10;self.assertFalse(self.check()['ready'])
    def test_declaration_not_access_token(self):self.assertIn('NOT_legal',self.check()['scope_check'])

class DataTests(unittest.TestCase):
    def write(self,rows):
        t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);p=Path(t.name)/'x.jsonl'
        p.write_text(''.join(json.dumps(x)+'\n' for x in rows));return p
    def test_hba1c_ifcc(self):
        d={'unit':'mmol/mol','matrix':'whole_blood','method_id':'test','standardization':'IFCC_traceable'}
        self.assertAlmostEqual(r.measurement(48,'hba1c',d),6.54304)
    def test_hba1c_wrong_standard(self):
        d={'unit':'mmol/mol','matrix':'whole_blood','method_id':'test','standardization':'unknown'}
        with self.assertRaises(ValueError):r.measurement(48,'hba1c',d)
    def test_matrices_not_equivalent(self):
        d=demo.dictionary()['measurements']['albumin'];d['matrix']='plasma'
        with self.assertRaises(ValueError):r.measurement(42,'albumin',d)
    def test_uacr_units(self):
        d=demo.dictionary()['measurements']['uacr'];d['unit']='mg/mmol'
        self.assertAlmostEqual(r.measurement(3,'uacr',d),26.52)
    def test_creatinine_units(self):
        d=demo.dictionary()['measurements']['creatinine'];d['unit']='umol/L'
        self.assertAlmostEqual(r.measurement(88.4,'creatinine',d),1)
    def test_measured_zero_not_missing(self):
        d=demo.dictionary()['measurements']['uacr']
        with self.assertRaises(ValueError):r.measurement(0,'uacr',d)
    def test_LOD_string_not_imputed(self):
        with self.assertRaises(ValueError):r.measurement('<3','uacr',demo.dictionary()['measurements']['uacr'])
    def test_duplicate_record_blocks(self):
        rows=demo.rows(3);rows[1]['record_id']=rows[0]['record_id']
        with self.assertRaises(ValueError):r.read_cohort(self.write(rows),demo.dictionary(),'iid')
    def test_new_sensitive_column_blocks(self):
        rows=demo.rows(3);rows[0]['name']='EXAMPLE_NAME'
        with self.assertRaises(ValueError):r.read_cohort(self.write(rows),demo.dictionary(),'iid')
    def test_common_domain_missing_retained(self):
        rows=r.read_cohort(self.write(demo.rows(10)),demo.dictionary(),'iid')
        self.assertEqual(len(rows),10);self.assertEqual(sum(x['complete'] for x in rows),9)
    def test_iid_not_survey_weights(self):
        rows=demo.rows(3);rows[0]['weight']=3
        with self.assertRaises(ValueError):r.read_cohort(self.write(rows),demo.dictionary(),'iid')
    def test_invalid_death_unknown_not_alive(self):
        rows=demo.rows(3);rows[0]['death']=None
        with self.assertRaises(ValueError):r.read_cohort(self.write(rows),demo.dictionary(),'iid')
    def test_zero_time_requires_adjudication(self):
        rows=demo.rows(3);rows[0]['followup_years']=0
        with self.assertRaises(ValueError):r.read_cohort(self.write(rows),demo.dictionary(),'iid')

class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.p=Path(cls.temp.name)/'demo';cls.config=demo.generate(cls.p)
        cls.rows=r.read_cohort(cls.p/'synthetic.jsonl',demo.dictionary(),'iid');cls.result=r.evaluate(cls.config,cls.rows)
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def test_pipeline_is_not_external_evidence(self):self.assertEqual(self.result['status'],'synthetic_pipeline_check')
    def test_country_preserved(self):self.assertEqual(self.result['country'],'DE');self.assertFalse(self.result['country_substituted_with_US'])
    def test_no_individual_outputs(self):
        s=json.dumps(self.result);self.assertNotIn('SYNTHETIC_0001',s);self.assertFalse(self.result['predictions_exported'])
    def test_no_coefficient_change(self):self.assertFalse(self.result['risk_coefficients_refitted'])
    def test_same_participants(self):
        self.assertEqual(self.result['metrics']['clinical3']['n'],self.result['metrics']['compact4']['n'])
    def test_early_censor_retained(self):self.assertEqual(self.result['metrics']['compact4']['early_censored_n'],1)
    def test_no_hidden_validation_statistics(self):self.assertFalse(self.result['cause_specific_predictions_exported'])
    def test_reproducible_bootstrap(self):
        other=r.evaluate(self.config,self.rows);self.assertEqual(other,self.result)
    def test_small_events_suppressed(self):
        rows=deepcopy(self.rows)
        for x in rows:x['death']=False;x['followup_years']=4
        self.assertEqual(r.evaluate(self.config,rows)['status'],'aggregate_output_suppressed')
    def test_public_UI_still_blocks_Russia(self):
        p=example('compact4');p['country']='RU';self.assertEqual(run(p)['status'],'blocked')
    def test_resampling_stratified(self):
        rows=[{'weight':1,'stratum':s,'cluster':c} for s in ['a','b'] for c in ['1','2']]
        w=r.resample_weights(rows,'stratified_psu',np.random.default_rng(2));self.assertEqual(w.sum(),4);self.assertEqual(set(w),{0,2})

class PrecisionTests(unittest.TestCase):
    def test_repeat_variance_hand(self):
        v=components(np.tile([9.,11.],(2,2,1)))
        self.assertEqual(v['variance_components']['repeatability'],2);self.assertEqual(v['within_lab_variance'],2)
    def test_between_day_hand(self):
        a=np.array([np.full((2,2),9),np.full((2,2),11)])
        v=components(a);self.assertEqual(v['variance_components']['between_day'],2)
    def test_wrong_shape(self):
        with self.assertRaises(ValueError):components(np.ones((2,2)))
    def test_variance_scaling(self):
        rng=np.random.default_rng(18);a=rng.normal(100,1,(5,2,2))
        self.assertAlmostEqual(components(2*a)['within_lab_variance'],4*components(a)['within_lab_variance'])
    def test_negative_truncation_exposed(self):
        v=components(np.tile([9.,11.],(2,2,1)));self.assertIn('between_run',v['truncated_components'])
    def doc(self):
        d={'schema':'precision18.1','data_origin':'synthetic','analyte':'SYNTHETIC_CONTROL','unit':'unit','matrix':'synthetic',
           'reference_method':'A','candidate_method':'B','reference_lot':'TEST','candidate_lot':'TEST','observations':[]}
        for day in range(1,4):
            for run in [1,2]:
                for rep in [1,2]:d['observations'].append({'level':'TEST','day':day,'run':run,'replicate':rep,'reference':20+day+rep/2,'candidate':30+day+rep})
        return d
    def test_no_automatic_lambda(self):
        a=analyse(self.doc());self.assertIsNone(a['automatic_Deming_lambda']);self.assertFalse(a['clinical_use_ready'])
    def test_no_missing_replicate(self):
        d=self.doc();d['observations'].pop()
        with self.assertRaises(ValueError):analyse(d)
    def test_no_duplicate_replicate(self):
        d=self.doc();d['observations'].append(d['observations'][0])
        with self.assertRaises(ValueError):analyse(d)

if __name__=='__main__':unittest.main()
