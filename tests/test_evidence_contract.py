"""Research-only fixtures. They are not biomedical observations."""
from copy import deepcopy
import json
import math
from pathlib import Path
import tempfile
import unittest

from research.mortality.evidence_contract import (audit_record, attach_increment,
    compare_contexts, digest, number_matches, project_estimates, source_review_gate)
from research.mortality.build_evidence_exchange import (build, literature_records,
    empirical_records, import_pr115_performance, measurement)

ROOT=Path(__file__).resolve().parents[1]/'research/mortality'


def fixture(metric='delta_AUC',value=.01):
    return {'id':'SYNTHETIC_TEST_ONLY','record_type':'prediction_performance',
      'context':{'endpoint':'all_cause_mortality','population':'test cohort',
        'country':'USA','sex':'female','age_range':[40,79],'horizon_years':10,
        'horizon_basis':'prediction_horizon','cohort_id':'C1','validation_set':'temporal_validation'},
      'measurement':{'id':'NT-proBNP','assay':'test_method','matrix':'plasma','units':'pg/mL'},
      'model':{'id':'extended','base_id':'base','component_kind':'single_marker','added_components':['NT-proBNP']},
      'estimate':{'metric':metric,'value':value,'interval':{'kind':'95% CI','lower':-.02,'upper':.03}},
      'source':{'kind':'test_fixture','review':'source_checked','locator':'fixture only','sha256':'0'*64}}


class NumericTests(unittest.TestCase):
    def test_sign_preserved(self):
        self.assertTrue(number_matches(-.01,'-0.01'))
        self.assertFalse(number_matches(.01,'-0.01'))
    def test_unicode_sign(self):
        self.assertTrue(number_matches(-.01,'−0.01'))
    def test_no_implicit_scale(self):
        self.assertFalse(number_matches(.004,'4'))
        self.assertTrue(number_matches(.004,'4',scale='per_thousand'))
    def test_percent_explicit(self):
        self.assertFalse(number_matches(.04,'4%'))
        self.assertTrue(number_matches(.04,'4%',scale='percent'))
    def test_scientific_notation(self):
        self.assertTrue(number_matches(.004,'4e-3'))
    def test_wrong_token_nonfinite_bool(self):
        for value,token in [(True,'1'),(math.nan,'nan'),(.04,'4-5'),(.04,'0.04 to 0.05'),(.04,'text 0.04')]:
            with self.subTest(value=value,token=token):self.assertFalse(number_matches(value,token))
    def test_scale_must_be_declared(self):
        with self.assertRaises(ValueError):number_matches(.04,'4',scale='automatic')


class MetricTests(unittest.TestCase):
    def test_negative_net_benefit_retained(self):
        r=fixture('net_benefit',.1-.4*.3/.7)
        r['estimate'].update(decision_threshold=.3,decision_action='test_action',net_benefit_definition='treated_unstandardized')
        self.assertEqual(audit_record(r)['errors'],[])
    def test_net_benefit_needs_action_and_threshold(self):
        a=audit_record(fixture('net_benefit',-.1))
        self.assertIn('missing_decision_threshold',a['gaps'])
        self.assertIn('missing_decision_action',a['gaps'])
    def test_nri_valid_beyond_one(self):
        for metric in ('continuous_NRI','categorical_NRI'):
            for value in (-2.,-1.2,1.2,2.):
                with self.subTest(metric=metric,value=value):self.assertEqual(audit_record(fixture(metric,value))['errors'],[])
    def test_nri_outside_bounds(self):
        self.assertIn('estimate_outside_theoretical_domain',audit_record(fixture('continuous_NRI',2.01))['errors'])
    def test_brier_not_pure_calibration(self):
        self.assertEqual(audit_record(fixture('Brier_score',.1))['evidence_kind'],'overall_prediction_error')
    def test_brier_change_negative_retained(self):
        a=audit_record(fixture('delta_Brier',-.01));self.assertEqual(a['errors'],[])
        self.assertEqual(a['evidence_kind'],'incremental_prediction_error')
    def test_iqr_not_ci_and_need_not_contain_mean(self):
        r=fixture('C_index',.9);r['estimate']['interval']={'kind':'IQR','lower':.7,'upper':.8}
        self.assertEqual(audit_record(r)['errors'],[]);self.assertEqual(r['estimate']['interval']['kind'],'IQR')
    def test_unconstrained_interval_not_clipped(self):
        r=fixture('C_index',.99);r['estimate']['interval']={'kind':'95% CI','lower':.9,'upper':1.03}
        self.assertEqual(audit_record(r)['errors'],[])
    def test_standard_error_is_separate_scalar(self):
        r=fixture('C_index',.8);r['estimate'].update(interval=None,standard_error=.02)
        self.assertEqual(audit_record(r)['errors'],[])
        r['estimate']['standard_error']=-.02;self.assertIn('invalid_standard_error',audit_record(r)['errors'])
    def test_zero_observed_expected_allowed(self):
        self.assertEqual(audit_record(fixture('observed_expected_ratio',0))['errors'],[])
    def test_malformed_inputs_reported(self):
        self.assertEqual(audit_record([])['errors'],['invalid_record_mapping'])
        r=fixture();r['estimate']=None;self.assertIn('invalid_mapping',audit_record(r)['errors'])
        r=fixture();r['estimate']['metric']=[];self.assertIn('invalid_metric_type',audit_record(r)['errors'])
        r=fixture();r['model']=['bad'];self.assertIn('invalid_model_mapping',audit_record(r)['errors'])
    def test_missing_metadata_does_not_discard(self):
        r=fixture();r['context']['country']=None
        self.assertTrue(audit_record(r)['retain']);self.assertFalse(audit_record(r)['comparison_ready'])


class ContextTests(unittest.TestCase):
    def test_identical_context(self):
        self.assertEqual(compare_contexts(fixture(),fixture())['status'],'same_context')
    def test_unknown_not_wildcard(self):
        r=fixture();r['context']['cohort_id']=None
        self.assertEqual(compare_contexts(r,r)['status'],'unresolved')
    def test_open_age_bound_unresolved(self):
        r=fixture();r['context']['age_range']=[40,None]
        self.assertEqual(compare_contexts(r,r)['status'],'unresolved')
    def test_each_context_axis_blocks_wrong_attachment(self):
        for key,value in [('horizon_years',1),('sex','male'),('country','DEU'),('age_range',[65,79]),('cohort_id','C2'),('validation_set','derivation'),('endpoint','incident_heart_failure')]:
            with self.subTest(key=key):
                r=fixture();r['context'][key]=value
                self.assertFalse(attach_increment(fixture(),r)['attach'])
    def test_invalid_increment_cannot_attach(self):
        r=fixture();r['estimate']['value']=2
        self.assertFalse(attach_increment(fixture(),r)['attach'])
    def test_malformed_record_cannot_attach(self):
        self.assertFalse(attach_increment([],fixture())['attach'])
    def test_correct_single_marker_can_attach(self):
        self.assertTrue(attach_increment(fixture(),fixture())['attach'])
    def test_joint_effect_never_attaches_to_one_component(self):
        r=fixture();r['model'].update(component_kind='marker_block',added_components=['NT-proBNP','cystatin_C'])
        self.assertFalse(attach_increment(fixture(),r)['attach'])
    def test_assay_transport_needs_evidence(self):
        r=fixture();r['measurement']['assay']='different_platform'
        self.assertFalse(attach_increment(fixture(),r)['attach'])
    def test_pending_source_not_attached(self):
        r=fixture();r['source']['review']='pending'
        self.assertFalse(attach_increment(fixture(),r)['attach'])
    def test_sex_effect_mapping_survives_projection(self):
        male,female=fixture('HR',1.2),fixture('HR',1.6)
        male['id']='synthetic_male';male['context']['sex']='male';female['id']='synthetic_female'
        rows=project_estimates([male,female])
        self.assertEqual({r['context']['sex']:r['estimate']['value'] for r in rows},{'male':1.2,'female':1.6})
    def test_projection_no_aliasing(self):
        r=fixture();out=project_estimates([r]);out[0]['context']['sex']='male'
        self.assertEqual(r['context']['sex'],'female')
    def test_followup_not_a_prediction_horizon(self):
        r=fixture();r['context']['horizon_basis']='median_followup'
        self.assertIn('followup_is_not_prediction_horizon',audit_record(r)['gaps'])
    def test_synthetic_never_evidence(self):
        r=fixture();r['source']['kind']='synthetic_control'
        a=audit_record(r);self.assertTrue(a['retain']);self.assertIn('synthetic_not_evidence',a['errors'])
    def test_no_clinical_upgrade(self):
        r=fixture();r['clinical_use_ready']=True
        self.assertFalse(audit_record(r)['clinical_use_ready'])
        self.assertFalse(audit_record(r)['usable_as_mortality_coefficient'])
    def test_country_and_endpoint_deployment_rejected(self):
        a=source_review_gate(fixture(),target_country='RUS',target_endpoint='noninfectious_natural_mortality')
        self.assertFalse(a['clinical_use_ready']);self.assertIn('country_transport_not_validated',a['reasons'])
        self.assertIn('endpoint_mismatch',a['reasons'])


class MigrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lit=literature_records(ROOT);cls.emp=empirical_records(ROOT)
    def test_all_legacy_estimates_retained(self):
        self.assertEqual(len(self.lit),58);self.assertEqual(len({r['id'] for r in self.lit}),58)
    def test_empirical_and_literature_not_mixed(self):
        self.assertEqual(len(self.emp),92)
        self.assertTrue(all(r['source']['kind']=='empirical_evaluation' and r['source']['pmid'] is None for r in self.emp))
    def test_prevent_iqr_preserved(self):
        rows=[r for r in self.lit if r['id'].startswith('literature:PREVENT:')]
        self.assertEqual(len(rows),6)
        self.assertTrue(all(r['estimate']['interval']['kind']=='IQR' for r in rows))
    def test_negated_mortality_not_misclassified(self):
        rows=[r for r in self.lit if 'HF_BIOMARKERS2025' in r['id']]
        self.assertTrue(all(r['context']['endpoint']=='incident_heart_failure' for r in rows))
    def test_meld_quarantine_retained(self):
        rows=[r for r in self.lit if 'MELD_GERMANY2025' in r['id']]
        self.assertEqual(len(rows),8)
        self.assertTrue(all('upstream_quality_quarantine' in audit_record(r)['errors'] for r in rows))
    def test_clinical_and_functional_not_laboratory_analytes(self):
        self.assertEqual(measurement('BODE_score')['kind'],'clinical_model_or_score')
        self.assertEqual(measurement('gait_speed_m_s')['kind'],'functional_measure')
        self.assertEqual(measurement('FEV1')['kind'],'physiological_measure')
        self.assertTrue(measurement('FIB4')['contains_age'])
    def test_posthoc_and_block_attribution_survive(self):
        r=next(r for r in self.emp if r['model'].get('base_id')=='M0c_clinical_only_posthoc')
        self.assertTrue(r['model']['post_hoc']);self.assertEqual(r['model']['component_kind'],'marker_block')
        r=next(r for r in self.emp if r['model'].get('base_id')=='M1_routine' and r['model']['id']=='M4_cystatinC')
        self.assertEqual(len(r['model']['added_components']),3)
    def test_frozen_conditional_interval_not_full_ci(self):
        r=next(r for r in self.emp if r['model'].get('base_id')=='M1_routine' and r['model']['id']=='M4_cystatinC' and r['estimate']['metric']=='delta_AUC')
        i=r['estimate']['interval'];self.assertEqual(i['kind'],'conditional_validation_percentiles')
        self.assertLess(i['lower'],0);self.assertGreater(i['upper'],0)
        self.assertIn('excludes_training',i['scope'])
    def test_no_false_historical_git_provenance_for_companion_table(self):
        r=next(r for r in self.emp if r['model'].get('base_id')=='M1_routine' and r['model']['id']=='M2_NTproBNP' and r['estimate']['metric']=='delta_AUC')
        source=r['estimate']['interval']['source']
        self.assertIsNone(source['commit']);self.assertIn('original_archive_sha256',source)
    def test_export_is_lossless_and_reproducible(self):
        with tempfile.TemporaryDirectory() as d:
            summary=build(ROOT,Path(d));payload=json.loads((Path(d)/'evidence_exchange.json').read_text())
            self.assertEqual(summary['records'],150);self.assertEqual(summary['retained'],150)
            self.assertFalse(summary['new_model_fit']);self.assertFalse(summary['upstream_snapshot_imported'])
            self.assertEqual(summary['records_with_blockers'],8)
            self.assertEqual(len(payload['records']),150)
    def test_pr115_adapter_keeps_unreviewed_and_synthetic_status(self):
        doc={'not_evidence':True,'rows':[{'id':'fixture','metric':'delta_C','value':-.01,'outcome_type':'all_cause_mortality',
          'sex':'female','horizon':'10 years','model':'M2','base_model':'M1','source':{'pmid':'1'},
          'checks':{'semantic':'supported'},'added_component_kind':'single_marker','added_component':'X'}]}
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'synthetic.json';path.write_text(json.dumps(doc));rows=import_pr115_performance(path)
        self.assertEqual(rows[0]['raw_record'],doc['rows'][0]);self.assertIsNone(rows[0]['context']['horizon_years'])
        self.assertIn('synthetic_not_evidence',audit_record(rows[0])['errors'])
        self.assertIn('source_verification_pending',audit_record(rows[0])['gaps'])
    def test_no_new_risk_coefficients(self):
        self.assertTrue(all(not r['usable_as_mortality_coefficient'] and not r['clinical_use_ready'] for r in self.lit+self.emp))

if __name__=='__main__':unittest.main()
