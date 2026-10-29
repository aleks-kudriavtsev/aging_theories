"""Synthetic metadata-only checks; no participant data, file IO or network."""
from copy import deepcopy
import unittest
from unittest.mock import patch
from research.mortality.freeze30 import qualification as q


def declared_complete_metadata():
    """Fictional declarations; never proof of an actual public cohort or rights."""
    c=q.blank_contract()
    c.update(cohort_id='FICTIONAL_EXTERNAL',cohort_name='Fictional independent cohort',release_id='release1',
        data_class='public_unrestricted',execution_context='ai_session')
    c['sources']=[dict(official_url='https://data.example.org/releases/release1/persons.csv',
        publisher_domain='example.org',release_id='release1',kind='person_level',direct_download=True,
        requires_account=False,requires_agreement=False,
        public_release_reference='https://data.example.org/releases/release1/description')]
    c['overlap']=dict(assessment='independent_participant_universe_documented',
        checked_against=list(q.REUSED_PROGRAMME_COHORTS),
        disjoint_evidence_urls=['https://data.example.org/protocol/participant-universe'])
    c['endpoint']=dict(kind='all_cause_death',baseline_origin='baseline_exam',target_age_min_years=40,
        target_age_max_years=75,horizon_years=8,baseline_visit='baseline_exam',
        linkage_same_person_verified=True,definition_reference='https://data.example.org/endpoint',
        followup_unit='months',conversion_to_years='divide_by_12',minimum_survivor_followup_years=10,
        early_censored_n=0,followup_support_verified=True)
    for name,row in c['raw_measurements'].items():
        row.update(available=True,source_variable='RAW_'+name,baseline_visit='baseline_exam',
            method='documented historical method',method_version='version1',
            definition_reference='https://data.example.org/dictionary/'+name,
            missing_and_substitution_policy='released missing codes excluded; no new imputation')
    scr=c['raw_measurements']['serum_creatinine']
    scr.update(calibration_basis='NHANESIII_original_CEP',transport=dict(kind='documented_same_historical_method',
        historical_method_id='NHANESIII_original_CEP_before_0.960x_minus0.184',
        evidence_url='https://data.example.org/assays/creatinine-comparison'))
    for name,row in c['routine_core'].items():
        row.update(available=True,baseline_visit='baseline_exam',derivation_reviewed=True,
            definition_reference='https://data.example.org/definitions/'+name,
            source_variables=[c['raw_measurements'][raw]['source_variable'] for raw in row['raw_inputs']])
    for name,row in c['proteins'].items():
        row.update(available=True,source_variable='PROTEIN_'+name,baseline_visit='baseline_exam',
            method='historical N Latex',method_version='historical documented version',
            quantification_policy='reported values; unknown individual QC/LOQ explicitly limited',
            source_reference='https://data.example.org/assays/'+name,
            transport=dict(kind='documented_same_historical_method',historical_method_id=q.HISTORICAL_ASSAYS[name],
                evidence_url='https://data.example.org/assays/'+name+'/comparison'))
    c['sampling']=dict(documented=True,design_reference='https://data.example.org/design',
        same_joint_domain_supported=True,design='complex_probability',weight_variable='WJOINT',
        stratum_variable='STRATUM',PSU_variable='PSU')
    c['protocol']=dict(sha256='a'*64,target_bundle_sha256='b'*64,
        immutable_lock_reference='https://data.example.org/protocol/commit123',
        frozen_at_utc='2026-10-01T00:00:00Z',fixed_before_outcomes=True,outcomes_accessed=False,
        first_outcome_access_at_utc=None)
    c['mapping_id']='operator_claimed_mapping'
    return c


class Qualification30Tests(unittest.TestCase):
    def audit(self,c=None):return q.audit(declared_complete_metadata() if c is None else c)
    def test_complete_declarations_cannot_unlock_empty_mapping_registry(self):
        r=self.audit()
        self.assertEqual(q.QUALIFIED_EXTERNAL_MAPPINGS,{})
        for category in ('access','overlap','endpoint','routine','assay','sampling','protocol'):
            self.assertEqual(r[category+'_blockers'],[],category)
        self.assertTrue(r['metadata_complete_except_mapping'])
        self.assertFalse(r['frozen_scoring_allowed'])
        self.assertEqual(r['mapping_blockers'],['no_reviewed_external_mapping_code_is_registered'])
        self.assertFalse(r['legal_authorization_verified'])
        self.assertFalse(r['independence_verified'])
        self.assertFalse(r['clinical_use_ready'])
        self.assertFalse(r['person_records_processed'])
    def test_gate_has_no_file_or_network_operations_and_does_not_mutate(self):
        c=declared_complete_metadata();before=deepcopy(c)
        with patch('builtins.open',side_effect=AssertionError('file opened')):
            r=q.audit(c)
        self.assertEqual(c,before);self.assertFalse(r['network_or_file_access_performed'])
    def test_blank_contract_is_incomplete_and_never_scores(self):
        r=q.audit(q.blank_contract());self.assertFalse(r['frozen_scoring_allowed'])
        self.assertFalse(r['metadata_complete_except_mapping'])
        self.assertTrue(r['routine_blockers']);self.assertTrue(r['assay_blockers'])
    def test_controlled_account_or_agreement_route_is_blocked(self):
        c=declared_complete_metadata();c['data_class']='controlled'
        self.assertIn('only_direct_public_unrestricted_data_supported',q.audit(c)['access_blockers'])
        for flag in ('requires_account','requires_agreement'):
            for value in (True,None,'false',0):
                with self.subTest(flag=flag,value=value):
                    c=declared_complete_metadata();c['sources'][0][flag]=value
                    self.assertTrue(q.audit(c)['access_blockers'])
    def test_unknown_access_does_not_become_public_by_url_presence(self):
        c=declared_complete_metadata();c['data_class']='unknown'
        self.assertTrue(q.audit(c)['access_blockers'])
    def test_hrs_cannot_be_relabelled_public_in_ai_or_actions(self):
        for context in ('ai_session','github_actions'):
            c=declared_complete_metadata();c.update(cohort_id='HRS_2016',execution_context=context)
            self.assertIn('HRS_person_records_prohibited_in_AI_or_unattended_cloud',q.audit(c)['access_blockers'])
        c=declared_complete_metadata();c['sources'][0].update(official_url='https://hrsdata.isr.umich.edu/data.csv',publisher_domain='umich.edu')
        self.assertIn('HRS_person_records_prohibited_in_AI_or_unattended_cloud',q.audit(c)['access_blockers'])
    def test_official_release_host_and_direct_file_declarations_required(self):
        for edit in ({'official_url':'https://evil.invalid/release/persons.csv'},
                     {'direct_download':False},{'kind':'codebook'},{'release_id':'other'},
                     {'official_url':'https://user:password@data.example.org/persons.csv'},
                     {'public_release_reference':None}):
            c=declared_complete_metadata();c['sources'][0].update(edit)
            self.assertTrue(q.audit(c)['access_blockers'])
    def test_all_nhanesIII_aliases_and_repeat_visits_rejected(self):
        for identity in ('NHANESIII_FIRST_EXAM','NHANESIII_SECOND_EXAM','NHANES_III','NHANESⅢ','NHANES3'):
            c=declared_complete_metadata();c['cohort_id']=identity
            self.assertIn('NHANESIII_first_second_or_renamed_release_already_used',q.audit(c)['overlap_blockers'])
    def test_used_continuous_nhanes_is_not_a_new_external_cohort(self):
        for identity in ('NHANES_1999_2000','NHANES_2017_2018'):
            c=declared_complete_metadata();c['cohort_id']=identity
            self.assertIn('NHANES1999_2018_programme_cycle_already_used',q.audit(c)['overlap_blockers'])
    def test_self_asserted_independence_needs_reference_and_programme_scope(self):
        c=declared_complete_metadata();c['overlap']['disjoint_evidence_urls']=[]
        self.assertTrue(q.audit(c)['overlap_blockers'])
        c=declared_complete_metadata();c['overlap']['checked_against']=['NHANESIII']
        self.assertTrue(q.audit(c)['overlap_blockers'])
    def test_canonical_time_origin_population_and_allcause_endpoint_required(self):
        for edit in ({'baseline_origin':'stored_serum_assay_year'},{'baseline_visit':'second_exam'},
                     {'kind':'CVD_composite'},{'horizon_years':3},{'target_age_min_years':60},
                     {'followup_unit':'days'},{'conversion_to_years':'guess'},{'linkage_same_person_verified':'true'}):
            c=declared_complete_metadata();c['endpoint'].update(edit)
            self.assertTrue(q.audit(c)['endpoint_blockers'])
    def test_short_or_early_censored_followup_needs_separate_implemented_design(self):
        for edit in ({'minimum_survivor_followup_years':7.99},{'early_censored_n':1},
                     {'early_censored_n':None},{'followup_support_verified':False}):
            c=declared_complete_metadata();c['endpoint'].update(edit)
            self.assertTrue(q.audit(c)['endpoint_blockers'])
    def test_missing_core_term_or_incompatible_semantics_cannot_default_normal(self):
        c=declared_complete_metadata();del c['routine_core']['log_uacr']
        self.assertIn('log_uacr:missing_engineered_term',q.audit(c)['routine_blockers'])
        c=declared_complete_metadata();c['routine_core']['current']['semantic_id']='any_tobacco'
        self.assertTrue(q.audit(c)['routine_blockers'])
    def test_raw_source_visit_unit_and_derived_mapping_are_checked(self):
        for edit in ({'available':False},{'unit':'mmol/L'},{'baseline_visit':'repeat_exam'},
                     {'missing_and_substitution_policy':None}):
            c=declared_complete_metadata();c['raw_measurements']['serum_creatinine'].update(edit)
            self.assertTrue(q.audit(c)['routine_blockers'])
        c=declared_complete_metadata();c['routine_core']['sbp']['source_variables']=['unrelated_column']
        self.assertIn('sbp:engineered_to_raw_source_variables_disagree',q.audit(c)['routine_blockers'])
    def test_modern_IDMS_creatinine_is_not_historical_CEP_by_boolean_claim(self):
        c=declared_complete_metadata();scr=c['raw_measurements']['serum_creatinine']
        scr.update(calibration_basis='modern_IDMS',transport_reviewed=True)
        self.assertTrue(any('modern_IDMS' in x for x in q.audit(c)['assay_blockers']))
    def test_missing_BTP_wrong_identity_matrix_or_RFU_is_blocked(self):
        c=declared_complete_metadata();del c['proteins']['BTP']
        self.assertIn('BTP:protein_missing',q.audit(c)['assay_blockers'])
        for edit in ({'uniprot':'TNFR1'},{'scale':'RFU'},{'unit':'NPX'},{'matrix':'plasma'},
                     {'baseline_visit':'second_exam'},{'quantification_policy':None}):
            c=declared_complete_metadata();c['proteins']['CST3'].update(edit)
            self.assertTrue(q.audit(c)['assay_blockers'])
    def test_assay_transport_cannot_be_approved_by_a_boolean(self):
        c=declared_complete_metadata();c['proteins']['CST3']['transport']={'reviewed':True}
        self.assertIn('CST3:historical_assay_transport_unresolved',q.audit(c)['assay_blockers'])
        c=declared_complete_metadata();c['proteins']['CST3']['transport'].update(kind='paired_bridge_to_historical_method',bridge_id='claimed')
        self.assertIn('CST3:paired_bridge_not_in_reviewed_external_mapping',q.audit(c)['assay_blockers'])
    def test_external_sampling_does_not_inherit_nhanes_subsample_weight(self):
        c=declared_complete_metadata();c['sampling']['weight_variable']='WTCYPEX6'
        self.assertIn('NHANES_specific_subsample_weight_cannot_be_external_default',q.audit(c)['sampling_blockers'])
        c=declared_complete_metadata();c['sampling']={'design':'unweighted_cohort','documented':True,
            'design_reference':'https://data.example.org/design','same_joint_domain_supported':True,
            'analysis_weight_rule':'unit_weights_received_cohort_estimand_only'}
        self.assertEqual(q.audit(c)['sampling_blockers'],[])
    def test_protocol_sha_lock_and_temporal_order_required(self):
        for edit in ({'sha256':'bad'},{'target_bundle_sha256':None},{'fixed_before_outcomes':'true'},
                     {'frozen_at_utc':'2026-10-01'},{'immutable_lock_reference':None},
                     {'outcomes_accessed':None},{'outcomes_accessed':False,'first_outcome_access_at_utc':'2026-10-02T00:00:00Z'},
                     {'outcomes_accessed':True,'first_outcome_access_at_utc':'2026-09-30T00:00:00Z'}):
            c=declared_complete_metadata();c['protocol'].update(edit)
            self.assertTrue(q.audit(c)['protocol_blockers'])
        c=declared_complete_metadata();c['protocol'].update(outcomes_accessed=True,first_outcome_access_at_utc='2026-10-02T00:00:00Z')
        self.assertEqual(q.audit(c)['protocol_blockers'],[])
    def test_untrusted_context_mapping_and_unit_shapes_fail_closed_without_typeerror(self):
        for value in ([],{},['ai_session'],{'claimed':'authorized'}):
            c=declared_complete_metadata();c['execution_context']=value
            self.assertTrue(q.audit(c)['access_blockers'])
            c=declared_complete_metadata();c['mapping_id']=value
            self.assertTrue(q.audit(c)['mapping_blockers'])
            c=declared_complete_metadata();c['endpoint']['followup_unit']=value
            self.assertTrue(q.audit(c)['endpoint_blockers'])
        c=declared_complete_metadata();c['overlap']['checked_against']=[{}]
        self.assertTrue(q.audit(c)['overlap_blockers'])
        c=declared_complete_metadata();c['routine_core']['sbp']['raw_inputs']=[{}]
        self.assertTrue(q.audit(c)['routine_blockers'])
    def test_nonfinite_endpoint_support_never_passes(self):
        for value in (float('nan'),float('inf'),float('-inf'),None,'10',True):
            c=declared_complete_metadata();c['endpoint']['minimum_survivor_followup_years']=value
            self.assertIn('eight_year_survivor_followup_not_established',q.audit(c)['endpoint_blockers'])

    def test_metadata_shape_and_accidental_participant_payload_rejected(self):
        with self.assertRaises(ValueError):q.audit([])
        for field in ('endpoint','routine_core','proteins','sampling'):
            c=declared_complete_metadata();c[field]=[]
            with self.assertRaises(ValueError):q.audit(c)
        c=declared_complete_metadata();c['participant_rows']=[{'id':'forbidden'}]
        with self.assertRaises(ValueError):q.audit(c)

if __name__=='__main__':unittest.main()
