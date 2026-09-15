"""Only invented presence records; no HRS microdata or trained models."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from research.mortality.hrs24.access import check_access, strict_json
from research.mortality.hrs24.feasibility import (
    ALL_FIELDS, CLINICAL, CORE13, PROTEINS, summarize_records, review_partition, run)


def synthetic_config():
    mappings = {k: {'reviewed': True, 'source_product': 'SYNTHETIC',
                   'source_variable': 'invented_' + k, 'source_locator': 'test-fixture',
                   'uniprot': PROTEINS.get(k), 'measured_form': 'NT-proBNP',
                   'assay_scale_and_QC_reviewed': True} for k in ALL_FIELDS}
    return {'data_class': 'synthetic_fixture', 'generated_without_hrs_microdata': True,
            'execution_context': 'ai_session', 'baseline_token': 'INVENTED_VISIT',
            'common_core_fields': ['albumin', 'creatinine'], 'mappings': mappings}


def authorized_declaration():
    # Invented declarations exercise rules; NOT an actual approval.
    d = synthetic_config()
    d.update(data_class='hrs_microdata', execution_context='authorized_human_local',
             approval_reference='TEST-NOT-AN-APPROVAL',
             user_dua_approved=True, institutional_signature_confirmed=True,
             secure_environment_reviewed=True, ai_file_access_disabled=True,
             third_party_transfer_disabled=True, research_only_use=True,
             no_consulting_or_licensing_obligation=True, source_mapping_reviewed=True,
             baseline_assignment_reviewed=True, requested_products=['TEST'],
             approved_products=['TEST'], commercial_affiliation=False,
             uses_ndi_or_cms=False, operator_country='RU', execution_country='RU')
    return d


def rows(n=12, omit=()):
    return [{'HHID': f'{900000+i:06d}', 'PN': '010', 'baseline_token': 'INVENTED_VISIT',
             'component': k, 'usable': '1'} for i in range(n)
            for k in sorted(ALL_FIELDS - set(omit))]


class AccessTests(unittest.TestCase):
    def test_public_metadata_allowed(self):
        self.assertEqual(check_access({'data_class':'public_documentation','contains_person_records':False})['status'], 'declared_conditions_met')
    def test_pseudonymized_still_denied_in_ai(self):
        d=authorized_declaration();d['execution_context']='ai_session'
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_github_blocked(self):
        d=authorized_declaration();d['execution_context']='github_actions'
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_public_app_blocked(self):
        d=authorized_declaration();d['execution_context']='public_web_app'
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_declared_not_legally_verified(self):
        r=check_access(authorized_declaration());self.assertEqual(r['status'],'declared_conditions_met')
        self.assertFalse(r['legal_authorization_verified']);self.assertFalse(r['environment_technically_verified'])
    def test_each_required_approval(self):
        for k in ('user_dua_approved','institutional_signature_confirmed','ai_file_access_disabled','baseline_assignment_reviewed'):
            d=authorized_declaration();d[k]=False
            with self.subTest(k=k):self.assertEqual(check_access(d)['status'],'blocked')
    def test_new_product_not_granted_by_old_dua(self):
        d=authorized_declaration();d['requested_products'].append('UNAPPROVED')
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_commercial_affiliation_disclosure(self):
        d=authorized_declaration();d['commercial_affiliation']=True
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_research_not_commercial(self):
        d=authorized_declaration();d['research_only_use']=False
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_remote_ndi_outside_US(self):
        d=authorized_declaration();d.update(uses_ndi_or_cms=True,execution_country='US',operator_country='RU')
        self.assertIn('ndi_cms_cannot_be_accessed_outside_US',check_access(d)['blockers'])
    def test_ai_generated_synthetic_not_real_data(self):
        self.assertEqual(check_access(synthetic_config())['status'],'declared_conditions_met')
    def test_derived_records_not_synthetic(self):
        d=synthetic_config();d['generated_without_hrs_microdata']=False
        self.assertEqual(check_access(d)['status'],'blocked')
    def test_aggregate_needs_review(self):
        self.assertEqual(check_access({'data_class':'reviewed_aggregate'})['status'],'blocked')
    def test_unknown_class(self):
        self.assertEqual(check_access({'data_class':'deidentified'})['status'],'blocked')
    def test_duplicate_json(self):
        with self.assertRaises(ValueError):strict_json('{"data_class":"a","data_class":"b"}')
    def test_nonfinite_json(self):
        with self.assertRaises(ValueError):strict_json('{"n":NaN}')


class FeasibilityTests(unittest.TestCase):
    def test_joint_same_people(self):
        r=summarize_records(rows(),synthetic_config())
        self.assertEqual(r['counts']['common_core_plus7_including_direct_TNF'],12)
        self.assertFalse(r['outcomes_read']);self.assertFalse(r['coefficients_fitted'])
    def test_missing_TNF_does_not_block_other6_or_become_normal(self):
        r=summarize_records(rows(omit=('TNF',)),synthetic_config())
        self.assertEqual(r['counts']['common_core_plus6_without_direct_TNF'],12)
        self.assertEqual(r['counts']['common_core_plus7_including_direct_TNF'],0)
    def test_belowlimit_not_usable(self):
        r=rows();r[0]['usable']='0';o=summarize_records(r,synthetic_config())
        self.assertEqual(o['counts']['common_core_plus7_including_direct_TNF'],11)
    def test_receptor_not_cytokine(self):
        c=synthetic_config();c['mappings']['TNF']['source_variable']='PTNFR1'
        with self.assertRaisesRegex(ValueError,'receptor'):summarize_records(rows(),c)
    def test_wrong_protein_identity(self):
        c=synthetic_config();c['mappings']['TNF']['uniprot']='P19438'
        with self.assertRaisesRegex(ValueError,'identity'):summarize_records(rows(),c)
    def test_NPPB_not_NTproBNP(self):
        c=synthetic_config();c['mappings']['NT_PROBNP']['measured_form']='NPPB'
        with self.assertRaises(ValueError):summarize_records(rows(),c)
    def test_other_wave_blocked(self):
        r=rows();r[0]['baseline_token']='INVENTED_LATER_VISIT'
        with self.assertRaisesRegex(ValueError,'baseline'):summarize_records(r,synthetic_config())
    def test_duplicate_not_averaged(self):
        r=rows();r.append(dict(r[0]))
        with self.assertRaisesRegex(ValueError,'duplicate'):summarize_records(r,synthetic_config())
    def test_leading_zeros_required(self):
        r=rows();r[0]['PN']='10'
        with self.assertRaisesRegex(ValueError,'PN'):summarize_records(r,synthetic_config())
    def test_unreviewed_mapping_blocked(self):
        c=synthetic_config();c['mappings']['GFAP']['reviewed']=False
        with self.assertRaises(ValueError):summarize_records(rows(),c)
    def test_reduced_core_explicit(self):
        c=synthetic_config();c['common_core_fields']=[]
        with self.assertRaises(ValueError):summarize_records(rows(),c)
    def test_results_not_auto_released(self):
        r=summarize_records(rows(),synthetic_config());self.assertFalse(r['approved_for_sharing'])
        self.assertFalse(r['existing_model_applicable']);self.assertFalse(r['participant_records_exported'])
        self.assertNotIn('900000',json.dumps(r))
    def test_small_cell_withholds_entire_partition(self):
        r=review_partition({'a':1,'b':20},reviewed=True)
        self.assertIsNone(r['counts']);self.assertIsNone(r['total'])
    def test_no_automatic_release_large_table(self):
        self.assertEqual(review_partition({'a':10,'b':20})['status'],'withheld')
    def test_reviewed_disjoint_table(self):
        self.assertEqual(review_partition({'a':10,'b':20},reviewed=True)['total'],30)
    def test_gate_before_source_open(self):
        with tempfile.TemporaryDirectory() as d:
            c=Path(d)/'cfg.json';c.write_text(json.dumps({'data_class':'hrs_microdata','execution_context':'ai_session'}))
            r=run(c,Path(d)/'DOES_NOT_EXIST.csv',Path(d)/'out.json')
            self.assertEqual(r['status'],'blocked');self.assertFalse((Path(d)/'out.json').exists())
    def test_wrong_schema_values_not_read(self):
        r=rows();r[0]['value']=12
        with self.assertRaisesRegex(ValueError,'schema'):summarize_records(r,synthetic_config())
    def test_input_not_modified(self):
        r=rows();c=synthetic_config();old=(deepcopy(r),deepcopy(c))
        summarize_records(r,c);self.assertEqual((r,c),old)

if __name__=='__main__':unittest.main()
