import copy
import unittest
from unittest.mock import patch

from research.mortality.cohort27.gate import (
    ASSAYS, AVAILABILITY_POLICY, audit, require_qualification,
)


def public_contract():
    return {
        'cohort_id': 'NHANESIII_FIRST_EXAM',
        'data_class': 'public_unrestricted',
        'execution_context': 'github_actions',
        'source_permission_evidence': 'https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx',
        'source_reviewed': True,
        'baseline_visit': 'first_exam',
        'endpoint': {
            'kind': 'all_cause_death', 'baseline_origin': 'exam', 'followup_unit': 'months',
            'followup_after_baseline_verified': True, 'linkage_same_person_verified': True,
            'horizon_years': 8,
        },
        'sampling': {'weight_variable': 'WTCYPEX6', 'subsample_weights_documented': True},
        'assays': {
            key: dict(identity, unit='mg/L', matrix='serum', scale='absolute_concentration',
                      baseline_visit='first_exam',
                      source_locator='https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.pdf',
                      method='documented immunonephelometry', quantification_policy=AVAILABILITY_POLICY)
            for key, identity in ASSAYS.items()
        },
        'frozen_model': {
            'required_features_mapping_verified': False,
            'clinical_semantics_verified': False,
            'assay_transport_reviewed': False,
            'baseline_model_sha256': 'a' * 64,
        },
    }


class PublicGateTests(unittest.TestCase):
    def test_public_qualification_not_model_validation(self):
        r = audit(public_contract())
        self.assertTrue(r['qualification_import_allowed'])
        self.assertFalse(r['frozen_model_evaluation_allowed'])
        self.assertFalse(r['quantified_or_QC_usable_measurements_established'])
        self.assertIn('frozen_model:clinical_semantics_verified', r['unresolved_metadata'])

    def test_assay_uncertainty_is_explicit_without_blocking_numeric_counts(self):
        r = require_qualification(public_contract())
        self.assertTrue(any('mg/dL' in note for note in r['notes']))
        self.assertEqual(r['qualification_scope'], 'numeric_availability_and_endpoint_counts_only')

    def test_gate_is_pure_metadata_and_does_not_mutate_contract(self):
        c = public_contract(); original = copy.deepcopy(c)
        with patch('builtins.open', side_effect=AssertionError('participant file opened')):
            audit(c)
        self.assertEqual(c, original)

    def test_strings_and_integers_cannot_bypass_boolean_review(self):
        for bad in ('true', 1, [], None):
            c = public_contract(); c['source_reviewed'] = bad
            self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_missing_or_wrong_followup_cannot_pass(self):
        for key, bad in [('kind', 'cardiovascular_event'), ('baseline_origin', 'interview'),
                         ('followup_unit', 'years'), ('followup_after_baseline_verified', 'true'),
                         ('linkage_same_person_verified', 1), ('horizon_years', True)]:
            with self.subTest(key=key):
                c = public_contract(); c['endpoint'][key] = bad
                self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_missing_endpoint_metadata_blocks_before_read(self):
        c = public_contract(); del c['endpoint']
        with self.assertRaisesRegex(ValueError, 'qualification_import_blocked'):
            require_qualification(c)

    def test_repeat_exam_not_interchangeable(self):
        c = public_contract(); c['baseline_visit'] = 'second_exam'
        self.assertFalse(audit(c)['qualification_import_allowed'])
        c = public_contract(); c['assays']['CST3']['baseline_visit'] = 'second_exam'
        self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_general_exam_weights_cannot_replace_subsample_weights(self):
        c = public_contract(); c['sampling']['weight_variable'] = 'WTPFEX6'
        self.assertFalse(audit(c)['qualification_import_allowed'])
        c = public_contract(); c['sampling']['subsample_weights_documented'] = 1
        self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_wrong_molecule_or_units_or_scale_block(self):
        for field, value in [('variable', 'PTNFR1'), ('uniprot', 'P01375'), ('unit', 'mg/dL'),
                             ('matrix', 'plasma'), ('scale', 'NPX')]:
            with self.subTest(field=field):
                c = public_contract(); c['assays']['B2M'][field] = value
                self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_missing_assay_or_method_or_policy_blocks(self):
        c = public_contract(); del c['assays']['BTP']
        self.assertFalse(audit(c)['qualification_import_allowed'])
        for field in ('method', 'quantification_policy'):
            c = public_contract(); c['assays']['CST3'][field] = 'unknown'
            self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_untrusted_urls_do_not_establish_public_release(self):
        bad_urls = ['https://wwwn.cdc.gov.evil.test/a', 'https://evil.test/cdc.gov',
                    'https://user:pass@wwwn.cdc.gov/a', 'http://wwwn.cdc.gov/a',
                    'https://wwwn.cdc.gov:444/a', 'https://wwwn.cdc.gov/a#fragment',
                    'https://wwwn.cdc.gov:broken/a']
        for url in bad_urls:
            c = public_contract(); c['source_permission_evidence'] = url
            self.assertFalse(audit(c)['qualification_import_allowed'], url)

    def test_model_assertions_cannot_authorize_evaluation(self):
        c = public_contract()
        c['frozen_model'].update(required_features_mapping_verified=True,
                                 clinical_semantics_verified=True, assay_transport_reviewed=True)
        r = audit(c)
        self.assertTrue(r['qualification_import_allowed'])
        self.assertFalse(r['frozen_model_evaluation_allowed'])
        self.assertIn('frozen_model:separate_evaluation_protocol_required', r['unresolved_metadata'])

    def test_bad_metadata_shapes_fail_closed(self):
        with self.assertRaises(ValueError): audit([])
        for key in ('endpoint', 'sampling', 'assays', 'frozen_model'):
            c = public_contract(); c[key] = []
            with self.assertRaisesRegex(ValueError, 'metadata_object'):
                audit(c)


class RestrictedAndSyntheticTests(unittest.TestCase):
    def test_hrs_cannot_be_relabelled_public_in_ai_or_cloud(self):
        for context in ('ai_session', 'github_actions'):
            c = public_contract(); c['cohort_id'] = 'HRS2016'; c['execution_context'] = context
            c['access_approved'] = True
            r = audit(c)
            self.assertFalse(r['qualification_import_allowed'])
            self.assertIn('hrs_person_records_prohibited_in_ai_or_unattended_cloud', r['blockers'])

    def test_controlled_data_blocked_even_with_approval_attestation(self):
        for approved in (False, True, 'true'):
            c = public_contract(); c['data_class'] = 'controlled'; c['access_approved'] = approved
            self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_unknown_cohort_context_and_data_class_block(self):
        for key, bad in [('cohort_id', 'OTHER_COHORT'), ('execution_context', 'public_web_app'),
                         ('data_class', 'public'), ('execution_context', [])]:
            c = public_contract(); c[key] = bad
            self.assertFalse(audit(c)['qualification_import_allowed'])

    def test_synthetic_fixture_never_counts_as_empirical_access(self):
        c = public_contract(); c['data_class'] = 'synthetic_fixture'
        c['generated_without_microdata'] = True
        r = audit(c)
        self.assertTrue(r['synthetic_fixture_import_allowed'])
        self.assertFalse(r['qualification_import_allowed'])
        self.assertFalse(r['frozen_model_evaluation_allowed'])
        with self.assertRaises(ValueError): require_qualification(c)

    def test_synthetic_origin_needs_exact_boolean(self):
        for value in (None, 'true', 1, False):
            c = public_contract(); c['data_class'] = 'synthetic_fixture'
            c['generated_without_microdata'] = value
            self.assertFalse(audit(c)['synthetic_fixture_import_allowed'])


if __name__ == '__main__':
    unittest.main()
