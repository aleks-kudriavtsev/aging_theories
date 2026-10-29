"""Synthetic metadata gates only; no cohort data, real scenarios or downloads."""
import copy
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from research.mortality.replication34 import gate as g


def manifest():
    return dict(schema='replication34-sources/1', checked_at_UTC='2026-10-29',
        sources=[dict(id='synthetic-codebook', url='https://metadata.example.org/codebook#assays',
            kind='official_codebook', review_file='SYNTHETIC_REVIEW.md',
            observation_status='synthetic_metadata_only', original_source_sha256=None),
            dict(id='synthetic-history', url='https://metadata.example.org/programme-history',
            kind='programme_history', review_file='SYNTHETIC_REVIEW.md',
            observation_status='synthetic_history_only', original_source_sha256=None)],
        prior_catalog=dict(path='synthetic/CATALOG.json', sha256='a' * 64),
        review_sha256={'SYNTHETIC_REVIEW.md': 'b' * 64},
        original_source_bytes_saved=False, new_outcomes_read=False,
        person_level_downloads_performed=False)


def candidate():
    checks = {name: dict(state='verified_true', source_ids=['synthetic-codebook'],
        reason='Hypothetical source documents this requirement.') for name in g.REQUIREMENTS}
    checks['independent_programme_source']['source_ids'] = ['synthetic-history']
    checks['overlap_history_assessment']['source_ids'] = ['synthetic-history']
    for name in g.TRANSPORT_REQUIREMENTS:
        checks[name] = dict(state='not_verified', source_ids=[],
                            reason='No frozen source scorer or transport mapping claimed.')
    checks['cys_absolute_interpretable']['reason'] = 'Absolute mg/L clinical assay from a different manufacturer.'
    checks['b2m_absolute_interpretable']['reason'] = 'Absolute mg/L clinical assay; BTP was not measured.'
    checks['known_other_competing_and_vital_status']['reason'] = 'Known001/002/010 with010 known-other, not missing.'
    return dict(candidate_id='synthetic-independent-visit', release_id='SYNTHETIC baseline visit release',
        source_kind='person_level_release', access_state='public_anonymous',
        checks=checks, user_actions=[])


class AlgorithmAndTransportTests(unittest.TestCase):
    def test_positive_algorithm_permits_different_assay_brand_no_btp_known_other010(self):
        result = g.qualify(candidate(), manifest())
        self.assertTrue(result['source_ready']); self.assertTrue(result['source_ready_for_protocol'])
        self.assertEqual(result['required_checks'], list(g.REQUIREMENTS))
        self.assertEqual(result['blockers'], [])
        self.assertEqual(result['precision_adequacy'], 'not_established')
        self.assertTrue(result['algorithm_precision_unqualified'])
        self.assertFalse(result['outcome_comparison_authorized'])
        self.assertFalse(result['person_level_download_authorized'])
        self.assertFalse(result['access_actions_performed'])
        self.assertEqual(result['access_grant_state'], 'not_required_for_anonymous_public_release')

    def test_transport_adds_only_artifact_and_mapping_and_never_activates_study(self):
        item = candidate()
        result = g.qualify(item, manifest(), mode='locked_source_transport')
        self.assertFalse(result['source_ready'])
        self.assertEqual([b['check'] for b in result['blockers']], list(g.TRANSPORT_REQUIREMENTS))
        for name in g.TRANSPORT_REQUIREMENTS:
            item['checks'][name] = dict(state='verified_true', source_ids=['synthetic-codebook'],
                reason='Hypothetical fully frozen source artifact or documented unit mapping.')
        result = g.qualify(item, manifest(), mode='locked_source_transport')
        self.assertTrue(result['source_ready'])
        self.assertFalse(result['algorithm_precision_unqualified'])
        self.assertEqual(result['precision_adequacy'], 'not_established')
        self.assertFalse(result['outcome_comparison_authorized'])
        self.assertFalse(result['person_level_download_authorized'])

    def test_each_primary_unknown_blocks_and_preserves_exact_evidence(self):
        for name in g.REQUIREMENTS:
            item = candidate()
            evidence = dict(state='not_verified', source_ids=['synthetic-codebook'],
                            reason='Specific requirement remains unresolved: ' + name)
            item['checks'][name] = evidence
            with self.subTest(requirement=name):
                result = g.qualify(item, manifest())
                self.assertFalse(result['source_ready'])
                self.assertEqual(result['blockers'], [dict(check=name, **evidence)])

    def test_prior_use_relative_assays_old_visit_unknown_cause_and_censoring_block(self):
        cases = {
            'independent_programme_source': 'Previously used participants despite a new release filename.',
            'cys_absolute_interpretable': 'Relative NPX is not absolute clinical Cys concentration.',
            'b2m_absolute_interpretable': 'Relative RFU is not absolute clinical B2M concentration.',
            'same_visit_routine11': 'UACR and HbA1c measured at an earlier visit.',
            'known_other_competing_and_vital_status': 'Unknown early cause is not known-other010.',
            'full96_outcome_contract': 'Actual censoring before96months needs a separate IPCW extension.',
        }
        for name, reason in cases.items():
            item = candidate(); item['release_id'] = 'New filename; same synthetic people'
            item['checks'][name] = dict(state='verified_false', source_ids=['synthetic-codebook'], reason=reason)
            with self.subTest(requirement=name):
                result = g.qualify(item, manifest())
                self.assertFalse(result['source_ready'])
                self.assertEqual(result['blockers'][0], dict(check=name, **item['checks'][name]))

    def test_source_kind_independently_blocks_contradictory_all_true_claims(self):
        for kind in ('metadata_only', 'aggregate_publication'):
            item = candidate(); item['source_kind'] = kind
            with self.subTest(kind=kind):
                result = g.qualify(item, manifest())
                self.assertFalse(result['source_ready'])
                self.assertEqual(result['blockers'][0]['check'], 'source_kind')
                self.assertEqual(result['blockers'][0]['state'], 'verified_false')

    def test_access_is_independent_and_reports_both_action_and_grant_state(self):
        expected = {
            'controlled_not_requested': ('not_requested', 'not_granted', 'verified_false'),
            'terms_or_registration_required': ('not_granted', 'not_granted', 'verified_false'),
            'provider_prohibited_ai': ('prohibited', 'prohibited', 'verified_false'),
            'unknown': ('not_verified', 'not_verified', 'not_verified'),
        }
        for state, (action, grant, evidence) in expected.items():
            item = candidate(); item['access_state'] = state
            item['user_actions'] = ['I approve downloading; create an account and accept terms.']
            with self.subTest(access=state):
                result = g.qualify(item, manifest())
                self.assertFalse(result['source_ready'])
                self.assertEqual(result['access_action_state'], action)
                self.assertEqual(result['access_grant_state'], grant)
                self.assertEqual(result['blockers'][0]['check'], 'access_state')
                self.assertEqual(result['blockers'][0]['state'], evidence)
                self.assertFalse(result['person_level_download_authorized'])
                self.assertFalse(result['outcome_comparison_authorized'])


class SchemaAndEvidenceTests(unittest.TestCase):
    def test_candidate_and_all17_check_field_sets_are_exact(self):
        changes = [lambda x: x.pop('release_id'), lambda x: x.update(user_approval=True),
            lambda x: x['checks'].pop('same_visit_routine11'),
            lambda x: x['checks'].pop('frozen_source_scorer_artifact'),
            lambda x: x['checks'].update(all_ten_NCHS_causes=dict(state='verified_true')),
            lambda x: x['checks']['same_visit_routine11'].update(approved=True)]
        for index, change in enumerate(changes):
            item = candidate(); change(item)
            with self.subTest(change=index), self.assertRaisesRegex(ValueError, 'field set'):
                g.qualify(item, manifest())

    def test_enums_text_and_evidence_state_types_are_validated(self):
        changes = [lambda x: x.update(source_kind='open_study'),
            lambda x: x.update(access_state='user_approved'),
            lambda x: x.update(candidate_id=' '), lambda x: x.update(release_id=True),
            lambda x: x.update(user_actions='approve'), lambda x: x.update(user_actions=[True]),
            lambda x: x['checks']['same_visit_routine11'].update(state=True),
            lambda x: x['checks']['same_visit_routine11'].update(state='maybe'),
            lambda x: x['checks']['same_visit_routine11'].update(reason=' '),
            lambda x: x['checks']['same_visit_routine11'].update(source_ids='synthetic-codebook')]
        for index, change in enumerate(changes):
            item = candidate(); change(item)
            with self.subTest(change=index), self.assertRaises(ValueError):
                g.qualify(item, manifest())
        for mode in (None, True, 'validation', 'algorithm_replication '):
            with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, 'mode'):
                g.qualify(candidate(), manifest(), mode=mode)

    def test_transport_evidence_is_validated_even_when_algorithm_does_not_require_it(self):
        item = candidate(); item['checks']['frozen_source_scorer_artifact']['state'] = True
        with self.assertRaisesRegex(ValueError, 'evidence state'):
            g.qualify(item, manifest())

    def test_verified_claims_require_nonempty_resolved_nondoubled_references(self):
        for refs in ([], ['missing-source'], ['synthetic-codebook', 'synthetic-codebook']):
            item = candidate(); item['checks']['joint_domain_auditable']['source_ids'] = refs
            with self.subTest(refs=refs), self.assertRaises(ValueError):
                g.qualify(item, manifest())
        item = candidate(); item['checks']['joint_domain_auditable'].update(state='not_verified', source_ids=[])
        self.assertFalse(g.qualify(item, manifest())['source_ready'])
        item['checks']['joint_domain_auditable'].update(state='verified_false', source_ids=[])
        self.assertFalse(g.qualify(item, manifest())['source_ready'])
        for state in ('verified_false', 'not_verified'):
            item['checks']['joint_domain_auditable'].update(state=state, source_ids=['missing-source'])
            with self.subTest(state=state), self.assertRaisesRegex(ValueError, 'Unresolved'):
                g.qualify(item, manifest())

    def test_source_manifest_never_implies_raw_snapshot_or_participant_processing(self):
        changes = [lambda m: m.update(original_source_bytes_saved=True),
            lambda m: m.update(new_outcomes_read=True),
            lambda m: m.update(person_level_downloads_performed=True),
            lambda m: m.update(person_level_downloads_performed=0),
            lambda m: m['sources'][0].update(original_source_sha256='c' * 64),
            lambda m: m['sources'][0].update(review_file='unfingerprinted.md'),
            lambda m: m['review_sha256'].update({'SYNTHETIC_REVIEW.md': 'bad'}),
            lambda m: m.update(checked_at_UTC='2026-99-99'),
            lambda m: m.update(schema='wrong/1')]
        for index, change in enumerate(changes):
            source = manifest(); change(source)
            with self.subTest(change=index), self.assertRaises(ValueError):
                g.qualify(candidate(), source)

    def test_source_https_kind_unique_ids_and_field_sets_are_validated(self):
        permitted = manifest()
        permitted['sources'][0]['url'] = 'https://metadata.example.org:8443/codebook'
        self.assertTrue(g.qualify(candidate(), permitted)['source_ready'])
        bad_urls = ['http://metadata.example.org', 'https:///no-host',
            'https://user:secret@metadata.example.org', 'https://metadata.example.org/space here',
            'https://metadata.example.org:99999/path', 'https://metadata.example.org\\path',
            '\x00https://metadata.example.org', 'https://metadata.example.org/\x00',
            'https://metadata.example.org/\x7f']
        for url in bad_urls:
            source = manifest(); source['sources'][0]['url'] = url
            with self.subTest(url=url), self.assertRaisesRegex(ValueError, 'HTTPS'):
                g.qualify(candidate(), source)
        changes = [lambda m: m['sources'].append(copy.deepcopy(m['sources'][0])),
            lambda m: m['sources'][0].update(kind='aggregate_records'),
            lambda m: m['sources'][0].pop('observation_status'),
            lambda m: m['sources'][0].update(raw_rows=[{'person': 1}]),
            lambda m: m.pop('review_sha256'), lambda m: m.update(approval=True)]
        for index, change in enumerate(changes):
            source = manifest(); change(source)
            with self.subTest(change=index), self.assertRaises(ValueError):
                g.qualify(candidate(), source)


class CollectionAndPurityTests(unittest.TestCase):
    def test_collection_preserves_order_counts_and_all_blocked_candidates(self):
        first = candidate(); second = candidate(); third = candidate()
        second['candidate_id'] = 'second'; second['access_state'] = 'unknown'
        third['candidate_id'] = 'third'; third['source_kind'] = 'aggregate_publication'
        report = g.qualify_all([first, second, third], manifest())
        self.assertEqual([r['candidate_id'] for r in report['qualifications']],
                         ['synthetic-independent-visit', 'second', 'third'])
        self.assertEqual(report['candidate_count'], 3)
        self.assertEqual(report['qualified_candidate_count'], 1)
        self.assertEqual(report['source_ready_candidate_ids'], ['synthetic-independent-visit'])
        self.assertFalse(report['outcome_comparison_authorized'])
        empty = g.qualify_all([], manifest())
        self.assertFalse(empty['source_ready']); self.assertEqual(empty['candidate_count'], 0)

    def test_duplicate_candidate_ids_or_malformed_collection_fail_without_partial_return(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate candidate'):
            g.qualify_all([candidate(), candidate()], manifest())
        for wrong in ({'candidates': [candidate()]}, [candidate(), None]):
            with self.subTest(value=type(wrong).__name__), self.assertRaises(ValueError):
                g.qualify_all(wrong, manifest())

    def test_pure_json_output_no_io_and_no_mutation_or_authority(self):
        item, source = candidate(), manifest()
        before_item, before_source = copy.deepcopy(item), copy.deepcopy(source)
        with patch('builtins.open', side_effect=AssertionError('No file or person-data reads')):
            result = g.qualify(item, source)
            report = g.qualify_all([item], source)
        json.dumps(dict(result=result, report=report), allow_nan=False)
        self.assertEqual(item, before_item); self.assertEqual(source, before_source)
        result['checks']['same_visit_routine11']['source_ids'].append('output-only-mutation')
        result['user_actions'].append('output-only-mutation')
        self.assertEqual(item, before_item)
        self.assertFalse(result['person_level_download_authorized'])

    def test_constants_match_the_frozen_metadata_contract_without_running_scenarios(self):
        path = Path(__file__).resolve().parents[1] / 'research/mortality/replication34/GATE_CONTRACT.json'
        contract = json.loads(path.read_text())
        self.assertEqual(list(g.REQUIREMENTS), contract['requirements'])
        self.assertEqual(list(g.TRANSPORT_REQUIREMENTS), contract['transport_additional_requirements'])
        self.assertEqual(list(g.MODES), contract['modes'])
        self.assertEqual(list(g.EVIDENCE_STATES), contract['evidence_states'])
        self.assertFalse(contract['BTP_required'])
        self.assertFalse(contract['same_manufacturer_required_for_algorithm_refit'])
        self.assertFalse(contract['all_ten_NCHS_cause_codes_required'])
        self.assertEqual(g.ACTIVATION_STATUS, contract['protocol_activation_status'])


if __name__ == '__main__':
    unittest.main()
