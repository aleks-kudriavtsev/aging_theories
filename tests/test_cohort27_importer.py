"""Synthetic parser/linkage/qualification tests; no participant data or network."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from research.mortality.cohort27 import importer, qualify, retrieve


def adult_row(**updates):
    values = dict(SEQN=1, HSSEX=2, HSAGEIR=45, HSAGEU=2, SDPPSU6=1, SDPSTRA6=12)
    values.update(updates)
    row = [' '] * 3346
    for key, (start, end) in importer.FIELDS.items():
        value = str(values[key]).rjust(end-start+1)
        if len(value) != end-start+1:
            raise AssertionError('Fixture field does not fit')
        row[start-1:end] = value
    return ''.join(row)


def sas_dictionary(**replacements):
    lines = ['INPUT']
    for key, (start, end) in importer.FIELDS.items():
        span = str(start) if start == end else f'{start}-{end}'
        lines.append(f'  {key} {replacements.get(key, span)}')
    return '\n'.join(lines) + '\n;\n'


def release_links(**updates):
    links = dict(adult='https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.dat',
                 cys='https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.xpt',
                 btp='https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.xpt')
    links.update(updates)
    return '\n'.join(f'<a href="{url}">data</a>' for url in links.values())


class AdultParserTests(unittest.TestCase):
    def read(self, text, dictionary=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'adult.dat').write_text(text, encoding='ascii')
            (root/'adult.sas').write_text(dictionary or sas_dictionary(), encoding='ascii')
            return importer.read_adult_checked(root/'adult.dat', root/'adult.sas')

    def test_two_independent_readers_agree_on_declared_positions(self):
        d, audit = self.read(adult_row()+'\n'+adult_row(SEQN=9, HSSEX=1, HSAGEIR=73, SDPPSU6=2)+'\n')
        np.testing.assert_array_equal(d.SEQN, [1., 9.])
        np.testing.assert_array_equal(d.HSAGEIR, [45., 73.])
        self.assertEqual(audit, dict(rows=2, fields=6, cell_comparisons=12, exact_match=True))

    def test_declared_dictionary_must_match_identifier_span(self):
        with self.assertRaisesRegex(ValueError, 'dictionary position'):
            self.read(adult_row(), sas_dictionary(SEQN='2-6'))

    def test_single_column_cannot_accept_a_wider_dictionary_span(self):
        with self.assertRaisesRegex(ValueError, 'dictionary position'):
            self.read(adult_row(), sas_dictionary(HSSEX='15-16'))

    def test_wrong_record_width_fails_before_interpretation(self):
        for row in (adult_row()[:-1], adult_row()+' '):
            with self.subTest(width=len(row)):
                with self.assertRaisesRegex(ValueError, 'record width'):
                    self.read(row)

    def test_non_numeric_values_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Non-numeric adult field'):
            self.read(adult_row(HSSEX='X'))

    def test_missing_design_is_preserved_as_missing_for_later_gate(self):
        d, _ = self.read(adult_row(SDPPSU6='.'))
        self.assertTrue(pd.isna(d.SDPPSU6.iloc[0]))

    def test_duplicate_and_missing_identifiers_rejected(self):
        for text in (adult_row()+'\n'+adult_row(), adult_row(SEQN='.')):
            with self.assertRaisesRegex(ValueError, 'identifiers'):
                self.read(text)

    def test_join_never_multiplies_duplicate_or_missing_ids(self):
        left = pd.DataFrame({'SEQN': [1., 2.], 'age': [45, 50]})
        for ids in ([1., 1.], [1., np.nan]):
            right = pd.DataFrame({'SEQN': ids, 'marker': [1., 2.]})
            with self.assertRaisesRegex(ValueError, 'identifiers'):
                importer.join_one_to_one(left, right)
        joined = importer.join_one_to_one(left, pd.DataFrame({'SEQN': [1.], 'marker': [5.]}))
        self.assertEqual(len(joined), 2)
        self.assertTrue(pd.isna(joined.marker.iloc[1]))


class ReleaseResolverTests(unittest.TestCase):
    def test_relative_links_and_filename_case_resolve_without_guessing(self):
        html = release_links(adult='/nchs/data/nhanes3/1a/ADULT.DAT',
                             cys='/NCHS/DATA/NHANES3/27A/SSCYSTAT.XPT')
        links = retrieve.resolve_links(html)
        self.assertEqual(set(links), set(retrieve.DYNAMIC))
        self.assertEqual(links['cystatin.xpt'], 'https://wwwn.cdc.gov/NCHS/DATA/NHANES3/27A/SSCYSTAT.XPT')

    def test_ftp_links_use_https_to_same_official_host(self):
        links = retrieve.resolve_links(release_links(adult='ftp://ftp.cdc.gov/pub/nhanes3/1a/adult.dat'))
        self.assertEqual(links['adult.dat'], 'https://ftp.cdc.gov/pub/nhanes3/1a/adult.dat')

    def test_ftp_scheme_case_is_normalized(self):
        links = retrieve.resolve_links(release_links(adult='FTP://ftp.cdc.gov/pub/nhanes3/1a/adult.dat'))
        self.assertEqual(links['adult.dat'], 'https://ftp.cdc.gov/pub/nhanes3/1a/adult.dat')

    def test_repeat_exam_37a_is_never_used_for_first_exam(self):
        for url in ('https://wwwn.cdc.gov/nchs/data/nhanes3/37a/SSCYSTAT.xpt',
                    'https://wwwn.cdc.gov/nchs/data/nhanes3/37a/SSNH3CYS.xpt'):
            with self.assertRaisesRegex(ValueError, 'cystatin.xpt'):
                retrieve.resolve_links(release_links(cys=url))

    def test_distinct_official_candidates_are_ambiguous(self):
        html = release_links()+'\n<a href="https://ftp.cdc.gov/pub/nhanes3/27a/sscystat.xpt">alternate</a>'
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            retrieve.resolve_links(html)

    def test_duplicate_identical_hrefs_do_not_create_false_ambiguity(self):
        html = release_links()+'\n<a href="https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.xpt">duplicate</a>'
        self.assertEqual(len(retrieve.resolve_links(html)), 3)

    def test_lookalike_and_credential_hosts_cannot_supply_release(self):
        bad = ['https://wwwn.cdc.gov.evil.test/nhanes3/27a/SSCYSTAT.xpt',
               'https://user:pass@wwwn.cdc.gov/nhanes3/27a/SSCYSTAT.xpt']
        for url in bad:
            with self.assertRaises(ValueError): retrieve.resolve_links(release_links(cys=url))
        self.assertFalse(retrieve.allowed_url('http://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.xpt'))


class QualificationSummaryTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame(dict(WTCYPEX6=[1., 2., 3., 4.], dead=[1, 1, 0, 0],
                                 t_exam=[96., 97., 95., 96.], SDPSTRA6=[1, 1, 2, 2], SDPPSU6=[1, 2, 1, 2]))

    def test_eight_year_boundary_and_early_censor_are_distinct(self):
        r = qualify.summarize_domain(self.frame())
        self.assertEqual(r['observed_deaths8'], 1)
        self.assertEqual(r['early_censored8'], 1)
        self.assertEqual(r['known_status8'], 3)
        self.assertEqual(r['minimum_survivor_followup_years'], 95/12)
        self.assertEqual(r['effective_weighted_n'], 100/30)
        self.assertEqual(r['strata_n'], 2)
        self.assertEqual(r['psu_n'], 4)

    def test_summary_exports_only_aggregates(self):
        d = self.frame(); d['SEQN'] = [100, 101, 102, 103]
        r = qualify.summarize_domain(d)
        self.assertNotIn('SEQN', r)
        self.assertFalse(any(isinstance(value, (list, pd.DataFrame, pd.Series)) for value in r.values()))

    def test_empty_common_domain_is_explicit_not_a_fake_performance_result(self):
        r = qualify.summarize_domain(self.frame().iloc[:0])
        self.assertEqual(r['n'], 0)
        self.assertEqual(r['observed_deaths8'], 0)
        self.assertIsNone(r['minimum_survivor_followup_years'])
        self.assertIsNone(r['effective_weighted_n'])
        self.assertNotIn('auc', r)
        self.assertNotIn('brier', r)


class GateOrderingTests(unittest.TestCase):
    def test_denied_contract_stops_before_manifest_or_participant_open(self):
        with patch.object(importer, 'require_qualification', side_effect=ValueError('contract_denied')), \
             patch.object(importer.base, 'sha', side_effect=AssertionError('source opened')), \
             patch.object(importer, 'read_adult_checked', side_effect=AssertionError('participant parsed')), \
             patch.object(importer.base, 'read_xpt_checked', side_effect=AssertionError('participant parsed')):
            with self.assertRaisesRegex(ValueError, 'contract_denied'):
                importer.load('/nonexistent/participant/source')

    def test_retrieval_gate_precedes_directory_creation_and_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'must_not_exist'
            with patch.object(retrieve, 'require_qualification', side_effect=ValueError('contract_denied')), \
                 patch.object(retrieve.urllib.request, 'urlopen', side_effect=AssertionError('network opened')):
                with self.assertRaisesRegex(ValueError, 'contract_denied'):
                    retrieve.retrieve(out)
            self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
