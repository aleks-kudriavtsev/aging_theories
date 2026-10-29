"""Synthetic/mocked source and domain checks; no public participant records."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from research.mortality.b2m31 import cohort, retrieve
from research.mortality.b2m31.features import CORE


INDEX_URL='https://wwwn.cdc.gov/catalogue.htm'
ASSAY_URL='https://wwwn.cdc.gov/nchs/data/nhanes/1999/SSCARD_A.xpt'
DEMO_URL='https://wwwn.cdc.gov/nchs/data/nhanes/2001/DEMO_B.xpt'


def small_sscard():
    frame=pd.DataFrame(dict(SEQN=[1.,2.,3.,4.,5.,6.],
        WTSSCB2Y=[100.,200.,300.,400.,500.,600.],
        WTSSCB4Y=[150.,250.,350.,450.,np.nan,np.nan],
        SSCYST=[1.,.16,8.5,np.nan,1.2,1.3],SSCYSTL=[0.,1.,2.,np.nan,0.,0.],
        SSB2M=[1.5,.5,24.,np.nan,1.6,1.7],SSB2ML=[0.,1.,2.,np.nan,0.,0.],
        SSPRIS=[1.,0.,1.,0.,1.,1.]))
    spec=dict(expected_rows_full_release=6,required_fields=list(frame),
        metadata_marginal_flag_counts={name:{'0':3,'1':1,'2':1,'missing':1} for name in ('SSCYSTL','SSB2ML')},
        metadata_marginal_pristine_counts={'SSPRIS':{'0':2,'1':4,'missing':0}},
        metadata_weight_counts={'WTSSCB2Y':{'positive':6,'zero':0,'missing':0},
                               'WTSSCB4Y':{'positive':4,'zero':0,'missing':2}})
    return frame,spec


class ReleaseLinkTests(unittest.TestCase):
    def test_unique_official_link_resolves_relative_and_preserves_case(self):
        url=retrieve.resolve_sscard('<a href="/nchs/data/nhanes/1999/SSCARD_A.XPT">XPT</a>',INDEX_URL)
        self.assertEqual(url,'https://wwwn.cdc.gov/nchs/data/nhanes/1999/SSCARD_A.XPT')
    def test_identical_duplicate_anchors_are_not_false_ambiguity(self):
        html=f'<a href="{ASSAY_URL}">XPT</a>'*2
        self.assertEqual(retrieve.resolve_sscard(html,INDEX_URL),ASSAY_URL)
    def test_distinct_official_links_are_ambiguous_and_never_guessed(self):
        html=f'<a href="{ASSAY_URL}">old</a><a href="https://ftp.cdc.gov/new/SSCARD_A.xpt">new</a>'
        with self.assertRaisesRegex(ValueError,'Unique official'):retrieve.resolve_sscard(html,INDEX_URL)
        for html in ('<a href="SSCARD_A.htm">documentation</a>','<p>No data anchor</p>'):
            with self.assertRaisesRegex(ValueError,'no guessed URL'):retrieve.resolve_sscard(html,INDEX_URL)
    def test_unsafe_or_lookalike_urls_cannot_supply_the_release(self):
        urls=['http://wwwn.cdc.gov/data/SSCARD_A.xpt','https://wwwn.cdc.gov.evil.invalid/SSCARD_A.xpt',
              'https://user:password@wwwn.cdc.gov/data/SSCARD_A.xpt','https://wwwn.cdc.gov:8443/data/SSCARD_A.xpt',
              'https://wwwn.cdc.gov/data/SSCARD_A.xpt#fragment','https://wwwn.cdc.gov:bad/data/SSCARD_A.xpt']
        for url in urls:
            with self.subTest(url=url):
                self.assertFalse(retrieve.allowed_url(url))
                with self.assertRaises(ValueError):retrieve.resolve_sscard(f'<a href="{url}">data</a>',INDEX_URL)
    def test_preflight_metadata_check_precedes_directory_and_network(self):
        with patch.object(retrieve,'check_study',side_effect=ValueError('locked metadata changed')):
            with patch.object(retrieve.urllib.request,'urlopen',side_effect=AssertionError('network')):
                with tempfile.TemporaryDirectory() as tmp:
                    out=Path(tmp)/'must_not_exist'
                    with self.assertRaisesRegex(ValueError,'locked metadata'):retrieve.retrieve(out)
                    self.assertFalse(out.exists())


class HistoricalPinTests(unittest.TestCase):
    def fixture(self,root):
        study=root/'study';study.mkdir();source=root/'sources';source.mkdir()
        history=root/'nhanes04';history.mkdir()
        (study/'SOURCE_SPEC.json').write_text(json.dumps({'sscard_index_url':INDEX_URL}))
        (study/'LOCK.json').write_text('{"synthetic_lock":true}\n')
        content={'index.htm':f'<html><a href="{ASSAY_URL}">XPT</a></html>'.encode(),
                 '2001_DEMO_B.xpt':b'HEADER RECORD synthetic DEMO fixture',
                 '1999_SSCARD_A.xpt':b'HEADER RECORD synthetic SSCARD fixture'}
        urls={'index.htm':INDEX_URL,'2001_DEMO_B.xpt':DEMO_URL,'1999_SSCARD_A.xpt':ASSAY_URL}
        rows=[]
        for name,value in content.items():
            (source/name).write_bytes(value)
            rows.append(dict(file=name,url=urls[name],resolved_url=urls[name],status='retrieved',
                             bytes=len(value),sha256=hashlib.sha256(value).hexdigest()))
        doc=dict(schema='b2m31-retrieval/1',plan_commit=cohort.PLAN_COMMIT,
                 lock_sha256=retrieve.sha(study/'LOCK.json'),files=rows)
        (source/'manifest.json').write_text(json.dumps(doc))
        (history/'SOURCE_DATA_SHA256.txt').write_text('\n'.join(f'{r["sha256"]} {r["file"]}' for r in rows if r['file'].endswith('.xpt'))+'\n')
        return study,source,doc,history
    def mocked(self,study):
        return (patch.object(cohort,'ROOT',study),patch.object(cohort,'check_study',return_value={}),
                patch.object(cohort,'static_urls',return_value={'2001_DEMO_B.xpt':DEMO_URL}))
    def test_full_manifest_and_historical_pins_are_verified_without_table_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            study,source,doc,_=self.fixture(Path(tmp));a,b,c=self.mocked(study)
            with a,b,c,patch.object(cohort,'read_xpt_checked',side_effect=AssertionError('table parsed')):
                observed,audit=cohort.verify_sources(source)
            self.assertEqual(observed,doc);self.assertEqual(audit['historical_data_hashes_verified'],2)
            self.assertEqual(audit['files_verified'],3)
    def test_tampered_bytes_fail_before_any_participant_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            study,source,_,_=self.fixture(Path(tmp));(source/'2001_DEMO_B.xpt').write_bytes(b'tampered')
            a,b,c=self.mocked(study)
            with a,b,c,patch.object(cohort,'read_xpt_checked',side_effect=AssertionError('table parsed')):
                with self.assertRaisesRegex(ValueError,'Source bytes changed'):cohort.load(source)
    def test_new_self_consistent_manifest_cannot_override_historical_data_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            study,source,doc,_=self.fixture(Path(tmp));new=b'HEADER RECORD changed synthetic release'
            (source/'2001_DEMO_B.xpt').write_bytes(new)
            row=next(r for r in doc['files'] if r['file']=='2001_DEMO_B.xpt')
            row.update(bytes=len(new),sha256=hashlib.sha256(new).hexdigest())
            (source/'manifest.json').write_text(json.dumps(doc));a,b,c=self.mocked(study)
            with a,b,c,patch.object(cohort,'read_xpt_checked',side_effect=AssertionError('table parsed')):
                with self.assertRaisesRegex(ValueError,'Historical public data release'):cohort.load(source)
    def test_missing_historical_pin_blocks_even_otherwise_matching_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            study,source,_,history=self.fixture(Path(tmp));(history/'SOURCE_DATA_SHA256.txt').write_text('')
            a,b,c=self.mocked(study)
            with a,b,c,patch.object(cohort,'read_xpt_checked',side_effect=AssertionError('table parsed')):
                with self.assertRaisesRegex(ValueError,'Historical public data release'):cohort.load(source)
    def test_protocol_rebinding_duplicate_or_wrong_release_is_rejected(self):
        changes=[lambda d:d.update(plan_commit='wrong'),lambda d:d['files'].append(deepcopy(d['files'][0])),
                 lambda d:d['files'][1].update(url='https://wwwn.cdc.gov/unapproved/DEMO_B.xpt'),
                 lambda d:d['files'][1].update(file='../DEMO_B.xpt')]
        for change in changes:
            with tempfile.TemporaryDirectory() as tmp:
                study,source,doc,_=self.fixture(Path(tmp));change(doc);(source/'manifest.json').write_text(json.dumps(doc))
                a,b,c=self.mocked(study)
                with a,b,c,patch.object(cohort,'read_xpt_checked',side_effect=AssertionError('table parsed')):
                    with self.assertRaises(ValueError):cohort.load(source)


class SixYearWeightTests(unittest.TestCase):
    def test_first_two_cycles_use_four_year_weight_third_uses_two_year(self):
        d=pd.DataFrame(dict(SDDSRVYR=[1.,2.,3.],WTSSCB4Y=[300.,600.,np.nan],WTSSCB2Y=[900.,1200.,1500.]),index=[9,4,12])
        w=cohort.six_year_weights(d)
        np.testing.assert_allclose(w,[200.,400.,500.]);self.assertEqual(w.index.tolist(),[9,4,12])
        self.assertNotEqual(w.iloc[0],d.WTSSCB2Y.iloc[0]/3)
    def test_missing_four_year_weight_in_third_cycle_does_not_exclude_it(self):
        d=pd.DataFrame(dict(SDDSRVYR=[3.],WTSSCB4Y=[np.nan],WTSSCB2Y=[900.]))
        self.assertAlmostEqual(cohort.six_year_weights(d).iloc[0],300.)
    def test_selected_nonpositive_nonfinite_weights_stay_missing(self):
        d=pd.DataFrame(dict(SDDSRVYR=[1,2,3,3],WTSSCB4Y=[0.,-1.,999.,999.],WTSSCB2Y=[999.,999.,np.nan,np.inf]))
        self.assertTrue(cohort.six_year_weights(d).isna().all())
    def test_unsupported_or_unknown_cycle_does_not_guess_a_formula(self):
        for value in (4,np.nan):
            d=pd.DataFrame(dict(SDDSRVYR=[value],WTSSCB4Y=[100.],WTSSCB2Y=[100.]))
            with self.assertRaisesRegex(ValueError,'Unknown survey cycle'):cohort.six_year_weights(d)


class AssayReleaseTests(unittest.TestCase):
    def test_generated_correct_release_satisfies_full_release_marginals(self):
        frame,spec=small_sscard();cohort.validate_sscard(frame,spec)
    def test_wrong_source_count_schema_and_duplicate_ids_fail(self):
        for edit in (lambda d:d.iloc[:-1],lambda d:d.drop(columns='SSB2ML'),lambda d:d.assign(SEQN=[1,1,3,4,5,6])):
            frame,spec=small_sscard()
            with self.assertRaises(ValueError):cohort.validate_sscard(edit(frame),spec)
    def test_wrong_flag_counts_unknown_codes_or_fake_epsilon_zero_fail(self):
        for value in (1.,3.,5.397605346934028e-79):
            frame,spec=small_sscard();frame.loc[0,'SSCYSTL']=value
            with self.assertRaises(ValueError):cohort.validate_sscard(frame,spec)
    def test_pristine_and_corrected_weight_marginals_cannot_be_silently_changed(self):
        for column,value in [('SSPRIS',3.),('WTSSCB2Y',0.),('WTSSCB4Y',np.inf)]:
            frame,spec=small_sscard();frame.loc[0,column]=value
            with self.assertRaises(ValueError):cohort.validate_sscard(frame,spec)
    def test_result_and_flag_missingness_must_be_aligned(self):
        frame,spec=small_sscard();frame.loc[0,'SSCYST']=np.nan
        with self.assertRaisesRegex(ValueError,'missingness'):cohort.validate_sscard(frame,spec)
    def test_nonpositive_or_nonfinite_released_result_cannot_pass_as_measurable(self):
        for value in (0.,-1.,np.inf):
            frame,spec=small_sscard();frame.loc[0,'SSB2M']=value
            with self.assertRaisesRegex(ValueError,'Invalid released protein'):cohort.validate_sscard(frame,spec)
    def test_primary_flags00_sensitivity01_and_no_flag2_in_either_arm_domain(self):
        flags=[(0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(np.nan,0),(0,0),(0,0)]
        d=pd.DataFrame({name:np.zeros(len(flags)) for name in CORE})
        d['SSCYST']=1.;d['SSB2M']=1.5
        d['SSCYSTL']=[p[0] for p in flags];d['SSB2ML']=[p[1] for p in flags]
        d.loc[7,'SSB2M']=0.;d.loc[8,'bmi']=np.nan
        primary,reported=cohort.domain_masks(d)
        np.testing.assert_array_equal(primary,[True,False,False,False,False,False,False,False,False])
        np.testing.assert_array_equal(reported,[True,True,True,True,False,False,False,False,False])
        self.assertFalse((primary&~reported).any())


class FullDesignInventoryTests(unittest.TestCase):
    def full(self):return pd.DataFrame(dict(SDMVSTRA=[1,1,1,14,14],SDMVPSU=[1,2,3,1,2]))
    def test_three_psu_stratum_is_preserved_and_not_assumed_two(self):
        full=self.full();eligible=pd.concat([full,full.iloc[[0,3]]],ignore_index=True)
        audit=cohort.validate_design_inventory(full,eligible)
        self.assertTrue(audit['full_DEMO_design_inventory_preserved'])
        self.assertEqual(audit['full_DEMO_masked_psu_n'],5)
        self.assertEqual(audit['eligible_psu_count_distribution'],{'2':1,'3':1})
    def test_dropped_full_demo_unit_is_a_blocker_even_if_two_psus_remain(self):
        full=self.full();eligible=full.loc[~((full.SDMVSTRA==1)&(full.SDMVPSU==3))]
        with self.assertRaisesRegex(ValueError,'removed a full-DEMO masked PSU'):cohort.validate_design_inventory(full,eligible)
    def test_singleton_blocks_even_when_inventory_exactly_matches(self):
        full=pd.DataFrame(dict(SDMVSTRA=[1,14,14],SDMVPSU=[1,1,2]))
        with self.assertRaisesRegex(ValueError,'Singleton'):cohort.validate_design_inventory(full,full.copy())
    def test_pair_identity_keeps_same_psu_code_in_different_strata_distinct(self):
        full=pd.DataFrame(dict(SDMVSTRA=[1,1,14,14],SDMVPSU=[1,2,1,2]))
        audit=cohort.validate_design_inventory(full,full.copy())
        self.assertEqual(audit['full_DEMO_masked_psu_n'],4)
        bad=full.copy();bad.loc[2,'SDMVSTRA']=1
        with self.assertRaises(ValueError):cohort.validate_design_inventory(full,bad)

if __name__=='__main__':unittest.main()
