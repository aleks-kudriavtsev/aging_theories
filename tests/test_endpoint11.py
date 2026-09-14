"""Artificial controls and optional data-integrity audit. Not clinical validation."""
from pathlib import Path
import json,math,os,unittest
from research.mortality.endpoint_contract11 import (PUBLIC_CAUSES,allowed_codes,classify_public,
 normalize_public_cause,observability,fixed_horizon_states,check_linkage_lineage,
 precision_planning,mec_weight_1999_2010,profile_transfer_gate)

class CodingTests(unittest.TestCase):
    def test_ten_plus_unknown(self):self.assertEqual(len(PUBLIC_CAUSES),11)
    def test_early_ten_groups(self):self.assertEqual(len(allowed_codes(1999)),10)
    def test_2013_ten_groups(self):self.assertEqual(len(allowed_codes(2013)),10)
    def test_late_three_groups(self):self.assertEqual(allowed_codes(2015),{'001','002','010'})
    def test_invalid_cycle(self):
        for y in (1998,2000,2019,True,'2011'):
            with self.subTest(y=y),self.assertRaises(ValueError):allowed_codes(y)
    def test_numeric_decoder(self):self.assertEqual(normalize_public_cause(5.),'005')
    def test_exact_string_decoder(self):self.assertEqual(normalize_public_cause(' 005 '),'005')
    def test_unknown_values(self):
        for value in (None,'','.',float('nan'),'UNK'):self.assertIsNone(normalize_public_cause(value))
    def test_unreviewed_values(self):
        for value in (True,5.5,math.inf,[],{},'5','I63','070'):
            with self.subTest(value=value),self.assertRaises(ValueError):normalize_public_cause(value)
    def test_alive(self):self.assertEqual(classify_public(0,None,survey_start=1999),'alive')
    def test_unknown_death_not_alive(self):self.assertEqual(classify_public(1,None,survey_start=1999),'UNK')
    def test_living_with_cause_refused(self):
        with self.assertRaises(ValueError):classify_public(0,'005',survey_start=2011)
    def test_ineligible_not_alive(self):
        with self.assertRaises(ValueError):classify_public(3,None,survey_start=2011)
    def test_unavailable_late_stroke_not_negative(self):
        with self.assertRaises(ValueError):classify_public(1,'005',survey_start=2015)
    def test_actual_public_group(self):self.assertEqual(classify_public(1,'005',survey_start=2009),'005')

class ObservabilityTests(unittest.TestCase):
    def test_ihd_not_all_heart(self):self.assertFalse(observability('ihd',survey_start=2009)['label_observable'])
    def test_stroke_group_exact(self):self.assertTrue(observability('cerebrovascular',survey_start=2009)['label_observable'])
    def test_ischemic_stroke_not_group(self):self.assertFalse(observability('ischemic_stroke',survey_start=2009)['label_observable'])
    def test_alzheimer_not_dementia(self):
        self.assertTrue(observability('alzheimer',survey_start=2009)['label_observable'])
        self.assertFalse(observability('dementia_F01_F03',survey_start=2009)['label_observable'])
    def test_renal_not_ckd(self):
        self.assertTrue(observability('renal_broad',survey_start=2009)['label_observable'])
        self.assertFalse(observability('ckd_N18',survey_start=2009)['label_observable'])
    def test_cancer_subtype_unobserved(self):self.assertFalse(observability('lung_cancer',survey_start=2009)['label_observable'])
    def test_no_noninfectious_rename(self):self.assertFalse(observability('noninfectious_natural',survey_start=2009)['label_observable'])
    def test_accidents_not_all_external(self):self.assertFalse(observability('all_external',survey_start=2009)['label_observable'])
    def test_unknown_target(self):
        with self.assertRaises(ValueError):observability('invented',survey_start=2009)
    def test_late_stroke_unobserved(self):self.assertEqual(observability('cerebrovascular',survey_start=2015)['reason'],'source_coarsened_for_cycle')
    def test_absent_assay_blocks(self):
        r=profile_transfer_gate(['egfr','crp'],['egfr'],target='cerebrovascular',survey_start=2011)
        self.assertFalse(r['same_model_replay_supported_by_metadata']);self.assertEqual(r['missing_assays'],['crp'])
    def test_same_names_not_clinical_ready(self):self.assertFalse(profile_transfer_gate(['egfr'],['egfr'],target='heart',survey_start=2009)['clinical_use_ready'])
    def test_metadata_shape_error(self):
        with self.assertRaises(ValueError):profile_transfer_gate('CRP',[],target='heart',survey_start=2009)
    def test_duplicate_assays(self):
        with self.assertRaises(ValueError):profile_transfer_gate(['CRP','CRP'],[],target='heart',survey_start=2009)

class FollowupTests(unittest.TestCase):
    def test_horizon_boundary(self):self.assertEqual(fixed_horizon_states([1,1,0],['005','001',None],[60,61,60],survey_start=2009,horizon_years=5),['005','alive','alive'])
    def test_death_unknown_at_horizon(self):self.assertEqual(fixed_horizon_states([1],[None],[12],survey_start=2009,horizon_years=1),['UNK'])
    def test_censoring_refused(self):
        with self.assertRaises(ValueError):fixed_horizon_states([0],[None],[59],survey_start=2009,horizon_years=5)
    def test_invalid_arrays(self):
        with self.assertRaises(ValueError):fixed_horizon_states([1],[],[12],survey_start=2009,horizon_years=1)
    def test_negative_followup(self):
        with self.assertRaises(ValueError):fixed_horizon_states([1],[None],[-1],survey_start=2009,horizon_years=1)
    def test_zero_month_death_kept(self):self.assertEqual(fixed_horizon_states([1],['005'],[0],survey_start=2009,horizon_years=1),['005'])

class LineageAndWeightTests(unittest.TestCase):
    def test_mixed_public_restricted_refused(self):
        with self.assertRaises(ValueError):check_linkage_lineage(['public_lmf2019','restricted_lmf2022'])
    def test_mixed_vintage_refused(self):
        with self.assertRaises(ValueError):check_linkage_lineage(['restricted_lmf2019','restricted_lmf2022'])
    def test_not_authorization(self):self.assertFalse(check_linkage_lineage(['restricted_lmf2022'])['data_access_authorized_by_this_check'])
    def test_early_weights_use_four_year(self):self.assertEqual(mec_weight_1999_2010(1999,300,600),200)
    def test_second_cycle_four_year(self):self.assertEqual(mec_weight_1999_2010(2001,300,600),200)
    def test_later_weights_two_year(self):self.assertEqual(mec_weight_1999_2010(2003,600),100)
    def test_missing_four_year_refused(self):
        with self.assertRaises(ValueError):mec_weight_1999_2010(1999,600)
    def test_count_precision(self):self.assertEqual(precision_planning(100)['poisson_relative_se_approx'],.1)
    def test_zero_count_not_precision(self):self.assertIsNone(precision_planning(0)['poisson_relative_se_approx'])
    def test_no_EPV_claim(self):self.assertFalse(precision_planning(30)['model_sample_size_criterion'])

class ActualDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        r=os.environ.get('ENDPOINT11_RESULT_DIR')
        if not r:raise unittest.SkipTest('Needs actual endpoint-feasibility aggregate output')
        import pandas as pd
        cls.path=Path(r);cls.r=json.loads((cls.path/'results.json').read_text())
        cls.e=pd.read_csv(cls.path/'event_support.csv',dtype={'period':str,'cause_code':str})
    def test_cohort_counts(self):self.assertEqual((self.r['n_eligible'],self.r['n_complete_profile']),(17629,15273))
    def test_complete_subsets_reproduce_history(self):self.assertEqual(sum(r['complete_profile_n'] for r in self.r['flow'] if r['cycle'] in [2005,2007]),5287)
    def test_reader_checks(self):self.assertEqual(self.r['independent_reader_cells_compared'],497280)
    def test_deaths_partition(self):
        for _,g in self.e.groupby(['period','domain','stratum','horizon_years']):self.assertEqual(g.cause_deaths.sum(),g.all_deaths.iloc[0])
    def test_complete_stroke_count(self):self.assertEqual(self.e.query("period=='1999_2010' and domain=='complete_profile' and stratum=='all' and horizon_years==8 and cause_code=='005'").cause_deaths.iloc[0],94)
    def test_no_validation_upgrade(self):self.assertFalse(self.r['new_independent_validation']);self.assertFalse(self.r['clinical_use_ready'])

if __name__=='__main__':unittest.main()
