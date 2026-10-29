"""Synthetic unit/definition tests; no participant records or model results."""
import unittest
import numpy as np
import pandas as pd
from research.mortality.routine29 import features as f


def source(**changes):
    values = dict(age_linear=[-1.], age_quadratic=[1.], male=[1.], current=[0.],
                  former=[0.], sbp=[0.], HSAGEIR=[50.], HSSEX=[1.],
                  CEP=[(1. + .184) / .960], GHP=[5.5], GHPMETH=[1.],
                  URP=[100.], UBP=[10.], BMPBMI=[27.], BMPWTFLG=[0.], BMPHTFLG=[0.])
    values.update(changes)
    return pd.DataFrame(values)


class RoutineFeaturesTests(unittest.TestCase):
    def test_serum_calibration_and_male_egfr_known_value(self):
        d = f.engineer(source())
        self.assertAlmostEqual(d.creatinine_calibrated.iloc[0], 1., places=14)
        self.assertAlmostEqual(d.egfr.iloc[0], 91.69147860983348, places=8)
        self.assertAlmostEqual(d.egfr_low.iloc[0], 4.)
        self.assertAlmostEqual(d.egfr_high.iloc[0], (91.69147860983348 - 60) / 30, places=8)
        self.assertTrue(f.primary_complete(d)[0])

    def test_female_egfr_and_no_race_term(self):
        d = f.engineer(source(HSSEX=[2.], male=[0.]))
        # Adult2021 published equation for woman50, IDMS-alignedScr1mg/dL.
        expected = 142 * (1 / .7) ** -1.2 * .9938 ** 50 * 1.012
        self.assertAlmostEqual(d.egfr.iloc[0], expected, places=12)
        other = source(HSSEX=[2.], male=[0.]);other['race'] = 'different'
        self.assertEqual(d.egfr.iloc[0], f.engineer(other).egfr.iloc[0])

    def test_serum_correction_never_applies_to_urine_creatinine(self):
        d = f.engineer(source(CEP=[1.], URP=[100.], UBP=[10.]))
        self.assertAlmostEqual(d.creatinine_calibrated.iloc[0], .776)
        self.assertEqual(d.URP.iloc[0], 100.)
        self.assertAlmostEqual(d.uacr.iloc[0], 10.)
        self.assertAlmostEqual(d.log_uacr.iloc[0], np.log2(10))

    def test_albumin_nondetect_released_fill_is_retained_and_flagged(self):
        d = f.engineer(source(UBP=[.4], URP=[100.]))
        self.assertAlmostEqual(d.uacr.iloc[0], .4)
        self.assertTrue(d.urine_albumin_reported_below_detection.iloc[0])
        self.assertTrue(f.primary_complete(d)[0])
        self.assertFalse(f.strict_metadata_mask(d)[0])

    def test_suspect_urine_creatinine_is_excluded_before_ratio(self):
        for value in (7.9, 9.99, 0., -1., np.nan, np.inf, 88888):
            with self.subTest(value=value):
                d = f.engineer(source(URP=[value]))
                self.assertTrue(pd.isna(d.uacr.iloc[0]))
                self.assertFalse(f.primary_complete(d)[0])
        d = f.engineer(source(URP=[10.], UBP=[1.]))
        self.assertAlmostEqual(d.uacr.iloc[0], 10.)

    def test_hba1c_affinity_is_primary_valid_but_not_strict_hplc(self):
        d = f.engineer(source(GHP=[6.2], GHPMETH=[4.]))
        self.assertEqual(d.hba1c.iloc[0], 6.2)
        self.assertTrue(f.primary_complete(d)[0])
        self.assertFalse(f.strict_metadata_mask(d)[0])
        for code in (0., 8., np.nan):
            d = f.engineer(source(GHPMETH=[code]))
            self.assertTrue(pd.isna(d.hba1c.iloc[0]))

    def test_source_missing_sentinels_removed_before_calibration_or_logs(self):
        for field, value, feature in [('CEP',8888,'egfr'),('GHP',8888,'hba1c'),
                                     ('URP',88888,'uacr'),('UBP',888888,'uacr'),
                                     ('BMPBMI',8888,'bmi')]:
            with self.subTest(field=field):
                d = f.engineer(source(**{field:[value]}))
                self.assertTrue(pd.isna(d[feature].iloc[0]))
                self.assertFalse(f.primary_complete(d)[0])

    def test_nonpositive_calibrated_creatinine_never_floored(self):
        d = f.engineer(source(CEP=[.1]))
        self.assertTrue(pd.isna(d.creatinine_calibrated.iloc[0]))
        self.assertTrue(pd.isna(d.egfr.iloc[0]))
        self.assertTrue(d.serum_creatinine_calibration_nonpositive.iloc[0])

    def test_positive_finite_values_not_clipped_by_frequency_ranges(self):
        d = f.engineer(source(CEP=[.2], GHP=[18.], BMPBMI=[90.], UBP=[20000.], URP=[900.]))
        self.assertEqual(d.hba1c.iloc[0],18.)
        self.assertAlmostEqual(d.bmi.iloc[0],(90-27)/5)
        self.assertAlmostEqual(d.creatinine_calibrated.iloc[0],.008)
        self.assertAlmostEqual(d.uacr.iloc[0],100*20000/900)
        self.assertTrue(f.primary_complete(d)[0])

    def test_substitution_flags_retained_in_main_and_restricted_only_in_sensitivity(self):
        for field, value in [('BMPWTFLG',3.),('BMPHTFLG',1.),('BMPHTFLG',np.nan)]:
            with self.subTest(field=field):
                d = f.engineer(source(**{field:[value]}))
                self.assertTrue(f.primary_complete(d)[0])
                self.assertFalse(f.strict_metadata_mask(d)[0])
        self.assertTrue(f.strict_metadata_mask(f.engineer(source()))[0])

    def test_invalid_sex_or_child_age_not_treated_as_adult_male(self):
        for changes in ({'HSSEX':[8.]},{'HSAGEIR':[17.]},{'HSAGEIR':[np.nan]}):
            d = f.engineer(source(**changes))
            self.assertTrue(pd.isna(d.egfr.iloc[0]))
            self.assertFalse(f.primary_complete(d)[0])

    def test_index_and_raw_source_values_preserved(self):
        raw = source();raw.index = [91];before = raw.copy(deep=True)
        d = f.engineer(raw)
        pd.testing.assert_frame_equal(raw,before)
        pd.testing.assert_frame_equal(d[before.columns],before)
        self.assertEqual(d.index.tolist(),[91])
        self.assertIsNot(d,raw)

    def test_core_has_one_creatinine_representation_and_no_protein_leakage(self):
        self.assertEqual(len(f.CORE),11)
        self.assertEqual(tuple(f.CORE[:6]),('age_linear','age_quadratic','male','current','former','sbp'))
        self.assertEqual(tuple(f.CORE[6:]),('bmi','hba1c','egfr_low','egfr_high','log_uacr'))
        for forbidden in ('CEP','creatinine_calibrated','SSCYSTAT','SSB2M','SSBTP','UBP','URP','diabetes_history'):
            self.assertNotIn(forbidden,f.CORE)

    def test_missing_required_source_field_fails(self):
        raw = source().drop(columns='URP')
        with self.assertRaisesRegex(ValueError,'Missing prespecified'):f.engineer(raw)

if __name__ == '__main__':unittest.main()
