"""Independent synthetic source/formula tests; no NHANES participant input."""
import unittest

import numpy as np
import pandas as pd

from research.mortality.b2m31 import features as f


def source(**changes):
    values = dict(cycle=[2001.], RIDAGEYR=[50.], RIAGENDR=[1.],
                  SMQ020=[1.], SMQ040=[1.], BPXSY1=[120.], BPXSY2=[130.],
                  BPXSY3=[140.], BMXBMI=[27.], LBXGH=[5.5],
                  raw_creatinine=[1.], URXUMA=[10.], URXUCR=[100.])
    values.update(changes)
    return pd.DataFrame(values)


def egfr2021(age, female, serum):
    """Independent scalar form of the published adult race-free equation."""
    kappa, alpha = (.7, -.241) if female else (.9, -.302)
    ratio = serum / kappa
    return (142 * min(ratio, 1) ** alpha * max(ratio, 1) ** -1.2
            * .9938 ** age * (1.012 if female else 1))


class ContinuousNHAFeaturesTests(unittest.TestCase):
    def test_cycle_specific_serum_standardization_matches_official_numeric_formula(self):
        raw = pd.concat([source(cycle=[cycle]) for cycle in (1999, 2001, 2003)], ignore_index=True)
        d = f.engineer(raw)
        np.testing.assert_allclose(d.creatinine_calibrated, [1.160, 1., 1.], atol=1e-15, rtol=0)
        np.testing.assert_array_equal(d.serum_creatinine_calibration_applied, [True, False, False])
        self.assertAlmostEqual(d.egfr.iloc[0], egfr2021(50, False, 1.160), places=12)
        self.assertAlmostEqual(d.egfr.iloc[1], 91.69147860983348, places=10)
        self.assertEqual(d.egfr.iloc[1], d.egfr.iloc[2])
        # Do not silently reuse the earlier NHANESIII serum correction.
        self.assertNotAlmostEqual(d.creatinine_calibrated.iloc[0], .960 - .184)

    def test_equivalent_calibrated_serum_produces_same_egfr_across_cycles(self):
        raw = pd.concat([source(cycle=[1999.], raw_creatinine=[(1 - .147) / 1.013]),
                         source(cycle=[2001.]), source(cycle=[2003.])], ignore_index=True)
        d = f.engineer(raw)
        np.testing.assert_allclose(d.egfr, [91.69147860983348] * 3, atol=1e-10, rtol=0)
        self.assertTrue(f.primary_complete(d).all())

    def test_female_equation_and_renal_piecewise_terms(self):
        d = f.engineer(source(RIAGENDR=[2.], RIDAGEYR=[70.], raw_creatinine=[2.]))
        expected = 142 * (2 / .7) ** -1.2 * .9938 ** 70 * 1.012
        self.assertAlmostEqual(d.egfr.iloc[0], expected, places=12)
        self.assertLess(expected, 60)
        self.assertAlmostEqual(d.egfr_low.iloc[0], expected / 15)
        self.assertEqual(d.egfr_high.iloc[0], 0.)
        high = f.engineer(source())
        self.assertEqual(high.egfr_low.iloc[0], 4.)
        self.assertAlmostEqual(high.egfr_high.iloc[0], (91.69147860983348 - 60) / 30)

    def test_all_core_terms_follow_registered_source_order_and_transforms(self):
        d = f.engineer(source(RIDAGEYR=[70.], RIAGENDR=[2.], SMQ040=[3.],
                              BPXSY1=[150.], BPXSY2=[170.], BPXSY3=[160.],
                              BMXBMI=[32.], LBXGH=[6.2], URXUMA=[3.], URXUCR=[50.]))
        renal = egfr2021(70, True, 1.)
        expected = [1., 1., 0., 0., 1., 1.5, 1., 6.2,
                    min(renal, 60) / 15, max(renal - 60, 0) / 30, np.log2(6.)]
        np.testing.assert_allclose(d[list(f.CORE)].iloc[0].to_numpy(float), expected, atol=1e-12, rtol=0)

    def test_below_kappa_equation_uses_sex_specific_exponent_without_race(self):
        for sex, female in ((1., False), (2., True)):
            raw = source(RIAGENDR=[sex], raw_creatinine=[.5])
            d = f.engineer(raw)
            self.assertAlmostEqual(d.egfr.iloc[0], egfr2021(50, female, .5), places=12)
            raw['RIDRETH1'] = 5.; raw['cystatin_c'] = 1000.
            self.assertEqual(d.egfr.iloc[0], f.engineer(raw).egfr.iloc[0])

    def test_smoking_never_current_some_days_and_former_definitions(self):
        raw = pd.concat([source(SMQ020=[2.], SMQ040=[np.nan]),
                         source(SMQ020=[1.], SMQ040=[1.]),
                         source(SMQ020=[1.], SMQ040=[2.]),
                         source(SMQ020=[1.], SMQ040=[3.])], ignore_index=True)
        d = f.engineer(raw)
        np.testing.assert_array_equal(d.current, [0., 1., 1., 0.])
        np.testing.assert_array_equal(d.former, [0., 0., 0., 1.])
        self.assertFalse(d.smoking_unknown.any()); self.assertTrue(f.primary_complete(d).all())

    def test_smoking_refusals_unknown_and_incomplete_skip_are_missing(self):
        for smoked, now in ((7., 1.), (9., 1.), (np.nan, 1.), (1., 7.),
                             (1., 9.), (1., np.nan), (1., 0.), (3., 3.)):
            with self.subTest(smoked=smoked, now=now):
                d = f.engineer(source(SMQ020=[smoked], SMQ040=[now]))
                self.assertTrue(pd.isna(d.current.iloc[0])); self.assertTrue(pd.isna(d.former.iloc[0]))
                self.assertTrue(d.smoking_unknown.iloc[0]); self.assertFalse(f.primary_complete(d)[0])

    def test_sbp_uses_only_first_three_positive_finite_research_readings(self):
        raw = source(BPXSY1=[0.], BPXSY2=[140.], BPXSY3=[160.])
        raw['BPXSY4'] = 1000.; raw['BPXSAR'] = 20.
        d = f.engineer(raw)
        self.assertEqual(d.sbp_measurement_count.iloc[0], 2)
        self.assertEqual(d.sbp.iloc[0], 1.)
        one = f.engineer(source(BPXSY1=[-1.], BPXSY2=[np.inf], BPXSY3=[135.]))
        self.assertEqual(one.sbp_measurement_count.iloc[0], 1)
        self.assertEqual(one.sbp.iloc[0], .25)
        none = f.engineer(source(BPXSY1=[0.], BPXSY2=[np.nan], BPXSY3=[-1.]))
        self.assertEqual(none.sbp_measurement_count.iloc[0], 0)
        self.assertTrue(none.routine_missing_sbp.iloc[0]); self.assertFalse(f.primary_complete(none)[0])

    def test_screening_age_bounds_and_wrong_sex_cannot_form_complete_core(self):
        for age in (39., 76., 50.5, np.nan, np.inf):
            with self.subTest(age=age):
                d = f.engineer(source(RIDAGEYR=[age]))
                self.assertTrue(pd.isna(d.age_linear.iloc[0])); self.assertTrue(pd.isna(d.egfr.iloc[0]))
                self.assertTrue(d.source_invalid_screening_age.iloc[0]); self.assertFalse(f.primary_complete(d)[0])
        for sex in (0., 3., 7., 9., np.nan, np.inf, True):
            with self.subTest(sex=sex):
                d = f.engineer(source(RIAGENDR=[sex]))
                self.assertTrue(pd.isna(d.male.iloc[0])); self.assertTrue(pd.isna(d.egfr.iloc[0]))
                self.assertTrue(d.source_invalid_sex.iloc[0]); self.assertFalse(f.primary_complete(d)[0])
        for age in (40., 75.):
            self.assertTrue(f.primary_complete(f.engineer(source(RIDAGEYR=[age])))[0])

    def test_urine_units_are_applied_once_and_small_positive_denominator_is_valid(self):
        d = f.engineer(source(cycle=[1999.], URXUMA=[10.], URXUCR=[2.]))
        self.assertEqual(d.URXUCR.iloc[0], 2.)
        self.assertEqual(d.uacr.iloc[0], 500.)
        self.assertEqual(d.log_uacr.iloc[0], np.log2(500.))
        self.assertTrue(f.primary_complete(d)[0])
        self.assertAlmostEqual(d.creatinine_calibrated.iloc[0], 1.160)

    def test_original_positive_hba1c_bmi_and_albumin_values_not_regressed_or_clipped(self):
        d = f.engineer(source(LBXGH=[18.], BMXBMI=[90.], URXUMA=[.1], URXUCR=[900.]))
        self.assertEqual(d.hba1c.iloc[0], 18.)
        self.assertEqual(d.bmi.iloc[0], (90 - 27) / 5)
        self.assertAlmostEqual(d.uacr.iloc[0], 100 * .1 / 900)
        self.assertTrue(f.primary_complete(d)[0])

    def test_invalid_lab_and_bmi_values_propagate_to_only_dependent_terms(self):
        for name, output in [('BMXBMI', 'bmi'), ('LBXGH', 'hba1c'),
                             ('raw_creatinine', 'egfr'), ('URXUMA', 'log_uacr'),
                             ('URXUCR', 'log_uacr')]:
            for value in (0., -1., np.nan, np.inf, -np.inf, True):
                with self.subTest(name=name, value=value):
                    d = f.engineer(source(**{name: [value]}))
                    self.assertTrue(pd.isna(d[output].iloc[0])); self.assertFalse(f.primary_complete(d)[0])
                    self.assertEqual(d.age_linear.iloc[0], -1.)
                    if name != 'LBXGH': self.assertEqual(d.hba1c.iloc[0], 5.5)

    def test_nonfinite_derived_urine_ratio_is_missing_without_floor_or_clipping(self):
        d = f.engineer(source(URXUMA=[1e308], URXUCR=[1.]))
        self.assertTrue(pd.isna(d.uacr.iloc[0])); self.assertTrue(pd.isna(d.log_uacr.iloc[0]))
        self.assertTrue(d.routine_missing_uacr.iloc[0]); self.assertFalse(f.primary_complete(d)[0])

    def test_unsupported_missing_cycle_and_missing_source_schema_fail(self):
        for cycle in (1998., 2000., 2002., 2005., np.nan, np.inf, True):
            with self.subTest(cycle=cycle), self.assertRaisesRegex(ValueError, 'cycle mapping'):
                f.engineer(source(cycle=[cycle]))
        for name in f.REQUIRED:
            with self.subTest(missing=name), self.assertRaisesRegex(ValueError, 'Missing prespecified'):
                f.engineer(source().drop(columns=name))
        for value in (None, [], {}):
            with self.assertRaisesRegex(ValueError, 'DataFrame'): f.engineer(value)

    def test_index_raw_sources_and_unrelated_columns_are_preserved(self):
        raw = source(); raw.index = [317]; raw['unrelated_source'] = 'synthetic'
        before = raw.copy(deep=True); d = f.engineer(raw)
        pd.testing.assert_frame_equal(raw, before)
        pd.testing.assert_frame_equal(d[before.columns], before)
        self.assertEqual(d.index.tolist(), [317]); self.assertIsNot(d, raw)

    def test_core_has_exact_eleven_ids_and_no_protein_or_diabetes_duplicate(self):
        self.assertEqual(f.CORE, ('age_linear', 'age_quadratic', 'male', 'current', 'former',
                                 'sbp', 'bmi', 'hba1c', 'egfr_low', 'egfr_high', 'log_uacr'))
        for name in ('raw_creatinine', 'creatinine_calibrated', 'egfr', 'URXUMA', 'URXUCR',
                     'LBXB2M', 'cystatin_c', 'diabetes_history', 'cycle'):
            self.assertNotIn(name, f.CORE)
        raw = source(); raw['LBXB2M'] = 1.; raw['cystatin_c'] = 1.
        original = f.engineer(raw)[list(f.CORE)]
        raw['LBXB2M'] = 1e9; raw['cystatin_c'] = np.nan; raw['diabetes_history'] = 1.
        pd.testing.assert_frame_equal(original, f.engineer(raw)[list(f.CORE)])


if __name__ == '__main__':
    unittest.main()
