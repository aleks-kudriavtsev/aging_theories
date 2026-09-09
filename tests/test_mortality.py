import math
import unittest
from research.mortality.engine import cumulative_incidence, egfr_ckd_epi_2021


class EngineTests(unittest.TestCase):
    def calc(self, ends=(20,), rows=None, horizon=10, target=("ncd",)):
        rows = [{"ncd": .02, "other": .01}] if rows is None else rows
        return cumulative_incidence(ends, rows, horizon, target, complete_cause_partition=True)

    def test_closed_form(self):
        r = self.calc()
        self.assertAlmostEqual(r["all_cause_death_risk"], 1-math.exp(-.3))
        self.assertAlmostEqual(r["target_death_risk"], (1-math.exp(-.3))*2/3)

    def test_competing_death_reduces_target_probability(self):
        r = self.calc()
        self.assertLess(r["target_death_risk"], 1-math.exp(-.2))

    def test_zero_hazards(self):
        r = self.calc(rows=[{"ncd": 0, "other": 0}])
        self.assertEqual(r["survival"], 1)
        self.assertEqual(r["all_cause_death_risk"], 0)

    def test_zero_horizon(self):
        self.assertEqual(self.calc(horizon=0)["target_death_risk"], 0)

    def test_piecewise_partial_interval(self):
        r = self.calc((2, 10), [{"ncd": .02, "other": .01}, {"ncd": .04, "other": .01}], 5)
        expected = -math.expm1(-.06)*2/3 + math.exp(-.06)*-math.expm1(-.15)*4/5
        self.assertAlmostEqual(r["target_death_risk"], expected)
        self.assertAlmostEqual(r["survival"], math.exp(-.21))

    def test_probability_partition_and_monotonicity(self):
        last = 0
        for h in (1, 5, 10, 20):
            r = self.calc(horizon=h)
            self.assertAlmostEqual(r["survival"] + sum(r["cause_specific_cif"].values()), 1)
            self.assertGreater(r["target_death_risk"], last)
            last = r["target_death_risk"]

    def test_tiny_hazard_precision(self):
        r = self.calc(rows=[{"ncd": 1e-20, "other": 0}], horizon=1)
        self.assertEqual(r["all_cause_death_risk"], 1e-20)

    def test_high_hazard(self):
        r = self.calc(rows=[{"ncd": 1000, "other": 1000}])
        self.assertEqual(r["all_cause_death_risk"], 1)
        self.assertEqual(r["target_death_risk"], .5)

    def test_requires_explicit_partition(self):
        with self.assertRaises(ValueError):
            cumulative_incidence([10], [{"ncd": .1}], 5, ["ncd"])

    def test_no_extrapolation(self):
        with self.assertRaises(ValueError):
            self.calc(horizon=21)

    def test_invalid_intervals(self):
        for ends in ((0,), (-1,), (2, 1), (1, 1), (float("nan"),)):
            with self.subTest(ends=ends), self.assertRaises(ValueError):
                self.calc(ends=ends)

    def test_invalid_hazards(self):
        for x in (-.01, float("nan"), float("inf"), True, "0.1"):
            with self.subTest(x=x), self.assertRaises(ValueError):
                self.calc(rows=[{"ncd": x, "other": .01}])

    def test_inconsistent_partition(self):
        with self.assertRaises(ValueError):
            self.calc((5, 10), [{"ncd": .1}, {"other": .1}])

    def test_invalid_targets(self):
        for t in (("absent",), ("ncd", "ncd"), "ncd"):
            with self.subTest(t=t), self.assertRaises(ValueError):
                self.calc(target=t)

    def test_egfr_unit_conversion(self):
        for sex in ("female", "male"):
            self.assertAlmostEqual(egfr_ckd_epi_2021(50, sex, 1), egfr_ckd_epi_2021(50, sex, 88.4, "umol/L"))

    def test_egfr_knot_reference(self):
        self.assertAlmostEqual(egfr_ckd_epi_2021(50, "male", .9), 142*.9938**50)
        self.assertAlmostEqual(egfr_ckd_epi_2021(50, "female", .7), 142*.9938**50*1.012)

    def test_egfr_monotonic(self):
        self.assertGreater(egfr_ckd_epi_2021(40, "male", 1), egfr_ckd_epi_2021(60, "male", 1))
        self.assertGreater(egfr_ckd_epi_2021(50, "female", .6), egfr_ckd_epi_2021(50, "female", 1))

    def test_egfr_rejects_invalid_inputs(self):
        for args in ((17, "male", 1), (50, "unknown", 1), (50, "male", 0), (50, "male", 1, "mmol/L"), (float("nan"), "male", 1), (50, "female", float("inf"))):
            with self.subTest(args=args), self.assertRaises(ValueError):
                egfr_ckd_epi_2021(*args)


if __name__ == "__main__":
    unittest.main()
