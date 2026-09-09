import json
import math
from pathlib import Path
import random
import unittest
from research.mortality.cohort_tools import fib4, observed_cif


class Fib4Tests(unittest.TestCase):
    def test_hand_calculation(self):
        self.assertAlmostEqual(fib4(50, 30, 25, 200)["fib4"], 1.5)

    def test_not_a_probability(self):
        result = fib4(50, 30, 25, 200)
        self.assertIsNone(result["mortality_probability"])
        self.assertEqual(result["status"], "numerical_index_only")

    def test_age_cautions(self):
        self.assertEqual(len(fib4(30, 30, 25, 200)["cautions"]), 2)
        self.assertEqual(len(fib4(70, 30, 25, 200)["cautions"]), 2)

    def test_invalid_inputs(self):
        for args in [(17,30,25,200),(50,0,25,200),(50,30,0,200),
                     (50,30,25,0),(True,30,25,200),(50,math.inf,25,200),
                     (50,30,math.nan,200),("50",30,25,200)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                fib4(*args)

    def test_acute_illness(self):
        for value in (True, 1, "false"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                fib4(50,30,25,200,acute_illness=value)


class CifTests(unittest.TestCase):
    def calc(self, times, events, h, causes=("target","other","unknown")):
        return observed_cif(times, events, h, causes=causes,
                            complete_cause_partition=True)

    def test_hand_calculation(self):
        r = self.calc([1,2,3,4],["target","other",None,"target"],4)
        self.assertAlmostEqual(r["cause_specific_cif"]["target"], .75)
        self.assertAlmostEqual(r["cause_specific_cif"]["other"], .25)
        self.assertAlmostEqual(r["survival"], 0)

    def test_ties_death_and_censor(self):
        r = self.calc([1,1,2,2],["target",None,"other",None],2)
        self.assertAlmostEqual(r["cause_specific_cif"]["target"], .25)
        self.assertAlmostEqual(r["cause_specific_cif"]["other"], .375)
        self.assertAlmostEqual(r["survival"], .375)

    def test_all_censored(self):
        r=self.calc([1,2],[None,None],2)
        self.assertEqual(r["survival"],1)
        self.assertEqual(r["all_cause_death_risk"],0)

    def test_zero_horizon(self):
        self.assertEqual(self.calc([1],["target"],0)["survival"],1)

    def test_unknown_is_event(self):
        r=self.calc([1],["unknown"],1)
        self.assertEqual(r["cause_specific_cif"]["unknown"],1)

    def test_partial_horizon(self):
        r=self.calc([1,2,3,4],["target","other",None,"target"],1.5)
        self.assertEqual(r["all_cause_death_risk"],.25)
        self.assertEqual(r["n_at_risk_just_before_horizon"],3)

    def test_partition_assertion_required(self):
        with self.assertRaises(ValueError):
            observed_cif([1],["target"],1,causes=["target"])

    def test_invalid_arrays(self):
        for times,events,h in [([],[],0),([1],[],1),([0],[None],0),
                              ([1],["undeclared"],1),([1],["target"],2),
                              ([1],["target"],-1),([math.nan],[None],0),
                              ([True],[None],0),([1],[1],1)]:
            with self.subTest(times=times,events=events,h=h), self.assertRaises(ValueError):
                self.calc(times,events,h)

    def test_invalid_cause_sets(self):
        for causes in ([],["target","target"],[""],"target"):
            with self.subTest(causes=causes), self.assertRaises(ValueError):
                self.calc([1],[None],1,causes)

    def test_permutation_and_probability_properties(self):
        rng=random.Random(42)
        for _ in range(100):
            records=[(rng.randint(1,10),rng.choice([None,"target","other","unknown"]))
                     for _ in range(30)]
            times,events=zip(*records)
            a=self.calc(times,events,max(times))
            rng.shuffle(records)
            times,events=zip(*records)
            b=self.calc(times,events,max(times))
            self.assertEqual(a["cause_specific_cif"],b["cause_specific_cif"])
            self.assertAlmostEqual(a["survival"]+sum(a["cause_specific_cif"].values()),1)
            previous={c:0 for c in a["cause_specific_cif"]}
            for row in a["event_time_curve"]:
                for cause,value in row["cause_specific_cif"].items():
                    self.assertGreaterEqual(value,previous[cause])
                    previous[cause]=value


class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[1]/"research/mortality/evidence_02.json"
        cls.data=json.loads(path.read_text(encoding="utf-8"))

    def test_study_ids_unique(self):
        ids=[s["id"] for s in self.data["studies"]]
        self.assertEqual(len(ids),len(set(ids)))
        self.assertEqual(len(ids),8)

    def test_intervals_and_endpoints(self):
        for study in self.data["studies"]:
            self.assertFalse(study["implemented_absolute_individual_risk"])
            for row in study["estimates"]:
                self.assertIn(row["endpoint_class"],("mortality","incidence","diagnosis"))
                ci=row.get("ci95")
                if ci is not None:
                    self.assertLessEqual(ci[0],row["value"])
                    self.assertLessEqual(row["value"],ci[1])

    def test_cdc_totals(self):
        rows=self.data["mortality_dataset"]["rows"]
        self.assertEqual(len(rows),12)
        for sex in ("all","male","female"):
            total=next(r["deaths"] for r in rows if r["sex"]==sex and r["age"]=="65+")
            self.assertEqual(total,sum(r["deaths"] for r in rows
                                      if r["sex"]==sex and r["age"]!="65+"))
        for age in ("65+","65-74","75-84","85+"):
            subset={r["sex"]:r["deaths"] for r in rows if r["age"]==age}
            self.assertEqual(subset["all"],subset["male"]+subset["female"])

    def test_rate_types_not_mixed(self):
        for row in self.data["mortality_dataset"]["rows"]:
            expected="age_standardized_2000_US" if row["age"]=="65+" else "age_specific"
            self.assertEqual(row["rate_type"],expected)


if __name__ == "__main__":
    unittest.main()
