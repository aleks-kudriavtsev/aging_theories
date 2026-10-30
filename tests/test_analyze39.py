import copy
import csv
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from research.mortality.period39 import analyze
from research.mortality.period39.plot import series
from research.mortality.period39.reference import calculate_reference


class ReplayContractTests(unittest.TestCase):
    def test_missing_freeze_blocks_before_source_load(self):
        with tempfile.TemporaryDirectory() as d, patch.object(analyze, "verify_freeze", side_effect=ValueError("no freeze")), patch.object(analyze, "load_bundle") as load:
            with self.assertRaisesRegex(ValueError, "no freeze"):
                analyze.run(Path(d)/"new")
            load.assert_not_called()
            self.assertFalse((Path(d)/"new").exists())

    def test_existing_output_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/"old").write_text("original")
            with self.assertRaisesRegex(ValueError, "no overwrite"):
                analyze.run(p)
            self.assertEqual((p/"old").read_text(), "original")

    def test_export_preserves_decimal_values_and_nonclinical_scope(self):
        # Manufactured constant-rate fixture, not a Russian probability.
        r = calculate_reference([.04]*5, [.039211]*5)
        row = dict(r, request=analyze.request_for(50, "F", 5), end_exact_birthday=55)
        rows = list(csv.DictReader(io.StringIO(analyze.csv_bytes([row]).decode())))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["CF_probability"], r["methods"]["constant_force"]["decimal"]["probability"])
        self.assertEqual(rows[0]["clinical_validation"], "False")
        self.assertEqual(rows[0]["current_calendar_forecast"], "False")
        self.assertEqual(rows[0]["output_scope"], "provisional_local_derived_period_reference")

    def test_plot_requires_exact_aggregate_grid_and_scope(self):
        rows = [{"request": analyze.request_for(a,s,h), "methods": {"constant_force": {"probability": .1}}}
                for s in ("F", "M") for a in range(40,81) for h in (5,10)]
        doc = dict(schema="period39-results/1", reference_rows=164, rows=rows,
                   clinical_validation=False, current_calendar_forecast=False,
                   source_year=2023, period_assumption="hold_2023_attained_age_schedule",
                   scope="provisional_local_derived_period_reference_under_declared_units")
        self.assertEqual(len(series(doc)),164)
        bad = copy.deepcopy(doc)
        bad["rows"].pop()
        with self.assertRaises(ValueError):
            series(bad)
        for field, value in (("country", "USA"), ("source_year", 2024), ("age_exact_years", 40.0),
                             ("period_assumption", "forecast"), ("purpose", "disease_free_baseline")):
            bad = copy.deepcopy(doc)
            bad["rows"][0]["request"][field] = value
            with self.assertRaises(ValueError):
                series(bad)
        bad = copy.deepcopy(doc)
        bad["source_year"] = 2024
        with self.assertRaises(ValueError):
            series(bad)
        bad = copy.deepcopy(doc)
        bad["clinical_validation"] = True
        with self.assertRaises(ValueError):
            series(bad)
        bad = copy.deepcopy(doc)
        bad["rows"][0]["methods"]["constant_force"]["probability"] = float("nan")
        with self.assertRaises(ValueError):
            series(bad)


if __name__ == "__main__":
    unittest.main()
