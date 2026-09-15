"""Exchange integrity and scientific-context checks; no patient-data tests."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from research.mortality import export_system_evidence as export
from research.mortality.workbench.runtime import engineer


class SystemEvidenceExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = export.build()

    def test_all_numbers_resolve_to_exact_frozen_source(self):
        bundle = self.bundle
        for row in bundle["performance_records"]:
            source = row["source"]
            raw = bundle["source_documents"][source["path"]]["raw_document"]
            self.assertEqual(sha256(raw.encode()).hexdigest(), source["sha256"])
            doc = json.loads(raw)
            self.assertEqual(row["value"], export.pointer(doc, source["json_pointer"]))
            if "interval" in row:
                for bound in ("lower", "upper"):
                    self.assertEqual(row["interval"][bound], export.pointer(doc,
                        row["interval"]["source_pointers"][bound]))

    def test_two_exact_source_documents_and_three_artifact_hashes(self):
        self.assertEqual(set(self.bundle["source_documents"]), set(export.SUMMARY_PATHS))
        for rel, expected in export.registry.ARTIFACTS.values():
            self.assertEqual(self.bundle["source_artifacts"]["research/mortality/"+rel]["sha256"], expected)
        self.assertEqual(self.bundle["source"]["commit"], export.SOURCE_COMMIT)

    def test_context_never_pools_horizons_or_state_partitions(self):
        rows = self.bundle["performance_records"]
        for row in rows:
            c = row["context"]
            self.assertEqual(c["country"], "US")
            self.assertEqual(c["age_range"], [40, 79])
            if row["source"]["path"] == export.COMPACT_PATH:
                self.assertEqual((c["horizon_years"], c["n"], c["deaths"]), (3, 2922, 93))
                self.assertIn(c["state_count"], (2, 4))
                self.assertEqual(c["censoring_method"], "IPCW_reverse_Kaplan_Meier")
            elif c["cohort_id"] == "NHANES2011_2012":
                self.assertEqual((c["horizon_years"], c["n"], c["deaths"]), (5, 2666, 130))
            else:
                self.assertEqual((c["horizon_years"], c["n"], c["deaths"]), (4, 3013, 113))
            if "multiclass" in row["metric"]:
                self.assertIn(c["state_count"], (4, 6))
                self.assertEqual(c["Brier_definition"], "sum_over_states_no_half_factor")
            else:
                self.assertEqual(c["state_count"], 2)
                self.assertEqual(c["endpoint"], "all_cause_mortality")

    def test_negative_and_cross_zero_intervals_not_clipped(self):
        contrasts = [r for r in self.bundle["performance_records"] if "interval" in r]
        self.assertTrue(any(r["interval"]["upper"] < 0 for r in contrasts))
        self.assertTrue(any(r["interval"]["lower"] < 0 < r["interval"]["upper"] for r in contrasts))
        for row in contrasts:
            self.assertEqual(row["interval"]["replicates"], 1000)
            self.assertEqual(row["component_attribution"], "model_or_block_only")

    def test_null_interval_is_retained_without_imputation(self):
        doc = {"value": -0.1, "lower": None, "upper": 0.2}
        r = export.record({"synthetic/test.json": json.dumps(doc).encode()}, "synthetic/test.json", "/value",
            "compact4", "delta_AUC", {}, interval_pointers=("/lower", "/upper"))
        self.assertIsNone(r["interval"]["lower"])
        self.assertEqual(r["interval"]["upper"], 0.2)

    def test_catalog_distinguishes_lab_predictors_from_complete_domain_inputs(self):
        models = {r["model_id"]: r for r in self.bundle["models"]}
        self.assertEqual(len(models), 5)
        clinical = models["clinical3"]
        self.assertEqual(clinical["laboratory_predictor_measurement_ids"], [])
        self.assertEqual(len(clinical["nonpredictive_domain_requirement_ids"]), 8)
        self.assertEqual(clinical["required_input_measurement_ids"], models["full8"]["required_input_measurement_ids"])
        compact = models["compact4"]
        self.assertEqual(compact["predictors"], 4)
        self.assertEqual(compact["determinations"], 5)
        mapping = compact["predictor_input_mappings"]
        self.assertTrue(mapping["creatinine"]["contains_age"])
        self.assertEqual(mapping["uacr"]["underlying_measurements"], ["urine_albumin", "urine_creatinine"])

    def test_no_coefficients_or_clinical_claims(self):
        b = self.bundle
        self.assertEqual(b["development_n"], 15280)
        self.assertFalse(b["safeguards"]["coefficients_exported"])
        self.assertFalse(b["external_validation"]["RU"])
        self.assertFalse(b["external_validation"]["DE"])
        self.assertTrue(b["external_validation"]["protocol_only"])
        for r in b["performance_records"]:
            self.assertFalse(r["clinical_use_ready"])
            self.assertFalse(r["usable_as_mortality_coefficient"])
        text = json.dumps(b)
        for prohibited in ('"coefficients"', '"baseline_hazard"', '"record_id"'):
            self.assertNotIn(prohibited, text)

    def test_mappings_point_to_frozen_features_with_correct_derivation(self):
        values = {"bmi": 26, "sbp": 125, "hba1c": 5.4, "creatinine": 0.9,
            "uacr": 8, "albumin": 43, "total_cholesterol": 195, "hdl": 55,
            "rdw": 12.8, "wbc": 6.5}
        clinical = {"smoking": "never", "bp_treatment": False, "diabetes_history": False,
            "cvd_history": False, "cancer_history": False}
        f = engineer(55, "female", clinical, values)
        older = engineer(70, "female", clinical, values)
        doubled = engineer(55, "female", clinical, dict(values, uacr=16))
        self.assertNotIn("creatinine", f)
        self.assertNotEqual((f["egfr_low"], f["egfr_high"]),
            (older["egfr_low"], older["egfr_high"]))
        self.assertAlmostEqual(doubled["log_uacr"] - f["log_uacr"], 1)
        for row in self.bundle["models"]:
            spec = export.registry.model_spec(row["model_id"])
            features = set(spec["preprocessing"]["features"])
            for mapping in row["predictor_input_mappings"].values():
                self.assertTrue(set(mapping["features"]) <= features)
                self.assertTrue(set(mapping["features"]) <= set(f))

    def test_wrong_commit_fails_before_export(self):
        for commit in ("latest", "1e5772ad", "f"*40):
            with self.assertRaises(ValueError):
                export.build(commit)

    def test_mutated_input_is_rejected(self):
        class FakeResult:
            stdout = b"different pinned bytes"
        with patch.object(export.subprocess, "run", return_value=FakeResult()):
            with self.assertRaisesRegex(ValueError, "differs from pinned"):
                export.build()

    def test_deterministic_stdlib_cli_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"exchange.json"
            cmd = [sys.executable, "-S", "-m", "research.mortality.export_system_evidence", "--out", str(path)]
            first = subprocess.run(cmd, cwd=export.REPO_ROOT, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(json.loads(path.read_text()), self.bundle)
            second = subprocess.run(cmd, cwd=export.REPO_ROOT, capture_output=True, text=True)
            self.assertNotEqual(second.returncode, 0)


if __name__ == "__main__":
    unittest.main()
