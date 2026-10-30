"""Manufactured fixtures only; no real164-cell period39 production is run."""
import copy
import csv
from dataclasses import FrozenInstanceError, replace
from decimal import Decimal, localcontext
import hashlib
import io
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from research.mortality.period39 import reference as ref


def qx_text(rate):
    with localcontext() as context:
        context.prec = 90
        return str((1 - (-Decimal(str(rate))).exp()).quantize(Decimal("0.000001")))


def fake_rows():
    return [["RUS", sex, "2023", str(age), "0.010000", qx_text("0.01"), "estimate"]
            for sex in ("F", "M") for age in range(40, 90)]


def csv_bytes(rows=None, header=ref.FIELDS):
    out = io.StringIO(newline="")
    writer = csv.writer(out)
    writer.writerow(header)
    writer.writerows(fake_rows() if rows is None else rows)
    return out.getvalue().encode()


def request(**kwargs):
    result = dict(schema="period39-request/1", model="RUS_WPP2023_period_reference", country="RUS",
        source_year=2023, source_sex="F", age_exact_years=40, horizon_years=5,
        age_origin="exact_birthday", period_assumption="hold_2023_attained_age_schedule",
        purpose="population_period_reference_only")
    result.update(kwargs)
    return result


def verified_fixture(bundle=None):
    """Manufactured rates with a mocked loader identity, never Russian evidence.

    Only tests call this helper. The matching metadata is explicitly patched
    through each positive adapter call; production never accepts these pins.
    """
    names = [ref.CSV_PATH, "fake_context_a", "fake_context_b", "fake_context_c"]
    hashes = {name: hashlib.sha256(name.encode()).hexdigest() for name in names}
    manifest = {"files": [{"path": name, "sha256": value} for name, value in hashes.items()],
                "source_commit": "manufactured-mocked-loader-version",
                "verified_unit_level": "manufactured test metadata; not upstream scaling"}
    bundle = ref.parse_rates(csv_bytes()) if bundle is None else bundle
    return replace(bundle, source_files_verified=True, source_hashes=hashes,
                   source_commit=manifest["source_commit"], source_unit_evidence=manifest["verified_unit_level"]), manifest


def score_fixture(r, bundle=None):
    """Exercise the strict adapter only with explicitly mocked fake provenance."""
    verified, manifest = verified_fixture(bundle)
    with patch.object(ref, "_load_metadata", return_value=({}, manifest)):
        return ref.score_reference(r, verified)


class PeriodReferenceTests(unittest.TestCase):
    def test_constant_rate_hand_equations(self):
        result = ref.calculate_reference([Decimal("0.01")] * 5, [Decimal(qx_text(".01"))] * 5)
        self.assertAlmostEqual(result["methods"]["constant_force"]["probability"], -math.expm1(-.05), places=15)
        expected = 1 - ((1 - .01 / (1 + .5 * .01)) ** 5)
        self.assertAlmostEqual(result["methods"]["UDD_ax_0_5"]["probability"], expected, places=15)
        self.assertEqual(result["years_count"], 5)
        self.assertEqual(result["precision"]["decimal_digits"], 85)

    def test_rounding_bounds_direct_equations_and_not_ci(self):
        result = ref.calculate_reference([".01"] * 2, [qx_text(".01")] * 2)
        for method in result["methods"].values():
            self.assertLess(method["rounding_bounds"]["lower"], method["probability"])
            self.assertLess(method["probability"], method["rounding_bounds"]["upper"])
        with localcontext() as context:
            context.prec = 95
            lo, hi = Decimal(".0099995"), Decimal(".0100005")
            expected_lo = 1 - (-2 * lo).exp()
            expected_hi = 1 - (-2 * hi).exp()
        got = result["methods"]["constant_force"]["decimal"]
        self.assertLess(abs(Decimal(got["lower"]) - expected_lo), Decimal("1e-83"))
        self.assertLess(abs(Decimal(got["upper"]) - expected_hi), Decimal("1e-83"))
        scored = score_fixture(request())
        self.assertFalse(scored["metadata"]["rounding_bounds_are_confidence_intervals"])

    def test_stored_qx_product_is_distinct_diagnostic(self):
        result = ref.calculate_reference([".01"] * 5, [qx_text(".01")] * 5)
        expected = 1 - (1 - float(qx_text(".01"))) ** 5
        self.assertAlmostEqual(result["diagnostic"]["stored_qx_product"], expected, places=15)
        self.assertNotEqual(result["diagnostic"]["delta_vs_constant_force"], 0)
        self.assertAlmostEqual(result["diagnostic"]["delta_vs_constant_force"],
            result["diagnostic"]["stored_qx_product"] - result["methods"]["constant_force"]["probability"], places=16)

    def test_rounding_lower_closes_at_zero(self):
        result = ref.calculate_reference([0], [0])
        for method in result["methods"].values():
            self.assertEqual(method["probability"], 0)
            self.assertEqual(method["rounding_bounds"]["lower"], 0)
            self.assertGreater(method["rounding_bounds"]["upper"], 0)

    def test_udd_upper_endpoint_two_and_reject_above(self):
        # m=1.9999995 is a pure arithmetic boundary fixture, not a saved6dec mx.
        # The actual rounded six-decimal boundary below2 is1.999999.
        got = ref.calculate_reference(["1.999999"], [qx_text("1.999999")])
        self.assertLess(got["methods"]["UDD_ax_0_5"]["rounding_bounds"]["upper"], 1)
        self.assertEqual(ref._udd([Decimal(2)]), 1)
        with self.assertRaises(ValueError):
            ref.calculate_reference([2], [qx_text(2)])

    def test_nonfinite_boolean_negative_and_precision_rejected(self):
        for bad in (True, None, math.nan, math.inf, Decimal("NaN"), "Infinity", "-.01", ".0100001", "1e1000000", object()):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                ref.calculate_reference([bad], [0])
        for bad in (True, math.nan, "-.01", "1.1", ".0100001"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                ref.calculate_reference([".01"], [bad])
        for rates, qs in (([], []), ([".01"], []), ((".01",), "0.00995")):
            with self.subTest(rates=rates, qs=qs), self.assertRaises(ValueError):
                ref.calculate_reference(rates, qs)

    def test_pure_parser_exact_target_shape_and_immutable_maps(self):
        bundle = ref.parse_rates(csv_bytes())
        self.assertEqual(len(bundle.rates), 100)
        self.assertEqual(bundle.rates["F", 40], Decimal(".01"))
        self.assertFalse(bundle.source_files_verified)
        with self.assertRaises(TypeError):
            bundle.rates["F", 40] = Decimal(".02")
        with self.assertRaises(FrozenInstanceError):
            bundle.source_files_verified = True
        rates, qxs = dict(bundle.rates), dict(bundle.stored_qx)
        clone = ref.ReferenceBundle(rates, qxs)
        rates["F", 40] = Decimal(".02")
        self.assertEqual(clone.rates["F", 40], Decimal(".01"))

    def test_non_target_countries_years_ages_not_used(self):
        rows = fake_rows()
        rows += [["USA", "B", "2023", "40", "bad", "bad", "projection"],
                 ["RUS", "F", "2022", "40", "bad", "bad", "estimate"],
                 ["RUS", "F", "2023", "90", "bad", "bad", "estimate"]]
        bundle = ref.parse_rates(csv_bytes(rows))
        self.assertEqual(len(bundle.rates), 100)
        self.assertNotIn(("F", 90), bundle.rates)

    def test_direct_bundle_rejects_float_age_keys(self):
        bundle = ref.parse_rates(csv_bytes())
        rates = dict(bundle.rates)
        rates["F", 40.0] = rates.pop(("F", 40))
        with self.assertRaisesRegex(ValueError, "integer age"):
            ref.ReferenceBundle(rates, bundle.stored_qx)

    def test_parser_missing_duplicate_unknown_sex_and_series_rejected(self):
        for mutation in ("missing", "duplicate", "sex", "series"):
            rows = fake_rows()
            if mutation == "missing":
                rows.pop()
            elif mutation == "duplicate":
                rows.append(rows[0])
            elif mutation == "sex":
                rows[0][1] = "B"
            else:
                rows[0][-1] = "projection"
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                ref.parse_rates(csv_bytes(rows))

    def test_parser_exact_header_width_and_integer_coordinates(self):
        with self.assertRaises(ValueError):
            ref.parse_rates(csv_bytes(header=tuple(reversed(ref.FIELDS))))
        for column, value in ((2, "2023.0"), (3, "40.0"), (3, "-1")):
            rows = fake_rows()
            rows[0][column] = value
            with self.subTest(column=column, value=value), self.assertRaises(ValueError):
                ref.parse_rates(csv_bytes(rows))
        rows = fake_rows()
        rows[0].append("extra")
        with self.assertRaises(ValueError):
            ref.parse_rates(csv_bytes(rows))

    def test_parser_qx_consistency_and_positive_six_decimal_rates(self):
        for column, value in ((4, "0"), (4, "NaN"), (4, "0.0100001"), (5, ".1"), (5, "1.1"), (5, "-0.1")):
            rows = fake_rows()
            rows[0][column] = value
            with self.subTest(column=column, value=value), self.assertRaises(ValueError):
                ref.parse_rates(csv_bytes(rows))

    def test_age80_h10_uses_80_through89_not90(self):
        bundle = ref.parse_rates(csv_bytes())
        scored = score_fixture(request(age_exact_years=80, horizon_years=10), bundle)
        self.assertEqual(scored["source_ages"], list(range(80, 90)))
        self.assertEqual(scored["end_exact_birthday"], 90)
        self.assertEqual(scored["years_count"], 10)

    def test_period_reference_chains_under_both_assumptions(self):
        bundle = ref.parse_rates(csv_bytes())
        a5 = score_fixture(request(), bundle)
        b5 = score_fixture(request(age_exact_years=45), bundle)
        a10 = score_fixture(request(horizon_years=10), bundle)
        for method in ("constant_force", "UDD_ax_0_5"):
            p, q = a5["methods"][method]["probability"], b5["methods"][method]["probability"]
            self.assertAlmostEqual(a10["methods"][method]["probability"], p + (1-p)*q, places=15)

    def test_strict_request_country_year_sex_age_horizon(self):
        bundle = ref.parse_rates(csv_bytes())
        changes = [dict(country="USA"), dict(source_year=2024), dict(source_year=2023.0),
                   dict(source_sex="B"), dict(source_sex=None), dict(age_exact_years=39),
                   dict(age_exact_years=81), dict(age_exact_years=40.5), dict(age_exact_years=True),
                   dict(horizon_years=1), dict(horizon_years=5.0), dict(horizon_years=True)]
        for change in changes:
            with self.subTest(change=change), self.assertRaises(ValueError):
                score_fixture(request(**change), bundle)

    def test_clocks_purpose_and_unknown_fields_not_silently_accepted(self):
        bundle = ref.parse_rates(csv_bytes())
        for change in (dict(age_origin="completed_age"), dict(period_assumption="calendar_forecast"),
                       dict(purpose="disease_free_baseline"), dict(model="clinical_personal_model"),
                       dict(biomarker_B2M=2), dict(disease_free=True), dict(horizon_months=60)):
            with self.subTest(change=change), self.assertRaises(ValueError):
                score_fixture(request(**change), bundle)
        r = request()
        r.pop("period_assumption")
        with self.assertRaises(ValueError):
            score_fixture(r, bundle)

    def test_metadata_retains_unit_and_clinical_limitations(self):
        scored = score_fixture(request())
        m = scored["metadata"]
        self.assertFalse(m["clinical_forecast"])
        self.assertFalse(m["clinical_validation"])
        self.assertFalse(m["official_raw_RDA_numeric_scaling_reverified"])
        self.assertFalse(m["official_raw_RDA_checksum_verified"])
        self.assertFalse(m["forecast_of_future_calendar_rates"])
        self.assertTrue(m["source_files_verified"])
        self.assertEqual(m["source_unit_evidence"], "manufactured test metadata; not upstream scaling")
        self.assertTrue(m["stored_qx_is_author_derived_not_official_UN_qx"])
        self.assertTrue(m["source_population_includes_diseased"])

    def test_unverified_numeric_bundle_cannot_produce_named_russian_score(self):
        bundle = ref.parse_rates(csv_bytes())
        with patch.object(ref, "_load_metadata") as metadata, patch.object(ref, "calculate_reference") as calculation:
            with self.assertRaisesRegex(ValueError, "verified pinned local source"):
                ref.score_reference(request(), bundle)
            calculation.assert_not_called()
            metadata.assert_not_called()

    def test_verified_flag_does_not_replace_pinned_provenance(self):
        bundle, manifest = verified_fixture()
        wrong = dict(bundle.source_hashes)
        first = next(iter(wrong))
        wrong[first] = "0" * 64
        alterations = [
            {"source_hashes": {}}, {"source_hashes": wrong},
            {"source_commit": "different-source-version"}, {"source_commit": None},
            {"source_unit_evidence": "different-or-upstream-verified-unit-claim"},
        ]
        for change in alterations:
            bad = replace(bundle, **change)
            with self.subTest(change=change), patch.object(ref, "_load_metadata", return_value=({}, manifest)), patch.object(ref, "calculate_reference") as calculation:
                with self.assertRaisesRegex(ValueError, "frozen manifest"):
                    ref.score_reference(request(), bad)
                calculation.assert_not_called()

    def source_fixture(self, directory):
        names = [ref.CSV_PATH, "context_a.txt", "context_b.txt", "context_c.txt"]
        entries = []
        for name in names:
            p = directory / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(csv_bytes() if name == ref.CSV_PATH else name.encode())
            entries.append({"path": name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
        manifest = {"source_root": str(directory), "files": entries, "source_commit": "manufactured-source",
                    "verified_unit_level": "manufactured fixture docs, not upstream scaling"}
        plan = {"source": {"csv": ref.CSV_PATH}}
        return plan, manifest

    def test_all_four_source_hashes_verified_before_any_csv_parse(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            plan, manifest = self.source_fixture(root)
            for item in manifest["files"]:
                p = root / item["path"]
                original = p.read_bytes()
                p.write_bytes(original + b"changed")
                with self.subTest(path=item["path"]), patch.object(ref, "_load_metadata", return_value=(plan, manifest)), patch.object(ref, "parse_rates") as parser:
                    with self.assertRaisesRegex(ValueError, "before CSV parsing"):
                        ref.load_bundle()
                    parser.assert_not_called()
                p.write_bytes(original)

    def test_successful_fake_load_is_immutable_and_source_pinned(self):
        with tempfile.TemporaryDirectory() as d:
            plan, manifest = self.source_fixture(Path(d))
            with patch.object(ref, "_load_metadata", return_value=(plan, manifest)):
                bundle = ref.load_bundle()
            self.assertTrue(bundle.source_files_verified)
            self.assertEqual(len(bundle.source_hashes), 4)
            with self.assertRaises(TypeError):
                bundle.source_hashes[ref.CSV_PATH] = "x" * 64

    def test_source_paths_missing_links_duplicates_and_counts_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _, manifest = self.source_fixture(root)
            for bad in (manifest["files"][:3], manifest["files"] + manifest["files"][:1]):
                with self.assertRaises(ValueError):
                    ref._verify_sources(root, bad)
            for name in ("../outside", "/tmp/outside"):
                entries = copy.deepcopy(manifest["files"])
                entries[0]["path"] = name
                with self.assertRaises(ValueError):
                    ref._verify_sources(root, entries)
            entries = copy.deepcopy(manifest["files"])
            entries[-1] = entries[0]
            with self.assertRaises(ValueError):
                ref._verify_sources(root, entries)

    def test_missing_and_linked_source_rejected_before_parse(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            plan, manifest = self.source_fixture(root)
            last = root / manifest["files"][-1]["path"]
            original = last.read_bytes()
            last.unlink()
            with patch.object(ref, "_load_metadata", return_value=(plan, manifest)), patch.object(ref, "parse_rates") as parser:
                with self.assertRaisesRegex(ValueError, "Missing, linked"):
                    ref.load_bundle()
                parser.assert_not_called()
            target = root / "target.txt"
            target.write_bytes(original)
            last.symlink_to(target)
            with patch.object(ref, "_load_metadata", return_value=(plan, manifest)), patch.object(ref, "parse_rates") as parser:
                with self.assertRaisesRegex(ValueError, "Missing, linked"):
                    ref.load_bundle()
                parser.assert_not_called()


if __name__ == "__main__":
    unittest.main()
