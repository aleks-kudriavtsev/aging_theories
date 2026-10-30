"""Provisional local-derived Russian2023 age/sex PERIOD mortality reference.

The pinned mx values are annual central rates according to the local extractor
and documentation. Original UN RDA numeric scaling has NOT been reverified.
Constant-force and UDD(a_x=.5) are explicitly separate actuarial assumptions.
This module never estimates biomarkers, a disease-free hazard, future calendar
rates, a personalised forecast, clinical calibration or statistical confidence.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import io
import json
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

HERE = Path(__file__).resolve().parent
PROTOCOL_COMMIT = "c0a570e0ac40de94da8a05b53f5e9d80aaf7aacf"
PLAN_SHA256 = "3c93c327c56a86a2015ad2d1f0e13d2041246fb8832655eccbb15d5ace864731"
MANIFEST_SHA256 = "13569348125ba78e44ec08168541c1d45825cd1a00d06518845c1b2182e11680"
CSV_PATH = "data/mortality/mx_qx_usa_deu_rus.csv"
FIELDS = ("country", "sex", "year", "age", "mx", "qx", "series")
KEYS = frozenset((sex, age) for sex in ("F", "M") for age in range(40, 90))
DIGITS = 85
HALFWIDTH = Decimal("0.0000005")
QX_TOLERANCE = Decimal("0.000001")
REQUEST_FIELDS = {
    "schema", "model", "country", "source_year", "source_sex", "age_exact_years",
    "horizon_years", "age_origin", "period_assumption", "purpose",
}


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _decimal(value, name):
    if type(value) not in (int, float, Decimal, str):
        raise ValueError("Explicit finite numeric value required: " + name)
    try:
        number = value if type(value) is Decimal else Decimal(str(value))
    except (InvalidOperation, ValueError) as error:
        raise ValueError("Invalid decimal number: " + name) from error
    if not number.is_finite():
        raise ValueError("Finite decimal number required: " + name)
    if number.as_tuple().exponent < -6:
        raise ValueError("At most six decimal places allowed in saved values: " + name)
    return number


def _rate(value, name, *, positive):
    number = _decimal(value, name)
    if number < 0 or (positive and number == 0):
        raise ValueError("Positive saved mx required" if positive else "Nonnegative mx required")
    if number > 2 or number + HALFWIDTH > 2:
        raise ValueError("UDD upper rounding rate exceeds2; would yield q>1; no clipping")
    return number


def _probability(value, name):
    number = _decimal(value, name)
    if not 0 <= number <= 1:
        raise ValueError("Stored qx must be in[0,1]")
    return number


def _check_qx(rate, qx):
    with localcontext() as context:
        context.prec = DIGITS
        if abs(qx - (1 - (-rate).exp())) > QX_TOLERANCE:
            raise ValueError("Stored author-qx differs from recomputed constant-force qx by more than1e-6")


@dataclass(frozen=True)
class ReferenceBundle:
    """Immutable100-cell input bundle; unverified fixtures remain labelled so.

    rates and stored_qx map (source_sex, attained_integer_age) to Decimal.
    The source hash/commit fields describe only already-checked local bytes;
    they never authenticate official upstream scaling or a clinical forecast.
    """
    rates: Mapping[tuple[str, int], Decimal]
    stored_qx: Mapping[tuple[str, int], Decimal]
    source_hashes: Mapping[str, str] = field(default_factory=dict)
    source_commit: str | None = None
    source_files_verified: bool = False
    source_unit_evidence: str = "unverified_or_manufactured_fixture_not_official_source"

    def __post_init__(self):
        rates, qxs = dict(self.rates), dict(self.stored_qx)
        if any(type(k) is not tuple or len(k) != 2 or type(k[0]) is not str or type(k[1]) is not int
               for k in [*rates, *qxs]):
            raise ValueError("Source keys must be exact(source_sex string, integer age) tuples")
        if set(rates) != KEYS or set(qxs) != KEYS:
            raise ValueError("Exactly100 source-sex/age keys F/M40..89 required")
        for key in KEYS:
            rates[key] = _rate(rates[key], "mx", positive=True)
            qxs[key] = _probability(qxs[key], "stored qx")
            _check_qx(rates[key], qxs[key])
        hashes = dict(self.source_hashes)
        if any(not isinstance(k, str) or not k or not isinstance(v, str) or not re.fullmatch(r"[0-9a-f]{64}", v) for k, v in hashes.items()):
            raise ValueError("Invalid source hash inventory")
        if type(self.source_files_verified) is not bool or not isinstance(self.source_unit_evidence, str) or not self.source_unit_evidence:
            raise ValueError("Explicit source-verification metadata required")
        if self.source_commit is not None and (not isinstance(self.source_commit, str) or not self.source_commit):
            raise ValueError("Invalid source commit metadata")
        object.__setattr__(self, "rates", MappingProxyType(rates))
        object.__setattr__(self, "stored_qx", MappingProxyType(qxs))
        object.__setattr__(self, "source_hashes", MappingProxyType(hashes))


def _integer(text, name):
    if not isinstance(text, str) or not re.fullmatch(r"0|[1-9][0-9]*", text):
        raise ValueError("Exact nonnegative integer CSV value required: " + name)
    return int(text)


def parse_rates(raw):
    """Pure target reader: whole CSV accepted, only RUS2023F/M40..89 retained.

    This parser does not verify a source checksum. Production load_bundle
    verifies ALL FOUR pinned source files before invoking it. No other country,
    year or age is used to create the100 target cells, and wrong series/sex in
    the target domain raises rather than silently shrinking that population.
    """
    if type(raw) is not bytes:
        raise ValueError("CSV bytes required")
    rates, qxs = {}, {}
    try:
        reader = csv.reader(io.StringIO(raw.decode("utf-8"), newline=""), strict=True)
        if tuple(next(reader, ())) != FIELDS:
            raise ValueError("Exact saved-rate CSV schema/order required")
        for row in reader:
            if len(row) != len(FIELDS):
                raise ValueError("Malformed saved-rate CSV row width")
            country, sex, year_text, age_text, mx_text, qx_text, series = row
            year, age = _integer(year_text, "year"), _integer(age_text, "age")
            if country != "RUS" or year != 2023 or not 40 <= age <= 89:
                continue
            if sex not in ("F", "M") or series != "estimate":
                raise ValueError("Target rows require source sexesF/M and estimate series")
            key = (sex, age)
            if key in rates:
                raise ValueError("Duplicate target source-sex/age key")
            rate = _rate(mx_text, "mx", positive=True)
            qx = _probability(qx_text, "qx")
            _check_qx(rate, qx)
            rates[key], qxs[key] = rate, qx
    except (UnicodeError, csv.Error) as error:
        raise ValueError("Invalid UTF-8 saved-rate CSV") from error
    return ReferenceBundle(rates, qxs)


def _strict_json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate frozen JSON field")
            result[key] = value
        return result
    return json.loads(raw, object_pairs_hook=pairs,
        parse_constant=lambda _: (_ for _ in ()).throw(ValueError("Nonfinite frozen JSON")))


def _load_metadata():
    plan_raw = (HERE / "STUDY_PLAN.json").read_bytes()
    manifest_raw = (HERE / "SOURCE_MANIFEST.json").read_bytes()
    if _sha(plan_raw) != PLAN_SHA256 or _sha(manifest_raw) != MANIFEST_SHA256:
        raise ValueError("Frozen period39 plan or source manifest changed")
    plan, manifest = _strict_json(plan_raw), _strict_json(manifest_raw)
    if (plan.get("schema") != "period39-study-plan/1"
            or manifest.get("schema") != "period39-source-manifest/1"
            or plan["source"]["manifest_sha256"] != MANIFEST_SHA256
            or plan["source"]["csv"] != CSV_PATH
            or manifest["official_raw_RDA_sha256"] is not None
            or manifest["qx_is_author_derived_not_official"] is not True
            or manifest["new_raw_source_downloads"] is not False):
        raise ValueError("Frozen source/provenance scope invalid")
    return plan, manifest


def _verify_sources(root, entries):
    if type(entries) is not list or len(entries) != 4:
        raise ValueError("Exactly four source-file pins required")
    root = Path(root)
    saved = {}
    for item in entries:
        if type(item) is not dict or set(item) != {"path", "sha256"}:
            raise ValueError("Source hash specification invalid")
        name, digest = item["path"], item["sha256"]
        if (not isinstance(name, str) or not name or Path(name).is_absolute()
                or ".." in Path(name).parts or name in saved
                or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest)):
            raise ValueError("Invalid or duplicate source path/hash")
        path = root / name
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root.resolve()):
            raise ValueError("Missing, linked or escaped source file")
        raw = path.read_bytes()
        if _sha(raw) != digest:
            raise ValueError("Source checksum mismatch before CSV parsing: " + name)
        saved[name] = raw
    return saved


def load_bundle():
    """Verify four local source files before parsing; never retrieves anything."""
    plan, manifest = _load_metadata()
    saved = _verify_sources(manifest["source_root"], manifest["files"])
    bundle = parse_rates(saved[plan["source"]["csv"]])
    return replace(bundle,
        source_hashes={item["path"]: item["sha256"] for item in manifest["files"]},
        source_commit=manifest["source_commit"], source_files_verified=True,
        source_unit_evidence=manifest["verified_unit_level"])


def _constant_force(rates):
    return 1 - (-sum(rates, Decimal(0))).exp()


def _udd(rates):
    survival = Decimal(1)
    for rate in rates:
        qx = rate / (1 + rate / 2)
        if not 0 <= qx <= 1:
            raise ValueError("UDD q outside[0,1]; no clipping")
        survival *= 1 - qx
    return 1 - survival


def _method(point, lower, upper):
    if not 0 <= lower <= point <= upper <= 1:
        raise ArithmeticError("Reference probability or rounding bounds invalid; no normalization")
    return {"probability": float(point),
            "rounding_bounds": {"lower": float(lower), "upper": float(upper)},
            "decimal": {"probability": str(point), "lower": str(lower), "upper": str(upper)}}


def calculate_reference(mx_values, storedqx_values):
    """Pure finite-array calculation under two predeclared actuarial assumptions.

    Six-decimal input mx is assumed rounded to nearest6decimal places. Bounds
    propagate ±0.5e-6 in EACH rate, lower-bounded by0, not statistical uncertainty.
    The author-derived stored qx product is a distinct arithmetic diagnostic.
    It is never called official UN qx. Arbitrary positive horizons are accepted
    here for mathematical fixtures; the strict public adapter accepts5/10only.
    """
    if type(mx_values) not in (list, tuple) or type(storedqx_values) not in (list, tuple) or not mx_values or len(mx_values) != len(storedqx_values):
        raise ValueError("Equal nonempty explicit mx/stored-qx sequences required")
    rates = tuple(_rate(x, "mx", positive=False) for x in mx_values)
    qxs = tuple(_probability(x, "stored qx") for x in storedqx_values)
    for rate, qx in zip(rates, qxs):
        _check_qx(rate, qx)
    with localcontext() as context:
        context.prec = DIGITS
        lower_rates = tuple(max(Decimal(0), rate - HALFWIDTH) for rate in rates)
        upper_rates = tuple(rate + HALFWIDTH for rate in rates)
        cf = _constant_force(rates)
        udd = _udd(rates)
        stored_survival = Decimal(1)
        for qx in qxs:
            stored_survival *= 1 - qx
        stored = 1 - stored_survival
        return {
            "methods": {
                "constant_force": _method(cf, _constant_force(lower_rates), _constant_force(upper_rates)),
                "UDD_ax_0_5": _method(udd, _udd(lower_rates), _udd(upper_rates)),
            },
            "diagnostic": {"stored_qx_product": float(stored), "delta_vs_constant_force": float(stored - cf),
                "decimal": {"stored_qx_product": str(stored), "delta_vs_constant_force": str(stored - cf)}},
            "precision": {"decimal_digits": DIGITS, "mx_rounding_halfwidth": "0.0000005"},
            "years_count": len(rates),
        }


def score_reference(request, bundle):
    """Strict exact-birthday population PERIOD adapter, never individual scoring."""
    if type(request) is not dict or set(request) != REQUEST_FIELDS:
        raise ValueError("Exact population-reference request schema required; no biomarkers/extra fields")
    constants = {
        "schema": "period39-request/1", "model": "RUS_WPP2023_period_reference", "country": "RUS",
        "age_origin": "exact_birthday", "period_assumption": "hold_2023_attained_age_schedule",
        "purpose": "population_period_reference_only",
    }
    if any(request[k] != v for k, v in constants.items()):
        raise ValueError("Population-only period/clock/country/model/purpose contract mismatch")
    if type(request["source_year"]) is not int or request["source_year"] != 2023:
        raise ValueError("Only the pinned source year2023 is supported")
    if type(request["source_sex"]) is not str or request["source_sex"] not in ("F", "M"):
        raise ValueError("Explicit source sexF/M required; no combined-sex averaging")
    age, horizon = request["age_exact_years"], request["horizon_years"]
    if type(age) is not int or not 40 <= age <= 80:
        raise ValueError("Exact integer birthday age40..80 required; no interpolation")
    if type(horizon) is not int or horizon not in (5, 10):
        raise ValueError("Explicit supported integer horizon5/10years required")
    if not isinstance(bundle, ReferenceBundle):
        raise ValueError("Immutable validated ReferenceBundle required")
    if bundle.source_files_verified is not True:
        raise ValueError("Named Russian period reference requires verified pinned local source files")
    # A validated numeric fixture is not automatically the named Russian source.
    # Match the immutable bundle's source identity to the frozen metadata too;
    # this reads only protocol/manifest bytes, not source CSV or new data.
    _, manifest = _load_metadata()
    expected_hashes = {item["path"]: item["sha256"] for item in manifest["files"]}
    if (len(expected_hashes) != 4 or len(manifest["files"]) != 4
            or bundle.source_commit != manifest["source_commit"]
            or dict(bundle.source_hashes) != expected_hashes
            or bundle.source_unit_evidence != manifest["verified_unit_level"]):
        raise ValueError("Bundle source version/hash/unit-evidence differs from the frozen manifest")
    ages = list(range(age, age + horizon))
    keys = [(request["source_sex"], a) for a in ages]
    if any(key not in bundle.rates or key not in bundle.stored_qx for key in keys):
        raise ValueError("Source attained-age support missing; no extrapolation")
    result = calculate_reference([bundle.rates[k] for k in keys], [bundle.stored_qx[k] for k in keys])
    result.update(schema="period39-reference-result/1", request=dict(request), source_ages=ages,
        end_exact_birthday=age + horizon,
        metadata={
            "research_only": True, "reference_status": "provisional_local_derived_population_period_reference",
            "clinical_forecast": False, "clinical_validation": False, "individual_calibration_performed": False,
            "source_population_includes_diseased": True, "disease_free_or_stage_specific_baseline": False,
            "forecast_of_future_calendar_rates": False, "future_schedule_held_at_source2023": True,
            "source_rates_declared_unit": "annual_deaths_per_person_year",
            "source_unit_evidence": bundle.source_unit_evidence,
            "official_raw_RDA_numeric_scaling_reverified": False,
            "official_raw_RDA_checksum_verified": False,
            "source_files_verified": bundle.source_files_verified,
            "source_commit": bundle.source_commit, "source_hashes": dict(bundle.source_hashes),
            "stored_qx_is_author_derived_not_official_UN_qx": True,
            "rounding_bounds_are_confidence_intervals": False,
            "rounding_assumption": "each_mx_rounded_to_nearest6decimals_with_independent_closed_halfwidth_0.5e-6",
            "assumption_difference_is_model_selection": False,
            "interpolation_or_extrapolation_performed": False, "sexes_averaged": False,
            "new_outcomes_or_participant_records_read": False,
        })
    return result
