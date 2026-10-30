"""Frozen 164-cell local-derived Russian period reference; no model fitting."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import platform

from .reference import load_bundle, score_reference

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROTOCOL = "c0a570e0ac40de94da8a05b53f5e9d80aaf7aacf"
PLAN_SHA256 = "3c93c327c56a86a2015ad2d1f0e13d2041246fb8832655eccbb15d5ace864731"
REQUIRED_FILES = {"research/mortality/period39/"+n for n in (
    "PROTOCOL_RU.md", "STUDY_PLAN.json", "SOURCE_MANIFEST.json", "METHOD_REVIEW.md",
    "SHFM_REVIEW.md", "SHFM_RIGHTS.md", "REVIEW_MANIFEST.json", "IMPLEMENTATION_NOTES_RU.md",
    "__init__.py", "reference.py", "analyze.py", "plot.py",
)} | {"tests/test_period39.py", "tests/test_analyze39.py"}
METHODS = ("constant_force", "UDD_ax_0_5")


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def verify_freeze():
    raw = (HERE/"CODE_MANIFEST.json").read_bytes()
    manifest = json.loads(raw)
    if (set(manifest) != {"schema", "protocol_commit", "files"}
            or manifest["schema"] != "period39-code-manifest/1" or manifest["protocol_commit"] != PROTOCOL):
        raise ValueError("Exact protocol manifest required")
    paths = [x["path"] for x in manifest["files"]]
    if set(paths) != REQUIRED_FILES or len(paths) != len(REQUIRED_FILES):
        raise ValueError("Exact code/source/test inventory required")
    for item in manifest["files"]:
        if set(item) != {"path", "sha256"} or sha((ROOT/item["path"]).read_bytes()) != item["sha256"]:
            raise ValueError("Frozen implementation byte mismatch")
    if sha((HERE/"STUDY_PLAN.json").read_bytes()) != PLAN_SHA256:
        raise ValueError("Frozen period-reference study changed")
    return sha(raw)


def request_for(age, sex, horizon):
    return dict(schema="period39-request/1", model="RUS_WPP2023_period_reference",
        country="RUS", source_year=2023, source_sex=sex, age_exact_years=age,
        horizon_years=horizon, age_origin="exact_birthday",
        period_assumption="hold_2023_attained_age_schedule", purpose="population_period_reference_only")


def csv_bytes(rows):
    fields = ["country", "source_year", "source_sex", "age_exact_years", "horizon_years",
              "end_exact_birthday", "CF_probability", "CF_rounding_lower", "CF_rounding_upper",
              "UDD_ax_0_5_probability", "UDD_rounding_lower", "UDD_rounding_upper",
              "UDD_minus_CF", "stored_author_qx_product", "stored_qx_minus_recomputed_CF",
              "output_scope", "clinical_validation", "current_calendar_forecast"]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        req, methods, diagnostic = row["request"], row["methods"], row["diagnostic"]
        cf, udd = (methods[k] for k in METHODS)
        writer.writerow(dict(country="RUS", source_year=2023, source_sex=req["source_sex"],
            age_exact_years=req["age_exact_years"], horizon_years=req["horizon_years"],
            end_exact_birthday=row["end_exact_birthday"], CF_probability=cf["decimal"]["probability"],
            CF_rounding_lower=cf["decimal"]["lower"], CF_rounding_upper=cf["decimal"]["upper"],
            UDD_ax_0_5_probability=udd["decimal"]["probability"], UDD_rounding_lower=udd["decimal"]["lower"],
            UDD_rounding_upper=udd["decimal"]["upper"], UDD_minus_CF=udd["probability"]-cf["probability"],
            stored_author_qx_product=diagnostic["decimal"]["stored_qx_product"],
            stored_qx_minus_recomputed_CF=diagnostic["decimal"]["delta_vs_constant_force"],
            output_scope="provisional_local_derived_period_reference", clinical_validation=False,
            current_calendar_forecast=False))
    return stream.getvalue().encode()


def run(out):
    out = Path(out)
    if out.exists():
        raise ValueError("New output directory required; no overwrite")
    code_hash = verify_freeze()
    plan = json.loads((HERE/"STUDY_PLAN.json").read_text())
    bundle = load_bundle()
    rows = [score_reference(request_for(age, sex, h), bundle) for sex in ("F", "M")
            for age in range(40, 81) for h in (5, 10)]
    if len(rows) != plan["output_plan"]["reference_grid_rows"]:
        raise ArithmeticError("Wrong reference-grid cardinality")
    keyed = {(r["request"]["source_sex"], r["request"]["age_exact_years"], r["request"]["horizon_years"]): r for r in rows}
    chain = []
    horizon_checks = 0
    for sex in ("F", "M"):
        for age in range(40, 81):
            for method in METHODS:
                p5 = keyed[sex, age, 5]["methods"][method]["probability"]
                p10 = keyed[sex, age, 10]["methods"][method]["probability"]
                if p10+1e-14 < p5:
                    raise ArithmeticError("Longer horizon probability cannot be smaller")
                horizon_checks += 1
        for age in range(40, 76):
            checks = {}
            for method in METHODS:
                a = keyed[sex, age, 5]["methods"][method]["probability"]
                b = keyed[sex, age+5, 5]["methods"][method]["probability"]
                total = keyed[sex, age, 10]["methods"][method]["probability"]
                error = abs(total-(a+(1-a)*b))
                if error > 2e-14:
                    raise ArithmeticError("Survival chaining identity failed")
                checks[method] = {"absolute_error": error, "passed": True}
            chain.append(dict(source_sex=sex, age_exact_years=age, methods=checks))
    max_width = {m: 0. for m in METHODS}
    for row in rows:
        h = row["request"]["horizon_years"]
        for method in METHODS:
            values = row["methods"][method]
            lo, hi = (values["rounding_bounds"][k] for k in ("lower", "upper"))
            p = values["probability"]
            if not 0 <= lo <= p <= hi <= 1 or hi-lo > h*1e-6*(1+1e-5):
                raise ArithmeticError("Quantization probability-domain/width check failed")
            max_width[method] = max(max_width[method], hi-lo)
    age_decreases = []
    for sex in ("F", "M"):
        for h in (5, 10):
            for age in range(40, 80):
                if keyed[sex, age+1, h]["methods"]["constant_force"]["probability"] < keyed[sex, age, h]["methods"]["constant_force"]["probability"]:
                    age_decreases.append(dict(source_sex=sex, age_from=age, horizon_years=h))
    maximum_difference = max(abs(r["methods"]["UDD_ax_0_5"]["probability"]-r["methods"]["constant_force"]["probability"]) for r in rows)
    maximum_qx_delta = max(abs(r["diagnostic"]["delta_vs_constant_force"]) for r in rows)
    selected = [r for r in rows if r["request"]["age_exact_years"] in plan["output_plan"]["report_selected_ages"]]
    output_csv = csv_bytes(rows)
    result = dict(schema="period39-results/1", protocol_commit=PROTOCOL, plan_sha256=PLAN_SHA256,
        code_manifest_sha256=code_hash, python=platform.python_version(),
        scope="provisional_local_derived_period_reference_under_declared_units",
        source_population="allcause_mixed_general_population_including_disease",
        source_year=2023, period_assumption="hold_2023_attained_age_schedule", source_rows=100,
        reference_rows=len(rows), rows=rows, selected_rows=selected, chaining_rows=chain,
        chaining_age_sex_rows=len(chain), chaining_scalar_method_checks=len(chain)*2,
        horizon_coherence_scalar_checks=horizon_checks,
        rounding_envelope_max_width=max_width, age_decreases_diagnostic_only=age_decreases,
        max_UDD_CF_absolute_difference=maximum_difference, max_stored_qx_CF_absolute_difference=maximum_qx_delta,
        reference_csv_sha256=sha(output_csv), rounding_bounds_are_confidence_intervals=False,
        model_assumption_difference_is_uncertainty_CI=False,
        disease_free_baseline=False, original_RDA_numeric_scale_reverified=False,
        current_calendar_forecast=False, biomarker_adjusted=False, clinical_validation=False,
        models_fitted=False, new_person_level_inputs=0,
        SHFM_updated_5year_reproduced=False)
    raw = (json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n").encode()
    if verify_freeze() != code_hash:
        raise ValueError("Code/source changed during calculation")
    out.mkdir(parents=True, exist_ok=False)
    (out/"RESULTS_SUMMARY.json").write_bytes(raw)
    (out/"REFERENCE_TABLE.csv").write_bytes(output_csv)
    print(json.dumps(dict(status="completed_local_period_reference", reference_rows=len(rows),
        result_sha256=sha(raw), csv_sha256=sha(output_csv), clinical_validation=False)))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    run(parser.parse_args().out)
