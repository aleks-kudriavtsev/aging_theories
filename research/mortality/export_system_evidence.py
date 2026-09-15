"""Export frozen model context to aging_biomarkers; never export coefficients.

python -S -m research.mortality.export_system_evidence --out exchange.json

Only two published aggregate summaries supply performance numbers. Each number
keeps its exact JSON pointer, source bytes, model identity and evaluation domain.
This is an exchange, not a new evaluation or clinical validation. The registered
clinical3 configuration is a coarse runtime view; transport12 also evaluated the
same clinical coefficient artifact with a six-state outcome. That historical
evaluation view is explicit and does not expand the runtime's supported inputs.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess

from .studio16 import registry
from .workbench.runtime import strict_json

SCHEMA_VERSION = "aging_system_evidence/1.0"
REPOSITORY = "aleks-kudriavtsev/aging_theories"
SOURCE_COMMIT = "1e5772ad813ef0d3f1a18f36c4bf34949972e039"
REPO_ROOT = Path(__file__).resolve().parents[2]
COMPACT_PATH = "research/mortality/compact14/RESULTS_SUMMARY.json"
TRANSPORT_PATH = "research/mortality/transport12/RESULTS_SUMMARY.json"
SUMMARY_PATHS = (COMPACT_PATH, TRANSPORT_PATH)
SOURCE_PATHS = SUMMARY_PATHS + (
    "research/mortality/studio16/registry.py",
    "research/mortality/workbench/runtime.py",
    "research/mortality/engine.py",
    "research/mortality/compact14/model_compact14.json",
    "research/mortality/transport12/model_bundle12.json",
    "research/mortality/transport12/model_sensitivity12.json",
    "research/mortality/external18/README.md",
    "research/mortality/external18/PROTOCOL_RU.md",
    "research/mortality/external18/COHORTS.json",
)
LAB_INPUTS = {
    "compact4": ["hba1c", "creatinine", "uacr", "albumin"],
    "six6": ["hba1c", "creatinine", "uacr", "albumin", "total_cholesterol", "hdl"],
    "full8": ["hba1c", "creatinine", "uacr", "albumin", "total_cholesterol", "hdl", "rdw", "wbc"],
    "clinical3": [],
}
INPUT_MAPPINGS = {
    "hba1c": {"features": ["hba1c"], "kind": "laboratory_measurement", "contains_age": False},
    "creatinine": {"features": ["egfr_low", "egfr_high"], "kind": "derived_measure_from_laboratory_input",
        "additional_inputs": ["age_years", "sex"], "derivation": "CKD_EPI_2021_creatinine_then_piecewise_scaling", "contains_age": True},
    "uacr": {"features": ["log_uacr"], "kind": "derived_laboratory_ratio",
        "underlying_measurements": ["urine_albumin", "urine_creatinine"], "derivation": "log2_of_UACR", "contains_age": False},
    "albumin": {"features": ["albumin"], "kind": "laboratory_analyte", "contains_age": False},
    "total_cholesterol": {"features": ["total_cholesterol"], "kind": "laboratory_measurement", "contains_age": False},
    "hdl": {"features": ["hdl"], "kind": "HDL_cholesterol_not_particle_count", "contains_age": False},
    "rdw": {"features": ["rdw"], "kind": "RDW_CV_not_RDW_SD", "contains_age": False},
    "wbc": {"features": ["log_wbc"], "kind": "cell_count", "derivation": "log2_of_WBC", "contains_age": False},
}


def pointer(document, path):
    """Read an RFC6901 pointer without guessing a missing field."""
    if not path.startswith("/"):
        raise ValueError("Absolute JSON pointer required")
    value = document
    for token in path[1:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        value = value[int(token)] if isinstance(value, list) else value[token]
    return value


def pinned_inputs(source_commit=SOURCE_COMMIT):
    """Allow new exporter files, but reject changed/missing scientific inputs."""
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit) or source_commit != SOURCE_COMMIT:
        raise ValueError("This exporter requires the reviewed full source commit " + SOURCE_COMMIT)
    result = {}
    for path in SOURCE_PATHS:
        try:
            historical = subprocess.run(["git", "show", source_commit + ":" + path],
                cwd=REPO_ROOT, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
            current = (REPO_ROOT / path).read_bytes()
        except (OSError, subprocess.CalledProcessError) as exc:
            raise ValueError("Pinned input unavailable: " + path) from exc
        if current != historical:
            raise ValueError("Scientific input differs from pinned source: " + path)
        result[path] = current
    return result


def models():
    result = registry.catalog()["models"]  # checks all three frozen artifact hashes
    for row in result:
        inputs = LAB_INPUTS["full8" if row["model_id"] == "full8_detailed" else row["model_id"]]
        row["laboratory_predictor_measurement_ids"] = deepcopy(inputs)
        row["predictor_input_mappings"] = {k: deepcopy(INPUT_MAPPINGS[k]) for k in inputs}
        row["required_input_measurement_ids"] = deepcopy(row["required"])
        row["nonpredictive_domain_requirement_ids"] = [
            k for k in row["required"] if k not in ("bmi", "sbp") and k not in inputs]
        row["measurement_mapping_is_evidence_of_added_value"] = False
    return result


def record(documents, path, value_pointer, model_id, metric, context, *,
           base_model_id=None, interval_pointers=None, interval_scope=None,
           replicates=None, evaluation_view=None):
    raw = documents[path]
    document = strict_json(raw)
    value = pointer(document, value_pointer)
    if type(value) not in (int, float):
        raise ValueError("Performance value is not a finite numeric estimate")
    result = {
        "id": "aging_theories:" + path.split("/")[-2] + ":" + value_pointer.lstrip("/"),
        "record_type": "empirical_model_evaluation", "model_id": model_id,
        "metric": metric, "value": value, "context": deepcopy(context),
        "component_attribution": "model_or_block_only",
        "source": {"repository": REPOSITORY, "commit": SOURCE_COMMIT, "path": path,
            "sha256": sha256(raw).hexdigest(), "json_pointer": value_pointer},
        "clinical_use_ready": False, "usable_as_mortality_coefficient": False,
    }
    if base_model_id is not None:
        result["base_model_id"] = base_model_id
    if evaluation_view is not None:
        result["evaluation_view"] = evaluation_view
    if interval_pointers is not None:
        lo, hi = (pointer(document, p) for p in interval_pointers)
        if lo is not None and hi is not None and lo > hi:
            raise ValueError("Reversed interval")
        result["interval"] = {"lower": lo, "upper": hi,
            "kind": "conditional_validation_percentiles", "scope": interval_scope,
            "replicates": replicates,
            "source_pointers": {"lower": interval_pointers[0], "upper": interval_pointers[1]}}
    return result


def context(*, cohort_id, horizon, n, deaths, states, censoring, endpoint=None):
    return {"endpoint": endpoint or ("all_cause_mortality" if states == 2 else
            "three_coarse_death_groups" if states == 4 else "five_competing_death_groups"),
        "horizon_years": horizon, "horizon_basis": "prediction_horizon",
        "cohort_id": cohort_id, "country": "US", "sex": "pooled", "age_range": [40, 79],
        "n": n, "deaths": deaths, "validation_kind": "temporal", "state_count": states,
        "Brier_definition": "binary_squared_error" if states == 2 else "sum_over_states_no_half_factor",
        "censoring_method": censoring, "population": "common_complete_full8_profile",
        "clinical_use_ready": False}


def compact_records(documents):
    doc = strict_json(documents[COMPACT_PATH]); rows = []
    ctx = dict(cohort_id=doc["evaluation_cycle"], horizon=3, n=doc["common_complete_domain_n"],
        deaths=doc["deaths3"], censoring="IPCW_reverse_Kaplan_Meier")
    binary = context(**ctx, states=2); multi = context(**ctx, states=4)
    binary["observed_IPCW_death_probability"] = doc["observed_IPCW_death_probability3"]
    multi["observed_IPCW_death_probability"] = doc["observed_IPCW_death_probability3"]
    metrics = (("four_state_brier", "multiclass_Brier_score", multi),
        ("all_cause_brier", "Brier_score", binary), ("auc", "fixed_horizon_AUC", binary),
        ("predicted", "predicted_death_proportion", binary),
        ("observed_expected", "observed_expected_ratio", binary),
        ("calibration_slope", "calibration_slope", binary),
        ("bias_pp", "prediction_bias_percentage_points", binary))
    for i, panel in enumerate(doc["panels"]):
        mid = "clinical3" if panel["panel"] == "clinical" else panel["panel"]
        for field, metric, c in metrics:
            rows.append(record(documents, COMPACT_PATH, f"/panels/{i}/{field}", mid, metric, c))
    contrasts = (
        ("primary_contrast", "difference", "delta_multiclass_Brier", multi, "clinical3", "lower", "upper"),
        ("all_cause_AUC_contrast", "difference", "delta_AUC", binary, "clinical3", "lower", "upper"),
        ("compact_minus_full8", "four_state_brier_difference", "delta_multiclass_Brier", multi, "full8", "interval/0", "interval/1"),
        ("compact_minus_full8", "auc_difference", "delta_AUC", binary, "full8", "auc_interval/0", "auc_interval/1"),
        ("compact_minus_six6", "four_state_brier_difference", "delta_multiclass_Brier", multi, "six6", "interval/0", "interval/1"),
    )
    for section, field, metric, c, base, lo, hi in contrasts:
        rows.append(record(documents, COMPACT_PATH, f"/{section}/{field}", "compact4", metric, c,
            base_model_id=base, interval_pointers=(f"/{section}/{lo}", f"/{section}/{hi}"),
            interval_scope=doc["interval_scope"], replicates=doc["bootstrap_replicates"]))
    return rows


def transport_records(documents):
    doc = strict_json(documents[TRANSPORT_PATH]); rows = []
    for i, result in enumerate(doc["results"]):
        ctx = dict(cohort_id=f"NHANES{result['cycle']}_{result['cycle']+1}",
            horizon=result["horizon"], n=result["n"], deaths=result["deaths"],
            censoring="complete_followup_through_prediction_horizon")
        binary = context(**ctx, states=2); multi = context(**ctx, states=6)
        metrics = (
            ("clinical_six_state_brier", "clinical3", "multiclass_Brier_score", multi),
            ("routine_six_state_brier", "full8_detailed", "multiclass_Brier_score", multi),
            ("clinical_all_cause_auc", "clinical3", "fixed_horizon_AUC", binary),
            ("routine_all_cause_auc", "full8_detailed", "fixed_horizon_AUC", binary),
            ("weighted_observed", "full8_detailed", "observed_death_proportion", binary),
            ("weighted_predicted", "full8_detailed", "predicted_death_proportion", binary),
            ("calibration_slope", "full8_detailed", "calibration_slope", binary),
        )
        for field, mid, metric, c in metrics:
            rows.append(record(documents, TRANSPORT_PATH, f"/results/{i}/{field}", mid, metric, c,
                evaluation_view="transport12_five_cause_evaluation_of_frozen_coefficients"))
        for field, metric, c, interval in (
            ("delta_six_state_brier", "delta_multiclass_Brier", multi, "delta_six_state_brier_interval"),
            ("delta_all_cause_auc", "delta_AUC", binary, "delta_auc_interval")):
            rows.append(record(documents, TRANSPORT_PATH, f"/results/{i}/{field}", "full8_detailed", metric, c,
                base_model_id="clinical3", interval_pointers=(f"/results/{i}/{interval}/0", f"/results/{i}/{interval}/1"),
                interval_scope=doc["interval_definition"], replicates=1000,
                evaluation_view="transport12_five_cause_evaluation_of_frozen_coefficients"))
    return rows


def build(source_commit=SOURCE_COMMIT):
    inputs = pinned_inputs(source_commit)
    compact = strict_json(inputs[COMPACT_PATH]); transport = strict_json(inputs[TRANSPORT_PATH])
    if compact["development_n"] != transport["development_n"] or compact["development_n"] != 15280:
        raise ValueError("Development population changed")
    source_artifacts = {path: {"sha256": sha256(raw).hexdigest()} for path, raw in inputs.items()}
    source_documents = {path: {"sha256": source_artifacts[path]["sha256"], "raw_document": inputs[path].decode("utf-8")}
        for path in SUMMARY_PATHS}
    rows = compact_records(inputs) + transport_records(inputs)
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate performance identity")
    return {"schema_version": SCHEMA_VERSION,
        "source": {"repository": REPOSITORY, "commit": source_commit}, "development_n": 15280,
        "models": models(), "performance_records": rows,
        "source_documents": source_documents, "source_artifacts": source_artifacts,
        "external_validation": {"RU": False, "DE": False, "protocol_only": True,
            "source_path": "research/mortality/external18/PROTOCOL_RU.md",
            "status": "protocol_and_software_checks_without_external_cohort_data",
            "clinical_use_ready": False},
        "safeguards": {"coefficients_exported": False, "participant_data_exported": False,
            "literature_effects_used_as_coefficients": False, "model_refit": False,
            "new_validation_performed": False, "individual_marker_added_value_inferred": False,
            "clinical_use_ready": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source-commit", default=SOURCE_COMMIT)
    args = parser.parse_args()
    doc = build(args.source_commit)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as stream:
        json.dump(doc, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"schema_version": SCHEMA_VERSION, "models": len(doc["models"]),
        "performance_records": len(doc["performance_records"]), "clinical_use_ready": False}))


if __name__ == "__main__":
    main()
