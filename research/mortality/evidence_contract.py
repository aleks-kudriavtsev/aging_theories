"""Typed research evidence exchange. This module never predicts patient risk.

Inspired by aging_biomarkers PR115, reviewed at commit 3c487a8. Literature,
synthetic controls, and empirical model evaluations remain different objects.
Missing metadata quarantine a use, not erase the record. Standard library only.
"""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import re
from typing import Any, Mapping

SCHEMA_VERSION = "0.5"
METRIC_KINDS = {
    "HR": "association", "RR": "association", "OR": "association",
    "sensitivity": "diagnostic_accuracy", "specificity": "diagnostic_accuracy",
    "PPV": "diagnostic_accuracy", "NPV": "diagnostic_accuracy",
    "C_index": "discrimination", "Harrell_C": "discrimination",
    "AUC": "discrimination", "fixed_horizon_AUC": "discrimination",
    "time_dependent_AUC": "discrimination",
    "delta_C": "incremental_discrimination", "delta_AUC": "incremental_discrimination",
    "continuous_NRI": "reclassification", "categorical_NRI": "reclassification",
    "IDI": "reclassification", "Brier_score": "overall_prediction_error",
    "delta_Brier": "incremental_prediction_error",
    "calibration_slope": "calibration", "calibration_in_the_large": "calibration",
    "observed_expected_ratio": "calibration", "net_benefit": "clinical_utility",
    "observed_death_proportion": "observed_frequency",
    "predicted_death_proportion": "predicted_frequency",
    "reported_relative_reduction": "association_not_a_verified_HR",
    "reported_relative_increase": "association_not_a_verified_HR",
}
DOMAINS = {
    **{m: (0., 1.) for m in ("C_index", "Harrell_C", "AUC", "fixed_horizon_AUC",
        "time_dependent_AUC", "Brier_score", "sensitivity", "specificity", "PPV", "NPV",
        "observed_death_proportion", "predicted_death_proportion")},
    "delta_C": (-1., 1.), "delta_AUC": (-1., 1.), "delta_Brier": (-1., 1.),
    "continuous_NRI": (-2., 2.), "categorical_NRI": (-2., 2.),
    # General IDI is a difference of two discrimination slopes. Each slope is
    # in [-1,1]; do not assume calibrated/nested models when importing a paper.
    "IDI": (-2., 2.),
}
INCREMENTS = {"delta_C", "delta_AUC", "delta_Brier", "IDI", "continuous_NRI", "categorical_NRI"}
INTERVAL_KINDS = {"95% CI", "90% CI", "99% CI", "credible_interval", "IQR", "range",
    "conditional_validation_percentiles", "not_reported", "unknown_reported_interval"}
CONTEXT_FIELDS = ("endpoint", "population", "country", "sex", "age_range", "horizon_years",
                  "horizon_basis", "cohort_id", "validation_set")
MORTALITY_ENDPOINTS = {"all_cause_mortality", "cardiovascular_mortality", "respiratory_mortality",
                       "liver_mortality", "noninfectious_natural_mortality"}


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def present(value: Any) -> bool:
    if value is None or value == "" or value == []:
        return False
    return not (isinstance(value, str) and value.lower() in {
        "unknown", "not_reported", "not reported", "unspecified", "none", "null"})


def finite(value: Any) -> bool:
    return type(value) in (int, float) and math.isfinite(value)


def number_matches(value: float, token: str, *, scale: str = "identity") -> bool:
    """Exact signed numeric witness, with an explicitly declared source scale.

    This checks transcription only, not whether the number refers to the right
    model, sex, endpoint or cohort. No searching for abs(value), no guessed *1000.
    Token must be the actual source substring selected for this estimate.
    """
    if not finite(value) or not isinstance(token, str):
        return False
    if scale not in {"identity", "percent", "per_thousand"}:
        raise ValueError("Source scale must be explicitly identity, percent or per_thousand")
    text = token.strip().replace("−", "-")
    suffix = "%" if scale == "percent" else "‰" if scale == "per_thousand" else ""
    if suffix and text.endswith(suffix):
        text = text[:-1].strip()
    if not re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", text):
        return False
    try:
        parsed = Decimal(text) / {"identity": Decimal(1), "percent": Decimal(100),
                                  "per_thousand": Decimal(1000)}[scale]
        return parsed == Decimal(str(value))
    except InvalidOperation:
        return False


def audit_record(record: Mapping[str, Any]) -> dict:
    """Return independent retention, comparability and clinical-use statuses.

    All input records remain retainable, including malformed/unsupported ones;
    errors bar analytic use. This is not a source-verification engine.
    """
    errors, gaps = [], []
    if not isinstance(record, Mapping):
        return {"retain": True, "errors": ["invalid_record_mapping"], "gaps": [],
                "comparison_ready": False, "clinical_use_ready": False,
                "usable_as_mortality_coefficient": False}
    context, estimate, source = (record.get(k, {}) for k in ("context", "estimate", "source"))
    if not all(isinstance(v, Mapping) for v in (context, estimate, source)):
        return {"retain": True, "errors": ["invalid_mapping"], "gaps": [],
                "comparison_ready": False, "clinical_use_ready": False,
                "usable_as_mortality_coefficient": False}
    metric, value = estimate.get("metric"), estimate.get("value")
    if not isinstance(metric, str):
        errors.append("invalid_metric_type")
        metric = "unresolved"
    if metric not in METRIC_KINDS:
        gaps.append("metric_definition_unresolved")
    if not finite(value):
        errors.append("nonfinite_estimate")
    elif metric in DOMAINS and not DOMAINS[metric][0] <= value <= DOMAINS[metric][1]:
        errors.append("estimate_outside_theoretical_domain")
    elif metric in {"HR", "RR", "OR"} and value <= 0:
        errors.append("effect_ratio_must_be_positive")
    elif metric == "observed_expected_ratio" and value < 0:
        errors.append("observed_expected_ratio_must_be_nonnegative")
    # CIs from unconstrained estimators may extend outside a theoretical domain;
    # IQR need not contain an independently reported mean. Preserve, do not clip.
    interval = estimate.get("interval")
    if interval is not None:
        if not isinstance(interval, Mapping):
            errors.append("invalid_interval_mapping")
        else:
            lo, hi = interval.get("lower"), interval.get("upper")
            if not finite(lo) or not finite(hi) or lo > hi:
                errors.append("invalid_interval_bounds")
            if not isinstance(interval.get("kind"), str) or interval.get("kind") not in INTERVAL_KINDS:
                errors.append("invalid_interval_kind")
            elif interval.get("kind") in {"not_reported", "unknown_reported_interval"}:
                gaps.append("uncertainty_kind_unresolved")
    se = estimate.get("standard_error")
    if se is not None and (not finite(se) or se < 0):
        errors.append("invalid_standard_error")
    if interval is None and se is None:
        gaps.append("uncertainty_not_reported")
    for key in CONTEXT_FIELDS:
        if not context_known(key, context.get(key)):
            gaps.append("missing_context_" + key)
    horizon = context.get("horizon_years")
    if horizon is not None and (not finite(horizon) or horizon <= 0):
        errors.append("invalid_horizon")
    if present(horizon) and context.get("horizon_basis") != "prediction_horizon":
        gaps.append("followup_is_not_prediction_horizon")
    model = record.get("model") or {}
    if not isinstance(model, Mapping):
        errors.append("invalid_model_mapping")
        model = {}
    if metric in INCREMENTS:
        for key in ("id", "base_id", "added_components"):
            if not present(model.get(key)):
                gaps.append("missing_comparison_" + key)
    if metric == "net_benefit":
        threshold = estimate.get("decision_threshold")
        if threshold is None:
            gaps.append("missing_decision_threshold")
        elif not finite(threshold) or not 0 < threshold < 1:
            errors.append("invalid_decision_threshold")
        if not present(estimate.get("decision_action")):
            gaps.append("missing_decision_action")
        if not present(estimate.get("net_benefit_definition")):
            gaps.append("missing_net_benefit_definition")
    if source.get("kind") == "synthetic_control" or source.get("synthetic") is True:
        errors.append("synthetic_not_evidence")
    if not present(source.get("locator")) or not present(source.get("sha256")):
        gaps.append("missing_source_provenance")
    if not isinstance(source.get("review"), str) or source.get("review") not in {"source_checked", "frozen_empirical_output"}:
        gaps.append("source_verification_pending")
    if source.get("quarantined"):
        errors.append("upstream_quality_quarantine")
    # A review/search-route count is not independent replication.
    if source.get("kind") == "literature_extraction" and not record.get("cohort_independence_adjudicated"):
        gaps.append("cohort_independence_not_adjudicated")
    if record.get("record_type") == "prediction_performance" and not present(model.get("id")):
        gaps.append("missing_model_id")
    return {"retain": True, "errors": sorted(set(errors)), "gaps": sorted(set(gaps)),
            "evidence_kind": METRIC_KINDS.get(metric, "unresolved"),
            "comparison_ready": not errors and not gaps,
            "clinical_use_ready": False, "usable_as_mortality_coefficient": False}


def context_known(key: str, value: Any) -> bool:
    if key == "age_range":
        return (isinstance(value, (list, tuple)) and len(value) == 2
                and all(finite(v) for v in value) and 0 <= value[0] <= value[1])
    return present(value)


def compare_contexts(left: Mapping, right: Mapping) -> dict:
    """Strict same-context check; unknown != wildcard. No averaging is done."""
    missing, different = [], []
    for key in CONTEXT_FIELDS:
        a, b = left.get("context", {}).get(key), right.get("context", {}).get(key)
        if not context_known(key, a) or not context_known(key, b):
            missing.append(key)
        elif a != b:
            different.append(key)
    return {"status": "incompatible" if different else "unresolved" if missing else "same_context",
            "different": different, "missing": missing}


def attach_increment(association: Mapping, performance: Mapping) -> dict:
    """A specific increment may annotate a matching context, never a block member.

    This is metadata attachment, not proof that the increment is reproducible.
    Exact assay/units are required for a single laboratory marker.
    """
    # A malformed estimate remains in the registry but cannot be attached.
    validation_errors = []
    for name, row in (("association", association), ("performance", performance)):
        validation_errors.extend(name + ":" + e for e in audit_record(row)["errors"])
    if validation_errors:
        return {"attach": False, "reasons": sorted(set(validation_errors))}
    result = compare_contexts(association, performance)
    reasons = list(result["different"]) + ["unknown_" + f for f in result["missing"]]
    model, a_measure, p_measure = (performance.get("model") or {},
        association.get("measurement") or {}, performance.get("measurement") or {})
    if performance.get("estimate", {}).get("metric") not in INCREMENTS:
        reasons.append("not_an_increment")
    if model.get("component_kind") != "single_marker":
        reasons.append("block_or_model_effect_not_single_marker")
    components = model.get("added_components")
    if not isinstance(components, list) or components != [a_measure.get("id")]:
        reasons.append("different_component")
    for key in ("assay", "units", "matrix"):
        if not present(a_measure.get(key)) or not present(p_measure.get(key)):
            reasons.append("unknown_measurement_" + key)
        elif a_measure[key] != p_measure[key]:
            reasons.append("different_measurement_" + key)
    if not present(model.get("base_id")):
        reasons.append("missing_base_model")
    for name, row in (("association", association), ("performance", performance)):
        source = row.get("source") or {}
        if not isinstance(source.get("review"), str) or source.get("review") not in {"source_checked", "frozen_empirical_output"}:
            reasons.append(name + "_source_not_reviewed")
        if source.get("quarantined") or source.get("synthetic") or source.get("kind") == "synthetic_control":
            reasons.append(name + "_source_quarantined")
    return {"attach": not reasons, "reasons": sorted(set(reasons))}


def source_review_gate(record: Mapping, *, target_country: str, target_endpoint: str) -> dict:
    """Never upgrade extraction/retrospective evidence into clinical deployment."""
    reasons = []
    if record.get("context", {}).get("endpoint") != target_endpoint:
        reasons.append("endpoint_mismatch")
    if record.get("context", {}).get("country") != target_country:
        reasons.append("country_transport_not_validated")
    reasons.extend(["prospective_utility_not_validated", "not_a_deployable_model_specification"])
    return {"clinical_use_ready": False, "reasons": sorted(set(reasons))}


def project_estimates(records: list[dict]) -> list[dict]:
    """No union-of-sex/age summary: context travels with every single estimate."""
    return [{"id": r["id"], "context": deepcopy(r["context"]),
             "measurement": deepcopy(r.get("measurement")),
             "model": deepcopy(r.get("model")), "estimate": deepcopy(r["estimate"]),
             "source": deepcopy(r["source"]), "audit": audit_record(r)} for r in records]
