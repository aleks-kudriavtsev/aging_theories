"""Research primitives, not a fitted individual mortality calculator.

FIB-4: Sterling et al., PMID 16729309; applicability: AASLD guidance,
DOI 10.1097/HEP.0000000000000323. AJ: https://pmc.ncbi.nlm.nih.gov/articles/PMC4558089/.
Only inception cohorts, right censoring, one terminal event/person are supported.
"""
from __future__ import annotations
import math
from collections import Counter, defaultdict
from collections.abc import Sequence
from numbers import Real


def _real(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def fib4(age_years: float, ast_u_l: float, alt_u_l: float,
         platelets_1e9_l: float, *, acute_illness: bool = False) -> dict:
    """Numerical FIB-4 only. No automatic diagnosis or mortality probability.

    Units are fixed in parameter names; e.g. platelets=200, NOT 200000.
    Adult restriction is a software scope, not evidence of validity at every age.
    """
    age, ast, alt, platelets = [_real(v, n) for v, n in zip(
        (age_years, ast_u_l, alt_u_l, platelets_1e9_l),
        ("age_years", "ast_u_l", "alt_u_l", "platelets_1e9_l"))]
    if age < 18 or min(ast, alt, platelets) <= 0:
        raise ValueError("Age >=18 and positive AST, ALT, platelets required")
    if not isinstance(acute_illness, bool) or acute_illness:
        raise ValueError("Do not use this FIB-4 pathway in acute illness")
    try:
        value = (age / platelets) * (ast / math.sqrt(alt))
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("Input outside numerical support") from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError("Input outside numerical support")
    cautions = ["Index is not a diagnosis or an absolute mortality risk"]
    if age < 35:
        cautions.append("Limited accuracy below 35; no automatic interpretation")
    if age >= 65:
        cautions.append("Older-age thresholds require the chosen clinical guideline")
    return {"fib4": value, "mortality_probability": None,
            "status": "numerical_index_only", "cautions": cautions}


def observed_cif(durations: Sequence[float], event_causes: Sequence[str | None],
                 horizon: float, *, causes: Sequence[str],
                 complete_cause_partition: bool = False) -> dict:
    """Aalen-Johansen cumulative incidence, not biomarker-adjusted prediction.

    Time units must be identical. None means right censoring, NOT unknown death.
    Unknown death must have its own explicit label in causes. Competing deaths
    must remain events. Assume independent censoring within the analysed stratum.
    Deaths and censors tied at a time are all at risk immediately before it.
    Baseline is time zero, no delayed entry, recurrent events or weights.
    Completeness is the caller's assertion, not a medical coding audit.
    No extrapolation, confidence intervals, or causal interpretation is provided.
    """
    if complete_cause_partition is not True:
        raise ValueError("A complete mutually exclusive cause partition is required")
    if isinstance(causes, (str, bytes)):
        raise ValueError("causes must be a sequence of labels")
    labels = tuple(causes)
    if (not labels or any(not isinstance(c, str) or not c.strip() for c in labels)
            or len(set(labels)) != len(labels)):
        raise ValueError("Cause labels must be nonempty unique strings")
    if isinstance(durations, (str, bytes)) or isinstance(event_causes, (str, bytes)):
        raise ValueError("Observation arrays must be sequences, not strings")
    times = [_real(t, "duration") for t in durations]
    events = list(event_causes)
    h = _real(horizon, "horizon")
    if not times or len(times) != len(events) or min(times) <= 0:
        raise ValueError("Nonempty equal-length arrays and durations >0 required")
    if h < 0 or h > max(times):
        raise ValueError("Horizon outside observed follow-up; no extrapolation")
    if any(e is not None and (not isinstance(e, str) or e not in labels) for e in events):
        raise ValueError("Each death must have a declared cause; only None is censoring")
    groups = defaultdict(list)
    for time, cause in zip(times, events):
        groups[time].append(cause)
    at_risk, survival = len(times), 1.0
    cif, curve = dict.fromkeys(labels, 0.0), []
    for time in sorted(groups):
        if time > h:
            break
        group = groups[time]
        deaths = Counter(e for e in group if e is not None)
        d_all = sum(deaths.values())
        if d_all:
            previous_survival = survival
            for label, count in deaths.items():
                cif[label] += previous_survival * count / at_risk
            survival *= 1.0 - d_all / at_risk
            curve.append({"time": time, "n_at_risk_before": at_risk,
                          "survival": survival, "cause_specific_cif": dict(cif)})
        at_risk -= len(group)
    if abs(survival + math.fsum(cif.values()) - 1.0) > 1e-10:
        raise ArithmeticError("Probability mass conservation failed")
    return {"horizon": h, "n": len(times), "n_at_risk_just_before_horizon":
            sum(t >= h for t in times), "survival": survival,
            "all_cause_death_risk": math.fsum(cif.values()),
            "cause_specific_cif": cif, "event_time_curve": curve,
            "uncertainty_intervals": None,
            "status": "unadjusted_cohort_estimate_not_individual_prediction"}
