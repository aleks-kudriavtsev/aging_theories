"""Research-only numerical primitives; NOT a trained mortality predictor.

No literature HR is used as a joint model coefficient. The caller must supply
all mutually exclusive death hazards, including competing/unknown causes.
Rates are per year; interval endpoints and horizons are years after a landmark.
"""
from __future__ import annotations
import math
from collections.abc import Mapping, Sequence
from numbers import Real


def _number(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def egfr_ckd_epi_2021(age: float, sex: str, creatinine: float,
                      unit: str = "mg/dL") -> float:
    """Adult eGFRcr in mL/min/1.73 m2, not mortality risk or a CKD diagnosis.

    Published sex variable: 'female' or 'male'; never silently impute it.
    Creatinine must be IDMS-standardized. Clinical applicability (e.g. unstable
    creatinine, pregnancy, extreme muscle mass) requires separate assessment.
    Source: NIDDK adult eGFR equations, reviewed May 2025.
    """
    age = _number(age, "age")
    scr = _number(creatinine, "creatinine")
    if age < 18 or scr <= 0:
        raise ValueError("Adult age >=18 and creatinine >0 are required")
    if sex not in ("female", "male"):
        raise ValueError("Published sex variable must be female or male")
    if unit == "umol/L":
        scr /= 88.4
    elif unit != "mg/dL":
        raise ValueError("unit must be mg/dL or umol/L")
    k, alpha, factor = (0.7, -0.241, 1.012) if sex == "female" else (0.9, -0.302, 1.0)
    ratio = scr / k
    result = 142 * min(ratio, 1.0)**alpha * max(ratio, 1.0)**(-1.2) * 0.9938**age * factor
    if not math.isfinite(result) or result <= 0:
        raise ValueError("Input outside numerical support")
    return result


def cumulative_incidence(interval_ends: Sequence[float],
                         hazards: Sequence[Mapping[str, float]],
                         horizon: float, target_causes: Sequence[str], *,
                         complete_cause_partition: bool = False) -> dict:
    """Exact CIF integration for piecewise-constant cause-specific hazards.

    The complete partition assertion is required; the routine cannot verify
    medical endpoint coverage. The result is mathematical integration, NOT
    external validation or a fitted prediction from a biomarker profile.
    No extrapolation past the last supplied interval is allowed.
    """
    if complete_cause_partition is not True:
        raise ValueError("Explicitly confirm a complete, nonoverlapping death-cause partition")
    ends = [_number(v, "interval end") for v in interval_ends]
    h = _number(horizon, "horizon")
    if not ends or len(ends) != len(hazards):
        raise ValueError("Provide one hazard mapping per nonempty interval")
    if ends[0] <= 0 or any(b <= a for a, b in zip(ends, ends[1:])):
        raise ValueError("Interval ends must be positive and strictly increasing")
    if not 0 <= h <= ends[-1]:
        raise ValueError("Horizon outside supported interval range")
    causes = tuple(hazards[0])
    if not causes or any(not isinstance(c, str) or not c.strip() for c in causes):
        raise ValueError("Cause names must be nonempty strings")
    cause_set = set(causes)
    rows = []
    for row in hazards:
        if set(row) != cause_set:
            raise ValueError("Every interval must contain the same complete cause set")
        normalized = {c: _number(row[c], "hazard") for c in causes}
        if any(v < 0 for v in normalized.values()):
            raise ValueError("Hazards must be nonnegative")
        try:
            total = math.fsum(normalized.values())
        except OverflowError as exc:
            raise ValueError("Total hazard exceeds numerical support") from exc
        if not math.isfinite(total):
            raise ValueError("Total hazard exceeds numerical support")
        rows.append((normalized, total))
    if isinstance(target_causes, (str, bytes)):
        raise ValueError("target_causes must be a sequence of cause names, not a string")
    target = tuple(target_causes)
    if any(not isinstance(c, str) for c in target):
        raise ValueError("Target names must be strings")
    if len(set(target)) != len(target) or not set(target) <= cause_set:
        raise ValueError("Target causes must be unique and present in the cause partition")
    cif = {c: 0.0 for c in causes}
    survival, previous, cumulative_hazard = 1.0, 0.0, 0.0
    for end, (row, total) in zip(ends, rows):
        dt = max(0.0, min(end, h) - previous)
        if total > 0 and dt > 0:
            exposure = total * dt
            event_mass = survival * -math.expm1(-exposure)
            for cause in causes:
                cif[cause] += event_mass * (row[cause] / total)
            survival *= math.exp(-exposure)
            cumulative_hazard += exposure
        previous = end
        if end >= h:
            break
    all_cause = -math.expm1(-cumulative_hazard)
    if abs(survival + math.fsum(cif.values()) - 1.0) > 1e-10:
        raise ArithmeticError("Probability mass conservation failed")
    return {"horizon_years": h, "survival": survival,
            "all_cause_death_risk": all_cause,
            "target_death_risk": math.fsum(cif[c] for c in target),
            "cause_specific_cif": cif,
            "status": "numerical_integration_only_not_clinically_validated"}
