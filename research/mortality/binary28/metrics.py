"""Aggregate weighted binary scores and diagnostic calibration of OOF risks.

Diagnostics never modify the evaluated predictions. Calibration models are
unpenalized weighted logistic likelihood fits; no confidence interval or
clinical calibration claim is produced by this module.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit
from ..validate_temporal08 import weighted_auc

CALIBRATION_CLIP = 1e-12


def _validated(y, p, w):
    y, p, w = [np.asarray(a, dtype=float) for a in (y, p, w)]
    if (y.ndim != 1 or not len(y) or p.shape != y.shape or w.shape != y.shape
            or not np.isfinite(np.r_[y, p, w]).all()
            or not np.isin(y, [0, 1]).all() or np.any((p < 0) | (p > 1))
            or np.any(w < 0) or not np.isfinite(w.sum()) or w.sum() <= 0):
        raise ValueError('Finite aligned binary outcomes, probabilities and nonnegative weights required')
    return y, p, w


def _calibration_fit(y, design, weights, initial, offset=None):
    """Fit diagnostics only; explicit finite likelihood/gradient checks."""
    w = (weights / weights.sum()) * len(weights)
    offset = np.zeros(len(y)) if offset is None else np.asarray(offset, float)

    def objective(beta):
        eta = offset + design @ beta
        if not np.isfinite(eta).all():
            return float('inf'), np.full(len(beta), np.nan)
        value = float(w @ (np.logaddexp(0., eta) - y * eta))
        gradient = design.T @ (w * (expit(eta) - y))
        return value, gradient

    result = minimize(objective, np.asarray(initial, float), jac=True,
                      method='BFGS', options={'gtol': 1e-8, 'maxiter': 500})
    value, gradient = objective(result.x)
    if (not np.isfinite(result.x).all() or not np.isfinite(value)
            or not np.isfinite(gradient).all()
            or np.max(np.abs(gradient)) > 1e-5):
        raise ArithmeticError('Calibration diagnostic fit not finite or converged')
    return result.x, float(np.max(np.abs(gradient)))


def metrics(y, p, w):
    """Score one unchanged OOF probability vector with survey/replicate weights."""
    y, p, w = _validated(y, p, w)
    positive = w > 0
    observed = float(w @ y / w.sum())
    predicted = float(w @ p / w.sum())
    two_classes = np.unique(y[positive]).size == 2
    result = dict(n=len(y), events=int(y.sum()), n_positive_weight=int(positive.sum()),
                  events_positive_weight=int(y[positive].sum()), weight_sum=float(w.sum()),
                  auc=weighted_auc(y, p, w) if two_classes else None,
                  brier=float(w @ ((y - p) ** 2) / w.sum()),
                  observed=observed, predicted=predicted,
                  bias_pp=100 * (predicted - observed),
                  observed_expected=observed / predicted if predicted > 0 else None,
                  calibration_in_large=None, calibration_intercept=None,
                  calibration_slope=None, calibration_status='undefined_single_class',
                  calibration_gradient_max=None,
                  calibration_probability_clip=CALIBRATION_CLIP,
                  calibration_clipped_positive_weight_n=int(((p <= CALIBRATION_CLIP)
                      | (p >= 1 - CALIBRATION_CLIP))[positive].sum()))
    if not two_classes:
        return result

    # Clipping applies only to this diagnostic logit, never to scored p.
    x = logit(np.clip(p, CALIBRATION_CLIP, 1 - CALIBRATION_CLIP))
    try:
        initial = float(logit(observed) - np.average(x, weights=w))
        offset_fit, offset_gradient = _calibration_fit(
            y, np.ones((len(y), 1)), w, [initial], offset=x)
        result['calibration_in_large'] = float(offset_fit[0])
        result['calibration_gradient_max'] = offset_gradient
    except (ValueError, ArithmeticError, FloatingPointError):
        result['calibration_status'] = 'undefined_offset_fit_nonconverged'
        return result

    if np.ptp(x[positive]) <= 1e-12:
        result['calibration_status'] = 'flat_prediction_slope_unidentified'
        return result
    controls, cases = x[positive & (y == 0)], x[positive & (y == 1)]
    # Exact or quasi separation makes finite intercept+slope MLE unidentified.
    if controls.max() <= cases.min() or cases.max() <= controls.min():
        result['calibration_status'] = 'separated_calibration_slope_unidentified'
        return result
    try:
        full_fit, full_gradient = _calibration_fit(
            y, np.column_stack([np.ones(len(y)), x]), w, [0., 1.])
        result['calibration_intercept'] = float(full_fit[0])
        result['calibration_slope'] = float(full_fit[1])
        result['calibration_gradient_max'] = max(offset_gradient, full_gradient)
        result['calibration_status'] = 'computed'
    except (ValueError, ArithmeticError, FloatingPointError):
        result['calibration_status'] = 'offset_only_full_fit_nonconverged'
    return result


def calibration_bins(y, p, w):
    """Up to five weighted-risk quintile bins; tied boundaries may collapse.

    Inverse weighted empirical-CDF boundaries use positive-weight observations.
    All observations sharing a risk remain together. Actual nonempty assigned
    bins are returned; zero-weight-only bins have undefined weighted summaries.
    Only aggregates are emitted, without participant identifiers or risks.
    """
    y, p, w = _validated(y, p, w)
    positive = w > 0
    values, inverse = np.unique(p[positive], return_inverse=True)
    masses = np.bincount(inverse, weights=w[positive], minlength=len(values))
    cumulative = np.cumsum(masses) / masses.sum()
    at = np.searchsorted(cumulative, [.2, .4, .6, .8], side='left')
    boundaries = np.unique(values[np.minimum(at, len(values) - 1)])
    assigned = np.searchsorted(boundaries, p, side='left')
    rows = []
    for number in np.unique(assigned):
        mask = assigned == number
        total = float(w[mask].sum())
        observed = float(w[mask] @ y[mask] / total) if total > 0 else None
        predicted = float(w[mask] @ p[mask] / total) if total > 0 else None
        rows.append(dict(bin=int(number + 1), n=int(mask.sum()), events=int(y[mask].sum()),
                         n_positive_weight=int((mask & positive).sum()), weight_sum=total,
                         observed=observed, predicted=predicted,
                         bias_pp=100 * (predicted - observed) if total > 0 else None))
    return rows
