"""Reported Cys/B2M nested comparisons with independent model registries.

Uses the preserved convex ridge numerical primitives. All preprocessing is
learned within the actual training split; earlier modules are never mutated.
Protein log2 input columns are ordered (Cys, B2M).
"""
from __future__ import annotations

import numpy as np

from ..binary28.model import (
    PENALTIES, _matrix, _rows, fit_ridge, predict_ridge, psu_partition,
)
from ..routine29.model import CORE_NAMES

PROTEIN_NAMES = ('log2_cystatin_c', 'log2_beta2_microglobulin')
MODEL_INDICES = {'routine': (), 'routine_b2m': (1,), 'routine_cys': (0,),
                 'routine_cys_b2m': (0, 1)}


def _features(core, proteins):
    core, proteins = _matrix(core), _matrix(proteins)
    if core.shape[1] != 11 or proteins.shape != (len(core), 2):
        raise ValueError('Eleven routine and two protein features required')
    return core, proteins


def _unscaled(core, proteins, specification):
    name = specification.get('model_name')
    if (name not in MODEL_INDICES or specification.get('core_indices') != list(range(11))
            or specification.get('protein_indices') != list(MODEL_INDICES[name])):
        raise ValueError('Incorrect registered b2m31 feature selection')
    return np.column_stack([core, proteins[:, MODEL_INDICES[name]]])


def fit_preprocessing(core, protein_log2, weights, model_name):
    core, proteins = _features(core, protein_log2)
    _, w = _rows(np.zeros(len(core)), weights, len(core))
    if model_name not in MODEL_INDICES:
        raise ValueError('Unregistered b2m31 model')
    spec = dict(model_name=model_name, core_indices=list(range(11)),
                protein_indices=list(MODEL_INDICES[model_name]))
    raw = _unscaled(core, proteins, spec)
    mean = np.average(raw, axis=0, weights=w)
    constant = np.ptp(raw[w > 0], axis=0) == 0
    mean[constant] = raw[w > 0][0, constant]
    scale = np.sqrt(np.average((raw - mean)**2, axis=0, weights=w))
    scale[constant] = 1.
    if not np.isfinite(mean).all() or not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ArithmeticError('Invalid training feature scaling')
    spec.update(mean=mean.tolist(), scale=scale.tolist(), constant=constant.tolist())
    return spec


def transform_features(core, protein_log2, specification):
    core, proteins = _features(core, protein_log2)
    raw = _unscaled(core, proteins, specification)
    mean, scale = np.asarray(specification['mean']), np.asarray(specification['scale'])
    if mean.shape != (raw.shape[1],) or scale.shape != mean.shape or np.any(scale <= 0):
        raise ValueError('Incorrect preprocessing layout')
    transformed = (raw - mean) / scale
    if not np.isfinite(transformed).all():
        raise ArithmeticError('Nonfinite transformed features')
    return transformed


def nested_crossfit(core_matrix, protein_log2, y, w, clusters, outer_folds,
                    model_name='routine'):
    core, proteins = _features(core_matrix, protein_log2)
    y, w = _rows(y, w, len(core))
    clusters, outer = np.asarray(clusters), np.asarray(outer_folds)
    if (clusters.shape != (len(core),) or outer.shape != (len(core),)
            or not np.issubdtype(outer.dtype, np.integer)
            or set(outer.tolist()) != set(range(5))
            or any(not isinstance(c, str) or not c for c in clusters.tolist())
            or model_name not in MODEL_INDICES):
        raise ValueError('Invalid PSU fold layout or b2m31 model')
    for cluster in set(clusters.tolist()):
        if len(set(outer[clusters == cluster].tolist())) != 1:
            raise ValueError('PSU crosses outer fold boundaries')
    predictions = np.full(len(core), np.nan)
    records = []
    for fold in range(5):
        training, scoring = outer != fold, outer == fold
        training_psu = set(clusters[training].tolist())
        scoring_psu = set(clusters[scoring].tolist())
        if training_psu & scoring_psu:
            raise ValueError('Outer PSU leakage')
        inner = psu_partition(clusters[training], f'b2m31:inner:{fold}', folds=3)
        train_indices = np.flatnonzero(training)
        heldout = {penalty: np.full(int(training.sum()), np.nan) for penalty in PENALTIES}
        for validation in range(3):
            inner_train = np.zeros(len(core), dtype=bool)
            inner_score = np.zeros(len(core), dtype=bool)
            inner_train[train_indices[inner != validation]] = True
            inner_score[train_indices[inner == validation]] = True
            if set(clusters[inner_train].tolist()) & set(clusters[inner_score].tolist()):
                raise ValueError('Inner PSU leakage')
            prep = fit_preprocessing(core[inner_train], proteins[inner_train], w[inner_train], model_name)
            ztrain = transform_features(core[inner_train], proteins[inner_train], prep)
            zscore = transform_features(core[inner_score], proteins[inner_score], prep)
            for penalty in PENALTIES:
                fitted = fit_ridge(ztrain, y[inner_train], w[inner_train], penalty)
                heldout[penalty][inner == validation] = predict_ridge(zscore, fitted)
        scores = {}
        for penalty, p in heldout.items():
            if not np.isfinite(p).all():
                raise ArithmeticError('Missing inner heldout prediction')
            scores[penalty] = float(np.average((y[training] - p)**2, weights=w[training]))
        best = min(scores.values())
        chosen = max(p for p in PENALTIES if scores[p] <= best + 1e-12)
        prep = fit_preprocessing(core[training], proteins[training], w[training], model_name)
        fitted = fit_ridge(transform_features(core[training], proteins[training], prep),
                           y[training], w[training], chosen)
        predictions[scoring] = predict_ridge(transform_features(core[scoring], proteins[scoring], prep), fitted)
        names = list(CORE_NAMES) + [PROTEIN_NAMES[j] for j in MODEL_INDICES[model_name]]
        records.append(dict(outer_fold=fold, model=model_name, selected_penalty=float(chosen),
            **{'lambda': float(chosen)},
            inner_brier_by_penalty={str(p): v for p, v in scores.items()},
            training_n=int(training.sum()), scoring_n=int(scoring.sum()),
            training_positive_weight_n=int((training & (w > 0)).sum()),
            training_events=int(y[training].sum()), scoring_events=int(y[scoring].sum()),
            training_positive_weight_events=int((training & (w > 0) & (y == 1)).sum()),
            training_weight_sum=float(w[training].sum()),
            training_psu_n=len(training_psu), scoring_psu_n=len(scoring_psu), PSU_overlap_n=0,
            core_feature_count=11, core_feature_indices=list(range(11)),
            protein_feature_indices=list(MODEL_INDICES[model_name]), feature_names=names,
            training_mean=prep['mean'], training_scale=prep['scale'], constant_features=prep['constant'],
            core_coefficients_standardized=fitted['coefficients'][:11],
            protein_coefficients_standardized=fitted['coefficients'][11:], intercept=fitted['intercept'],
            gradient_max=fitted['gradient_max'], iterations=fitted['iterations'], objective=fitted['objective'],
            normalized_weight_sum=fitted['normalized_weight_sum'], convergence_status=fitted['convergence_status'],
            protein_transformation='reported_log2'))
    if not np.isfinite(predictions).all() or np.any((predictions < 0) | (predictions > 1)):
        raise ArithmeticError('Invalid outer heldout probabilities')
    return predictions, records
