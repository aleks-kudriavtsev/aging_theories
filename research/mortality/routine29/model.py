"""New routine-core comparison using preserved binary28 numerical primitives.

The legacy module and its registries are never mutated. All routine/core and
protein transforms are learned separately inside each actual training split.
"""
from __future__ import annotations
import numpy as np
from ..binary28.model import (
    PENALTIES, PROTEIN_NAMES, _matrix, _rows, _weighted_support,
    logistic_objective, fit_ridge, predict_ridge, psu_partition,
)

CORE_NAMES = ('age_centered', 'age_centered_squared', 'male', 'smoking_current',
              'smoking_former', 'systolic_bp', 'bmi', 'hba1c', 'egfr_low',
              'egfr_high', 'log2_uacr')
MODEL_INDICES = {'simple': (), 'simple_all3': (0, 1, 2), 'routine': (),
                 'routine_all3': (0, 1, 2), 'rank_all3': (0, 1, 2),
                 'clip_all3': (0, 1, 2)}
MODEL_CORE_INDICES = {name: tuple(range(6 if name.startswith('simple') else 11))
                      for name in MODEL_INDICES}


def _features(core, proteins):
    core, proteins = _matrix(core), _matrix(proteins)
    if core.shape[1] != 11 or proteins.shape != (len(core), 3):
        raise ValueError('Eleven routine and three protein features required')
    return core, proteins


def fit_preprocessing(core, protein_log2, weights, model_name):
    core, proteins = _features(core, protein_log2)
    _, w = _rows(np.zeros(len(core)), weights, len(core))
    if model_name not in MODEL_INDICES:
        raise ValueError('Unregistered routine29 model')
    transforms = []
    for j in MODEL_INDICES[model_name]:
        if model_name == 'rank_all3':
            values, _, cdf = _weighted_support(proteins[:, j], w)
            transforms.append({'column': j, 'kind': 'right_continuous_weighted_ecdf',
                               'values': values.tolist(), 'cdf': cdf.tolist()})
        elif model_name == 'clip_all3':
            values, mass, _ = _weighted_support(proteins[:, j], w)
            midpoint = (np.cumsum(mass) - .5 * mass) / mass.sum()
            lower, upper = np.interp([.01, .99], midpoint, values)
            transforms.append({'column': j, 'kind': 'weighted_midpoint_tail_clip',
                               'lower': float(lower), 'upper': float(upper)})
        else:
            transforms.append({'column': j, 'kind': 'identity_log2'})
    specification = {'model_name': model_name,
                     'core_indices': list(MODEL_CORE_INDICES[model_name]),
                     'protein_transforms': transforms}
    raw = _unscaled(core, proteins, specification)
    mean = np.average(raw, axis=0, weights=w)
    constant = np.ptp(raw[w > 0], axis=0) == 0
    mean[constant] = raw[w > 0][0, constant]
    scale = np.sqrt(np.average((raw - mean)**2, axis=0, weights=w))
    scale[constant] = 1.
    if not np.isfinite(mean).all() or not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ArithmeticError('Invalid training feature scaling')
    specification.update(mean=mean.tolist(), scale=scale.tolist(), constant=constant.tolist())
    return specification


def _unscaled(core, proteins, specification):
    name = specification.get('model_name')
    if name not in MODEL_CORE_INDICES or specification.get('core_indices') != list(MODEL_CORE_INDICES[name]):
        raise ValueError('Wrong registered core feature selection')
    columns = [core[:, MODEL_CORE_INDICES[name]]]
    for transform in specification['protein_transforms']:
        values = proteins[:, transform['column']]
        if transform['kind'] == 'right_continuous_weighted_ecdf':
            support, cdf = np.asarray(transform['values']), np.asarray(transform['cdf'])
            index = np.searchsorted(support, values, side='right') - 1
            transformed = np.zeros(len(values)); valid = index >= 0
            transformed[valid] = cdf[index[valid]]
        elif transform['kind'] == 'weighted_midpoint_tail_clip':
            transformed = np.clip(values, transform['lower'], transform['upper'])
        elif transform['kind'] == 'identity_log2':
            transformed = values
        else:
            raise ValueError('Unknown protein transformation')
        columns.append(transformed[:, None])
    return np.column_stack(columns)


def transform_features(core, protein_log2, specification):
    core, proteins = _features(core, protein_log2)
    raw = _unscaled(core, proteins, specification)
    mean, scale = np.asarray(specification['mean']), np.asarray(specification['scale'])
    if mean.shape != (raw.shape[1],) or scale.shape != mean.shape or np.any(scale <= 0):
        raise ValueError('Incorrect preprocessing layout')
    result = (raw - mean) / scale
    if not np.isfinite(result).all():
        raise ArithmeticError('Nonfinite transformed features')
    return result


def _fit_predict(core, proteins, y, weights, training, scoring, name, penalty):
    preprocessing = fit_preprocessing(core[training], proteins[training], weights[training], name)
    fit = fit_ridge(transform_features(core[training], proteins[training], preprocessing),
                    y[training], weights[training], penalty)
    prediction = predict_ridge(transform_features(core[scoring], proteins[scoring], preprocessing), fit)
    return prediction, fit, preprocessing


def nested_crossfit(core_matrix, protein_log2, y, w, clusters, outer_folds, model_name='routine'):
    core, proteins = _features(core_matrix, protein_log2)
    y, w = _rows(y, w, len(core))
    clusters, outer = np.asarray(clusters), np.asarray(outer_folds)
    if (clusters.shape != (len(core),) or outer.shape != (len(core),)
            or not np.issubdtype(outer.dtype, np.integer)
            or set(outer.tolist()) != set(range(5))
            or any(not isinstance(c, str) or not c for c in clusters.tolist())
            or model_name not in MODEL_INDICES):
        raise ValueError('Invalid PSU fold layout or routine29 model')
    for cluster in set(clusters.tolist()):
        if len(set(outer[clusters == cluster].tolist())) != 1:
            raise ValueError('PSU crosses outer fold boundaries')
    predictions = np.full(len(core), np.nan); records = []
    for fold in range(5):
        training, scoring = outer != fold, outer == fold
        training_psu, scoring_psu = set(clusters[training].tolist()), set(clusters[scoring].tolist())
        if training_psu & scoring_psu:
            raise ValueError('Outer PSU leakage')
        inner = psu_partition(clusters[training], f'routine29:inner:{fold}', folds=3)
        train_indices = np.flatnonzero(training)
        heldout = {penalty: np.full(int(training.sum()), np.nan) for penalty in PENALTIES}
        for validation in range(3):
            inner_train = np.zeros(len(core), dtype=bool); inner_score = np.zeros(len(core), dtype=bool)
            inner_train[train_indices[inner != validation]] = True
            inner_score[train_indices[inner == validation]] = True
            if set(clusters[inner_train].tolist()) & set(clusters[inner_score].tolist()):
                raise ValueError('Inner PSU leakage')
            preprocessing = fit_preprocessing(core[inner_train], proteins[inner_train], w[inner_train], model_name)
            ztrain = transform_features(core[inner_train], proteins[inner_train], preprocessing)
            zscore = transform_features(core[inner_score], proteins[inner_score], preprocessing)
            for penalty in PENALTIES:
                fit = fit_ridge(ztrain, y[inner_train], w[inner_train], penalty)
                heldout[penalty][inner == validation] = predict_ridge(zscore, fit)
        scores = {}
        for penalty, prediction in heldout.items():
            if not np.isfinite(prediction).all():
                raise ArithmeticError('Missing inner heldout prediction')
            scores[penalty] = float(np.average((y[training] - prediction)**2, weights=w[training]))
        best = min(scores.values())
        chosen = max(p for p in PENALTIES if scores[p] <= best + 1e-12)
        prediction, fit, preprocessing = _fit_predict(core, proteins, y, w, training, scoring, model_name, chosen)
        predictions[scoring] = prediction
        transforms = []
        for specification in preprocessing['protein_transforms']:
            diagnostic = {k: v for k, v in specification.items() if k not in ('values', 'cdf')}
            if 'values' in specification:
                diagnostic['training_support_count'] = len(specification['values'])
            transforms.append(diagnostic)
        core_indices = MODEL_CORE_INDICES[model_name]; count = len(core_indices)
        feature_names = [CORE_NAMES[j] for j in core_indices] + [PROTEIN_NAMES[j] for j in MODEL_INDICES[model_name]]
        records.append(dict(outer_fold=fold, model=model_name, selected_penalty=chosen, **{'lambda': chosen},
            inner_brier_by_penalty={str(p): score for p, score in scores.items()},
            training_n=int(training.sum()), scoring_n=int(scoring.sum()),
            training_positive_weight_n=int((training & (w > 0)).sum()),
            training_events=int(y[training].sum()), scoring_events=int(y[scoring].sum()),
            training_positive_weight_events=int(((y == 1) & training & (w > 0)).sum()),
            training_weight_sum=float(w[training].sum()),
            training_psu_n=len(training_psu), scoring_psu_n=len(scoring_psu), PSU_overlap_n=0,
            core_feature_count=count, core_feature_indices=list(core_indices), feature_names=feature_names,
            training_mean=preprocessing['mean'], training_scale=preprocessing['scale'],
            constant_features=preprocessing['constant'], protein_transform_parameters=transforms,
            core_coefficients_standardized=fit['coefficients'][:count],
            protein_coefficients_standardized=fit['coefficients'][count:], intercept=fit['intercept'],
            gradient_max=fit['gradient_max'], iterations=fit['iterations'], objective=fit['objective'],
            normalized_weight_sum=fit['normalized_weight_sum'], convergence_status=fit['convergence_status'],
            protein_transformation='right_continuous_weighted_ecdf' if model_name == 'rank_all3'
                else 'weighted_midpoint_tail_clip_1_99' if model_name == 'clip_all3' else 'reported_log2'))
    if not np.isfinite(predictions).all() or np.any((predictions < 0) | (predictions > 1)):
        raise ArithmeticError('Invalid outer heldout probabilities')
    return predictions, records
