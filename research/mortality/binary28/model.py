"""Nested whole-PSU assessment of new weighted binary logistic models.

This module does not load participant files or use a historical mortality engine.
All feature processing and ridge selection are fitted within each training split.
"""
from __future__ import annotations
import hashlib
import numpy as np
from scipy.special import expit

PENALTIES = (1.0, 10.0, 100.0)
CORE_NAMES = ('age_centered', 'age_centered_squared', 'male', 'smoking_current',
              'smoking_former', 'systolic_bp')
PROTEIN_NAMES = ('log2_cystatin_c', 'log2_beta2_microglobulin', 'log2_beta_trace_protein')
MODEL_INDICES = {'core': (), 'all3': (0, 1, 2), 'cys': (0,), 'b2m': (1,),
                 'btp': (2,), 'rank_all3': (0, 1, 2), 'clip_all3': (0, 1, 2)}


def _matrix(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or not len(x) or not np.isfinite(x).all():
        raise ValueError('Finite nonempty feature matrix required')
    return x


def _rows(y, w, n):
    y, w = np.asarray(y, float), np.asarray(w, float)
    if (y.shape != (n,) or w.shape != (n,) or not np.isfinite(y).all()
            or not np.isfinite(w).all() or not np.isin(y, [0, 1]).all()
            or np.any(w < 0) or not np.isfinite(w.sum()) or w.sum() <= 0):
        raise ValueError('Invalid binary outcomes or weights')
    return y, w


def logistic_objective(beta, design, y, normalized_weights, penalty):
    """Finite weighted Bernoulli negative log likelihood, gradient and Hessian.

    The first design column is the intercept and is not penalized.
    The caller supplies weights normalized to the number of training records.
    """
    beta = np.asarray(beta, float)
    x, y, w = np.asarray(design, float), np.asarray(y, float), np.asarray(normalized_weights, float)
    eta = x @ beta
    p = expit(eta)
    penalized = beta.copy(); penalized[0] = 0.
    objective = float(w @ (np.logaddexp(0., eta) - y * eta)
                      + .5 * penalty * (penalized @ penalized))
    gradient = x.T @ (w * (p - y)) + penalty * penalized
    hessian = x.T @ ((w * p * (1. - p))[:, None] * x)
    hessian.flat[::len(beta)+1] += np.r_[0., np.full(len(beta)-1, penalty)]
    if not np.isfinite(objective) or not np.isfinite(gradient).all() or not np.isfinite(hessian).all():
        raise ArithmeticError('Nonfinite logistic objective, gradient or Hessian')
    return objective, gradient, hessian


def fit_ridge(features, y, weights, penalty):
    """Newton fit on already transformed features; intercept is unpenalized."""
    z = _matrix(features); y, w = _rows(y, weights, len(z))
    if type(penalty) not in (float, int) or not np.isfinite(penalty) or penalty <= 0:
        raise ValueError('Positive fixed ridge penalty required')
    positive = w > 0
    if np.unique(y[positive]).size != 2:
        raise ValueError('Training needs positive-weight events and non-events')
    w = w * (len(z) / w.sum())
    x = np.column_stack([np.ones(len(z)), z])
    observed = float(w @ y / w.sum())
    beta = np.zeros(x.shape[1]); beta[0] = np.log(observed / (1. - observed))
    objective, gradient, hessian = logistic_objective(beta, x, y, w, penalty)
    iterations = 0
    tiny_step = False
    for iteration in range(100):
        if np.max(np.abs(gradient)) <= 1e-6:
            break
        try:
            direction = np.linalg.solve(hessian, gradient)
        except np.linalg.LinAlgError as error:
            raise ArithmeticError('Singular logistic Hessian') from error
        directional = float(gradient @ direction)
        if not np.isfinite(direction).all() or not np.isfinite(directional) or directional <= 0:
            raise ArithmeticError('Invalid logistic Newton direction')
        if np.max(np.abs(direction)) <= 1e-10 and np.max(np.abs(gradient)) <= 1e-4:
            tiny_step = True
            break
        step = 1.
        accepted = False
        for _ in range(40):
            candidate = beta - step * direction
            if not np.isfinite(candidate).all():
                step *= .5; continue
            new_objective, new_gradient, new_hessian = logistic_objective(candidate, x, y, w, penalty)
            tolerance = 1e-12 * max(1., abs(objective))
            if new_objective <= objective - 1e-4 * step * directional + tolerance:
                beta, objective, gradient, hessian = candidate, new_objective, new_gradient, new_hessian
                accepted = True
                break
            step *= .5
        if not accepted:
            raise ArithmeticError('Logistic line search did not converge')
        iterations = iteration + 1
    if (not np.isfinite(beta).all() or not np.isfinite(objective)
            or not np.isfinite(gradient).all()
            or not (np.max(np.abs(gradient)) <= 1e-6 or (tiny_step and np.max(np.abs(gradient)) <= 1e-4))):
        raise ArithmeticError('Logistic gradient convergence criterion failed')
    return {'intercept': float(beta[0]), 'coefficients': beta[1:].tolist(),
            'gradient_max': float(np.max(np.abs(gradient))),
            'objective': objective, 'iterations': iterations, 'penalty': float(penalty),
            'convergence_status': 'tiny_step_gradient_1e-4' if tiny_step else 'gradient_1e-6',
            'normalized_weight_sum': float(w.sum())}


def predict_ridge(features, fit):
    x = _matrix(features)
    coefficients = np.asarray(fit['coefficients'], float)
    if coefficients.shape != (x.shape[1],) or not np.isfinite(coefficients).all():
        raise ValueError('Invalid fitted coefficient layout')
    if not np.isfinite(fit['intercept']):
        raise ValueError('Nonfinite fitted intercept')
    eta = fit['intercept'] + x @ coefficients
    p = expit(eta)
    if not np.isfinite(p).all():
        raise ArithmeticError('Nonfinite fitted probability')
    return p


def _weighted_support(values, weights):
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    positive = weights > 0
    if not positive.any():
        raise ValueError('Positive training weight required')
    unique, inverse = np.unique(values[positive], return_inverse=True)
    mass = np.bincount(inverse, weights=weights[positive], minlength=len(unique))
    return unique, mass, np.cumsum(mass) / mass.sum()


def fit_preprocessing(core, protein_log2, weights, model_name):
    core, proteins = _matrix(core), _matrix(protein_log2)
    if core.shape[1] != 6 or proteins.shape != (len(core), 3):
        raise ValueError('Six core and three protein features required')
    _, w = _rows(np.zeros(len(core)), weights, len(core))
    if model_name not in MODEL_INDICES:
        raise ValueError('Unregistered model')
    transforms = []
    for j in MODEL_INDICES[model_name]:
        if model_name == 'rank_all3':
            values, mass, cdf = _weighted_support(proteins[:, j], w)
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
    preliminary = {'model_name': model_name, 'protein_transforms': transforms}
    raw = _unscaled(core, proteins, preliminary)
    mean = np.average(raw, axis=0, weights=w)
    constant = np.ptp(raw[w > 0], axis=0) == 0
    mean[constant] = raw[w > 0][0, constant]
    scale = np.sqrt(np.average((raw - mean)**2, axis=0, weights=w))
    scale[constant] = 1.
    if not np.isfinite(mean).all() or not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ArithmeticError('Invalid training feature scaling')
    preliminary.update(mean=mean.tolist(), scale=scale.tolist(), constant=constant.tolist())
    return preliminary


def _unscaled(core, proteins, specification):
    columns = [core]
    for transform in specification['protein_transforms']:
        values = proteins[:, transform['column']]
        if transform['kind'] == 'right_continuous_weighted_ecdf':
            support = np.asarray(transform['values']); cdf = np.asarray(transform['cdf'])
            index = np.searchsorted(support, values, side='right') - 1
            transformed = np.zeros(len(values))
            valid = index >= 0
            transformed[valid] = cdf[index[valid]]
        elif transform['kind'] == 'weighted_midpoint_tail_clip':
            transformed = np.clip(values, transform['lower'], transform['upper'])
        elif transform['kind'] == 'identity_log2':
            transformed = values
        else:
            raise ValueError('Unknown feature transformation')
        columns.append(transformed[:, None])
    return np.column_stack(columns)


def transform_features(core, protein_log2, specification):
    core, proteins = _matrix(core), _matrix(protein_log2)
    if core.shape[1] != 6 or proteins.shape != (len(core), 3):
        raise ValueError('Six core and three protein features required')
    raw = _unscaled(core, proteins, specification)
    result = (raw - np.asarray(specification['mean'])) / np.asarray(specification['scale'])
    if not np.isfinite(result).all():
        raise ArithmeticError('Nonfinite transformed features')
    return result


def psu_partition(clusters, salt, folds=3):
    clusters = np.asarray(clusters)
    if (clusters.ndim != 1 or not len(clusters) or not isinstance(salt, str)
            or any(not isinstance(c, str) or not c for c in clusters.tolist())
            or type(folds) is not int or folds < 2):
        raise ValueError('Nonempty string PSU identities and valid fold count required')
    keys = sorted(set(clusters.tolist()),
                  key=lambda key: hashlib.sha256((salt + ':' + key).encode()).hexdigest())
    if len(keys) < folds:
        raise ValueError('Insufficient distinct PSUs')
    mapping = {key: rank % folds for rank, key in enumerate(keys)}
    return np.asarray([mapping[key] for key in clusters.tolist()], dtype=int)


def _fit_predict(core, proteins, y, weights, training, scoring, name, penalty):
    preprocessing = fit_preprocessing(core[training], proteins[training], weights[training], name)
    ztrain = transform_features(core[training], proteins[training], preprocessing)
    fit = fit_ridge(ztrain, y[training], weights[training], penalty)
    zscore = transform_features(core[scoring], proteins[scoring], preprocessing)
    return predict_ridge(zscore, fit), fit, preprocessing


def nested_crossfit(core_matrix, protein_log2, y, w, clusters, outer_folds, model_name='core'):
    core, proteins = _matrix(core_matrix), _matrix(protein_log2)
    if core.shape[1] != 6 or proteins.shape != (len(core), 3):
        raise ValueError('Six core and three protein features required')
    y, w = _rows(y, w, len(core))
    clusters, outer = np.asarray(clusters), np.asarray(outer_folds)
    if (clusters.shape != (len(core),) or outer.shape != (len(core),)
            or not np.issubdtype(outer.dtype, np.integer) or set(outer.tolist()) != set(range(5))
            or any(not isinstance(c, str) or not c for c in clusters.tolist())
            or model_name not in MODEL_INDICES):
        raise ValueError('Invalid PSU fold layout or model')
    for cluster in set(clusters.tolist()):
        if len(set(outer[clusters == cluster].tolist())) != 1:
            raise ValueError('PSU crosses outer fold boundaries')
    predictions = np.full(len(core), np.nan)
    records = []
    for fold in range(5):
        training, scoring = outer != fold, outer == fold
        training_psu = set(clusters[training].tolist()); scoring_psu = set(clusters[scoring].tolist())
        if training_psu & scoring_psu:
            raise ValueError('Outer PSU leakage')
        inner = psu_partition(clusters[training], f'binary28:inner:{fold}', folds=3)
        train_indices = np.flatnonzero(training)
        inner_predictions = {penalty: np.full(int(training.sum()), np.nan) for penalty in PENALTIES}
        for validation in range(3):
            inner_train = np.zeros(len(core), dtype=bool); inner_score = np.zeros(len(core), dtype=bool)
            inner_train[train_indices[inner != validation]] = True
            inner_score[train_indices[inner == validation]] = True
            if set(clusters[inner_train].tolist()) & set(clusters[inner_score].tolist()):
                raise ValueError('Inner PSU leakage')
            # Feature transformations are independent of ridge lambda. Cache
            # them within this training split without sharing across folds.
            preprocessing = fit_preprocessing(core[inner_train], proteins[inner_train], w[inner_train], model_name)
            ztrain = transform_features(core[inner_train], proteins[inner_train], preprocessing)
            zscore = transform_features(core[inner_score], proteins[inner_score], preprocessing)
            for penalty in PENALTIES:
                fit = fit_ridge(ztrain, y[inner_train], w[inner_train], penalty)
                inner_predictions[penalty][inner == validation] = predict_ridge(zscore, fit)
        scores = {}
        for penalty, inner_prediction in inner_predictions.items():
            if not np.isfinite(inner_prediction).all():
                raise ArithmeticError('Missing inner heldout prediction')
            scores[penalty] = float(np.average((y[training] - inner_prediction)**2, weights=w[training]))
        minimum = min(scores.values())
        chosen = max(p for p in PENALTIES if scores[p] <= minimum + 1e-12)
        prediction, fit, preprocessing = _fit_predict(core, proteins, y, w, training, scoring, model_name, chosen)
        predictions[scoring] = prediction
        coefficients = fit['coefficients']
        transform_diagnostics = []
        for specification in preprocessing['protein_transforms']:
            diagnostic = {k: v for k, v in specification.items() if k not in ('values', 'cdf')}
            if 'values' in specification:
                diagnostic['training_support_count'] = len(specification['values'])
            transform_diagnostics.append(diagnostic)
        names = list(CORE_NAMES) + [PROTEIN_NAMES[j] for j in MODEL_INDICES[model_name]]
        records.append(dict(outer_fold=fold, model=model_name, selected_penalty=chosen, **{'lambda': chosen},
            inner_brier_by_penalty={str(p): v for p, v in scores.items()},
            training_n=int(training.sum()), scoring_n=int(scoring.sum()),
            training_positive_weight_n=int((training & (w > 0)).sum()),
            training_events=int(y[training].sum()), scoring_events=int(y[scoring].sum()),
            training_positive_weight_events=int(((y == 1) & training & (w > 0)).sum()),
            training_weight_sum=float(w[training].sum()),
            training_psu_n=len(training_psu), scoring_psu_n=len(scoring_psu),
            PSU_overlap_n=0, feature_names=names, training_mean=preprocessing['mean'],
            training_scale=preprocessing['scale'], constant_features=preprocessing['constant'],
            protein_transform_parameters=transform_diagnostics,
            core_coefficients_standardized=coefficients[:6], protein_coefficients_standardized=coefficients[6:],
            intercept=fit['intercept'], gradient_max=fit['gradient_max'], iterations=fit['iterations'],
            objective=fit['objective'], normalized_weight_sum=fit['normalized_weight_sum'],
            convergence_status=fit['convergence_status'],
            protein_transformation='right_continuous_weighted_ecdf' if model_name == 'rank_all3'
                else 'weighted_midpoint_tail_clip_1_99' if model_name == 'clip_all3' else 'reported_log2'))
    if not np.isfinite(predictions).all() or np.any((predictions < 0) | (predictions > 1)):
        raise ArithmeticError('Invalid outer heldout probabilities')
    return predictions, records
