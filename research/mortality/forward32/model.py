"""Pure early-training/late-scoring ridge fits; no late outcomes or weights.

The caller qualifies the1999/2001 training and2003 scoring domains and passes
log2 proteins ordered (Cys, B2M). All preprocessing and tuning use early rows.
This module performs no file I/O and never emits participant-level records.
"""
from __future__ import annotations

import numpy as np

from ..b2m31 import model as numerical

ARMS = ('routine', 'routine_b2m', 'routine_cys', 'routine_cys_b2m')
INNER_SALT = 'forward32:inner:v1'
SCHEMA = 'forward32-fit/1'


def _clusters(values, n):
    values = np.asarray(values)
    if values.shape != (n,) or any(not isinstance(c, str) or not c for c in values.tolist()):
        raise ValueError('Aligned nonempty string PSU identities required')
    return values


def _inner_support(y, w, clusters, folds):
    rows = []
    for fold in range(3):
        training, scoring = folds != fold, folds == fold
        train_psu = set(clusters[training].tolist())
        score_psu = set(clusters[scoring].tolist())
        if train_psu & score_psu:
            raise ValueError('PSU crosses inner training/scoring boundary')
        rows.append(dict(fold=fold, training_n=int(training.sum()), scoring_n=int(scoring.sum()),
            training_events=int(y[training].sum()), scoring_events=int(y[scoring].sum()),
            training_positive_weight_n=int((training & (w > 0)).sum()),
            scoring_positive_weight_n=int((scoring & (w > 0)).sum()),
            training_positive_weight_events=int((training & (w > 0) & (y == 1)).sum()),
            scoring_positive_weight_events=int((scoring & (w > 0) & (y == 1)).sum()),
            training_weight_sum=float(w[training].sum()), scoring_weight_sum=float(w[scoring].sum()),
            training_psu_n=len(train_psu), scoring_psu_n=len(score_psu), PSU_overlap_n=0))
    return rows


def _fit_arm(core, proteins, y, w, folds, arm):
    heldout = {penalty: np.full(len(y), np.nan) for penalty in numerical.PENALTIES}
    diagnostics = []
    for fold in range(3):
        training, scoring = folds != fold, folds == fold
        prep = numerical.fit_preprocessing(core[training], proteins[training], w[training], arm)
        ztrain = numerical.transform_features(core[training], proteins[training], prep)
        zscore = numerical.transform_features(core[scoring], proteins[scoring], prep)
        for penalty in numerical.PENALTIES:
            fitted = numerical.fit_ridge(ztrain, y[training], w[training], penalty)
            heldout[penalty][scoring] = numerical.predict_ridge(zscore, fitted)
            diagnostics.append(dict(fold=fold, penalty=float(penalty),
                gradient_max=fitted['gradient_max'], iterations=fitted['iterations'],
                convergence_status=fitted['convergence_status'],
                normalized_weight_sum=fitted['normalized_weight_sum']))
    scores = {}
    for penalty, p in heldout.items():
        if not np.isfinite(p).all():
            raise ArithmeticError('Missing early inner-heldout probability')
        scores[penalty] = float(np.average((y - p)**2, weights=w))
    minimum = min(scores.values())
    chosen = max(p for p in numerical.PENALTIES if scores[p] <= minimum + 1e-12)
    prep = numerical.fit_preprocessing(core, proteins, w, arm)
    fitted = numerical.fit_ridge(numerical.transform_features(core, proteins, prep), y, w, chosen)
    return dict(preprocessing=prep, fit=fitted,
        feature_names=list(numerical.CORE_NAMES) + [numerical.PROTEIN_NAMES[i] for i in numerical.MODEL_INDICES[arm]],
        inner_brier_by_penalty={str(p): v for p, v in scores.items()}, selected_penalty=float(chosen),
        inner_fit_diagnostics=diagnostics)


def predict_forward(core_late, protein_late, artifact):
    """Score with fixed training parameters only; return arrays in memory.

    The caller retains the qualified split/provenance. This helper does not
    qualify an external assay or accept outcomes, weights, calibration or fit.
    """
    core, proteins = numerical._features(core_late, protein_late)
    if (not isinstance(artifact, dict) or artifact.get('schema') != SCHEMA
            or not isinstance(artifact.get('models'), dict) or set(artifact['models']) != set(ARMS)):
        raise ValueError('Complete forward32 training artifact required')
    predictions = {}
    for arm in ARMS:
        row = artifact['models'][arm]
        if (not isinstance(row, dict) or not isinstance(row.get('preprocessing'), dict)
                or row['preprocessing'].get('model_name') != arm or not isinstance(row.get('fit'), dict)
                or row.get('selected_penalty') != row['fit'].get('penalty')):
            raise ValueError('Invalid forward32 model parameters')
        transformed = numerical.transform_features(core, proteins, row['preprocessing'])
        predictions[arm] = numerical.predict_ridge(transformed, row['fit'])
    return predictions


def _validated_split(core_train, protein_train, y_train, w_train, clusters_train,
                     core_late, protein_late, clusters_late):
    core, proteins = numerical._features(core_train, protein_train)
    late_core, late_proteins = numerical._features(core_late, protein_late)
    y, w = numerical._rows(y_train, w_train, len(core))
    train_clusters = _clusters(clusters_train, len(core))
    late_clusters = _clusters(clusters_late, len(late_core))
    training_psu, scoring_psu = set(train_clusters.tolist()), set(late_clusters.tolist())
    if training_psu & scoring_psu:
        raise ValueError('Training and late domain share PSU identities')
    if np.unique(y[w > 0]).size != 2:
        raise ValueError('Early training needs positive-weight events and nonevents')
    folds = numerical.psu_partition(train_clusters, INNER_SALT, folds=3)
    return core, proteins, y, w, train_clusters, late_core, late_proteins, late_clusters, folds


def fit_training_arm(core_train, protein_train, y_train, w_train, clusters_train,
                     core_late, protein_late, clusters_late, arm):
    """Fit one arm with the identical forward guards for paired refit bootstrap.

    No late outcomes/weights enter this API. An arm-local fitting failure is
    raised for the caller to record without discarding another defined pair.
    """
    if arm not in ARMS:
        raise ValueError('Unregistered forward32 model')
    core, proteins, y, w, _, late_core, late_proteins, _, folds = _validated_split(
        core_train, protein_train, y_train, w_train, clusters_train,
        core_late, protein_late, clusters_late)
    row = _fit_arm(core, proteins, y, w, folds, arm)
    p = numerical.predict_ridge(numerical.transform_features(late_core, late_proteins,
                                                           row['preprocessing']), row['fit'])
    return p, row


def fit_forward(core_train, protein_train, y_train, w_train, clusters_train,
                core_late, protein_late, clusters_late):
    """Tune/refit on early data; score the disjoint late domain in memory.

    Returns ``(late_predictions, artifact)``. No late outcomes or late weights
    are accepted. The caller must not publish returned participant vectors.
    """
    core, proteins, y, w, train_clusters, late_core, late_proteins, late_clusters, folds = _validated_split(
        core_train, protein_train, y_train, w_train, clusters_train,
        core_late, protein_late, clusters_late)
    training_psu, scoring_psu = set(train_clusters.tolist()), set(late_clusters.tolist())
    support = _inner_support(y, w, train_clusters, folds)
    models = {arm: _fit_arm(core, proteins, y, w, folds, arm) for arm in ARMS}
    artifact = dict(schema=SCHEMA, models=models, inner_folds_aggregate=support,
        training=dict(n=len(y), events8=int(y.sum()), positive_weight_n=int((w > 0).sum()),
            positive_weight_events8=int(((w > 0) & (y == 1)).sum()),
            weight_sum=float(w.sum()), psu_n=len(training_psu), weight_name='WTSSCB4Y',
            cycles=[1999, 2001], inner_folds=3, inner_salt=INNER_SALT,
            weight_normalization='sum_to_raw_training_n'),
        late_domain=dict(n=len(late_core), psu_n=len(scoring_psu), cycle=2003,
            evaluation_weight_name='WTSSCB2Y', training_psu_overlap_n=0,
            outcomes_used_for_fitting=False, weights_used_for_fitting=False),
        fit_scope='early_training_only', selection_scope='early_three_fold_PSU_inner_Brier',
        no_outer_cross_validation=True, late_feature_preprocessing_fitted=False,
        independent_external_validation=False, participant_records_exported=False)
    return predict_forward(late_core, late_proteins, artifact), artifact
