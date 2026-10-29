"""Separate exploratory eight-year binary CIF comparisons, never summed.

Both causes and every arm share one fixed domain, partitions and each full
eligible-frame PSU draw. Only point parameters and aggregate scores leave RAM.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import scipy

from ..b2m31 import model
from ..b2m31.analyze import (
    ARMS, CONTRASTS, METRICS, assign_outer, _survey_design,
    design_bootstrap_weights, flatten_scores, intervals,
)
from ..b2m31.features import CORE
from ..binary28.metrics import metrics, calibration_bins

ROOT = Path(__file__).parent
CAUSES = ('002', '001')
PROTEINS = ('SSCYST', 'SSB2M')
REPLICATES = 500
SEED = 20261105
PLAN_COMMIT = 'd772b874198a9f7739225bbce5d1b6b6dfcd0c08'


def _bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf-8')


def _sha(value):
    return hashlib.sha256(value).hexdigest()


def _weight_sha(weights):
    """Aggregate reproducibility fingerprint, never an exported weight vector."""
    return _sha(np.asarray(weights, dtype='<f8').tobytes(order='C'))


def _validate_arrays(core, proteins, targets, weights, clusters, folds, eligible, common):
    core, proteins = model._features(core, proteins)
    if not isinstance(targets, dict) or set(targets) != set(CAUSES):
        raise ValueError('Exactly the two locked cause targets required')
    checked = {}
    for cause in CAUSES:
        checked[cause], w = model._rows(targets[cause], weights, len(core))
        if np.unique(checked[cause][w > 0]).size != 2:
            raise ValueError('Each qualified cause requires events and known noncases')
    if np.any(checked['002'] + checked['001'] > 1):
        raise ValueError('Underlying-cause event indicators must be disjoint')
    clusters, folds, common = np.asarray(clusters), np.asarray(folds), np.asarray(common)
    if (clusters.shape != (len(core),) or any(not isinstance(c, str) or not c for c in clusters.tolist())
            or folds.shape != clusters.shape or not np.issubdtype(folds.dtype, np.integer)
            or set(folds.tolist()) != set(range(5))
            or common.shape != (len(eligible),) or common.dtype != bool):
        raise ValueError('Aligned whole-PSU folds and boolean common domain required')
    _survey_design(eligible)
    if (int(common.sum()) != len(core)
            or not np.array_equal(eligible.weight.to_numpy(float)[common], w)
            or not np.array_equal(eligible.cluster.to_numpy()[common], clusters)
            or set(clusters) != set(eligible.cluster)
            or not np.array_equal(assign_outer(eligible.cluster.to_numpy())[common], folds)):
        raise ValueError('Shared weights, rows, eligible inventory or fixed31 folds changed')
    for cluster in set(clusters):
        if len(set(folds[clusters == cluster])) != 1:
            raise ValueError('PSU crosses an outer fold')
    return core, proteins, checked, w, clusters, folds, common


def _fit_score(core, proteins, y, w, clusters, folds, cause):
    points, statuses, parameters, predictions = {}, {}, [], {}
    for arm in ARMS:
        try:
            p, rows = model.nested_crossfit(core, proteins, y, w, clusters, folds, arm)
        except (ValueError, ArithmeticError, RuntimeError) as error:
            statuses[arm] = dict(status='undefined', error='fit|' + type(error).__name__ + ': ' + str(error))
            continue
        parameters.extend(dict(cause=cause, **row) for row in rows)
        try:
            point = metrics(y, p, w)
        except (ValueError, ArithmeticError, RuntimeError) as error:
            statuses[arm] = dict(status='undefined', error='metric|' + type(error).__name__ + ': ' + str(error))
            continue
        predictions[arm], points[arm], statuses[arm] = p, point, dict(status='computed')
    return points, statuses, parameters, predictions


def _check_shared_preprocessing(parameters):
    """Outcome changes cannot change transforms of identical X/w/folds."""
    previous = {}
    invariant = ('training_mean', 'training_scale', 'constant_features', 'feature_names',
                 'core_feature_indices', 'protein_feature_indices', 'training_n', 'scoring_n',
                 'training_positive_weight_n', 'training_weight_sum', 'training_psu_n',
                 'scoring_psu_n', 'PSU_overlap_n')
    for row in parameters:
        key = row['model'], row['outer_fold']
        values = {field: row[field] for field in invariant}
        if key in previous and values != previous[key]:
            raise ArithmeticError('Identical X/weights/folds have cause-dependent preprocessing')
        previous[key] = values


def _resampling_summary(eligible_weights, weights, targets, clusters):
    positive = weights > 0
    result = dict(eligible_weight_sum=float(eligible_weights.sum()),
        eligible_weight_squared_sum=float(eligible_weights @ eligible_weights),
        eligible_positive_weight_rows=int((eligible_weights > 0).sum()),
        common_weight_sum=float(weights.sum()), common_weight_squared_sum=float(weights @ weights),
        common_positive_weight_rows=int(positive.sum()), common_positive_weight_PSUs=len(set(clusters[positive])),
        shared_common_weight_sha256=_weight_sha(weights),
        shared_eligible_weight_sha256=_weight_sha(eligible_weights))
    for cause in CAUSES:
        y = targets[cause]
        result[cause + '|event_weight_sum'] = float(weights @ y)
        result[cause + '|positive_weight_event_rows'] = int((positive & (y == 1)).sum())
        result[cause + '|positive_weight_noncase_rows'] = int((positive & (y == 0)).sum())
    return result


def _failure_summary(draws):
    model_failures, contrast_failures, errors, calibration = {}, {}, {}, {}
    for arm in ARMS:
        model_failures[arm] = sum(row[arm + '|status'] != 'computed' for row in draws)
        errors[arm] = dict(Counter(row[arm + '|error'] for row in draws if row[arm + '|error'] is not None))
        calibration[arm] = dict(Counter(row[arm + '|calibration_status'] or 'predictive_fit_or_metric_unavailable'
                                       for row in draws if row[arm + '|calibration_status'] != 'computed'))
    for contrast in CONTRASTS:
        contrast_failures[contrast] = sum(row[contrast + '|status'] != 'computed' for row in draws)
    return dict(model_failure_counts=model_failures, contrast_failure_counts=contrast_failures,
                fit_or_metric_error_reasons=errors, calibration_undefined_reasons=calibration)


def analyze_arrays(core, proteins, targets, weights, clusters, folds, eligible, common,
                   replicates=REPLICATES, seed=SEED):
    """Testable numerical orchestration; inputs are already qualified arrays.

    Returns only aggregate dictionaries/rows. It never loads data or writes
    participant outcomes, probabilities, bootstrap coefficients or raw weights.
    """
    if type(replicates) is not int or replicates <= 0 or type(seed) is not int:
        raise ValueError('Positive fixed replicate count and integer seed required')
    core, proteins, targets, w, clusters, folds, common = _validate_arrays(
        core, proteins, targets, weights, clusters, folds, eligible, common)
    point, statuses, contrasts, parameters, bins = {}, {}, {}, [], []
    for cause in CAUSES:
        points_c, statuses_c, rows, predictions = _fit_score(core, proteins, targets[cause], w, clusters, folds, cause)
        point[cause], statuses[cause] = points_c, statuses_c
        contrasts[cause] = flatten_scores(points_c, statuses_c)
        parameters.extend(rows)
        for arm, p in predictions.items():
            bins.extend(dict(cause=cause, model=arm, **row) for row in calibration_bins(targets[cause], p, w))
    _check_shared_preprocessing(parameters)
    rng = np.random.default_rng(seed)
    wide, paired = [], []
    cause_draws = {cause: [] for cause in CAUSES}
    selection = {cause: {arm: {str(p): 0 for p in (1., 10., 100.)} for arm in ARMS} for cause in CAUSES}
    for rep in range(replicates):
        full_weights = design_bootstrap_weights(eligible, rng)
        wr = full_weights[common]
        row = dict(replicate=rep, **_resampling_summary(full_weights, wr, targets, clusters))
        replicate_parameters = []
        for cause in CAUSES:
            points_c, status_c, rows, _ = _fit_score(core, proteins, targets[cause], wr, clusters, folds, cause)
            replicate_parameters.extend(rows)
            flattened = flatten_scores(points_c, status_c)
            cause_draws[cause].append(dict(replicate=rep, **flattened))
            row.update({cause + '|' + name: value for name, value in flattened.items()})
            for fit in rows:
                selection[cause][fit['model']][str(float(fit['lambda']))] += 1
            for contrast in CONTRASTS:
                brier, auc = flattened[contrast + '|delta_brier'], flattened[contrast + '|delta_auc']
                paired.append(dict(replicate=rep, cause=cause, contrast=contrast,
                    status=flattened[contrast + '|status'], delta_brier=brier, delta_auc=auc,
                    brier_status='computed' if brier is not None else 'undefined',
                    auc_status='computed' if auc is not None else 'undefined'))
        _check_shared_preprocessing(replicate_parameters)
        wide.append(row)
        if (rep + 1) % 50 == 0:
            print('CA33B_SHARED_FULL_REFITS', rep + 1, flush=True)
    bounds = {cause: intervals(cause_draws[cause], replicates) for cause in CAUSES}
    for cause in CAUSES:
        for bound in bounds[cause].values():
            bound.update(marginal_exploratory=True, simultaneous_familywise_coverage=False,
                         endpoint_selection_uncertainty_included=False)
    result = dict(causes=list(CAUSES), n=len(w),
        events8={cause: int(targets[cause].sum()) for cause in CAUSES},
        controls8={cause: int((targets[cause] == 0).sum()) for cause in CAUSES},
        point=point, point_statuses=statuses, contrasts=contrasts, intervals=bounds,
        bootstrap_replicates=replicates, bootstrap_seed=seed,
        bootstrap_failures={cause: _failure_summary(cause_draws[cause]) for cause in CAUSES},
        bootstrap_penalty_selection_counts=selection,
        shared_preprocessing_verified_for_available_matched_fits=True,
        shared_preprocessing_verification_scope='available_same_arm_outer_fold_pairs_only; unavailable_pairs_after_fit_failure_not_verified',
        one_eligible_PSU_draw_per_replicate_shared_across_causes_and_arms=True,
        retained_zero_weight_rows_and_fixed_folds=True, bootstrap_redraw=False,
        marginal_exploratory_intervals_not_familywise_confirmation=True,
        between_cause_effect_comparisons_performed=False,
        independent_binary_probabilities_summed=False,
        joint_integration_and_Fine_Gray='deferred_not_implemented_or_run')
    return result, parameters, bins, wide, paired


def _contract():
    contract = json.loads((ROOT / 'CONTRACT.json').read_text())
    required = dict(causes=list(CAUSES), arms=list(ARMS),
        contrasts={key: list(value) for key, value in CONTRASTS.items()}, horizon_months=96,
        bootstrap_replicates=REPLICATES, bootstrap_seed=SEED,
        bootstrap_shared_across_causes_and_arms=True, bootstrap_refit_all_preprocessing_tuning_and_fits=True,
        bootstrap_redraw=False, interval_quantiles=[.025, .975], independent_probability_summing=False,
        no_between_cause_effect_contrast=True)
    if any(contract.get(key) != value for key, value in required.items()):
        raise ValueError('Cause/arm/bootstrap contract differs from the fixed analysis')
    return contract


def _check_gate(audit):
    equality, gate = audit.get('scientific_gate_equality', {}), audit.get('reproduced33A_gate', {})
    if (equality.get('exact_structure_and_counts') is not True
            or equality.get('expected_proceed_codes') != list(CAUSES)
            or gate.get('proceed_codes') != list(CAUSES)
            or gate.get('gate', {}).get('source_endpoint_qualified') is not True
            or gate.get('gate', {}).get('later_comparisons_allowed') is not True
            or any(gate.get('causes', {}).get(cause, {}).get('proceed_to_later_protocol') is not True for cause in CAUSES)):
        raise ValueError('Complete reproduced33A gate and frozen-count identity required before fitting')


def run(source, out):
    from .cohort import load, cause_targets

    out = Path(out)
    if out.exists():
        raise ValueError('New analysis output directory required')
    if PLAN_COMMIT == 'PENDING_BEFORE_EMPIRICAL_RUN':
        raise ValueError('Plan commit must be fixed before comparative execution')
    contract = _contract()
    eligible, common, audit, manifest = load(Path(source))
    _check_gate(audit)
    targets = cause_targets(eligible, common)
    d = eligible.loc[common].copy()
    if (len(d) != contract['common_n'] or len(eligible) != contract['eligible_n']
            or any(int(np.asarray(targets[cause]).sum()) != contract['cause_events8'][cause] for cause in CAUSES)
            or any(len(d) - int(np.asarray(targets[cause]).sum()) != contract['controls8'][cause] for cause in CAUSES)):
        raise ValueError('Locked common-domain/cause counts changed before fitting')
    core = d[list(CORE)].to_numpy(float)
    raw_proteins = d[list(PROTEINS)].to_numpy(float)
    if not np.isfinite(raw_proteins).all() or (raw_proteins <= 0).any():
        raise ValueError('Qualified positive reported protein values required')
    proteins = np.log2(raw_proteins)
    folds = assign_outer(eligible.cluster.to_numpy())[np.asarray(common)]
    result, parameters, bins, wide, paired = analyze_arrays(core, proteins, targets,
        d.weight.to_numpy(float), d.cluster.to_numpy(), folds, eligible, common,
        replicates=REPLICATES, seed=SEED)
    result.update(analysis='exploratory_separate_recorded_cause_binary_CIF_at8years_reused_NHANES1999_2004',
        plan_commit=PLAN_COMMIT, horizon_months=96, eligible_n=len(eligible),
        comparative_model_fitting_attempted=True, comparative_models_fitted=bool(parameters),
        routine_features=list(CORE), protein_features=list(PROTEINS), source_audit=audit,
        auc_controls=contract['auc_controls'], weight='pooled_corrected_six_year_SSCARD_weight',
        plan_sha256=_sha((ROOT / 'ANALYSIS_PLAN.md').read_bytes()),
        contract_sha256=_sha((ROOT / 'CONTRACT.json').read_bytes()),
        lock_sha256=_sha((ROOT / 'LOCK.json').read_bytes()),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)): _sha(p.read_bytes())
                                  for p in sorted(ROOT.parent.rglob('*.py'))},
        versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        independent_external_validation=False, clinical_use_ready=False, participant_records_exported=False,
        original_models_and_lab_report_v3_presentation_changed=False)
    payloads = {'FOLD_PARAMETERS.json': _bytes(parameters), 'CALIBRATION_BINS.json': _bytes(bins),
        'SOURCE_MANIFEST.json': _bytes(manifest),
        'bootstrap_aggregates.csv': pd.DataFrame(wide).to_csv(index=False).encode('utf-8'),
        'paired_bootstrap.csv': pd.DataFrame(paired).to_csv(index=False).encode('utf-8')}
    result['outputs_sha256'] = {name: _sha(value) for name, value in payloads.items()}
    result['source_manifest_export_sha256'] = result['outputs_sha256']['SOURCE_MANIFEST.json']
    payloads['RESULTS_SUMMARY.json'] = _bytes(result)
    # Serialize everything before creating the new output directory.
    out.mkdir(parents=True)
    for name, content in payloads.items():
        (out / name).write_bytes(content)
        print('CA33B_FILE_BEGIN ' + name, flush=True)
        print(content.decode('utf-8'), end='', flush=True)
        print('CA33B_FILE_END ' + name, flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.source, args.out)
