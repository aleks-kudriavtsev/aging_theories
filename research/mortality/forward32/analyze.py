"""Forward-calendar assessment with two separately labelled uncertainty modes.

Early fitting APIs receive no late outcomes or weights. Only aggregates and
new early-training parameter artifacts are exported; predictions stay in memory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import scipy

from ..b2m31.analyze import (
    ARMS, CONTRASTS, METRICS, _survey_design, design_bootstrap_weights,
    flatten_scores, intervals,
)
from ..b2m31.features import CORE
from ..binary28.metrics import metrics, calibration_bins
from . import model

ROOT = Path(__file__).parent
PROTEINS = ('SSCYST', 'SSB2M')
MODES = ('fixed_training', 'joint_refit')
PLAN_COMMIT = '377778c'
REPLICATES = 500
LATE_SEED = 20261103
TRAIN_SEED = 20261104


def _bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf-8')


def _digest(value):
    return hashlib.sha256(value).hexdigest()


def _configuration():
    lock = json.loads((ROOT / 'LOCK.json').read_text())
    required = dict(replicates_per_mode=REPLICATES, late_seed=LATE_SEED, train_seed=TRAIN_SEED,
                    primary_uncertainty_mode='joint_refit', inner_salt=model.INNER_SALT,
                    inner_folds=3, penalties=[1., 10., 100.])
    if any(lock.get(key) != value for key, value in required.items()):
        raise ValueError('Forward uncertainty/training configuration differs from locked plan')
    return required


def _provenance():
    return dict(plan_commit=PLAN_COMMIT,
        plan_sha256=_digest((ROOT / 'ANALYSIS_PLAN.md').read_bytes()),
        lock_sha256=_digest((ROOT / 'LOCK.json').read_bytes()),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)): _digest(p.read_bytes())
                                  for p in sorted(ROOT.parent.rglob('*.py'))},
        versions=dict(python=platform.python_version(), numpy=np.__version__,
                      pandas=pd.__version__, scipy=scipy.__version__))


def _validate_schema(eligible, common, audit, manifest):
    required = set(CORE) | set(PROTEINS) | {'partition', 'cycle', 'weight', 'stratum', 'cluster',
                                          'dead', 't_exam', 'WTSSCB4Y', 'WTSSCB2Y'}
    if not isinstance(eligible, pd.DataFrame) or not required <= set(eligible.columns):
        raise ValueError('Qualified temporal source schema required')
    common = np.asarray(common)
    if common.shape != (len(eligible),) or common.dtype != bool:
        raise ValueError('Aligned primary common-domain boolean mask required')
    if not set(eligible.partition.unique()) <= {'train', 'late'}:
        raise ValueError('Unknown temporal partition label')
    early = eligible.partition.eq('train')
    if (not eligible.loc[early, 'cycle'].isin([1999, 2001]).all()
            or not eligible.loc[~early, 'cycle'].eq(2003).all()):
        raise ValueError('Partition labels do not match the fixed calendar split')
    if (not isinstance(audit, dict) or not isinstance(audit.get('gate'), dict)
            or type(audit['gate'].get('temporal_evaluation_allowed')) is not bool
            or not isinstance(audit['gate'].get('blockers'), list)
            or any(not isinstance(item, str) for item in audit['gate']['blockers'])):
        raise ValueError('Explicit temporal feasibility gate required')
    if (not isinstance(manifest, dict) or manifest.get('schema') != 'forward32-source-manifest/1'
            or not isinstance(manifest.get('intake'), dict) or not isinstance(manifest.get('source'), dict)
            or not isinstance(manifest['source'].get('files'), list)):
        raise ValueError('Complete forward retrieval and source manifest required')
    return common


def _prepare_period(eligible, common, part):
    selected = eligible.partition.eq(part).to_numpy()
    base = eligible.loc[selected].copy()
    mask = common[selected]
    native = base.WTSSCB4Y if part == 'train' else base.WTSSCB2Y
    if not np.array_equal(base.weight.to_numpy(float), native.to_numpy(float)):
        raise ValueError('Native-period weights were not assigned before analysis')
    _survey_design(base)
    if (not np.isfinite(base.weight).all() or (base.weight <= 0).any()
            or not base.dead.isin([0, 1]).all() or not np.isfinite(base.t_exam).all()
            or (base.t_exam < 0).any() or (base.dead.eq(0) & base.t_exam.lt(96)).any()):
        raise ValueError('Qualified native weights and known eight-year endpoint required')
    d = base.loc[mask].copy()
    if not len(d) or set(d.cluster) != set(base.cluster):
        raise ValueError('Common domain must preserve the eligible PSU inventory')
    raw_proteins = d[list(PROTEINS)].to_numpy(float)
    if not np.isfinite(raw_proteins).all() or (raw_proteins <= 0).any():
        raise ValueError('Positive primary reported proteins required')
    core, proteins = model.numerical._features(d[list(CORE)].to_numpy(float), np.log2(raw_proteins))
    y = (d.dead.eq(1) & d.t_exam.le(96)).to_numpy(int)
    w = d.weight.to_numpy(float)
    if np.unique(y[w > 0]).size != 2:
        raise ValueError('Allowed temporal domain must contain both outcome classes')
    return dict(eligible=base, common=mask, core=core, protein=proteins, y=y,
                weight=w, clusters=d.cluster.to_numpy(), frame=d)


def _check_inner_mapping(train, late, audit):
    # The cohort gate requires identical eligible/common PSU sets. Thus hashing
    # common early rows is exactly hashing the eligible precomplete inventory.
    values = model._validated_split(train['core'], train['protein'], train['y'], train['weight'],
        train['clusters'], late['core'], late['protein'], late['clusters'])
    folds = values[-1]
    mapping = {key: int(folds[np.flatnonzero(train['clusters'] == key)[0]])
               for key in sorted(set(train['clusters']))}
    if mapping != audit.get('eligible_derived_early_inner_fold_mapping'):
        raise ValueError('Training folds differ from the eligible-derived fixed mapping')
    return folds


def _fit_arms(train, late, weights):
    predictions, rows, statuses = {}, {}, {}
    for arm in ARMS:
        try:
            p, row = model.fit_training_arm(train['core'], train['protein'], train['y'], weights,
                train['clusters'], late['core'], late['protein'], late['clusters'], arm)
        except (ValueError, ArithmeticError, RuntimeError) as error:
            statuses[arm] = dict(status='undefined', error='training|' + type(error).__name__ + ': ' + str(error))
            continue
        predictions[arm], rows[arm], statuses[arm] = p, row, dict(status='computed')
    return predictions, rows, statuses


def _score_predictions(y, predictions, weights, model_statuses):
    points, statuses = {}, {}
    for arm in ARMS:
        if arm not in predictions:
            statuses[arm] = dict(model_statuses[arm])
            continue
        try:
            points[arm] = metrics(y, predictions[arm], weights)
        except (ValueError, ArithmeticError, RuntimeError) as error:
            statuses[arm] = dict(status='undefined', error='evaluation|' + type(error).__name__ + ': ' + str(error))
            continue
        statuses[arm] = dict(status='computed')
    return points, statuses


def _point_artifact(train, late, rows, statuses, folds, provenance):
    y, w = train['y'], train['weight']
    return dict(schema=model.SCHEMA, models=rows, model_statuses=statuses,
        inner_folds_aggregate=model._inner_support(y, w, train['clusters'], folds),
        training=dict(n=len(y), events8=int(y.sum()), positive_weight_n=int((w > 0).sum()),
            positive_weight_events8=int(((w > 0) & (y == 1)).sum()), weight_sum=float(w.sum()),
            psu_n=len(set(train['clusters'])), weight_name='WTSSCB4Y', cycles=[1999, 2001],
            inner_folds=3, inner_salt=model.INNER_SALT, weight_normalization='sum_to_raw_training_n'),
        late_domain=dict(n=len(late['core']), psu_n=len(set(late['clusters'])), cycle=2003,
            evaluation_weight_name='WTSSCB2Y', training_psu_overlap_n=0,
            outcomes_used_for_fitting=False, weights_used_for_fitting=False),
        fit_scope='early_training_only', selection_scope='early_three_fold_PSU_inner_Brier',
        no_outer_cross_validation=True, late_feature_preprocessing_fitted=False,
        inner_scores_are_training_diagnostics_not_validation=True,
        independent_external_validation=False, participant_records_exported=False, **provenance)


def _weight_summary(weights, y, clusters, prefix):
    positive = weights > 0
    return {prefix + 'weight_sum': float(weights.sum()),
        prefix + 'weight_squared_sum': float(weights @ weights),
        prefix + 'positive_weight_n': int(positive.sum()),
        prefix + 'positive_weight_events8': int((positive & (y == 1)).sum()),
        prefix + 'positive_weight_nonevents8': int((positive & (y == 0)).sum()),
        prefix + 'event_weight_sum': float(weights @ y),
        prefix + 'positive_weight_psu_n': len(set(clusters[positive]))}


def bootstrap_modes(train, late, point_predictions, point_statuses, replicates=REPLICATES,
                    late_seed=LATE_SEED, train_seed=TRAIN_SEED):
    """Retain every paired draw; no late outcomes/weights enter arm fitting."""
    if type(replicates) is not int or replicates <= 0 or late_seed == train_seed:
        raise ValueError('Positive replicate count and independent RNG seeds required')
    late_rng, train_rng = np.random.default_rng(late_seed), np.random.default_rng(train_seed)
    draws = {mode: [] for mode in MODES}
    selection = {arm: {str(p): 0 for p in (1., 10., 100.)} for arm in ARMS}
    for rep in range(replicates):
        # The full qualified period frame is sampled BEFORE its fixed complete
        # mask. Zero-weight common rows retain their original cluster/fold IDs.
        late_full = design_bootstrap_weights(late['eligible'], late_rng)
        late_weights = late_full[late['common']]
        early_full = design_bootstrap_weights(train['eligible'], train_rng)
        early_weights = early_full[train['common']]
        late_summary = _weight_summary(late_weights, late['y'], late['clusters'], 'late_')
        late_summary.update(late_eligible_weight_sum=float(late_full.sum()),
                            late_eligible_weight_squared_sum=float(late_full @ late_full))
        fixed_points, fixed_status = _score_predictions(late['y'], point_predictions, late_weights, point_statuses)
        fixed = dict(mode='fixed_training', replicate=rep, early_resample_used_for_fitting=False,
            **late_summary, **_weight_summary(train['weight'], train['y'], train['clusters'], 'early_'),
            **flatten_scores(fixed_points, fixed_status))
        predictions, rows, fitted_status = _fit_arms(train, late, early_weights)
        joint_points, joint_status = _score_predictions(late['y'], predictions, late_weights, fitted_status)
        joint = dict(mode='joint_refit', replicate=rep, early_resample_used_for_fitting=True,
            **late_summary, **_weight_summary(early_weights, train['y'], train['clusters'], 'early_'),
            **flatten_scores(joint_points, joint_status))
        for arm, row in rows.items():
            selection[arm][str(float(row['selected_penalty']))] += 1
        draws['fixed_training'].append(fixed)
        draws['joint_refit'].append(joint)
        if (rep + 1) % 50 == 0:
            print('FORWARD32_PAIRED_DRAWS', rep + 1, flush=True)
    return draws, selection


def _mode_summary(draws, total, mode):
    values = intervals(draws, total)
    scope = ('late_PSU_resampling_with_fixed_early_training' if mode == 'fixed_training' else
             'independent_early_and_late_PSU_resampling_with_full_early_preprocessing_tuning_refitting')
    for row in values.values():
        row.update(scope=scope, uncertainty_mode=mode, early_training_refitted=mode == 'joint_refit')
    return dict(intervals=values, attempted_draws=total,
        model_failure_counts={arm: sum(d[arm + '|status'] != 'computed' for d in draws) for arm in ARMS},
        contrast_failure_counts={label: sum(d[label + '|status'] != 'computed' for d in draws) for label in CONTRASTS},
        scope=scope, is_primary_uncertainty=mode == 'joint_refit')


def _base_result(audit, provenance):
    return dict(analysis='exploratory_reused_NHANES1999_2004_fixed_forward_calendar_B2M',
        horizon_years=8, horizon_months=96, gate=audit['gate'], source_audit=audit,
        train_cycles=[1999, 2001], late_cycles=[2003], train_weight='WTSSCB4Y', late_weight='WTSSCB2Y',
        routine_features=list(CORE), protein_features=list(PROTEINS),
        primary_contrast='routine_b2m_minus_routine', secondary_contrast='routine_cys_b2m_minus_routine_cys',
        primary_score='late_survey_weighted_eight_year_Brier', primary_uncertainty_mode='joint_refit',
        same_late_draw_all_arms_and_both_modes=True, uncertainty_modes_pooled=False,
        same_common_domain_all_four_models_within_each_period=True,
        late_recalibration_performed=False, independent_external_validation=False,
        clinical_use_ready=False, participant_records_exported=False,
        previous_frozen_coefficients_used=False, **provenance)


def _export(out, result, manifest, artifacts=None, bins=None, csvs=None):
    payloads = {'SOURCE_MANIFEST.json': _bytes(manifest)}
    if artifacts is not None:
        payloads['MODEL_ARTIFACTS.json'] = _bytes(artifacts)
        payloads['CALIBRATION_BINS.json'] = _bytes(bins)
    for filename, text in (csvs or {}).items():
        payloads[filename] = text.encode('utf-8')
    result['outputs_sha256'] = {name: _digest(value) for name, value in payloads.items()}
    result['source_manifest_export_sha256'] = result['outputs_sha256']['SOURCE_MANIFEST.json']
    payloads['RESULTS.json'] = _bytes(result)
    out.mkdir(parents=True)
    for filename, value in payloads.items():
        (out / filename).write_bytes(value)
    blocks = [('RESULTS', result), ('MANIFEST', manifest)]
    if artifacts is not None:
        blocks.extend([('MODEL_ARTIFACTS', artifacts), ('CALIBRATION_BINS', bins)])
    for label, value in blocks:
        print('FORWARD32_' + label + '_BEGIN')
        print(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False))
        print('FORWARD32_' + label + '_END')
    for filename, text in (csvs or {}).items():
        label = filename.removesuffix('.csv').upper()
        print('FORWARD32_' + label + '_CSV_BEGIN')
        print(text, end='')
        print('FORWARD32_' + label + '_CSV_END')


def run(source, out):
    from .cohort import load

    out = Path(out)
    if out.exists():
        raise ValueError('New forward output directory required')
    if PLAN_COMMIT == 'PENDING_BEFORE_EMPIRICAL_RUN':
        raise ValueError('Protocol commit must be pinned before empirical execution')
    config = _configuration()
    eligible, common, audit, manifest = load(Path(source))
    common = _validate_schema(eligible, common, audit, manifest)
    provenance = _provenance()
    result = _base_result(audit, provenance)
    if not audit['gate']['temporal_evaluation_allowed']:
        result.update(status='blocked_feasibility_no_model_fits', model_fits_attempted=0,
                      uncertainty_draws_attempted=0, blocked_reasons=audit['gate']['blockers'])
        _export(out, result, manifest)
        return result
    if audit['gate']['blockers']:
        raise ValueError('Allowed temporal gate contains unresolved blockers')
    if set(eligible.partition.unique()) != {'train', 'late'}:
        raise ValueError('Allowed temporal gate requires both fixed partitions')
    train, late = _prepare_period(eligible, common, 'train'), _prepare_period(eligible, common, 'late')
    folds = _check_inner_mapping(train, late, audit)
    # All structural/domain checks precede the arm-local numerical failure loop.
    predictions, fitted_rows, model_statuses = _fit_arms(train, late, train['weight'])
    artifact = _point_artifact(train, late, fitted_rows, model_statuses, folds, provenance)
    point, statuses = _score_predictions(late['y'], predictions, late['weight'], model_statuses)
    point_scores = flatten_scores(point, statuses)
    bins = []
    for arm, p in predictions.items():
        bins.extend(dict(model=arm, analysis_domain='late_point', **row)
                    for row in calibration_bins(late['y'], p, late['weight']))
    print('FORWARD32_POINT_COMPUTED', json.dumps(dict(train_n=len(train['y']), late_n=len(late['y']),
        primary_delta_brier=point_scores['primary|delta_brier'],
        secondary_delta_brier=point_scores['secondary|delta_brier'])), flush=True)
    draws, selection = bootstrap_modes(train, late, predictions, model_statuses,
        replicates=config['replicates_per_mode'], late_seed=config['late_seed'], train_seed=config['train_seed'])
    paired = []
    for mode in MODES:
        for row in draws[mode]:
            paired.append({key: row[key] for key in ['mode', 'replicate']
                + [label + '|' + key for label in CONTRASTS for key in ('status', 'delta_brier', 'delta_auc')]})
    csvs = {mode + '_bootstrap.csv': pd.DataFrame(draws[mode]).to_csv(index=False) for mode in MODES}
    csvs['paired_bootstrap.csv'] = pd.DataFrame(paired).to_csv(index=False)
    result.update(status='computed' if all(s['status'] == 'computed' for s in statuses.values())
                        else 'computed_with_point_arm_failures',
        model_fits_attempted=4, train_n=len(train['y']), train_events8=int(train['y'].sum()),
        late_n=len(late['y']), late_events8=int(late['y'].sum()),
        point=point, point_model_statuses=model_statuses, point_statuses=statuses, contrasts=point_scores,
        bootstrap_replicates_per_mode=config['replicates_per_mode'], late_seed=config['late_seed'], train_seed=config['train_seed'],
        uncertainty={mode: _mode_summary(draws[mode], config['replicates_per_mode'], mode) for mode in MODES},
        joint_refit_penalty_selection_counts=selection,
        training_inner_fold_aggregates_match_eligible_mapping=True,
        bootstrap_training_inner_partition_fixed=True,
        fixed_training_excludes_training_uncertainty=True,
        assay_storage_calibration_selection_history_uncertainty_included=False)
    _export(out, result, manifest, artifact, bins, csvs)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.source, args.out)
