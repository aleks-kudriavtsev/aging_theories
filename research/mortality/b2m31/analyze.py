"""Exploratory Cys/B2M comparisons in the reused continuous NHANES cohort.

Only aggregate outputs are exported. Every bootstrap replicate repeats all
training preprocessing, tuning and fitting with the original PSU partitions.
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

from ..binary28.metrics import metrics, calibration_bins
from .model import nested_crossfit

ROOT = Path(__file__).parent
ARMS = ('routine', 'routine_b2m', 'routine_cys', 'routine_cys_b2m')
PROTEINS = ('SSCYST', 'SSB2M')
CONTRASTS = {'primary': ('routine', 'routine_b2m'),
             'secondary': ('routine_cys', 'routine_cys_b2m')}
METRICS = ('auc', 'brier', 'observed', 'predicted', 'bias_pp', 'observed_expected',
           'calibration_in_large', 'calibration_intercept', 'calibration_slope')
REPLICATES = 500
SEED = 20261102
PLAN_COMMIT = '6774fbe'


def _bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf-8')


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assign_outer(clusters):
    clusters = np.asarray(clusters)
    if clusters.ndim != 1 or any(not isinstance(c, str) or not c for c in clusters.tolist()):
        raise ValueError('Nonempty string PSU identities required')
    keys = sorted(set(clusters), key=lambda c: hashlib.sha256(('b2m31:outer:' + c).encode()).hexdigest())
    if len(keys) < 15:
        raise ValueError('Insufficient PSU support for five-by-three nested design')
    table = {key: i % 5 for i, key in enumerate(keys)}
    return np.array([table[c] for c in clusters], dtype=int)


def _survey_design(d):
    """Validate the precomplete resampling frame before any model evaluation."""
    required = {'weight', 'stratum', 'cluster'}
    if not required <= set(d.columns) or not len(d):
        raise ValueError('Survey resampling frame missing required fields')
    w = d.weight.to_numpy(float)
    if not np.isfinite(w).all() or (w < 0).any() or not np.isfinite(w.sum()) or w.sum() <= 0:
        raise ValueError('Invalid qualified survey weights')
    if d[['stratum', 'cluster']].isna().any().any():
        raise ValueError('Missing survey design identity')
    if any(not isinstance(c, str) or not c for c in d.cluster.tolist()):
        raise ValueError('Nonempty string PSU identities required')
    if (d.groupby('cluster', dropna=False).stratum.nunique() != 1).any():
        raise ValueError('PSU identity crosses survey strata')
    groups = []
    for stratum in sorted(d.stratum.unique(), key=str):
        indices = np.flatnonzero(d.stratum.to_numpy() == stratum)
        units = np.unique(d.cluster.to_numpy()[indices])
        if len(units) < 2:
            raise ValueError('Singleton stratum requires a predeclared variance strategy')
        groups.append((indices, units))
    return w, groups


def design_bootstrap_weights(d, rng):
    """Rescaled within-stratum whole-PSU bootstrap of the qualified frame."""
    w, groups = _survey_design(d)
    multiplier = np.zeros(len(d))
    cluster = d.cluster.to_numpy()
    for indices, units in groups:
        m = len(units)
        counts = Counter(rng.choice(units, size=m - 1, replace=True))
        for unit, count in counts.items():
            multiplier[indices[cluster[indices] == unit]] = count * m / (m - 1)
    return w * multiplier


def flatten_scores(points, statuses):
    """Keep arm/contrast status separate from metric-local calibration failures."""
    row = {}
    for arm in ARMS:
        row[arm + '|status'] = statuses[arm]['status']
        row[arm + '|error'] = statuses[arm].get('error')
        point = points.get(arm)
        for name in METRICS:
            row[arm + '|' + name] = None if point is None else point[name]
        row[arm + '|calibration_status'] = None if point is None else point['calibration_status']
    for label, (baseline, expanded) in CONTRASTS.items():
        defined = baseline in points and expanded in points
        row[label + '|status'] = 'computed' if defined else 'undefined'
        for name in ('brier', 'auc'):
            a = row[expanded + '|' + name]
            b = row[baseline + '|' + name]
            row[label + '|delta_' + name] = None if a is None or b is None else a - b
    return row


def intervals(draws, total):
    keys = [arm + '|' + metric for arm in ARMS for metric in METRICS]
    keys += [label + '|delta_' + metric for label in CONTRASTS for metric in ('brier', 'auc')]
    result = {}
    if len(draws) != total:
        raise ValueError('Every attempted bootstrap replicate must be retained')
    for key in keys:
        values = [d[key] for d in draws if d.get(key) is not None and np.isfinite(d[key])]
        lower, upper = map(float, np.quantile(values, [.025, .975])) if values else (None, None)
        result[key] = dict(lower=lower, upper=upper, valid=len(values), undefined=total - len(values),
            scope='nested_preprocessing_tuning_refitting_with_fixed_partition',
            defined_draws_only=len(values) != total)
    return result


def fit_models(core, proteins, y, w, clusters, folds, domain, tolerate_fit_failure=False):
    predictions, points, statuses, parameters = {}, {}, {}, []
    for arm in ARMS:
        try:
            p, rows = nested_crossfit(core, proteins, y, w, clusters, folds, arm)
            point = metrics(y, p, w)
        except (ValueError, ArithmeticError, RuntimeError) as error:
            if not tolerate_fit_failure:
                raise
            statuses[arm] = dict(status='undefined', error=type(error).__name__ + ': ' + str(error))
            continue
        predictions[arm] = p
        points[arm] = point
        statuses[arm] = dict(status='computed')
        parameters.extend(dict(analysis_domain=domain, **row) for row in rows)
    return predictions, points, statuses, parameters


def _domain_arrays(d, core_names):
    core = d[list(core_names)].to_numpy(float)
    raw_proteins = d[list(PROTEINS)].to_numpy(float)
    if core.shape != (len(d), 11) or not np.isfinite(core).all():
        raise ValueError('Eleven finite qualified routine features required')
    if not np.isfinite(raw_proteins).all() or (raw_proteins <= 0).any():
        raise ValueError('Positive finite qualified reported proteins required')
    y = (d.dead.eq(1) & d.t_exam.le(96)).to_numpy(int)
    w = d.weight.to_numpy(float)
    if not len(y) or not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
        raise ValueError('Nonempty positive-weight qualified common domain required')
    return core, np.log2(raw_proteins), y, w, d.cluster.to_numpy()


def run(source, out):
    # Lazy imports keep pure aggregate/statistical helpers testable without data.
    from .cohort import load
    from .features import CORE

    source, out = Path(source), Path(out)
    if out.exists():
        raise ValueError('New analysis output directory required')
    if PLAN_COMMIT == 'PENDING_BEFORE_EMPIRICAL_RUN':
        raise ValueError('Analysis plan commit must be pinned before empirical execution')
    eligible, common, reported, audit, manifest = load(source)
    common, reported = np.asarray(common), np.asarray(reported)
    if (common.shape != (len(eligible),) or reported.shape != common.shape
            or common.dtype != bool or reported.dtype != bool or (common & ~reported).any()):
        raise ValueError('Aligned primary/reported common-domain boolean masks required')
    if not isinstance(manifest, dict) or not isinstance(manifest.get('files'), list):
        raise ValueError('Complete source retrieval manifest required')
    _survey_design(eligible)
    if (not eligible.dead.isin([0, 1]).all() or not np.isfinite(eligible.t_exam).all()
            or (eligible.t_exam < 0).any()):
        raise ValueError('Known binary mortality and finite nonnegative exam followup required')
    if (eligible.dead.eq(0) & eligible.t_exam.lt(96)).any():
        raise ValueError('Eight-year binary outcome has early censoring')
    folds_all = assign_outer(eligible.cluster.to_numpy())
    d = eligible.loc[common].copy()
    if not len(d):
        raise ValueError('Empty primary common domain')
    core, proteins, y, w, clusters = _domain_arrays(d, CORE)
    # Validate both prespecified data domains before comparisons; schema/source
    # errors must not become a tolerated sensitivity fitting failure.
    s = eligible.loc[reported].copy()
    cs, ps, ys, ws, clusters_s = _domain_arrays(s, CORE)
    folds = folds_all[common]
    predictions, point, statuses, parameters = fit_models(core, proteins, y, w, clusters, folds, 'primary_common')
    scores = flatten_scores(point, statuses)
    bins = []
    for arm, p in predictions.items():
        bins.extend(dict(model=arm, analysis_domain='primary_common', **row) for row in calibration_bins(y, p, w))
    print('B2M31_POINT_COMPUTED', json.dumps(dict(n=len(d), events8=int(y.sum()),
        primary_delta_brier=scores['primary|delta_brier'], secondary_delta_brier=scores['secondary|delta_brier'])), flush=True)
    rng = np.random.default_rng(SEED)
    draws = []
    selection = {arm: {str(p): 0 for p in (1., 10., 100.)} for arm in ARMS}
    for rep in range(REPLICATES):
        wr = design_bootstrap_weights(eligible, rng)[common]
        _, points_r, status_r, rows = fit_models(core, proteins, y, wr, clusters, folds,
                                                'bootstrap', tolerate_fit_failure=True)
        row = dict(replicate=rep, **flatten_scores(points_r, status_r))
        draws.append(row)
        for fit in rows:
            selection[fit['model']][str(float(fit['lambda']))] += 1
        if (rep + 1) % 50 == 0:
            print('B2M31_NESTED_REFITS', rep + 1, flush=True)

    sensitivity = dict(status='undefined', n=int(reported.sum()),
        interpretation='reported_values_flags01_different_conditional_domain_within_domain_pairs_only',
        bootstrap_intervals_computed=False)
    try:
        preds_s, points_s, status_s, params_s = fit_models(cs, ps, ys, ws, clusters_s,
                                                        folds_all[reported], 'reported_values_sensitivity')
        parameters.extend(params_s)
        sensitivity.update(status='computed', events8=int(ys.sum()), point=points_s,
                           contrasts=flatten_scores(points_s, status_s))
        for arm, p in preds_s.items():
            bins.extend(dict(model=arm, analysis_domain='reported_values_sensitivity', **row)
                        for row in calibration_bins(ys, p, ws))
    except (ValueError, ArithmeticError, RuntimeError) as error:
        sensitivity['error'] = type(error).__name__ + ': ' + str(error)

    wide_csv = pd.DataFrame(draws).to_csv(index=False)
    narrow = {}
    for label in CONTRASTS:
        rows = [dict(replicate=row['replicate'], status=row[label + '|status'],
                     delta_brier=row[label + '|delta_brier'], delta_auc=row[label + '|delta_auc']) for row in draws]
        narrow[label] = pd.DataFrame(rows).to_csv(index=False)
    source_bytes = _bytes(manifest)
    result = dict(analysis='exploratory_reused_NHANES1999_2004_new_B2M_binary_models', plan_commit=PLAN_COMMIT,
        horizon_years=8, horizon_months=96, n=len(d), events8=int(y.sum()), eligible_n=len(eligible),
        common_domain='all11_routine_plus_positive_Cys_B2M_both_flag0',
        weight='source_qualified_six_year_weight', effective_n=float(w.sum()**2 / (w @ w)),
        strata_n=int(d.stratum.nunique()), psu_n=len(set(clusters)),
        common_weighted_coverage=float(np.average(common, weights=eligible.weight)),
        point=point, contrasts=scores,
        primary_contrast='routine_b2m_minus_routine', secondary_contrast='routine_cys_b2m_minus_routine_cys',
        primary_score='survey_weighted_eight_year_Brier', intervals=intervals(draws, REPLICATES),
        bootstrap_replicates=REPLICATES, bootstrap_seed=SEED,
        bootstrap_model_failure_counts={arm: sum(r[arm + '|status'] != 'computed' for r in draws) for arm in ARMS},
        bootstrap_contrast_failure_counts={label: sum(r[label + '|status'] != 'computed' for r in draws) for label in CONTRASTS},
        bootstrap_penalty_selection_counts=selection,
        bootstrap_partition_fixed=True, full_nested_preprocessing_tuning_and_refitting=True,
        reported_values_sensitivity=sensitivity, source_audit=audit,
        routine_features=list(CORE), protein_features=list(PROTEINS),
        same_common_domain_all_four_models=True, no_single_protein_winner_search=True,
        independent_external_validation=False, clinical_use_ready=False, participant_records_exported=False,
        previous_frozen_coefficients_used=False,
        model_family='eight_year_survey_weighted_ridge_logistic_not_cause_specific_engine',
        assay_calibration_and_prior_selection_uncertainty_included=False,
        source_manifest_export_sha256=hashlib.sha256(source_bytes).hexdigest(),
        bootstrap_aggregates_csv_sha256=hashlib.sha256(wide_csv.encode()).hexdigest(),
        primary_bootstrap_csv_sha256=hashlib.sha256(narrow['primary'].encode()).hexdigest(),
        secondary_bootstrap_csv_sha256=hashlib.sha256(narrow['secondary'].encode()).hexdigest(),
        plan_sha256=_sha(ROOT / 'ANALYSIS_PLAN.md'), lock_sha256=_sha(ROOT / 'LOCK.json'),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)): _sha(p) for p in sorted(ROOT.parent.rglob('*.py'))},
        versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__))
    out.mkdir(parents=True)
    for filename, value in [('RESULTS.json', result), ('fold_parameters.json', parameters), ('calibration_bins.json', bins)]:
        (out / filename).write_bytes(_bytes(value))
    (out / 'SOURCE_MANIFEST.json').write_bytes(source_bytes)
    (out / 'bootstrap_aggregates.csv').write_text(wide_csv)
    for label, text in narrow.items():
        (out / (label + '_bootstrap.csv')).write_text(text)
    for label, value in [('RESULTS', result), ('MANIFEST', manifest), ('PARAMETERS', parameters), ('CALIBRATION_BINS', bins)]:
        print('B2M31_' + label + '_BEGIN')
        print(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False))
        print('B2M31_' + label + '_END')
    for label, text in [('BOOTSTRAP', wide_csv), ('PRIMARY_BOOTSTRAP', narrow['primary']),
                        ('SECONDARY_BOOTSTRAP', narrow['secondary'])]:
        print('B2M31_' + label + '_CSV_BEGIN')
        print(text, end='')
        print('B2M31_' + label + '_CSV_END')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.source, args.out)
