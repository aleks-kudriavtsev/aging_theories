"""Four-year urate update on a wider NHANES MEC domain; aggregate outputs only."""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import sklearn
from ..nfl22 import analyze as base
from ..joint23.robustness import crossfit, contrast

ROOT = Path(__file__).parent
SEED = 20261029
PLAN_COMMIT = '90663b1'


def join_component(left, right, ids):
    if right.SEQN.isna().any() or right.SEQN.duplicated().any() or not set(right.SEQN) <= ids:
        raise ValueError('Invalid component IDs')
    return left.merge(right[['SEQN'] + [c for c in right if c not in left]],
                      on='SEQN', how='left', validate='1:1')


def prepare_domain(d):
    eligible = (d.elig.eq(1) & d.dead.isin([0, 1]) & d.t_exam.ge(0)
                & np.isfinite(d.t_exam) & d.RIAGENDR.isin([1, 2])
                & d.WTMEC2YR.gt(0) & np.isfinite(d.WTMEC2YR))
    d = d.loc[eligible].copy().reset_index(drop=True)
    if d.empty or d[['SDMVSTRA', 'SDMVPSU']].isna().any().any():
        raise ValueError('Empty domain or missing survey design')
    if ((d.dead == 0) & (d.t_exam < 48)).any():
        raise ValueError('Early censoring requires new protocol')
    d['weight'] = d.WTMEC2YR
    d['stratum'] = '2013:' + d.SDMVSTRA.astype(int).astype(str)
    d['cluster'] = d.SDMVSTRA.astype(int).astype(str) + ':' + d.SDMVPSU.astype(int).astype(str)
    d['time_years'] = d.t_exam.clip(lower=.5) / 12
    return d


def cohort(source):
    source = Path(source)
    audit = base.check_sources(source)
    demo = base.read_xpt_checked(source / '2013_DEMO_H.xpt')
    if demo.SEQN.isna().any() or demo.SEQN.duplicated().any():
        raise ValueError('Invalid demographic IDs')
    ids = set(demo.SEQN)
    cols = ['SEQN', 'RIDAGEYR', 'RIAGENDR', 'SDMVSTRA', 'SDMVPSU', 'WTMEC2YR']
    d = demo.loc[demo.RIDAGEYR.between(40, 75), cols].copy()
    age_n = len(d)
    for c in base.COMPONENTS[1:]:
        d = join_component(d, base.read_xpt_checked(source / f'2013_{c}_H.xpt'), ids)
    mort = base.independent_mortality_reader(source / '2013_mortality.dat')
    if set(mort.SEQN) != ids:
        raise ValueError('Mortality ID mismatch')
    d = prepare_domain(join_component(d, mort, ids))
    d = base.engineer_expanded(d)
    d['log_urate'] = np.log2(d.LBXSUA.where(d.LBXSUA > 0))
    audit.update(survey_n=len(demo), age40_75_n=age_n, eligible_mec_n=len(d),
                 reader_audit=base.audit_mortality_reader(source / '2013_mortality.dat'))
    return d, audit


def reproduce_baseline(source, model):
    e, _ = base.cohort(Path(source))
    common = np.isfinite(e[model['preprocessing']['features']]).all(axis=1).to_numpy()
    d = e.loc[common].copy()
    y = ((d.dead == 1) & (d.time_years <= 4)).to_numpy(int)
    if len(d) != 1270 or y.sum() != 37:
        raise ValueError('NfL22 sample mismatch')
    p, _ = crossfit(y, base.baseline_hazard(d, model, np.minimum(d.time_years, 4)),
                   base.baseline_hazard(d, model, 4), d.weight.to_numpy(),
                   np.log2(d.SSSNFL.to_numpy()), base.assign_folds(e.cluster)[common])
    result = contrast(y, p, d.weight.to_numpy())
    targets = dict(calibration_brier=.02307227, nfl_brier=.02339697,
                   calibration_auc=.78016395, nfl_auc=.77144764)
    if any(abs(result[k] - v) > 1e-8 for k, v in targets.items()):
        raise ValueError('NfL22 metric mismatch')
    return dict(n=len(d), events4=int(y.sum()), **result)


def urate_crossfit(y, hs, h4, w, x, folds):
    # Reuse the audited single-covariate likelihood, not its NfL field labels.
    p, _ = crossfit(y, hs, h4, w, x, folds)
    return {'calibration': p['calibration'], 'urate': p['nfl']}


def score(y, p, w):
    rows = {k: base.binary_metrics(y, v, w) for k, v in p.items()}
    a, b = rows['urate'], rows['calibration']
    return dict(metrics=rows, delta_brier=a['brier'] - b['brier'],
                delta_auc=None if a['auc'] is None or b['auc'] is None else a['auc'] - b['auc'])


def summarize(draws, total):
    result = {}
    for key in ('delta_brier', 'delta_auc'):
        values = [r[key] for r in draws if r.get(key) is not None and np.isfinite(r[key])]
        limits = [float(v) for v in np.quantile(values, [.025, .975])] if values else [None, None]
        result[key] = dict(lower=limits[0], upper=limits[1], valid=len(values),
                           undefined=total - len(values),
                           scope='fixed_baseline_partition_update_refitted',
                           defined_draws_only=len(values) != total)
    return result


def run(source, out):
    source, out = Path(source), Path(out)
    if out.exists():
        raise ValueError('New output directory required')
    model = base.load_model()
    reproduction = reproduce_baseline(source, model)
    e, audit = cohort(source)
    features = model['preprocessing']['features']
    complete_base = np.isfinite(e[features]).all(axis=1).to_numpy()
    common = complete_base & np.isfinite(e.log_urate).to_numpy()
    d = e.loc[common].copy()
    folds = base.assign_folds(e.cluster)[common]
    y = ((d.dead == 1) & (d.time_years <= 4)).to_numpy(int)
    hs = base.baseline_hazard(d, model, np.minimum(d.time_years, 4))
    h4 = base.baseline_hazard(d, model, 4)
    w, x = d.weight.to_numpy(), d.log_urate.to_numpy()
    p = urate_crossfit(y, hs, h4, w, x, folds)
    p['frozen'] = -np.expm1(-h4)
    point = score(y, p, w)
    rng = np.random.default_rng(SEED)
    draws = []
    for rep in range(1000):
        wr = base.design_bootstrap_weights(e, rng)[common]
        try:
            q = urate_crossfit(y, hs, h4, wr, x, folds)
            s = score(y, q, wr)
            row = dict(replicate=rep, status='computed', delta_brier=s['delta_brier'], delta_auc=s['delta_auc'])
        except (ValueError, ArithmeticError, RuntimeError) as exc:
            row = dict(replicate=rep, status='undefined', error=type(exc).__name__ + ': ' + str(exc))
        draws.append(row)
    result = dict(analysis='exploratory_urate_MEC_not_external_validation', n=len(d),
                  events4=int(y.sum()), eligible_n=len(e), complete_baseline_n=int(complete_base.sum()),
                  missing_baseline_n=int((~complete_base).sum()),
                  urate_excluded_after_baseline_n=int((complete_base & ~common).sum()),
                  weight='WTMEC2YR', horizon_years=4, outcome_states=2,
                  baseline_cause_components=len(model['causes']), point=point,
                  bootstrap=summarize(draws, 1000), bootstrap_refit_failures=sum(r['status']=='undefined' for r in draws),
                  seed=SEED, replicates=1000, prior_baseline_reproduction=reproduction,
                  source_audit=audit, baseline_sha256=base.MODEL_SHA,
                  plan_commit=PLAN_COMMIT, plan_sha256=base.sha(ROOT/'ANALYSIS_PLAN.md'),
                  code_sha256=base.sha(__file__),
                  versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                                scipy=scipy.__version__, sklearn=sklearn.__version__),
                  clinical_use_ready=False, participant_records_exported=False)
    out.mkdir(parents=True)
    base.dump(out/'RESULTS.json', result)
    pd.DataFrame(draws).to_csv(out/'bootstrap_aggregates.csv', index=False)
    (out/'SOURCE_MANIFEST.json').write_bytes((source/'manifest.json').read_bytes())
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.out), indent=2, allow_nan=False))
