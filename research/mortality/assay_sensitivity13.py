"""Post-validation assay sensitivity on already examined cohorts. No refitting.

Fixed plan: assay13/ANALYSIS_PLAN.md. Artificial analytical perturbations are NOT
observed platform biases, interventions, confidence intervals or allowable error.
The code only exports aggregates and does not infer correction factors from death.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .transport_no_crp12 import (load_evaluation, engineer_no_crp, predict, labels,
    multiclass_brier, GROUPS)
from .validate_temporal08 import design_bootstrap_weights, weighted_auc
from .replay_transport12 import MODEL_SHA, MODEL_PATH

PLAN_COMMIT = '95458703478ff686ad29b0ae20bcd4430e682309'
SOURCE_MANIFEST_SHA = '0fd86c5dca2b55ec9a9a869f513ac9ae72ccb4c46020fce1f22ef618b12803d3'
BASELINE_REFERENCE_SHA = '0ea29545c98545506a758edc590c26817f681b505d6f7ab8e360c2179c2ea930'
SECONDARY_SHA = '7c7dc918f042ff44f141a9ad3b49e6573db7f9732c5eb9d93e29b09b0686fa24'
RAW = {'total_cholesterol': ('LBXTC', 1.), 'hdl': ('LBDHDD', 1.),
       'hba1c': ('LBXGH', 1.), 'creatinine': ('LBXSCR', 1.),
       'uacr': ('URXUMA', None), 'albumin': ('LBXSAL', 10.),
       'rdw': ('LBXRDW', 1.), 'wbc': ('LBXWBCSI', 1.)}
SENTINELS = ('rdw_add_+0.5', 'albumin_mul_0.95', 'hba1c_add_+0.2', 'cbc_joint_positive')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scenarios():
    rows = [{'id': 'baseline', 'changes': {}}]
    for marker in RAW:
        for factor in (.9, .95, 1.05, 1.1):
            rows.append({'id': f'{marker}_mul_{factor:g}',
                         'changes': {marker: {'factor': factor, 'offset': 0.}}})
    for marker, offsets in [('rdw', (-.5, -.3, .3, .5)),
                            ('hba1c', (-.4, -.2, .2, .4)), ('albumin', (-2., 2.))]:
        for offset in offsets:
            rows.append({'id': f'{marker}_add_{offset:+g}',
                         'changes': {marker: {'factor': 1., 'offset': offset}}})
    for name, off, factor in [('positive', .3, 1.05), ('negative', -.3, .95)]:
        rows.append({'id': 'cbc_joint_' + name, 'changes': {
            'rdw': {'factor': 1., 'offset': off}, 'wbc': {'factor': factor, 'offset': 0.}}})
    return rows


def perturb(frame, changes):
    """Alter raw assay result, then recalculate ALL derived predictors."""
    if not isinstance(changes, dict) or set(changes) - set(RAW):
        raise ValueError('Unsupported measurement perturbation')
    d = frame.copy()
    for marker, spec in changes.items():
        if not isinstance(spec, dict) or set(spec) != {'factor', 'offset'}:
            raise ValueError('Explicit factor and canonical offset required')
        factor, offset = spec['factor'], spec['offset']
        if any(type(v) not in (int, float) or not np.isfinite(v) for v in (factor, offset)) or factor <= 0:
            raise ValueError('Invalid perturbation magnitude')
        raw, scale = RAW[marker]
        if marker == 'uacr':
            if offset != 0: raise ValueError('Only proportional UACR perturbations are specified')
            d[raw] = d[raw] * factor  # URXUCR unchanged; UACR recomputed from both inputs.
        else:
            d[raw] = d[raw] * factor + offset / scale
        if not np.isfinite(d[raw]).all() or (d[raw] <= 0).any():
            raise ValueError('Scenario outside positive measurement domain; no row dropping')
    return engineer_no_crp(d)


def weighted_quantile(values, weights, probability):
    v, w = np.asarray(values, float), np.asarray(weights, float)
    if v.ndim != 1 or v.shape != w.shape or not len(v) or not 0 <= probability <= 1:
        raise ValueError('Invalid quantile arrays')
    if not np.isfinite(v).all() or not np.isfinite(w).all() or np.any(w < 0) or w.sum() <= 0:
        raise ValueError('Invalid quantile support')
    order = np.argsort(v, kind='mergesort'); v, w = v[order], w[order]
    positive = w > 0; v, w = v[positive], w[positive]
    index = min(np.searchsorted(np.cumsum(w), probability * w.sum(), side='left'), len(v)-1)
    return float(v[index])


def clipping(frame, model):
    pp = model['preprocessing']; x = frame[pp['features']].to_numpy(float)
    return np.any((x < np.array(pp['lower'])) | (x > np.array(pp['upper'])), axis=1)


def metric_row(target, prediction, baseline, weights):
    y = np.asarray(target) > 0; w = np.asarray(weights, float)
    r, b = 1 - prediction[:, 0], 1 - baseline[:, 0]; delta = r - b
    return {'n': len(y), 'events': int(y.sum()),
        'six_state_brier': multiclass_brier(target, prediction, w),
        'all_cause_auc': weighted_auc(y, r, w),
        'all_cause_brier': float(np.average((y-r)**2, weights=w)),
        'observed_weighted': float(np.average(y, weights=w)),
        'mean_risk': float(np.average(r, weights=w)),
        'mean_shift_pp': float(np.average(delta, weights=w)*100),
        'mean_absolute_shift_pp': float(np.average(np.abs(delta), weights=w)*100),
        'p95_absolute_shift_pp': weighted_quantile(np.abs(delta), w, .95)*100,
        'maximum_absolute_shift_pp': float(np.max(np.abs(delta))*100),
        'fraction_shift_gt_1pp': float(np.average(np.abs(delta)>0.01, weights=w))}


def conditional_draws(target, pred, base, w):
    r, b = 1-pred[:, 0], 1-base[:, 0]
    return {'delta_six_state_brier': multiclass_brier(target, pred, w)-multiclass_brier(target, base, w),
            'delta_auc': weighted_auc(target>0, r, w)-weighted_auc(target>0, b, w),
            'mean_shift_pp': float(np.average(r-b, weights=w)*100)}


def main(source, output, model_dir=MODEL_PATH.parent):
    source, output, model_dir = map(Path, (source, output, model_dir))
    if output.exists(): raise ValueError('New output directory required')
    primary, secondary = model_dir/'model_bundle12.json', model_dir/'model_sensitivity12.json'
    if sha(primary) != MODEL_SHA or sha(secondary) != SECONDARY_SHA: raise ValueError('Frozen model changed')
    manifest = json.loads((source/'manifest.json').read_text(encoding='utf-8'))
    if sha(source/'manifest.json') != SOURCE_MANIFEST_SHA:
        raise ValueError('Not the previously examined source version')
    for row in manifest['files']:
        if row['status'] != 'retrieved' or Path(row['file']).name != row['file'] or sha(source/row['file']) != row['sha256']:
            raise ValueError('Changed or unavailable source')
    models = {'routine_no_crp': json.loads(primary.read_text())['models']['routine_no_crp'],
              'no_cbc_sensitivity': json.loads(secondary.read_text())['model']}
    reference_path = Path(__file__).parent/'assay13/BASELINE_REFERENCE.json'
    if sha(reference_path) != BASELINE_REFERENCE_SHA: raise ValueError('Changed baseline reference')
    original = json.loads(reference_path.read_text(encoding='utf-8'))['rows']
    rows, causes, subgroups, intervals, audit = [], [], [], [], []
    grid = scenarios(); max_mass_error = 0.; snapshot_seen = set()
    for year, horizon, expected_n in [(2011, 5, 2666), (2013, 4, 3013)]:
        d, source_audit = load_evaluation(source, year, MODEL_SHA)
        if len(d) != expected_n or set(d.SEQN) & snapshot_seen: raise ValueError('Cohort changed or duplicated')
        snapshot_seen |= set(d.SEQN)
        target, w = labels(d, horizon), d.weight.to_numpy(float)
        bases = {name: predict(d, m, horizon) for name, m in models.items()}
        for name, p in bases.items():
            old = next(r for r in original if r['cycle'] == year and r['model'] == name)
            if abs(multiclass_brier(target, p, w)-old['six_state_brier']) > 1e-12:
                raise ValueError('Baseline not reproduced')
        predictions = {}
        for scenario in grid:
            changed = perturb(d, scenario['changes'])
            for name, model in models.items():
                p = predict(changed, model, horizon); base = bases[name]
                mass_error = float(np.max(np.abs(p.sum(axis=1)-1)))
                max_mass_error = max(max_mass_error, mass_error)
                if mass_error > 1e-10: raise ArithmeticError('Incoherent competing probabilities')
                if name == 'no_cbc_sensitivity' and set(scenario['changes']) <= {'rdw','wbc'}:
                    if not np.array_equal(p, base): raise ArithmeticError('Unused assay changed no-CBC model')
                item = {'cycle': year, 'horizon': horizon, 'model': name, 'scenario': scenario['id'],
                        **metric_row(target, p, base, w),
                        'any_feature_clipped_fraction': float(np.average(clipping(changed, model), weights=w))}
                rows.append(item)
                for k, cause in enumerate(GROUPS, 1):
                    causes.append({'cycle':year, 'horizon':horizon, 'model':name, 'scenario':scenario['id'],
                      'cause':cause, 'mean_risk':float(np.average(p[:, k],weights=w)),
                      'mean_shift_pp':float(np.average(p[:, k]-base[:, k],weights=w)*100)})
                if scenario['id'] in SENTINELS:
                    predictions[name, scenario['id']] = p
                    for group, mask in [('female',d.RIAGENDR.eq(2)), ('male',d.RIAGENDR.eq(1)),
                                       ('age40_59',d.RIDAGEYR.lt(60)), ('age60_79',d.RIDAGEYR.ge(60))]:
                        mask=mask.to_numpy();subgroups.append({'cycle':year,'horizon':horizon,'model':name,
                          'scenario':scenario['id'],'subgroup':group,**metric_row(target[mask],p[mask],base[mask],w[mask])})
        print('SCENARIOS_COMPLETED',year,len(grid),flush=True)
        rng=np.random.default_rng(20260913+year); draws={key:[] for key in predictions}
        for _ in range(1000):
            wr=design_bootstrap_weights(d,rng)
            for key,pred in predictions.items():draws[key].append(conditional_draws(target,pred,bases[key[0]],wr))
        for (name,scenario),draw in draws.items():
            for metric in draw[0]:
                lo,hi=np.quantile([r[metric] for r in draw],[.025,.975])
                intervals.append({'cycle':year,'horizon':horizon,'model':name,'scenario':scenario,'metric':metric,
                  'lower':float(lo),'upper':float(hi),'replicates':1000,
                  'interval_kind':'conditional_PSU_percentiles_not_assay_uncertainty'})
        audit.append(source_audit)
    output.mkdir(parents=True)
    for filename,data in [('metrics.csv',rows),('cause_shifts.csv',causes),('subgroup_sensitivity.csv',subgroups),('intervals.csv',intervals)]:
        pd.DataFrame(data).to_csv(output/filename,index=False)
    summary={'version':'0.13_assay_stress','analysis_plan_commit':PLAN_COMMIT,
       'analysis_type':'post_validation_simulated_measurement_sensitivity_on_real_profiles',
       'new_independent_validation':False,'new_model_fit':False,'actual_paired_assay_data':False,
       'clinical_use_ready':False,'risk_change_threshold_pp':1.,'threshold_is_clinical':False,
       'scenario_count':len(grid),'models':list(models),'profiles_reused':sum(a['n_complete'] for a in audit),
       'metric_rows':len(rows),'cause_rows':len(causes),'subgroup_rows':len(subgroups),
       'interval_rows':len(intervals),'source_files_verified':len(manifest['files']),
       'maximum_probability_mass_error':max_mass_error,'primary_sha256':MODEL_SHA,
       'secondary_sha256':SECONDARY_SHA,'source_manifest_sha256':sha(source/'manifest.json'),
       'code_sha256':sha(__file__),'bootstrap_replicates':1000,
       'excludes_uncertainty':'model_fitting_and_real_assay_parameter_estimation','raw_records_exported':False,
       'source_audit':audit}
    for filename,doc in [('RESULTS_SUMMARY.json',summary),('SCENARIOS.json',grid)]:
        (output/filename).write_text(json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(main(args.source,args.out),indent=2))
