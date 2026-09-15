"""Post-result audit of PSU-key serialization, not a new validation.

Both outcome-blind allocations are reported, never selected on performance.
The original implementation and its original result file remain untouched.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from . import analyze as a


def assess(source,out):
    source,out=Path(source),Path(out)
    if out.exists():raise ValueError('A new result directory is required')
    model=a.load_model();eligible,audit=a.cohort(source)
    mask=np.isfinite(eligible[model['preprocessing']['features']]).all(axis=1).to_numpy()
    d=eligible[mask].reset_index(drop=True);w=d.weight.to_numpy();x=np.log2(d.SSSNFL.to_numpy())
    y=((d.dead==1)&(d.time_years<=4)).to_numpy(int)
    H=a.baseline_hazard(d,model,4);Hstop=a.baseline_hazard(d,model,np.minimum(d.time_years,4))
    schemes={'unprefixed_PSU_key':a.assign_folds(eligible.cluster)[mask],
             'prefixed_PSU_key':a.assign_folds('2013:'+eligible.cluster)[mask]}
    predictions={};rows=[];intervals=[]
    for allocation,folds in schemes.items():
        preds={k:np.zeros(len(d)) for k in ('calibration_only','calibration_plus_NfL')}
        for fold in range(5):
            train=folds!=fold;test=~train
            for panel,include in [('calibration_only',False),('calibration_plus_NfL',True)]:
                fit=a.fit_update(y[train],Hstop[train],w[train],x[train],include)
                preds[panel][test]=a.predict_update(H[test],x[test],fit)
        predictions[allocation]=preds
        for panel,p in preds.items():
            rows.append({'allocation':allocation,'panel':panel,**a.binary_metrics(y,p,w)})
    rng=np.random.default_rng(a.SEED);draws=[]
    for _ in range(1000):
        wr=a.design_bootstrap_weights(eligible,rng)[mask];row={}
        for allocation,preds in predictions.items():
            metrics={panel:a.binary_metrics(y,p,wr) for panel,p in preds.items()}
            for metric in ('brier','auc'):
                u,v=metrics['calibration_plus_NfL'][metric],metrics['calibration_only'][metric]
                row[allocation+'|'+metric]=u-v if u is not None and v is not None else None
        draws.append(row)
    draws=pd.DataFrame(draws)
    for key in draws:
        allocation,metric=key.split('|');valid=draws[key].dropna()
        bounds=np.quantile(valid,[.025,.975]).tolist() if len(valid) else [None,None]
        interval={'allocation':allocation,'metric':metric,'lower':bounds[0],'upper':bounds[1],
                  'valid_draws':len(valid),'undefined_draws':1000-len(valid)}
        intervals.append(interval)
    summary={'date':'2026-09-15','design':'post_result_fold_sensitivity_not_new_validation',
        'n':len(d),'deaths4':int(y.sum()),'model_sha256':a.MODEL_SHA,
        'changed_fold_participants':int(np.sum(schemes['unprefixed_PSU_key']!=schemes['prefixed_PSU_key'])),
        'hash_key_difference':'nfl22:<stratum>:<PSU> versus nfl22:2013:<stratum>:<PSU>',
        'metrics':rows,'intervals':intervals,'bootstrap_replicates':1000,
        'interval_scope':'conditional_on_fixed_cross_fitted_predictions_excludes_all_fitting_and_fold_selection_uncertainty',
        'automatic_preferred_result':None,'clinical_use_ready':False,'individual_rows_exported':False,
        'conclusion':'NfL raises point-estimate Brier in both allocations; incremental AUC changes sign. Do not select the favourable split.'}
    out.mkdir(parents=True)
    (out/'FOLD_SENSITIVITY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    pd.DataFrame(rows).to_csv(out/'fold_metrics.csv',index=False)
    pd.DataFrame(intervals).to_csv(out/'fold_intervals.csv',index=False)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();print(json.dumps(assess(args.source,args.out),indent=2))
