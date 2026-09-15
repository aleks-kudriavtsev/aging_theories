"""Post-hoc design-bootstrap refitting of NfL22, not new validation.

The baseline, ridge penalty, old folds and observed population remain fixed.
No individual records or predictions are exported.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from ..nfl22.analyze import (cohort,load_model,baseline_hazard,assign_folds,
    fit_update,predict_update,sha,MODEL_SHA,design_bootstrap_weights,binary_metrics)
SEED=20260923
PLAN_COMMIT='f4e6d73a4dd2d68d5fc69d85dae87e55e39b7110'

def partition(clusters,salt):
    """Only cluster identity enters this deterministic partition."""
    if not isinstance(salt,str):raise ValueError('String salt required')
    keys=sorted(set(clusters),key=lambda x:hashlib.sha256((salt+':'+str(x)).encode()).hexdigest())
    if len(keys)<5:raise ValueError('Need at least five clusters')
    mapping={key:i%5 for i,key in enumerate(keys)}
    return np.array([mapping[k] for k in clusters],dtype=int)

def crossfit(y,hstop,horizon,w,x,folds):
    arrays=[np.asarray(v) for v in (y,hstop,horizon,w,x,folds)]
    if arrays[0].ndim!=1 or not len(arrays[0]) or any(a.shape!=arrays[0].shape for a in arrays):
        raise ValueError('Equal nonempty one-dimensional arrays required')
    y,hstop,horizon,w,x,folds=arrays
    if not np.isfinite(np.concatenate([a.astype(float) for a in arrays])).all():raise ValueError('Nonfinite data')
    if not np.isin(folds,range(5)).all() or len(set(folds))!=5:raise ValueError('Five folds required')
    p={name:np.full(len(y),np.nan) for name in ('calibration','nfl')};fit_rows=[]
    for fold in range(5):
        tr=folds!=fold;te=~tr
        for name,use in [('calibration',False),('nfl',True)]:
            fit=fit_update(y[tr],hstop[tr],w[tr],x[tr],use)
            p[name][te]=predict_update(horizon[te],x[te],fit)
            fit_rows.append({'fold':fold,'model':name,'training_events':int(y[tr].sum()),
                'positive_weight_training_events':int(((y==1)&(w>0)&tr).sum()),**fit})
    if any(not np.isfinite(a).all() or np.any((a<0)|(a>1)) for a in p.values()):raise ArithmeticError('Invalid OOF probabilities')
    return p,fit_rows

def contrast(y,p,w):
    a,b=[binary_metrics(y,p[name],w) for name in ('nfl','calibration')]
    return {'delta_brier':a['brier']-b['brier'],
        'delta_auc':a['auc']-b['auc'] if a['auc'] is not None and b['auc'] is not None else None,
        'calibration_brier':b['brier'],'nfl_brier':a['brier'],'calibration_auc':b['auc'],'nfl_auc':a['auc']}

def run(source,out,replicates=1000,split_count=20):
    out=Path(out)
    if out.exists():raise ValueError('New output directory required')
    if replicates!=1000 or split_count!=20:raise ValueError('Declared1000draws/20partitions required')
    m=load_model();eligible,audit=cohort(Path(source))
    common=np.isfinite(eligible[m['preprocessing']['features']]).all(axis=1).to_numpy()
    d=eligible.loc[common].copy().reset_index(drop=True);y=((d.dead==1)&(d.time_years<=4)).to_numpy(int)
    if len(d)!=1270 or y.sum()!=37:raise ValueError('Prior complete domain changed')
    w=d.weight.to_numpy(float);x=np.log2(d.SSSNFL.to_numpy(float))
    hstop=baseline_hazard(d,m,np.minimum(d.time_years,4));h4=baseline_hazard(d,m,4)
    folds=assign_folds(eligible.cluster)[common];p,original_fits=crossfit(y,hstop,h4,w,x,folds)
    original=contrast(y,p,w)
    if abs(original['delta_brier']-.00032470)>1e-8 or abs(original['delta_auc']+.0087163)>1e-7:raise ValueError('Old contrast mismatch')
    print('ORIGINAL_REPRODUCED',original,flush=True)
    rng=np.random.default_rng(SEED);draws=[]
    for rep in range(replicates):
        wr=design_bootstrap_weights(eligible,rng)[common];row={'replicate':rep,'refit_status':'computed'}
        row.update({'fixed_'+k:v for k,v in contrast(y,p,wr).items()})
        try:
            q,fr=crossfit(y,hstop,h4,wr,x,folds)
            row.update({'refit_'+k:v for k,v in contrast(y,q,wr).items()})
            row['negative_beta_folds']=sum(r['log2_nfl_coefficient']<0 for r in fr if r['model']=='nfl')
        except (ValueError,ArithmeticError,RuntimeError) as error:
            row['refit_status']='undefined';row['error']=type(error).__name__+': '+str(error)
        draws.append(row)
        if (rep+1)%100==0:print('REFITS',rep+1,flush=True)
    frame=pd.DataFrame(draws);intervals=[]
    for mode in ('fixed','refit'):
        for metric in ('delta_brier','delta_auc'):
            vals=frame[mode+'_'+metric].dropna().to_numpy();lo,hi=np.quantile(vals,[.025,.975]) if len(vals) else (None,None)
            intervals.append({'approach':mode,'metric':metric,'point_estimate':original[metric],
                'lower':lo,'upper':hi,'valid_replicates':len(vals),'undefined_replicates':replicates-len(vals),
                'interpretation':'conditional_fixed_baseline_and_partition' if len(vals)==replicates else 'defined_draws_only_not_nominal_coverage',
                'update_refitted':mode=='refit'})
    splits=[]
    for s in range(split_count):
        fs=partition(eligible.cluster,f'joint23:split:{s}')[common]
        try:
            q,fr=crossfit(y,hstop,h4,w,x,fs);row={'split':s,'status':'computed',**contrast(y,q,w)}
            betas=[r['log2_nfl_coefficient'] for r in fr if r['model']=='nfl'];row.update(nfl_beta_min=min(betas),nfl_beta_max=max(betas))
        except (ValueError,ArithmeticError,RuntimeError) as error:row={'split':s,'status':'undefined','error':str(error)}
        splits.append(row)
    sf=pd.DataFrame(splits)
    result={'date':'2026-09-15','plan_commit':PLAN_COMMIT,'analysis':'post_hoc_robustness_not_new_validation',
        'n':len(d),'events4':int(y.sum()),'baseline_sha256':MODEL_SHA,'source_audit':audit,
        'old_source_reused':True,'point':original,'bootstrap_replicates':replicates,'seed':SEED,
        'refit_failures':int((frame.refit_status!='computed').sum()),
        'additional_partitions':split_count,'partition_failures':int((sf.status!='computed').sum()),
        'partitions_with_brier_improvement':int((sf.delta_brier<0).sum()),
        'partitions_with_auc_improvement':int((sf.delta_auc>0).sum()),
        'partition_delta_brier_range':[float(sf.delta_brier.min()),float(sf.delta_brier.max())],
        'partition_delta_auc_range':[float(sf.delta_auc.min()),float(sf.delta_auc.max())],
        'baseline_fitting_uncertainty_included':False,'partition_uncertainty_in_interval':False,
        'public_model_updated':False,'clinical_use_ready':False,'participant_records_exported':False,'code_sha256':sha(__file__)}
    out.mkdir(parents=True)
    for name,t in [('intervals',pd.DataFrame(intervals)),('bootstrap_aggregates',frame),('partition_sensitivity',sf),('original_fold_parameters',pd.DataFrame(original_fits))]:t.to_csv(out/(name+'.csv'),index=False)
    (out/'RESULTS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(pd.DataFrame(intervals).to_string(index=False),flush=True);print(json.dumps(result,indent=2),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
