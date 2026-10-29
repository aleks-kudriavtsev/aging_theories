"""Prespecified nested whole-PSU binary model assessment with paired refitting."""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from ..nfl22.analyze import sha,dump
from ..validate_temporal08 import design_bootstrap_weights
from .cohort import load,CORE,PROTEINS,ROOT
from .model import nested_crossfit
from .metrics import metrics,calibration_bins
SEED=20261031
REPLICATES=500
ARMS=('core','all3','cys','b2m','btp','rank_all3','clip_all3')


def assign_outer(clusters):
    keys=sorted(set(clusters),key=lambda x:hashlib.sha256(('binary28:outer:'+str(x)).encode()).hexdigest())
    if len(keys)<15:raise ValueError('Insufficient PSU support for nested design')
    mapping={key:i%5 for i,key in enumerate(keys)}
    return np.array([mapping[c] for c in clusters],int)


def flatten_pair(y,p,w):
    out={}
    for name in ('core','all3'):
        m=metrics(y,p[name],w)
        for key in ('auc','brier','observed','predicted','bias_pp','observed_expected','calibration_in_large','calibration_intercept','calibration_slope'):
            out[name+'|'+key]=m[key]
        out[name+'|calibration_undefined']=int(m['calibration_status']!='computed')
    for key in ('auc','brier'):
        a,b=out['all3|'+key],out['core|'+key]
        out['delta_'+key]=a-b if a is not None and b is not None else None
    return out


def intervals(draws,total):
    keys=sorted({'delta_brier','delta_auc'}|{k for d in draws for k in d if k not in ('replicate','status','error') and not k.endswith('calibration_undefined')})
    result={}
    for key in keys:
        vals=[d[key] for d in draws if d.get(key) is not None and np.isfinite(d[key])]
        lo,hi=map(float,np.quantile(vals,[.025,.975])) if vals else (None,None)
        result[key]=dict(lower=lo,upper=hi,valid=len(vals),undefined=total-len(vals),
            scope='nested_preprocessing_tuning_refitting_with_fixed_partition',defined_draws_only=len(vals)!=total)
    return result


def run(source,out):
    out=Path(out)
    if out.exists():raise ValueError('New output directory required')
    eligible,common,audit,manifest=load(source);d=eligible.loc[common].copy()
    core=d[list(CORE)].to_numpy(float);protein=np.log2(d[list(PROTEINS)].to_numpy(float))
    y=(d.dead.eq(1)&d.t_exam.le(96)).to_numpy(int);w=d.weight.to_numpy(float);clusters=d.cluster.to_numpy()
    folds=assign_outer(eligible.cluster)[common]
    predictions={};point={};parameters=[];bins=[]
    for name in ARMS:
        p,rows=nested_crossfit(core,protein,y,w,clusters,folds,name)
        predictions[name]=p;point[name]=metrics(y,p,w);parameters.extend(rows)
        for row in calibration_bins(y,p,w):bins.append(dict(model=name,**row))
    primary=flatten_pair(y,predictions,w)
    print('BINARY28_POINT_COMPUTED',json.dumps({'n':len(d),'events8':int(y.sum()),'delta_brier':primary['delta_brier'],'delta_auc':primary['delta_auc']}),flush=True)
    rng=np.random.default_rng(SEED);draws=[];selection_counts={name:{str(lam):0 for lam in (1.,10.,100.)} for name in ('core','all3')}
    for rep in range(REPLICATES):
        wr=design_bootstrap_weights(eligible,rng)[common]
        try:
            p={};fit_rows={}
            for name in ('core','all3'):
                p[name],rows=nested_crossfit(core,protein,y,wr,clusters,folds,name)
                fit_rows[name]=rows
            row=dict(replicate=rep,status='computed',**flatten_pair(y,p,wr))
            for name,rows in fit_rows.items():
                for fit in rows:selection_counts[name][str(float(fit['lambda']))]+=1
        except (ValueError,ArithmeticError,RuntimeError) as exc:
            row=dict(replicate=rep,status='undefined',error=type(exc).__name__+': '+str(exc))
        draws.append(row)
        if (rep+1)%50==0:print('BINARY28_NESTED_REFITS',rep+1,flush=True)
    result=dict(analysis='new_reported_value_binary_logistic_internal_nested_PSU_CV',cohort='NHANESIII1988_1994_first_exam',
        plan_commit='d81972a',horizon_years=8,n=len(d),events8=int(y.sum()),weight='WTCYPEX6',
        effective_n=float(w.sum()**2/(w@w)),strata_n=int(d.stratum.nunique()),psu_n=len(set(clusters)),
        point=point,primary_contrast='all3_minus_core',primary=primary,intervals=intervals(draws,REPLICATES),
        bootstrap_replicates=REPLICATES,bootstrap_seed=SEED,bootstrap_failures=sum(r['status']=='undefined' for r in draws),
        bootstrap_penalty_selection_counts=selection_counts,source_audit=audit,
        clinical_core=list(CORE),protein_features=list(PROTEINS),
        exploratory_arms=list(ARMS[2:]),exploratory_intervals_computed=False,
        model_family='eight_year_survey_weighted_ridge_logistic_not_cause_specific_engine',
        comparison_scope='reported_values_beyond_age_sex_smoking_SBP_not_full_renal_or_clinical_panel',
        all_three_reported_values_are_not_QC_LOQ_validated=True,independent_external_validation=False,
        no_legacy_model_coefficients_used=True,clinical_use_ready=False,participant_records_exported=False,
        plan_sha256=sha(ROOT/'ANALYSIS_PLAN.md'),lock_sha256=sha(ROOT/'LOCK.json'),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)):sha(p) for p in sorted(ROOT.parent.rglob('*.py'))},
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
    out.mkdir(parents=True)
    dump(out/'RESULTS.json',result);dump(out/'SOURCE_MANIFEST.json',manifest);dump(out/'fold_parameters.json',parameters)
    pd.DataFrame(draws).to_csv(out/'bootstrap_aggregates.csv',index=False);pd.DataFrame(bins).to_csv(out/'calibration_bins.csv',index=False)
    print('BINARY28_RESULTS_BEGIN');print(json.dumps(result,indent=2,allow_nan=False));print('BINARY28_RESULTS_END')
    print('BINARY28_MANIFEST_BEGIN');print(json.dumps(manifest,indent=2,allow_nan=False));print('BINARY28_MANIFEST_END')
    print('BINARY28_PARAMETERS_BEGIN');print(json.dumps(parameters,indent=2,allow_nan=False));print('BINARY28_PARAMETERS_END')
    print('BINARY28_CALIBRATION_BINS_BEGIN');print(json.dumps(bins,indent=2,allow_nan=False));print('BINARY28_CALIBRATION_BINS_END')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
