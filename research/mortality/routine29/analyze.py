"""Exploratory renal/routine comparison on a reused cohort, fully nested refits."""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from ..binary28.analyze import assign_outer,intervals
from ..binary28.metrics import metrics,calibration_bins
from ..binary28.cohort import PROTEINS
from ..validate_temporal08 import design_bootstrap_weights
from ..nfl22.analyze import sha,dump
from .cohort import load,ROOT
from .features import CORE
from .model import nested_crossfit
ARMS=('simple','simple_all3','routine','routine_all3','rank_all3','clip_all3')
PRIMARY=('routine','routine_all3')
REPLICATES=500
SEED=20261101


def pair_scores(y,p,w):
    row={}
    for name in PRIMARY:
        m=metrics(y,p[name],w)
        for key in ('auc','brier','observed','predicted','bias_pp','observed_expected','calibration_in_large','calibration_intercept','calibration_slope'):
            row[name+'|'+key]=m[key]
        row[name+'|calibration_undefined']=int(m['calibration_status']!='computed')
    for key in ('auc','brier'):
        a,b=row[PRIMARY[1]+'|'+key],row[PRIMARY[0]+'|'+key]
        row['delta_'+key]=a-b if a is not None and b is not None else None
    return row


def partition(clusters,index):
    keys=sorted(set(clusters),key=lambda x:hashlib.sha256((f'routine29:partition:{index}:'+str(x)).encode()).hexdigest())
    if len(keys)<15:raise ValueError('Insufficient PSU support')
    table={key:i%5 for i,key in enumerate(keys)}
    return np.array([table[c] for c in clusters],int)


def fit_pair(core,protein,y,w,clusters,folds,domain):
    predictions={};parameters=[]
    for name in PRIMARY:
        predictions[name],rows=nested_crossfit(core,protein,y,w,clusters,folds,name)
        parameters.extend([dict(analysis_domain=domain,**r) for r in rows])
    return predictions,parameters


def run(source,out):
    out=Path(out)
    if out.exists():raise ValueError('New output directory required')
    eligible,common,strict,audit,manifest=load(source);d=eligible.loc[common].copy()
    core=d[list(CORE)].to_numpy(float);protein=np.log2(d[list(PROTEINS)].to_numpy(float))
    y=(d.dead.eq(1)&d.t_exam.le(96)).to_numpy(int);w=d.weight.to_numpy(float);clusters=d.cluster.to_numpy()
    folds_all=assign_outer(eligible.cluster);folds=folds_all[common]
    if len(d)==0 or np.unique(y).size!=2:raise ValueError('Empty/no-event comparison domain')
    predictions={};point={};parameters=[];bins=[]
    for name in ARMS:
        p,rows=nested_crossfit(core,protein,y,w,clusters,folds,name)
        predictions[name]=p;point[name]=metrics(y,p,w);parameters.extend([dict(analysis_domain='main',**r) for r in rows])
        for b in calibration_bins(y,p,w):bins.append(dict(model=name,analysis_domain='main',**b))
    primary=pair_scores(y,predictions,w)
    rng=np.random.default_rng(SEED);draws=[];selection={name:{str(v):0 for v in (1.,10.,100.)} for name in PRIMARY}
    print('ROUTINE29_POINT_COMPUTED',json.dumps({'n':len(d),'events8':int(y.sum()),'delta_brier':primary['delta_brier'],'delta_auc':primary['delta_auc']}),flush=True)
    for rep in range(REPLICATES):
        wr=design_bootstrap_weights(eligible,rng)[common]
        try:
            p,rows=fit_pair(core,protein,y,wr,clusters,folds,'bootstrap')
            row=dict(replicate=rep,status='computed',**pair_scores(y,p,wr))
            for fit in rows:selection[fit['model']][str(float(fit['lambda']))]+=1
        except (ValueError,ArithmeticError,RuntimeError) as error:
            row=dict(replicate=rep,status='undefined',error=type(error).__name__+': '+str(error))
        draws.append(row)
        if (rep+1)%50==0:print('ROUTINE29_NESTED_REFITS',rep+1,flush=True)
    partitions=[]
    for index in range(5):
        f=partition(eligible.cluster,index)[common]
        try:
            p,rows=fit_pair(core,protein,y,w,clusters,f,f'partition_{index}')
            parameters.extend(rows);partitions.append(dict(partition=index,status='computed',**pair_scores(y,p,w)))
        except (ValueError,ArithmeticError,RuntimeError) as error:
            partitions.append(dict(partition=index,status='undefined',error=type(error).__name__+': '+str(error)))
    strictr=dict(status='undefined',n=int(strict.sum()),interpretation='different_conditional_domain_within_domain_pair_only')
    try:
        s=eligible.loc[strict].copy();ys=(s.dead.eq(1)&s.t_exam.le(96)).to_numpy(int);ws=s.weight.to_numpy(float)
        p,rows=fit_pair(s[list(CORE)].to_numpy(float),np.log2(s[list(PROTEINS)].to_numpy(float)),ys,ws,s.cluster.to_numpy(),folds_all[strict],'strict_metadata')
        parameters.extend(rows);strictr.update(status='computed',events8=int(ys.sum()),point={k:metrics(ys,v,ws) for k,v in p.items()},primary=pair_scores(ys,p,ws))
    except (ValueError,ArithmeticError,RuntimeError) as error:strictr['error']=type(error).__name__+': '+str(error)
    primary_csv=pd.DataFrame(draws).reindex(columns=['replicate','status','delta_brier','delta_auc']).to_csv(index=False)
    valid=[r for r in partitions if r['status']=='computed']
    stress={metric:dict(minimum=min([r[metric] for r in valid if r.get(metric) is not None],default=None),
                            maximum=max([r[metric] for r in valid if r.get(metric) is not None],default=None)) for metric in ('delta_brier','delta_auc')}
    result=dict(analysis='exploratory_reused_NHANESIII_new_routine_binary_model',plan_commit='7800458',horizon_years=8,
        n=len(d),events8=int(y.sum()),weight='WTCYPEX6',effective_n=float(w.sum()**2/(w@w)),
        strata_n=int(d.stratum.nunique()),psu_n=len(set(clusters)),point=point,primary_contrast='routine_all3_minus_routine',primary=primary,
        primary_intervals=intervals(draws,REPLICATES),bootstrap_replicates=REPLICATES,bootstrap_seed=SEED,
        bootstrap_failures=sum(r['status']=='undefined' for r in draws),bootstrap_penalty_selection_counts=selection,
        partition_sensitivity=partitions,partition_sensitivity_range=stress,partition_failures=5-len(valid),
        partition_sensitivity_pooled_into_CI=False,strict_metadata_sensitivity=strictr,
        common_weighted_coverage=float(np.average(common,weights=eligible.weight)),source_audit=audit,
        routine_features=list(CORE),protein_features=list(PROTEINS),
        serum_creatinine_equation='0.960*originalCEP-0.184mg/dL_serum_only',eGFR='estimated_2021_CKDEPI_creatinine_only_not_measured_GFR',
        assay_and_creatinine_calibration_uncertainty_included=False,reported_protein_QC_LOQ_validated=False,
        independent_external_validation=False,clinical_use_ready=False,participant_records_exported=False,
        same_common_domain_all_primary_and_simple_comparators=True,no_single_protein_winner_search=True,
        primary_bootstrap_csv_sha256=hashlib.sha256(primary_csv.encode()).hexdigest(),
        plan_sha256=sha(ROOT/'ANALYSIS_PLAN.md'),lock_sha256=sha(ROOT/'LOCK.json'),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)):sha(p) for p in sorted(ROOT.parent.rglob('*.py'))},
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
    out.mkdir(parents=True);dump(out/'RESULTS.json',result);dump(out/'SOURCE_MANIFEST.json',manifest);dump(out/'fold_parameters.json',parameters)
    pd.DataFrame(draws).to_csv(out/'bootstrap_aggregates.csv',index=False);(out/'primary_bootstrap.csv').write_text(primary_csv)
    dump(out/'calibration_bins.json',bins)
    for label,value in [('RESULTS',result),('MANIFEST',manifest),('PARAMETERS',parameters),('CALIBRATION_BINS',bins)]:
        print('ROUTINE29_'+label+'_BEGIN');print(json.dumps(value,indent=2,allow_nan=False));print('ROUTINE29_'+label+'_END')
    print('ROUTINE29_PRIMARY_BOOTSTRAP_CSV_BEGIN');print(primary_csv,end='');print('ROUTINE29_PRIMARY_BOOTSTRAP_CSV_END')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
