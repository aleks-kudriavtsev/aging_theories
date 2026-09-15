"""NfL incremental information with cluster cross-fitting, fixed baseline.

ANALYSIS_LOCK.md and its provenance correction must accompany these results.
No individual records are exported. This is not a temporal/external validation
of the newly learned NfL coefficient. Bootstrap intervals condition on OOF fits.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from ..transport_no_crp12 import engineer_no_crp, COMPONENTS
from ..xpt_checked import read_xpt_checked
from ..validate_temporal08 import independent_mortality_reader, audit_mortality_reader, design_bootstrap_weights, diagnostics
from ..competing_risks09 import transform, binary_metrics

ROOT = Path(__file__).parent
MODEL_PATH = ROOT/'frozen_expanded13.json'
MODEL_SHA = '87a8623b9fd5cb793ce2da5311a17a56b7d94251da5df335de3d406e437a55da'
PARENT_SHA = 'ed37a1b0e8be5eaa40556eb0c169166fb0f06a3d52dfd5ae6882d608788562ec'
HORIZONS = (1,4)
SEED = 20260922

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def check_sources(root):
    rows=json.loads((root/'manifest.json').read_text())['files'];seen=set()
    required={f'2013_{c}_H.xpt' for c in (*COMPONENTS,'SSSNFL')}|{'2013_mortality.dat'}
    valid=[]
    for row in rows:
        name=row['file']
        if name in seen or Path(name).name!=name:raise ValueError('Duplicate/unsafe input')
        seen.add(name)
        if row['status']!='retrieved':
            if name in required:raise ValueError('Required source failed: '+name)
            continue
        path=root/name
        if not path.is_file() or path.stat().st_size!=row['bytes'] or sha(path)!=row['sha256']:
            raise ValueError('Source version changed: '+name)
        valid.append(name)
    if not required<=set(valid):raise ValueError('Missing required input')
    return {'files_verified':len(valid),'scientific_files':len(required),'source_manifest_sha256':sha(root/'manifest.json')}

def load_model(path=MODEL_PATH):
    if sha(path)!=MODEL_SHA:raise ValueError('Wrong frozen baseline')
    m=json.loads(Path(path).read_text())
    if m['interval_ends_years']!=[5.,8.] or set(m['causes'])!=set(map(str,range(1,11))):
        raise ValueError('Invalid baseline model layout')
    return m

def engineer_expanded(d):
    d=engineer_no_crp(d)
    for key,col in [('log_alp','LBXSAPSI'),('log_ggt','LBXSGTSI'),('log_bilirubin','LBXSTB'),('log_platelets','LBXPLTSI')]:
        d[key]=np.log2(d[col].where(d[col]>0))
    ratio=d.LBXNEPCT.where(d.LBXNEPCT>0)/d.LBXLYPCT.where(d.LBXLYPCT>0)
    d['log_nlr']=np.log2(ratio)
    return d

def cohort(source):
    source=Path(source);audit=check_sources(source)
    demo=read_xpt_checked(source/'2013_DEMO_H.xpt')
    if demo.SEQN.isna().any() or demo.SEQN.duplicated().any():raise ValueError('Invalid demographic IDs')
    ids=set(demo.SEQN)
    d=demo.loc[demo.RIDAGEYR.between(40,75),['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU']].copy()
    age_n=len(d)
    for c in (*COMPONENTS[1:],'SSSNFL'):
        raw=read_xpt_checked(source/f'2013_{c}_H.xpt')
        if raw.SEQN.isna().any() or raw.SEQN.duplicated().any() or not set(raw.SEQN)<=ids:
            raise ValueError('Invalid covariate identifiers: '+c)
        if c=='SSSNFL':
            nfl_audit={'source_rows':len(raw),'source_measured':int(raw.SSSNFL.notna().sum()),
                       'below_LLOQ':int(raw.SSNFLL.eq(1).sum()),'above_ULOQ':int(raw.SSNFLH.eq(1).sum())}
            if not raw.SSNFLL.dropna().isin([0,1]).all() or not raw.SSNFLH.dropna().isin([0,1]).all():raise ValueError('Unknown NfL limit flag')
            if (raw.SSNFLL.eq(1)&raw.SSNFLH.eq(1)).any():raise ValueError('Both NfL limits flagged')
        d=d.merge(raw[['SEQN']+[x for x in raw if x not in d]],on='SEQN',how='left',validate='1:1')
    mort=independent_mortality_reader(source/'2013_mortality.dat')
    if set(mort.SEQN)!=ids:raise ValueError('Mortality IDs differ')
    d=d.merge(mort,on='SEQN',how='left',validate='1:1')
    eligible=(d.elig.eq(1)&d.dead.isin([0,1])&d.t_exam.notna()&d.t_exam.ge(0)&d.RIAGENDR.isin([1,2])
        &d.WTSSNH2Y.gt(0)&np.isfinite(d.WTSSNH2Y)&d.SSSNFL.gt(0)&np.isfinite(d.SSSNFL))
    d=d.loc[eligible].copy().reset_index(drop=True)
    if d[['SDMVSTRA','SDMVPSU','SSNFLL','SSNFLH']].isna().any().any():raise ValueError('Missing design/flags')
    d=engineer_expanded(d);d['weight']=d.WTSSNH2Y
    d['stratum']='2013:'+d.SDMVSTRA.astype(int).astype(str)
    d['cluster']=d.SDMVSTRA.astype(int).astype(str)+':'+d.SDMVPSU.astype(int).astype(str)
    d['time_years']=d.t_exam.clip(lower=.5)/12
    if ((d.dead==0)&(d.time_years<4)).any():raise ValueError('Early censoring: protocol amendment needed BEFORE metrics')
    audit.update(survey_n=len(demo),age40_75_n=age_n,eligible_nfl_n=len(d),nfl_source=nfl_audit,
        minimum_survivor_followup=float(d.loc[d.dead==0,'time_years'].min()),
        zero_month_deaths=int((d.dead.eq(1)&d.t_exam.eq(0)).sum()),
        strata=int(d.stratum.nunique()),clusters=int(d.cluster.nunique()),
        reader_audit=audit_mortality_reader(source/'2013_mortality.dat'))
    return d,audit

def assign_folds(clusters):
    # Deterministic cluster-only allocation. Neither outcome nor analyte is read.
    keys=sorted(set(clusters),key=lambda s:hashlib.sha256(('nfl22:'+str(s)).encode()).hexdigest())
    if len(keys)<5:raise ValueError('At least five clusters required')
    mapping={key:i%5 for i,key in enumerate(keys)}
    return np.array([mapping[x] for x in clusters],int)

def baseline_hazard(d,m,times):
    z=transform(d,m['preprocessing']);t=np.broadcast_to(np.asarray(times,float),(len(d),))
    if not np.isfinite(t).all() or np.any((t<0)|(t>8)):raise ValueError('Unsupported times')
    dt=np.column_stack([np.minimum(t,5),np.maximum(t-5,0)])
    H=np.zeros(len(d))
    for c in m['causes'].values():
        rates=np.asarray(c['baseline_rates']);lp=z@np.asarray(c['coefficients'])
        H+=np.exp(lp)*(dt@rates)
    if not np.isfinite(H).all() or np.any(H<0):raise ValueError('Invalid baseline hazard')
    return H

def fit_update(y,h,w,log_nfl,with_nfl):
    y,h,w,x=map(lambda a:np.asarray(a,float),(y,h,w,log_nfl))
    if y.ndim!=1 or len(y)==0 or any(a.shape!=y.shape for a in (h,w,x)):
        raise ValueError('Invalid update shape')
    if not np.isfinite(np.r_[y,h,w,x]).all() or not np.isin(y,[0,1]).all() or np.any(h<=0) or np.any(w<0) or w.sum()<=0 or (w@y)<=0:
        raise ValueError('Invalid update/event support')
    w=w/w.mean();center=float(np.average(x,weights=w));xc=x-center
    a0=float(np.log((w@y)/(w@h)))
    if not with_nfl:return {'intercept':a0,'log2_nfl_coefficient':0.,'log2_nfl_center':center,'gradient_max':0.,'ridge':0.}
    X=np.column_stack([np.ones(len(y)),xc])
    def fun(b):
        lp=X@b;mu=h*np.exp(lp)
        loss=float(w@(mu-y*lp)+.5*b[1]**2)
        grad=X.T@(w*(mu-y))+np.array([0.,b[1]])
        return loss,grad
    opt=minimize(fun,[a0,0.],jac=True,method='BFGS',options={'maxiter':500,'gtol':1e-7})
    gradient=float(np.max(np.abs(fun(opt.x)[1])))
    if not np.isfinite(opt.x).all() or gradient>1e-4:raise ValueError('NfL update did not converge')
    return {'intercept':float(opt.x[0]),'log2_nfl_coefficient':float(opt.x[1]),'log2_nfl_center':center,'gradient_max':gradient,'ridge':1.}

def predict_update(h,x,fit):
    lp=fit['intercept']+fit['log2_nfl_coefficient']*(np.asarray(x)-fit['log2_nfl_center'])
    return -np.expm1(-np.asarray(h)*np.exp(lp))

def run(source,out):
    source,out=Path(source),Path(out)
    if out.exists():raise ValueError('New result directory required')
    m=load_model();parent=ROOT.parent/'integrated20/model_bundle20.json'
    parent_checked=False
    if parent.exists():
        if sha(parent)!=PARENT_SHA or json.loads(parent.read_text())['models']['expanded13']!=m:raise ValueError('Parent/extract disagreement')
        parent_checked=True
    eligible,audit=cohort(source)
    common=np.isfinite(eligible[m['preprocessing']['features']]).all(axis=1).to_numpy()
    d=eligible.loc[common].copy().reset_index(drop=True);fold_all=assign_folds(eligible.cluster)
    folds=fold_all[common];w=d.weight.to_numpy(float);x=np.log2(d.SSSNFL.to_numpy(float))
    events=((d.dead==1)&(d.time_years<=4)).to_numpy(int)
    Hstop=baseline_hazard(d,m,np.minimum(d.time_years,4));H={h:baseline_hazard(d,m,h) for h in HORIZONS}
    preds={(name,h):np.zeros(len(d)) for name in ['frozen','calibration_oof','nfl_oof'] for h in HORIZONS}
    for h in HORIZONS:preds['frozen',h]=-np.expm1(-H[h])
    fold_records=[]
    for fold in range(5):
        train=folds!=fold;test=~train
        if set(d.loc[train,'cluster'])&set(d.loc[test,'cluster']):raise AssertionError('Cluster leakage')
        for name,use_nfl in [('calibration_oof',False),('nfl_oof',True)]:
            fit=fit_update(events[train],Hstop[train],w[train],x[train],use_nfl)
            fold_records.append(dict(fold=fold,model=name,train_n=int(train.sum()),test_n=int(test.sum()),
                train_events=int(events[train].sum()),test_events=int(events[test].sum()),**fit))
            for h in HORIZONS:preds[name,h][test]=predict_update(H[h][test],x[test],fit)
    masks={'all':np.ones(len(d),bool),'female':d.RIAGENDR.eq(2).to_numpy(),'male':d.RIAGENDR.eq(1).to_numpy(),
           'age40_59':d.RIDAGEYR.lt(60).to_numpy(),'age60_75':d.RIDAGEYR.ge(60).to_numpy(),
           'quantified_only':d.SSNFLL.eq(0).to_numpy()&d.SSNFLH.eq(0).to_numpy()}
    metrics=[]
    for h in HORIZONS:
        y=((d.dead==1)&(d.time_years<=h)).to_numpy(int)
        for name in ['frozen','calibration_oof','nfl_oof']:
            p=preds[name,h]
            for group,mask in masks.items():
                r=binary_metrics(y[mask],p[mask],w[mask]);r.update(model=name,horizon=h,subgroup=group)
                try:r.update(diagnostics(y[mask],p[mask],w[mask]));r['diagnostic_status']='computed'
                except ValueError:r['diagnostic_status']='nonconverged_not_refitted'
                metrics.append(r)
    rng=np.random.default_rng(SEED);draws=[]
    for rep in range(1000):
        wr=design_bootstrap_weights(eligible,rng)[common]
        rows={name:binary_metrics(events,preds[name,4],wr) for name in ['frozen','calibration_oof','nfl_oof']}
        record={}
        for name,r in rows.items():
            for metric in ['auc','brier','bias_pp']:
                record[name+'|'+metric]=r[metric]
        for left,right in [('nfl_oof','calibration_oof'),('nfl_oof','frozen')]:
            for metric in ['auc','brier','bias_pp']:
                a,b=rows[left][metric],rows[right][metric]
                record[left+'_minus_'+right+'|'+metric]=a-b if a is not None and b is not None else None
        draws.append(record)
    boot=pd.DataFrame(draws);intervals=[]
    for name in boot:
        values=boot[name].dropna();low,high=np.quantile(values,[.025,.975]) if len(values) else [None,None]
        intervals.append({'comparison_metric':name,'lower':low,'upper':high,'valid_replicates':len(values),
            'undefined_replicates':1000-len(values),'scope':'conditional_fixed_cross_fitted_predictions_excludes_all_fitting'})
    fullfit=fit_update(events,Hstop,w,x,True)
    fullfit.update(hr_per_doubling=float(np.exp(fullfit['log2_nfl_coefficient'])),n=len(d),events4=int(events.sum()),
                   role='full_sample_association_not_heldout_performance',clinical_use_ready=False)
    audit.update(common_complete_n=len(d),events4=int(events.sum()),events1=int(((d.dead==1)&(d.time_years<=1)).sum()),
        complete_below_LLOQ=int(d.SSNFLL.eq(1).sum()),complete_above_ULOQ=int(d.SSNFLH.eq(1).sum()),
        model_sha256=MODEL_SHA,parent_model_sha256=PARENT_SHA,parent_structural_equality_verified=parent_checked,
        analysis='cluster_cross_fitted_incremental_NfL_not_external_validation',fold_allocation='outcome_blind_PSU_sha256_round_robin5',
        unique_participant_ids=int(d.SEQN.nunique()),bootstrap_replicates=1000,bootstrap_seed=SEED,
        clinical_use_ready=False,raw_participant_records_exported=False,source_assay='Siemens AE chemiluminescence, SSSNFL_H',
        censoring_at4=False,full_sample_association=fullfit,
        quantified_sensitivity='restrict_same_OOF_predictions_no_refit_descriptive')
    missing=[{'feature':f,'missing_n':int((~np.isfinite(eligible[f])).sum())} for f in m['preprocessing']['features']]
    out.mkdir(parents=True)
    for name,rows in [('metrics',metrics),('intervals',intervals),('folds',fold_records),('missingness',missing)]:pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    dump(out/'RESULTS.json',audit);dump(out/'nfl_update_investigational.json',fullfit)
    (out/'SOURCE_MANIFEST.json').write_bytes((source/'manifest.json').read_bytes())
    print(json.dumps(audit,indent=2));print(pd.DataFrame(metrics).query('subgroup=="all"').to_string(index=False))
    return audit

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
