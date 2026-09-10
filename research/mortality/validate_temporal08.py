"""Locked participant-disjoint temporal validation. Aggregates only are exported.

See validation08/ANALYSIS_LOCK.md (Git commit f61b86c before input retrieval).
No model fitting/selection occurs on the new evaluation cycles. The sex-specific
calibration factors use the previously examined 2003-2004 subset only.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.optimize import minimize, brentq

from .nhanes_benchmark import Preprocessor, read_mortality, yes_no, evaluate
from .nhanes_cohort06 import make_cohort_checked
from .xpt_checked import read_xpt_checked

MODEL_SHA = '7b287f3f4778b2127f2226f79d2e503e1e39087306a172e5f524a561964b07d0'
MODELS = ('M0_age_sex', 'M0c_clinical_only_posthoc', 'M1_routine')
HORIZONS = (1, 5, 10)
COMPONENTS = ('DEMO','BPX','BMX','SMQ','MCQ','BPQ','DIQ','BIOPRO',
              'TCHOL','HDL','GHB','ALB_CR','CBC','CRP')
LOCK_COMMIT = 'f61b86c6fc00a43ae2490cf001b06c680271af63'
SEED = 20260910


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_new_sources(source: Path) -> dict:
    manifest = json.loads((source/'manifest.json').read_text())
    if manifest.get('analysis_lock_commit') != LOCK_COMMIT:
        raise ValueError('Source manifest does not identify the locked analysis')
    required = {f'{year}_{comp}_{suffix}.xpt' for year,suffix in [(2005,'D'),(2007,'E')]
                for comp in COMPONENTS} | {'2005_mortality.dat','2007_mortality.dat'}
    seen, checked, failures = set(), 0, []
    for row in manifest['files']:
        name = row['file']
        if Path(name).name != name or name in seen:
            raise ValueError('Invalid or duplicate source filename')
        seen.add(name)
        if row['status'] != 'retrieved':
            failures.append({'file':name,'error':row.get('error')})
            if name in required: raise ValueError('Missing required scientific data: '+name)
            continue
        path = source/name
        if not path.is_file() or path.stat().st_size != row['bytes'] or sha(path) != row['sha256']:
            raise ValueError('Changed source: '+name)
        checked += 1
    if required-seen: raise ValueError('Manifest omits required data')
    return {'files_verified':checked,'required_data_files':len(required),
            'documentation_retrieval_failures':failures,
            'manifest_sha256':sha(source/'manifest.json')}


def independent_mortality_reader(path: Path) -> pd.DataFrame:
    """Manual 1-based positions, independently parsed without pandas.read_fwf.

    Uses the public CDC file specification, not inference from outcome patterns.
    Preserves missingness and all declared cause categories. No cause recoding.
    """
    specification = [('SEQN',1,6),('elig',15,15),('dead',16,16),('cause',17,19),
                     ('death_diabetes',20,20),('death_hypertension',21,21),
                     ('t_int',43,45),('t_exam',46,48)]
    records=[]
    for line in path.read_text(encoding='ascii').splitlines():
        line=line.ljust(48)
        row={}
        for label,start,end in specification:
            token=line[start-1:end].strip()
            if token in {'','.'}: row[label]=np.nan
            elif token.isdigit(): row[label]=float(int(token))
            else: raise ValueError(f'Unexpected token at {start}:{end}')
        records.append(row)
    frame=pd.DataFrame(records)
    if frame.SEQN.isna().any() or frame.SEQN.duplicated().any():
        raise ValueError('Invalid linkage identifiers')
    return frame


def audit_mortality_reader(path: Path) -> dict:
    a=read_mortality(path).sort_values('SEQN').reset_index(drop=True)
    b=independent_mortality_reader(path).sort_values('SEQN').reset_index(drop=True)
    pd.testing.assert_frame_equal(a.astype(float),b.astype(float),check_dtype=False,check_exact=True)
    return {'file':path.name,'rows':len(a),'fields':len(a.columns),
            'cell_comparisons':int(a.size),'exact_match':True,'sha256':sha(path)}


def engineer_new(frame: pd.DataFrame, cycle: int) -> pd.DataFrame:
    d=frame.copy()
    raw=d.LBXSCR.where(d.LBXSCR>0)
    # CDC 2005-2006 equation. 2007-2008 file already applies its crossover correction.
    scr=(-.016+.978*raw) if cycle==2005 else raw
    scr=scr.where(scr>0)
    female=d.RIAGENDR.eq(2)
    r=scr/np.where(female,.7,.9)
    egfr=142*np.minimum(r,1)**np.where(female,-.241,-.302)*np.maximum(r,1)**-1.2*.9938**d.RIDAGEYR*np.where(female,1.012,1.)
    d['creatinine_standardized']=scr;d['egfr']=egfr
    d['age']=(d.RIDAGEYR-60)/10
    a=d.RIDAGEYR
    d['age_rcs']=(np.maximum(a-40,0)**3-(39/19)*np.maximum(a-60,0)**3+(20/19)*np.maximum(a-79,0)**3)/(39**2)
    d['male']=d.RIAGENDR.map({1.:1.,2.:0.});d['age_male']=d.age*d.male
    never=d.SMQ020.eq(2);smoker=d.SMQ020.eq(1)&d.SMQ040.isin([1,2,3])
    d['smoker_current']=np.where(never,0,np.where(smoker,d.SMQ040.isin([1,2]).astype(float),np.nan))
    d['smoker_former']=np.where(never,0,np.where(smoker,d.SMQ040.eq(3).astype(float),np.nan))
    d['bmi']=(d.BMXBMI-27)/5;d['bmi_sq']=d.bmi**2
    d['sbp']=(d[[f'BPXSY{i}' for i in [1,2,3]]].where(lambda x:x>0).mean(axis=1)-130)/20
    d['bp_treatment']=yes_no(d.BPQ050A)
    d.loc[d.BPQ020.eq(2)|d.BPQ040A.eq(2),'bp_treatment']=0.
    d['diabetes_history']=d.DIQ010.map({1.:1.,2.:0.,3.:0.})
    cvd=d[['MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']]
    d['cvd_history']=np.where(cvd.eq(1).any(axis=1),1.,np.where(cvd.eq(2).all(axis=1),0.,np.nan))
    d['cancer_history']=yes_no(d.MCQ220)
    d['total_cholesterol']=d.LBXTC;d['hdl']=d.LBDHDD
    d['hba1c']=d.LBXGH  # Withdrawn crossover equation is NOT applied.
    d['albumin']=d.LBXSAL*10
    d['egfr_low']=np.minimum(egfr,60)/15;d['egfr_high']=np.maximum(egfr-60,0)/30
    d['uacr']=100*d.URXUMA/d.URXUCR.where(d.URXUCR>0)
    d['log_uacr']=np.log2(d.uacr.where(d.uacr>0));d['rdw']=d.LBXRDW
    d['log_wbc']=np.log2(d.LBXWBCSI.where(d.LBXWBCSI>0))
    d['log_crp']=np.log2((d.LBXCRP*10).where(d.LBXCRP>0))
    return d


def new_cohort(source: Path) -> tuple[pd.DataFrame,list, list]:
    frames=[];flow=[];readers=[]
    for year,suffix in [(2005,'D'),(2007,'E')]:
        demo=read_xpt_checked(source/f'{year}_DEMO_{suffix}.xpt')
        survey_n=len(demo);ids=set(demo.SEQN)
        d=demo.loc[demo.RIDAGEYR.between(40,79),
                   ['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']].copy()
        age_n=len(d)
        for comp in COMPONENTS[1:]:
            data=read_xpt_checked(source/f'{year}_{comp}_{suffix}.xpt')
            if not set(data.SEQN).issubset(ids): raise ValueError('Foreign covariate participant IDs')
            columns=['SEQN']+[c for c in data if c not in d]
            d=d.merge(data[columns],on='SEQN',how='left',validate='1:1')
        readers.append(audit_mortality_reader(source/f'{year}_mortality.dat'))
        mort=independent_mortality_reader(source/f'{year}_mortality.dat')
        if set(mort.SEQN)!=ids:raise ValueError('Demographic/mortality identifiers differ')
        d=d.merge(mort,on='SEQN',how='left',validate='1:1')
        linked=d.elig.eq(1)&d.t_exam.notna()&d.t_exam.ge(0)&d.dead.isin([0,1])
        linked_n=int(linked.sum());d=d[linked].copy()
        valid=d.WTMEC2YR.gt(0)&np.isfinite(d.WTMEC2YR)&d.RIAGENDR.isin([1,2])
        d=d[valid].copy()
        if d[['SDMVSTRA','SDMVPSU']].isna().any().any():raise ValueError('Missing survey design field')
        d=engineer_new(d,year);d['cycle']=year;d['weight']=d.WTMEC2YR/2
        d['time_years']=d.t_exam.clip(lower=.5)/12
        d['stratum']=d.cycle.astype(str)+':'+d.SDMVSTRA.astype(str)
        if ((d.dead==0)&(d.time_years<10)).any():
            raise ValueError('Insufficient survivor follow-up for locked binary analysis')
        flow.append({'cycle':year,'survey_n':survey_n,'age40_79_n':age_n,
                     'eligible_linked_n':linked_n,'analysis_n':len(d),
                     'excluded_linkage_or_followup_n':age_n-linked_n,
                     'excluded_weight_or_sex_n':linked_n-len(d),
                     'zero_month_events':int(((d.t_exam==0)&(d.dead==1)).sum()),
                     'minimum_survivor_followup_years':float(d.loc[d.dead==0,'time_years'].min())})
        frames.append(d)
    result=pd.concat(frames,ignore_index=True)
    if result.SEQN.duplicated().any():raise ValueError('Repeated participant across new cycles')
    return result,flow,readers


def frozen_design(frame: pd.DataFrame, bundle: dict, name: str) -> np.ndarray:
    model=bundle['models'][name]
    pp=Preprocessor();pp.features=list(model['features']);pp.stats=bundle['preprocessing']
    return pp.transform(frame)


def cumulative_hazard(z: np.ndarray,model: dict,times) -> np.ndarray:
    z=np.asarray(z,float);t=np.broadcast_to(np.asarray(times,float),(len(z),))
    if not np.isfinite(z).all() or not np.isfinite(t).all() or np.any((t<0)|(t>10)):
        raise ValueError('Invalid cumulative hazard input')
    ends=np.asarray(model['interval_ends_years'],float);starts=np.r_[0,ends[:-1]]
    dt=np.maximum(0,np.minimum(t[:,None],ends)-starts)
    lp=z@np.asarray(model['coefficients'],float)
    return (dt@np.exp(model['baseline_log_hazards']))*np.exp(lp)


def fit_hazard_factor(events: np.ndarray,hazard: np.ndarray,weights: np.ndarray) -> float:
    y,h,w=map(lambda x:np.asarray(x,float),(events,hazard,weights))
    if not (y.ndim==h.ndim==w.ndim==1 and len(y)==len(h)==len(w) and len(y)>0):
        raise ValueError('Invalid calibration array shape')
    if not np.isfinite(np.r_[y,h,w]).all() or not np.isin(y,[0,1]).all() or np.any(h<0) or np.any(w<=0):
        raise ValueError('Invalid calibration input')
    numerator=float(w@y);denominator=float(w@h)
    if numerator<=0 or denominator<=0:raise ValueError('Insufficient calibration support')
    return numerator/denominator


def diagnostics(y,p,w) -> dict:
    y,p,w=map(lambda a:np.asarray(a,float),(y,p,w));w=w/w.mean()
    if np.unique(y).size<2:return {'calibration_slope':None,'calibration_in_the_large_logit':None,'diagnostic_gradient':None}
    x=logit(np.clip(p,1e-12,1-1e-12));design=np.column_stack([np.ones(len(x)),x])
    def fun(beta):
        lp=design@beta;err=w*(expit(lp)-y)
        return float(np.sum(w*(np.logaddexp(0,lp)-y*lp))),design.T@err
    fit=minimize(fun,[0.,1.],jac=True,method='BFGS',options={'gtol':1e-5,'maxiter':1000})
    grad=float(np.max(np.abs(fun(fit.x)[1])))
    if grad>1e-3:raise ValueError('Calibration diagnostic regression has not converged')
    cil=brentq(lambda a:float(w@(expit(x+a)-y)),-40,40)
    return {'calibration_intercept':float(fit.x[0]),'calibration_slope':float(fit.x[1]),
            'calibration_in_the_large_logit':float(cil),'diagnostic_gradient':grad}


def weighted_auc(y,p,w) -> float:
    # Independent rank-sum implementation, handles ties and zero resample weights.
    order=np.argsort(p,kind='mergesort');s=np.asarray(p)[order];y=np.asarray(y)[order];w=np.asarray(w)[order]
    if not np.isin(y,[0,1]).all() or np.any(w<0):raise ValueError('Invalid AUC input')
    starts=np.r_[0,np.flatnonzero(np.diff(s))+1]
    pos=np.add.reduceat(w*y,starts);neg=np.add.reduceat(w*(1-y),starts)
    if pos.sum()<=0 or neg.sum()<=0:raise ValueError('AUC requires both outcomes')
    return float(np.sum(pos*(np.cumsum(neg)-neg+.5*neg))/(pos.sum()*neg.sum()))


def design_bootstrap_weights(d: pd.DataFrame,rng) -> np.ndarray:
    mult=np.zeros(len(d));strata=d['stratum'].to_numpy();psu=d['SDMVPSU'].to_numpy()
    for st in np.unique(strata):
        indices=np.flatnonzero(strata==st);units=np.unique(psu[indices]);m=len(units)
        if m<2:raise ValueError('Singleton stratum requires a declared variance strategy')
        counts=Counter(rng.choice(units,size=m-1,replace=True))
        for unit,count in counts.items():mult[indices[psu[indices]==unit]]=count*m/(m-1)
    return d.weight.to_numpy(float)*mult


def write_json(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(old_source: Path,new_source: Path,model_path: Path,out: Path,replicates: int=1000) -> dict:
    if replicates!=1000: raise ValueError('Locked release requires 1000 paired PSU replicates')
    if out.exists():raise ValueError('Use a new output directory; do not overwrite a scientific run')
    if sha(model_path)!=MODEL_SHA:raise ValueError('Frozen model checksum mismatch')
    source_audit=verify_new_sources(new_source)
    bundle=json.loads(model_path.read_text())
    old,old_flow=make_cohort_checked(old_source)
    train=old[old.cycle<2003].copy();cal=old[old.cycle==2003].copy()
    # Old source data fingerprint set, not just a self-reported manifest.
    fingerprint_file=Path(__file__).parent/'nhanes04/SOURCE_DATA_SHA256.txt'
    old_hashes=0
    for line in fingerprint_file.read_text().splitlines():
        if not line.strip():continue
        expected,name=line.split(maxsplit=1);name=name.strip().lstrip('*')
        if sha(old_source/name)!=expected:raise ValueError('Historical data version mismatch: '+name)
        old_hashes+=1
    factors={};factor_rows=[]
    for name in MODELS:
        z=frozen_design(cal,bundle,name);m=bundle['models'][name]
        h=cumulative_hazard(z,m,np.minimum(cal.time_years,10))
        e=((cal.dead==1)&(cal.time_years<=10)).to_numpy(float)
        factors[name]={}
        for sex,code in [('female',2),('male',1)]:
            idx=cal.RIAGENDR.eq(code).to_numpy();w=cal.weight.to_numpy(float)[idx]
            factor=fit_hazard_factor(e[idx],h[idx],w);factors[name][sex]=factor
            factor_rows.append({'model':name,'sex':sex,'factor':factor,'n':int(idx.sum()),'events_10y':int(e[idx].sum())})
    # New outcomes are first processed below; all design and calibration choices above are fixed.
    d,flow,readers=new_cohort(new_source)
    if set(old.SEQN)&set(d.SEQN):raise ValueError('Evaluation/development participant overlap')
    features=bundle['models']['M1_routine']['features']
    complete=np.isfinite(d[features].to_numpy(float)).all(axis=1)
    records=[];predictions={};data_quality=[];bins=[]
    for f in features:
        st=bundle['preprocessing'][f];x=d[f].to_numpy(float);valid=np.isfinite(x)
        data_quality.append({'feature':f,'missing_n':int((~valid).sum()),'missing_fraction':float((~valid).mean()),
              'clipped_n':int((valid&((x<st['lower'])|(x>st['upper']))).sum()),
              'missing_indicator_trained':bool(st['missing_indicator'])})
    subgroups={'all':np.ones(len(d),bool),'female':d.RIAGENDR.eq(2).to_numpy(),
       'male':d.RIAGENDR.eq(1).to_numpy(),'age40_59':d.RIDAGEYR.lt(60).to_numpy(),
       'age60_79':d.RIDAGEYR.ge(60).to_numpy(),'cycle2005':d.cycle.eq(2005).to_numpy(),
       'cycle2007':d.cycle.eq(2007).to_numpy(),'complete_case':complete}
    for name in MODELS:
        model=bundle['models'][name];z=frozen_design(d,bundle,name)
        sf=np.where(d.RIAGENDR.eq(2),factors[name]['female'],factors[name]['male'])
        for horizon in HORIZONS:
            e=((d.dead==1)&(d.time_years<=horizon)).to_numpy(float)
            h=cumulative_hazard(z,model,horizon)
            for version,factor in [('frozen',1.),('sex_recalibrated',sf)]:
                p=-np.expm1(-h*factor);predictions[(name,version,horizon)]=p
                for group,mask in subgroups.items():
                    w=d.weight.to_numpy(float)[mask];result=evaluate(e[mask],p[mask],w)
                    result.update(diagnostics(e[mask],p[mask],w))
                    result.update(model=name,version=version,horizon_years=horizon,subgroup=group)
                    records.append(result)
                    if horizon==10 and group=='all':
                        if abs(weighted_auc(e,p,d.weight)-result['auc'])>1e-12:
                            raise ValueError('Independent AUC verification failed')
                if horizon==10:
                    cuts=np.quantile(p,np.linspace(0,1,11));assignment=np.digitize(p,cuts[1:-1],right=True)
                    for b in range(10):
                        mask=assignment==b
                        bins.append({'model':name,'version':version,'bin':b+1,**evaluate(e[mask],p[mask],d.weight.to_numpy()[mask])})
    print('PRIMARY_EVALUATION_COMPUTED',flush=True)
    summary_metrics=pd.DataFrame(records)
    print(summary_metrics.query("horizon_years==10 and subgroup=='all'")[['model','version','n','events','auc','brier','weighted_observed','weighted_predicted']].to_string(index=False),flush=True)
    rng=np.random.default_rng(SEED);draws=[]
    y=((d.dead==1)&(d.time_years<=10)).to_numpy(float)
    for rep in range(replicates):
        w=design_bootstrap_weights(d,rng);row={'replicate':rep}
        for version in ('frozen','sex_recalibrated'):
            for name in MODELS:
                p=predictions[(name,version,10)];prefix=f'{name}|{version}'
                row[prefix+'|auc']=weighted_auc(y,p,w)
                row[prefix+'|brier']=float(np.average((y-p)**2,weights=w))
            for metric in ('auc','brier'):
                row[f'delta_M1_M0c|{version}|{metric}']=row[f'M1_routine|{version}|{metric}']-row[f'M0c_clinical_only_posthoc|{version}|{metric}']
        draws.append(row)
    draws=pd.DataFrame(draws);intervals=[]
    for key in draws.columns[1:]:
        name,version,metric=key.split('|');lo,hi=np.quantile(draws[key],[.025,.975])
        intervals.append({'model_or_contrast':name,'version':version,'metric':metric,'horizon_years':10,
         'lower':float(lo),'upper':float(hi),'replicates':replicates,'interval_kind':'conditional_PSU_percentiles',
         'excludes_uncertainty':'model_training_and_calibration_factor_estimation'})
    results={'date':'2026-09-10','analysis_lock_commit':LOCK_COMMIT,'frozen_model_sha256':MODEL_SHA,
      'status':'participant_disjoint_temporal_validation','new_model_coefficients_fitted':False,
      'clinical_use_ready':False,'endpoint':'all_cause_mortality','country':'USA',
      'development_n':len(train),'calibration_n':len(cal),'evaluation_n':len(d),
      'evaluation_events_10y':int(y.sum()),'evaluation_female_n':int(d.RIAGENDR.eq(2).sum()),
      'complete_case_n':int(complete.sum()),'complete_case_events_10y':int(y[complete].sum()),
      'participant_overlap':0,'historical_data_hashes_verified':old_hashes,
      'new_source_audit':source_audit,'flow':flow,'old_flow':old_flow,'reader_audit':readers,
      'calibration_factors':factors,'bootstrap_replicates':replicates,'bootstrap_seed':SEED,
      'evaluation_strata':int(d.stratum.nunique()),
      'raw_participant_records_exported':False,'locked_models':list(MODELS)}
    out.mkdir(parents=True)
    summary_metrics.to_csv(out/'metrics.csv',index=False)
    pd.DataFrame(intervals).to_csv(out/'intervals.csv',index=False)
    pd.DataFrame(data_quality).to_csv(out/'feature_quality.csv',index=False)
    pd.DataFrame(factor_rows).to_csv(out/'calibration_factors.csv',index=False)
    pd.DataFrame(bins).to_csv(out/'calibration_bins.csv',index=False)
    write_json(out/'results.json',results)
    write_json(out/'recalibration_spec.json',{'version':'0.8_research_temporal_validation',
      'parent_model_sha256':MODEL_SHA,'factors':factors,'factor_estimation_cohort':'NHANES2003_2004_cardio_subset',
      'evaluation_cohorts':['NHANES2005_2006','NHANES2007_2008'],'endpoint':'all_cause_mortality',
      'age_range':[40,79],'horizons_years':list(HORIZONS),'clinical_use_ready':False,
      'not_automatically_activated_in_workbench':True})
    (out/'SOURCE_MANIFEST.json').write_bytes((new_source/'manifest.json').read_bytes())
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old-source',type=Path,required=True);p.add_argument('--new-source',type=Path,required=True)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();print(json.dumps(main(args.old_source,args.new_source,args.model,args.out),indent=2))
