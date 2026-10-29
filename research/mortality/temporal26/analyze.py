"""One-year assessment on a new-to-checked-history cycle; frozen models only."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import sklearn
from ..nfl22 import analyze as base
from ..panel25.analyze import join_component
from ..ipcw14 import outcome_weights, binary_scores, aalen_johansen_reference
from .retrieve import source_urls

ROOT=Path(__file__).parent
MORTALITY=ROOT.parent
SEED=20261030


def load_models():
    lock=json.loads((ROOT/'LOCK.json').read_text())
    if base.sha(ROOT/'ANALYSIS_PLAN.md')!=lock['plan_sha256']:
        raise ValueError('Analysis plan changed')
    for name,digest in lock['models'].items():
        if base.sha(MORTALITY/name)!=digest: raise ValueError('Frozen model changed')
    expanded=base.load_model()
    compact=json.loads((MORTALITY/'compact14/model_compact14.json').read_text())['model']
    return {'expanded13':expanded,'compact4':compact},lock


def verify_sources(source):
    source=Path(source); doc=json.loads((source/'manifest.json').read_text())
    if doc.get('analysis_lock_commit')!='5fbb68c' or doc.get('lock_sha256')!=base.sha(ROOT/'LOCK.json'):
        raise ValueError('Sources not tied to locked protocol')
    expected=source_urls(); seen=set(); stable={}
    for row in doc['files']:
        name=row['file']
        if name in seen or name not in expected or Path(name).name!=name: raise ValueError('Unexpected/duplicate source')
        seen.add(name)
        if row.get('url')!=expected[name] or row.get('status')!='retrieved': raise ValueError('Wrong or missing official source')
        path=source/name
        if not path.is_file() or path.stat().st_size!=row['bytes'] or base.sha(path)!=row['sha256']:
            raise ValueError('Source changed: '+name)
        stable[name]=row['sha256']
    if seen!=set(expected): raise ValueError('Missing required source')
    stable_digest=hashlib.sha256(json.dumps(stable,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return doc,stable,stable_digest


def bridge(frame):
    """CDC published backward equations, Cobas6000 to DxC660i units."""
    d=frame.copy()
    for col in ('LBXSAL','LBXSAPSI','LBXSCR','LBXSGTSI'):
        d[col]=d[col].where(d[col]>0)
    d['LBXSAL']=1.044*d.LBXSAL+.01128
    d['LBXSAPSI']=10**(1.001*np.log10(d.LBXSAPSI)-.04294)
    d['LBXSCR']=1.051*d.LBXSCR-.06945
    d['LBXSGTSI']=.8042*d.LBXSGTSI+2.363
    for col in ('LBXSAL','LBXSAPSI','LBXSCR','LBXSGTSI'):
        d[col]=d[col].where((d[col]>0)&np.isfinite(d[col]))
    return d


def cohort(source):
    source=Path(source); doc,stable,digest=verify_sources(source)
    demo=base.read_xpt_checked(source/'2017_DEMO_J.xpt')
    if demo.SEQN.isna().any() or demo.SEQN.duplicated().any(): raise ValueError('Invalid demographic IDs')
    # Official2017–2018 SEQN range; every older survey ID is outside this range.
    if not demo.SEQN.between(93703,102956).all(): raise ValueError('Unexpected cycle ID support')
    ids=set(demo.SEQN)
    d=demo.loc[demo.RIDAGEYR.between(40,75),['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']].copy()
    age_n=len(d)
    for c in base.COMPONENTS[1:]:
        d=join_component(d,base.read_xpt_checked(source/f'2017_{c}_J.xpt'),ids)
    mort=base.independent_mortality_reader(source/'2017_mortality.dat')
    if set(mort.SEQN)!=ids: raise ValueError('Mortality/demo ID mismatch')
    d=join_component(d,mort,ids)
    valid=(d.elig.eq(1)&d.dead.isin([0,1])&d.t_exam.ge(0)&np.isfinite(d.t_exam)
           &d.RIAGENDR.isin([1,2])&d.WTMEC2YR.gt(0)&np.isfinite(d.WTMEC2YR))
    d=d.loc[valid].copy().reset_index(drop=True)
    if d.empty or d[['SDMVSTRA','SDMVPSU']].isna().any().any(): raise ValueError('Empty/missing design')
    d['weight']=d.WTMEC2YR;d['stratum']='2017:'+d.SDMVSTRA.astype(int).astype(str)
    d['time_years']=d.t_exam.clip(lower=.5)/12
    audit=dict(survey_n=len(demo),age40_75_n=age_n,eligible_n=len(d),
               reader=base.audit_mortality_reader(source/'2017_mortality.dat'),
               source_files=len(stable),source_fingerprints=stable,source_content_sha256=digest,
               source_manifest_sha256=base.sha(source/'manifest.json'),
               zero_month_deaths=int((d.dead.eq(1)&d.t_exam.eq(0)).sum()),
               minimum_survivor_followup_years=float(d.loc[d.dead.eq(0),'time_years'].min()),
               unknown_cause_deaths=int((d.dead.eq(1)&d.cause.isna()).sum()))
    return d,audit,doc


def total_risk(frame,model,horizon=1):
    ends=np.asarray(model['interval_ends_years'],float)
    if horizon<=0 or horizon>ends[-1]: raise ValueError('Unsupported horizon')
    starts=np.r_[0,ends[:-1]];dt=np.maximum(0,np.minimum(horizon,ends)-starts)
    z=base.transform(frame,model['preprocessing']);hazard=np.zeros(len(frame))
    for cause in model['causes'].values():
        rates=np.asarray(cause['baseline_rates'],float) if 'baseline_rates' in cause else np.exp(cause['baseline_log_hazards'])
        hazard+=np.exp(z@np.asarray(cause['coefficients']))*(rates@dt)
    p=-np.expm1(-hazard)
    if not np.isfinite(p).all() or np.any((p<0)|(p>1)): raise ValueError('Invalid frozen risk')
    return p


def assess(d,predictions,weights):
    w=np.asarray(weights,float)
    o=outcome_weights(d.dead,d.cause,d.time_years,w,1)
    target=o['target']
    reference=aalen_johansen_reference(d.dead,d.cause,d.time_years,w,1)
    empirical=np.bincount(target,weights=o['effective_weights'],minlength=4)/w.sum()
    if np.max(np.abs(reference-empirical))>1e-10: raise ArithmeticError('IPCW/Aalen-Johansen disagreement')
    rows={name:binary_scores(target,p,w,o) for name,p in predictions.items()}
    a,b=rows['expanded13_bridge'],rows['compact4_bridge']
    return dict(metrics=rows,delta_brier=a['brier']-b['brier'],
                delta_auc=None if a['auc'] is None or b['auc'] is None else a['auc']-b['auc'],
                censoring_G=o['G_horizon_left'],early_censored_n=o['early_censored_n'],
                AJ_check_max_abs=float(np.max(np.abs(reference-empirical))))


def summarize(draws,total):
    keys=sorted({'delta_brier','delta_auc'} | {k for row in draws for k in row if k not in ('replicate','status','error')})
    result={}
    for key in keys:
        vals=[r[key] for r in draws if r.get(key) is not None and np.isfinite(r[key])]
        lo,hi=map(float,np.quantile(vals,[.025,.975])) if vals else (None,None)
        result[key]=dict(lower=lo,upper=hi,valid=len(vals),undefined=total-len(vals),
                         scope='fixed_models_bridge_censoring_refitted',defined_draws_only=len(vals)!=total)
    return result


def file_provenance():
    # Hash the complete scientific Python subtree to cover imported readers and helpers.
    return {str(p.relative_to(MORTALITY)):base.sha(p) for p in sorted(MORTALITY.rglob('*.py'))}


def run(source,out):
    source,out=Path(source),Path(out)
    if out.exists(): raise ValueError('New output directory required')
    models,lock=load_models();eligible,audit,manifest=cohort(source)
    corrected=base.engineer_expanded(bridge(eligible));raw=base.engineer_expanded(eligible)
    good={name:np.isfinite(corrected[m['preprocessing']['features']]).all(axis=1).to_numpy() for name,m in models.items()}
    common=good['expanded13']&good['compact4']
    # Raw sensitivity must keep exactly these participants, with no silent drop.
    features=sorted(set().union(*(m['preprocessing']['features'] for m in models.values())))
    if not np.isfinite(raw.loc[common,features]).all().all(): raise ValueError('Raw sensitivity needs same common participants')
    d=corrected.loc[common].copy();r=raw.loc[common].copy();w=d.weight.to_numpy(float)
    if d.empty: raise ValueError('Empty common domain')
    predictions={name+'_'+variant:total_risk(frame,m) for name,m in models.items() for variant,frame in [('bridge',d),('raw',r)]}
    point=assess(d,predictions,w)
    draws=[];rng=np.random.default_rng(SEED)
    for rep in range(1000):
        try:
            wr=base.design_bootstrap_weights(eligible,rng)[common]
            score=assess(d,predictions,wr)
            row=dict(replicate=rep,status='computed',delta_brier=score['delta_brier'],delta_auc=score['delta_auc'])
            for name,m in score['metrics'].items():
                for k in ('brier','auc','observed','predicted','bias_pp'): row[name+'|'+k]=m[k]
        except (ValueError,ArithmeticError,RuntimeError) as exc:
            row=dict(replicate=rep,status='undefined',error=type(exc).__name__+': '+str(exc))
        draws.append(row)
    quality={name:{f:int(((d[f]<m['preprocessing']['lower'][j])|(d[f]>m['preprocessing']['upper'][j])).sum())
                   for j,f in enumerate(m['preprocessing']['features'])} for name,m in models.items()}
    result=dict(analysis='locked_new_cycle_temporal_frozen_panels_one_year',cycle='2017-2018',horizon_years=1,
        n=len(d),events1=int((d.dead.eq(1)&d.time_years.le(1)).sum()),
        coverage={name:dict(complete_n=int(mask.sum()),weighted_fraction=float(np.average(mask,weights=eligible.weight))) for name,mask in good.items()},
        common_domain_n=int(common.sum()),excluded_n=int((~common).sum()),
        effective_n=float(w.sum()**2/(w@w)),strata_n=int(d.stratum.nunique()),
        psu_n=int(d[['stratum','SDMVPSU']].drop_duplicates().shape[0]),
        missingness={f:int((~np.isfinite(corrected[f])).sum()) for f in features},clipped_n=quality,
        point=point,intervals=summarize(draws,1000),bootstrap_replicates=1000,seed=SEED,
        bootstrap_failures=sum(r['status']=='undefined' for r in draws),source_audit=audit,
        models=lock['models'],plan_sha256=lock['plan_sha256'],lock_sha256=base.sha(ROOT/'LOCK.json'),
        analysis_plan_commit='5fbb68c',scientific_python_sha256=file_provenance(),
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__),
        model_fitting_uncertainty_included=False,bridge_estimation_uncertainty_included=False,
        clinical_use_ready=False,participant_records_exported=False)
    out.mkdir(parents=True)
    base.dump(out/'RESULTS.json',result);base.dump(out/'SOURCE_MANIFEST.json',manifest)
    pd.DataFrame(draws).to_csv(out/'bootstrap_aggregates.csv',index=False)
    print('TEMPORAL26_RESULTS_BEGIN');print(json.dumps(result,indent=2,allow_nan=False));print('TEMPORAL26_RESULTS_END')
    print('TEMPORAL26_MANIFEST_BEGIN');print(json.dumps(manifest,indent=2));print('TEMPORAL26_MANIFEST_END')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
