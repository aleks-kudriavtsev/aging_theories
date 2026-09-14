"""Five-cause historical-US model, locked before independent 2011-2014 input retrieval.

Train and freeze first; evaluate later without refitting. Only aggregates/models
are exported. Complete-profile domain is identical for both panel comparisons.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from sklearn.metrics import roc_auc_score
from .audit_endpoint_cohorts11 import cohort, FEATURES as OLD_FEATURES
from .competing_risks09 import preprocess_fit, transform, integrate_hazards, multiclass_brier, binary_metrics
from .validate_temporal08 import (independent_mortality_reader, audit_mortality_reader,
    design_bootstrap_weights, diagnostics, yes_no)
from .xpt_checked import read_xpt_checked

LOCK = '55279f23e10b9318ceddfa266d8d404e717006ee'
GROUPS = ('heart_diseases','malignant_neoplasms','chronic_lower_respiratory',
          'cerebrovascular','other_or_unknown')
ENDS = np.array([5.,8.])
CLINICAL = list(OLD_FEATURES[:13])
ROUTINE = [f for f in OLD_FEATURES if f != 'log_crp']
PANELS = {'clinical':CLINICAL,'routine_no_crp':ROUTINE}
FOLDS = ((1999,2001),(2003,2005),(2007,2009))
PENALTIES = (1.,10.,100.)
SEED = 20260912
COMPONENTS = ('DEMO','BPX','BMX','SMQ','MCQ','BPQ','DIQ','BIOPRO','TCHOL','HDL','GHB','ALB_CR','CBC')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def classify(dead, cause):
    d,c=np.asarray(dead,float),np.asarray(cause,float)
    if d.ndim!=1 or c.shape!=d.shape or not np.isin(d,[0,1]).all() or np.isinf(c).any():
        raise ValueError('Invalid vital-status/cause shape or values')
    known=np.isfinite(c)
    if not np.isin(c[known],range(1,11)).all() or np.any((d==0)&known):
        raise ValueError('Unsupported cause code or cause assigned to survivor')
    result=np.full(len(d),4,int)
    for code,label in [(1,0),(2,1),(3,2),(5,3)]:result[c==code]=label
    result[d==0]=-1
    return result


def exposure(time, group):
    t,g=np.asarray(time,float),np.asarray(group)
    if t.ndim!=1 or g.shape!=t.shape or not len(t) or not np.isfinite(t).all() or (t<=0).any() or not np.isin(g,range(-1,5)).all():
        raise ValueError('Invalid observation time or terminal event')
    dt=np.maximum(0,np.minimum(t[:,None],ENDS)-np.r_[0.,ENDS[:-1]])
    ev=np.zeros((len(t),2,5));idx=np.flatnonzero((g>=0)&(t<=8))
    ev[idx,np.searchsorted(ENDS,t[idx],side='left'),g[idx]]=1.
    return dt,ev


def fit_cause(z, dt, events, weights, penalty):
    z,dt,y,w=[np.asarray(v,float) for v in (z,dt,events,weights)]
    if (z.ndim!=2 or dt.shape!=(len(z),2) or y.shape!=dt.shape or w.shape!=(len(z),)
        or not len(z) or not all(np.isfinite(v).all() for v in (z,dt,y,w))
        or np.any(w<=0) or np.any(dt<0) or not np.isin(y,[0,1]).all()
        or np.any(y.sum(axis=1)>1) or np.any((y>0)&(dt<=0))
        or type(penalty) not in (int,float) or not np.isfinite(penalty) or penalty<0):
        raise ValueError('Invalid cause likelihood inputs')
    w=w/w.mean();counts=(w[:,None]*y).sum(axis=0)
    if np.any(counts<=0):raise ValueError('Empty event cell: no pseudo-events permitted')
    logmass=np.full_like(dt,-np.inf);pos=dt>0;logmass[pos]=np.log((w[:,None]*dt)[pos])
    ex=z.T@(w*y.sum(axis=1))
    def objective(beta):
        terms=logmass+(z@beta)[:,None];totals=logsumexp(terms,axis=0)
        grad=z.T@(np.exp(terms-totals)@counts)-ex+penalty*beta
        return float(counts@totals-ex@beta+.5*penalty*(beta@beta)),grad
    fit=minimize(objective,np.zeros(z.shape[1]),jac=True,method='BFGS',options={'gtol':1e-6,'maxiter':2000})
    val,grad=objective(fit.x);gmax=float(np.max(np.abs(grad)))
    if not np.isfinite(val) or gmax>1e-4:raise ValueError(f'Optimization failed: gradient {gmax}')
    baselines=np.log(counts)-logsumexp(logmass+(z@fit.x)[:,None],axis=0)
    return {'coefficients':fit.x.tolist(),'baseline_log_hazards':baselines.tolist(),
        'raw_events':y.sum(axis=0).astype(int).tolist(),'gradient_max':gmax,'iterations':int(fit.nit)}


def fit_panel(d, features, penalty):
    pp=preprocess_fit(d,features);z=transform(d,pp)
    dt,ev=exposure(d.time_years,classify(d.dead,d.cause))
    return {'preprocessing':pp,'causes':{g:fit_cause(z,dt,ev[:,:,k],d.weight,penalty) for k,g in enumerate(GROUPS)},
        'interval_ends_years':ENDS.tolist(),'penalty':penalty,'training_n':len(d)}


def predict(d, model, horizon):
    z=transform(d,model['preprocessing'])
    rates=np.stack([np.exp((z@np.asarray(model['causes'][g]['coefficients']))[:,None]+
            np.asarray(model['causes'][g]['baseline_log_hazards'])[None,:]) for g in GROUPS],axis=2)
    return integrate_hazards(rates,horizon,model['interval_ends_years'])


def labels(d, horizon):
    if not np.isfinite(horizon) or not 0<horizon<=8:raise ValueError('Unsupported horizon')
    if ((d.dead==0)&(d.time_years<horizon)).any():raise ValueError('Early censoring; fixed-horizon outcome not identified')
    g=classify(d.dead,d.cause)
    return np.where((g>=0)&(d.time_years<=horizon),g+1,0)


def train(old_source, source08, source09, out):
    if out.exists():raise ValueError('New output directory required')
    roots={y:old_source if y<2005 else source08 if y<2009 else source09 for y in (1999,2001,2003,2005,2007,2009)}
    d,flow,used,reader=cohort(roots)
    expected=json.loads((Path(__file__).parent/'transport12/expected_development_sources.json').read_text())
    for row in expected:
        if used.get(row['file'],{}).get('sha256')!=row['sha256']:raise ValueError('Historical data version changed')
    included=np.isfinite(d[ROUTINE]).all(axis=1);d=d[included].copy().reset_index(drop=True)
    d['time_years']=d.t_exam.clip(lower=.5)/12;d['weight']=d.pooled_weight
    models={};cv=[]
    for panel,features in PANELS.items():
        scores={}
        for penalty in PENALTIES:
            fold_scores=[]
            for heldout in FOLDS:
                idx=d.cycle.isin(heldout);m=fit_panel(d[~idx],features,penalty)
                score=multiclass_brier(labels(d[idx],5),predict(d[idx],m,5),d.loc[idx,'weight'])
                cv.append({'panel':panel,'penalty':penalty,'heldout_cycles':str(heldout),'n':int(idx.sum()),'six_state_brier':score})
                fold_scores.append(score)
            scores[penalty]=float(np.mean(fold_scores));print('CV',panel,penalty,scores[penalty],flush=True)
        selected=min(PENALTIES,key=lambda p:(scores[p],-p))
        models[panel]=fit_panel(d,features,selected)
        print('FROZEN',panel,selected,len(d),flush=True)
    bundle={'version':'0.12_five_cause_no_crp','analysis_lock_commit':LOCK,'base_commit':'85f9df19ea87bbe1cffd94d008af2a9acf84167d',
      'groups':list(GROUPS),'models':models,'development_cycles':[1999,2001,2003,2005,2007,2009],
      'development_n':len(d),'age_range':[40,79],'country':'historical_USA','clinical_use_ready':False,
      'noninfectious_risk_identified':False,'imputation':False,'CRP_used':False,
      'evaluated_horizons_planned':{'2011':[1,5],'2013':[1,4]},'fitted_at_utc':datetime.now(timezone.utc).isoformat(),
      'fit_code_sha256':sha(Path(__file__)),'evaluation_data_seen':False}
    counts=[]
    for horizon in (5,8):
        target=labels(d,horizon)
        counts += [{'horizon':horizon,'cause':g,'events':int((target==k+1).sum())} for k,g in enumerate(GROUPS)]
    out.mkdir(parents=True);dump(out/'model_bundle12.json',bundle)
    pd.DataFrame(cv).to_csv(out/'training_cv.csv',index=False)
    pd.DataFrame(counts).to_csv(out/'training_events.csv',index=False)
    dump(out/'development_audit.json',{'n_eligible':int(len(included)),'n_no_crp_complete':len(d),
      'data_files_verified':len(used),'flow_legacy_crp_counts':flow,'reader_audit':reader,
      'participant_overlap_between_cycles':0,'model_sha256':sha(out/'model_bundle12.json')})
    dump(out/'development_sources.json',list(used.values()))
    # Local-only cache for disjointness verification. Never packaged or committed.
    np.save(out/'LOCAL_ONLY_DEVELOPMENT_IDS.npy',d.SEQN.to_numpy())
    print('MODEL_SHA256',sha(out/'model_bundle12.json'),flush=True)
    return bundle


def engineer_no_crp(d):
    d=d.copy();female=d.RIAGENDR.eq(2);r=d.LBXSCR.where(d.LBXSCR>0)/np.where(female,.7,.9)
    egfr=142*np.minimum(r,1)**np.where(female,-.241,-.302)*np.maximum(r,1)**-1.2*.9938**d.RIDAGEYR*np.where(female,1.012,1.)
    a=d.RIDAGEYR;d['age']=(a-60)/10;d['age_rcs']=(np.maximum(a-40,0)**3-(39/19)*np.maximum(a-60,0)**3+(20/19)*np.maximum(a-79,0)**3)/(39**2)
    d['male']=d.RIAGENDR.map({1.:1.,2.:0.});d['age_male']=d.age*d.male
    never=d.SMQ020.eq(2);smoker=d.SMQ020.eq(1)&d.SMQ040.isin([1,2,3])
    d['smoker_current']=np.where(never,0,np.where(smoker,d.SMQ040.isin([1,2]).astype(float),np.nan))
    d['smoker_former']=np.where(never,0,np.where(smoker,d.SMQ040.eq(3).astype(float),np.nan))
    d['bmi']=(d.BMXBMI-27)/5;d['bmi_sq']=d.bmi**2
    d['sbp']=(d[[f'BPXSY{i}' for i in (1,2,3)]].where(lambda x:x>0).mean(axis=1)-130)/20
    d['bp_treatment']=yes_no(d.BPQ050A);d.loc[d.BPQ020.eq(2)|d.BPQ040A.eq(2),'bp_treatment']=0.
    d['diabetes_history']=d.DIQ010.map({1.:1.,2.:0.,3.:0.})
    cvd=d[['MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']]
    d['cvd_history']=np.where(cvd.eq(1).any(axis=1),1.,np.where(cvd.eq(2).all(axis=1),0.,np.nan))
    d['cancer_history']=yes_no(d.MCQ220)
    d['total_cholesterol']=d.LBXTC;d['hdl']=d.LBDHDD;d['hba1c']=d.LBXGH
    d['albumin']=d.LBXSAL*10;d['egfr_low']=np.minimum(egfr,60)/15;d['egfr_high']=np.maximum(egfr-60,0)/30
    uacr=100*d.URXUMA/d.URXUCR.where(d.URXUCR>0)
    d['log_uacr']=np.log2(uacr.where(uacr>0));d['rdw']=d.LBXRDW
    d['log_wbc']=np.log2(d.LBXWBCSI.where(d.LBXWBCSI>0))
    return d


def load_evaluation(source, year, expected_model_sha):
    if year not in (2011,2013):raise ValueError('Unplanned evaluation cycle')
    doc=json.loads((source/'manifest.json').read_text())
    if doc.get('analysis_lock_commit')!=LOCK or doc.get('model_sha256_before_retrieval')!=expected_model_sha:
        raise ValueError('Retrieval is not tied to frozen model and protocol')
    manifests={r['file']:r for r in doc['files']};used=[]
    def read(name):
        row=manifests.get(name,{})
        if row.get('status')!='retrieved' or sha(source/name)!=row.get('sha256'):raise ValueError('Missing/changed evaluation input: '+name)
        used.append(name)
        return read_xpt_checked(source/name) if name.endswith('.xpt') else independent_mortality_reader(source/name)
    suffix='G' if year==2011 else 'H'
    demo=read(f'{year}_DEMO_{suffix}.xpt');ids=set(demo.SEQN)
    if demo.SEQN.duplicated().any():raise ValueError('Duplicate demographic ID')
    cols=['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']
    d=demo.loc[demo.RIDAGEYR.between(40,79),cols].copy();n_age=len(d)
    for comp in COMPONENTS[1:]:
        raw=read(f'{year}_{comp}_{suffix}.xpt')
        if not set(raw.SEQN)<=ids:raise ValueError('Out-of-cycle predictor IDs')
        d=d.merge(raw[['SEQN']+[c for c in raw if c not in d]],on='SEQN',how='left',validate='1:1')
    mort=read(f'{year}_mortality.dat')
    if set(mort.SEQN)!=ids:raise ValueError('Mortality IDs differ')
    d=d.merge(mort,on='SEQN',how='left',validate='1:1')
    idx=d.elig.eq(1)&d.dead.isin([0,1])&d.t_exam.notna()&d.t_exam.ge(0)&d.WTMEC2YR.gt(0)&np.isfinite(d.WTMEC2YR)&d.RIAGENDR.isin([1,2])
    d=engineer_no_crp(d[idx].copy());eligible=len(d);complete=np.isfinite(d[ROUTINE]).all(axis=1)
    missing={f:int((~np.isfinite(d[f])).sum()) for f in ROUTINE}
    d=d[complete].copy().reset_index(drop=True);d['cycle']=year;d['time_years']=d.t_exam.clip(lower=.5)/12
    d['weight']=d.WTMEC2YR;d['stratum']=str(year)+':'+d.SDMVSTRA.astype(str)
    if not len(d) or d[['SDMVSTRA','SDMVPSU']].isna().any().any():raise ValueError('Empty sample/missing design')
    return d,{'year':year,'n_survey':len(demo),'n_age40_79':n_age,'n_eligible':eligible,'n_complete':len(d),
      'missing_feature_counts':missing,'minimum_survivor_followup':float(d.loc[d.dead==0,'time_years'].min()),
      'source_data_files_verified':len(used),'zero_month_deaths':int(((d.dead==1)&(d.t_exam==0)).sum()),
      'reader_audit':audit_mortality_reader(source/f'{year}_mortality.dat')}


def evaluate_new(source, trained, out, expected_model_sha):
    if out.exists():raise ValueError('New evaluation output directory required')
    modelpath=trained/'model_bundle12.json'
    if sha(modelpath)!=expected_model_sha:raise ValueError('Frozen model changed')
    bundle=json.loads(modelpath.read_text());old_ids=set(np.load(trained/'LOCAL_ONLY_DEVELOPMENT_IDS.npy'))
    rows=[];intervals=[];flow=[];seen=set();quality=[];bins=[];mass=0.;counts=[]
    for year,primary_h in [(2011,5),(2013,4)]:
        d,audit=load_evaluation(source,year,expected_model_sha);ids=set(d.SEQN)
        if ids&old_ids or ids&seen:raise ValueError('Development/evaluation participant overlap')
        seen|=ids;flow.append(audit);w=d.weight.to_numpy();preds={}
        masks={'all':np.ones(len(d),bool),'female':d.RIAGENDR.eq(2).to_numpy(),'male':d.RIAGENDR.eq(1).to_numpy(),
               'age40_59':d.RIDAGEYR.lt(60).to_numpy(),'age60_79':d.RIDAGEYR.ge(60).to_numpy()}
        for sex in ('female','male'):
            for age in ('age40_59','age60_79'):masks[sex+'_'+age]=masks[sex]&masks[age]
        for horizon in (1,primary_h):
            target=labels(d,horizon)
            counts += [{'year':year,'horizon':horizon,'cause':g,'events':int((target==k+1).sum())} for k,g in enumerate(GROUPS)]
            for panel,model in bundle['models'].items():
                p=predict(d,model,horizon);preds[panel,horizon]=p;mass=max(mass,float(abs(p.sum(axis=1)-1).max()))
                for subgroup,mask in masks.items():
                    if not mask.any():continue
                    rows.append({'year':year,'horizon':horizon,'panel':panel,'subgroup':subgroup,'outcome':'six_state',
                        'n':int(mask.sum()),'events':int((target[mask]>0).sum()),'brier':multiclass_brier(target[mask],p[mask],w[mask])})
                    for k,name in list(enumerate(GROUPS,1))+[(0,'all_cause')]:
                        y=target==k if k else target>0;risk=p[:,k] if k else 1-p[:,0]
                        met=binary_metrics(y[mask],risk[mask],w[mask])
                        row={'year':year,'horizon':horizon,'panel':panel,'subgroup':subgroup,'outcome':name,**met}
                        if subgroup=='all' and horizon==primary_h:
                            try:row.update(diagnostics(y,risk,w));row['calibration_diagnostic_status']='computed'
                            except ValueError:row['calibration_diagnostic_status']='not_converged_no_recalibration'
                            if met['auc'] is not None and abs(met['auc']-roc_auc_score(y,risk,sample_weight=w))>1e-12:raise ArithmeticError('AUC mismatch')
                            cut=np.quantile(risk,np.linspace(0,1,9));assign=np.digitize(risk,cut[1:-1],right=True)
                            for b in range(8):
                                s=assign==b
                                if s.any():bins.append({'year':year,'horizon':horizon,'panel':panel,'outcome':name,'bin':b+1,**binary_metrics(y[s],risk[s],w[s])})
                        rows.append(row)
        for panel,model in bundle['models'].items():
            if (preds[panel,primary_h][:,1:]<preds[panel,1][:,1:]-1e-12).any():raise ArithmeticError('Nonmonotone CIF')
            pp=model['preprocessing']
            for j,f in enumerate(pp['features']):quality.append({'year':year,'panel':panel,'feature':f,'clipped_n':int(((d[f]<pp['lower'][j])|(d[f]>pp['upper'][j])).sum())})
        target=labels(d,primary_h);draws=[];rng=np.random.default_rng(SEED+year-2011)
        for rep in range(1000):
            wr=design_bootstrap_weights(d,rng);r={}
            for panel in PANELS:
                p=preds[panel,primary_h]
                r[panel+'|six_state|brier']=multiclass_brier(target,p,wr)
                for k,name in list(enumerate(GROUPS,1))+[(0,'all_cause')]:
                    y=target==k if k else target>0;risk=p[:,k] if k else 1-p[:,0]
                    m=binary_metrics(y,risk,wr)
                    for metric in ('brier','auc','bias_pp'):r[panel+'|'+name+'|'+metric]=m[metric]
            for outcome in ['six_state',*GROUPS,'all_cause']:
                for metric in (['brier'] if outcome=='six_state' else ['brier','auc','bias_pp']):
                    a,b=r['routine_no_crp|'+outcome+'|'+metric],r['clinical|'+outcome+'|'+metric]
                    r['delta_routine_clinical|'+outcome+'|'+metric]=a-b if a is not None and b is not None else None
            draws.append(r)
        bootstrap=pd.DataFrame(draws)
        for key in bootstrap:
            panel,outcome,metric=key.split('|');x=bootstrap[key].dropna();lo,hi=np.quantile(x,[.025,.975]) if len(x) else (None,None)
            intervals.append({'year':year,'horizon':primary_h,'panel_or_contrast':panel,'outcome':outcome,'metric':metric,
                'lower':lo,'upper':hi,'valid_replicates':len(x),'undefined_replicates':1000-len(x),
                'interval_kind':'conditional_1000_paired_rescaled_PSU_percentiles'})
        print('EVALUATION_COMPLETE',year,len(d),int((target>0).sum()),flush=True)
    out.mkdir(parents=True);pd.DataFrame(rows).to_csv(out/'metrics.csv',index=False)
    for name,data in [('intervals.csv',intervals),('feature_transport.csv',quality),('calibration_bins.csv',bins),('cause_counts.csv',counts)]:pd.DataFrame(data).to_csv(out/name,index=False)
    result={'date':'2026-09-10','analysis_lock_commit':LOCK,'model_sha256':expected_model_sha,'flow':flow,'development_n':bundle['development_n'],
       'participant_overlap':0,'primary':{'cycle':2011,'horizon':5,'metric':'six_state_brier'},'secondary':{'cycle':2013,'horizon':4},
       'bootstrap_replicates':1000,'training_uncertainty_included':False,'evaluated_on_new_cycles':True,
       'fit_on_evaluation':False,'clinical_use_ready':False,'noninfectious_risk_identified':False,'CRP_used':False,
       'maximum_probability_mass_error':mass,'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
       'participant_rows_exported':False,'default_workbench_changed':False,'source_manifest_sha256':sha(source/'manifest.json')}
    dump(out/'results.json',result);dump(out/'source_manifest.json',json.loads((source/'manifest.json').read_text()))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    t=sub.add_parser('train')
    for field in ('old-source','source08','source09','out'):t.add_argument('--'+field,type=Path,required=True)
    e=sub.add_parser('evaluate')
    for field in ('source','trained','out'):e.add_argument('--'+field,type=Path,required=True)
    e.add_argument('--expected-model-sha',required=True);a=p.parse_args()
    if a.command=='train':train(a.old_source,a.source08,a.source09,a.out)
    else:evaluate_new(a.source,a.trained,a.out,a.expected_model_sha)
