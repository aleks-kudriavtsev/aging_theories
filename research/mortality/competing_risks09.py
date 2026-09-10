"""Coherent cause-specific mortality models; research only, no patient deployment.

Development NHANES2005-2008; new test2009-2010. ANALYSIS_LOCK.md defines the
three mutually exclusive groups and training-only penalty selection. This file
never drops a death because its cause is unknown and never allocates by WHO shares.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from sklearn.metrics import roc_auc_score
from .validate_temporal08 import (new_cohort, engineer_new, independent_mortality_reader,
    audit_mortality_reader, design_bootstrap_weights, weighted_auc, verify_new_sources,
    COMPONENTS, MODEL_SHA)
from .xpt_checked import read_xpt_checked

LOCK = '99e8c150afc34d2c8c53077e7263547628638c8b'
GROUPS = ('heart_diseases', 'malignant_neoplasms', 'other_or_unknown')
ENDS = np.array([1., 5., 8.])
HORIZONS = (1., 5., 8.)
PENALTIES = (1., 10., 100.)
SEED = 20260911


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2,
                                   allow_nan=False)+'\n', encoding='utf-8')


def classify(dead, causes):
    dead, causes = np.asarray(dead, float), np.asarray(causes, float)
    if dead.ndim != 1 or causes.shape != dead.shape or not np.isin(dead, [0, 1]).all():
        raise ValueError('Invalid vital status/shape')
    if np.isinf(causes).any(): raise ValueError('Infinite cause code is not a missing cause')
    known = np.isfinite(causes)
    if not np.isin(causes[known], range(1, 11)).all():
        raise ValueError('Unexpected cause code, review source definition')
    if np.any((dead == 0) & known):
        raise ValueError('Cause of death assigned to a survivor')
    # -1 is right censoring/alive; 2 explicitly includes deaths of unknown cause.
    return np.where(dead == 0, -1, np.where(causes == 1, 0,
                    np.where(causes == 2, 1, 2))).astype(int)


def exposure_design(time, event_group):
    time = np.asarray(time, float); event_group = np.asarray(event_group)
    if (time.ndim != 1 or event_group.shape != time.shape or not len(time)
            or not np.isfinite(time).all() or np.any(time <= 0)
            or not np.isin(event_group, [-1, 0, 1, 2]).all()):
        raise ValueError('Invalid time or event group')
    starts = np.r_[0, ENDS[:-1]]
    exposure = np.maximum(0, np.minimum(time[:, None], ENDS)-starts)
    event_interval = np.searchsorted(ENDS, time, side='left')
    events = np.zeros((len(time), len(ENDS), len(GROUPS)))
    idx = np.flatnonzero((event_group >= 0) & (time <= ENDS[-1]))
    events[idx, event_interval[idx], event_group[idx]] = 1
    return exposure, events


def preprocess_fit(frame, features):
    x = frame[features].to_numpy(float)
    if not len(x) or not np.isfinite(x).all():
        raise ValueError('Complete finite training predictors required')
    lower, upper = np.quantile(x, [.005, .995], axis=0)
    clipped = np.clip(x, lower, upper)
    mean = clipped.mean(axis=0); sd = clipped.std(axis=0)
    sd = np.where(sd < 1e-10, 1., sd)
    return {'features':list(features), 'lower':lower.tolist(), 'upper':upper.tolist(),
            'mean':mean.tolist(), 'sd':sd.tolist(),
            'method':'training_only_unweighted_winsor005_995_population_sd',
            'imputation':False}


def transform(frame, preprocessing):
    x = frame[preprocessing['features']].to_numpy(float)
    if not np.isfinite(x).all():
        raise ValueError('Missing predictor; no imputation in this model')
    return ((np.clip(x, preprocessing['lower'], preprocessing['upper']) -
             np.asarray(preprocessing['mean'])) / np.asarray(preprocessing['sd']))


def fit_one_cause(z, exposure, events, weights, penalty):
    """Profile unpenalized interval baselines out of the exact Poisson likelihood.

    Objective differs from negative log likelihood by data-only constants.
    Profiling avoids unstable intercept iterations; it does not add pseudo-events.
    """
    z, exposure, events, weights = [np.asarray(v, float) for v in (z, exposure, events, weights)]
    if (z.ndim != 2 or exposure.shape != (len(z), len(ENDS))
            or events.shape != exposure.shape or weights.shape != (len(z),)
            or not len(z) or not all(np.isfinite(v).all() for v in (z, exposure, events, weights))
            or np.any(weights <= 0) or np.any(exposure < 0)
            or not np.isin(events, [0., 1.]).all()
            or np.any(events.sum(axis=1) > 1) or np.any((events > 0) & (exposure <= 0))
            or isinstance(penalty, (bool, np.bool_)) or not np.isscalar(penalty)
            or not np.isfinite(penalty) or penalty < 0):
        raise ValueError('Invalid likelihood inputs')
    w = weights/weights.mean()
    counts = (w[:, None]*events).sum(axis=0)
    if np.any(counts <= 0):
        raise ValueError('No events in a cause/interval cell; locked model cannot be fitted')
    log_mass = np.full_like(exposure, -np.inf)
    positive = exposure > 0
    log_mass[positive] = np.log((w[:, None]*exposure)[positive])
    event_x = z.T @ (w*events.sum(axis=1))
    def objective(beta):
        log_terms = log_mass + (z @ beta)[:, None]
        totals = logsumexp(log_terms, axis=0)
        allocation = np.exp(log_terms-totals)
        gradient = z.T@(allocation@counts) - event_x + penalty*beta
        value = counts@totals - event_x@beta + .5*penalty*(beta@beta)
        return float(value), gradient
    result = minimize(objective, np.zeros(z.shape[1]), jac=True, method='BFGS',
                      options={'gtol':1e-6, 'maxiter':2000})
    value, gradient = objective(result.x)
    max_gradient = float(np.max(np.abs(gradient)))
    if not np.isfinite(value) or max_gradient > 1e-4:
        raise ValueError(f'Optimization did not converge: gradient={max_gradient}; {result.message}')
    baseline = np.log(counts)-logsumexp(log_mass+(z@result.x)[:,None], axis=0)
    return {'coefficients':result.x.tolist(), 'baseline_log_hazards':baseline.tolist(),
            'weighted_events_mean1':counts.tolist(), 'raw_events':events.sum(axis=0).astype(int).tolist(),
            'gradient_max':max_gradient, 'iterations':int(result.nit), 'objective':value}


def fit_panel(frame, features, penalty):
    pp = preprocess_fit(frame, features); z = transform(frame, pp)
    group = classify(frame.dead, frame.cause)
    exposure, events = exposure_design(frame.time_years, group)
    models = [fit_one_cause(z, exposure, events[:,:,k], frame.weight, penalty) for k in range(3)]
    return {'preprocessing':pp, 'causes':dict(zip(GROUPS, models)),
            'interval_ends_years':ENDS.tolist(), 'penalty':float(penalty),
            'training_n':len(frame), 'method':'survey_weighted_ridge_piecewise_exponential'}


def integrate_hazards(rates, horizon, ends=ENDS):
    rates = np.asarray(rates, float); ends = np.asarray(ends, float)
    if (rates.ndim != 3 or rates.shape[1] != len(ends) or rates.shape[2] < 1
            or ends.ndim != 1 or not len(ends) or not np.isfinite(ends).all()
            or np.any(np.diff(np.r_[0., ends]) <= 0) or not np.isfinite(rates).all()
            or np.any(rates < 0) or isinstance(horizon, bool)
            or not np.isscalar(horizon) or not np.isfinite(horizon)
            or not 0 <= horizon <= ends[-1]):
        raise ValueError('Invalid hazard array/horizon; no extrapolation')
    survival = np.ones(rates.shape[0]); cif = np.zeros((rates.shape[0], rates.shape[2]))
    for j, (lower, upper) in enumerate(zip(np.r_[0., ends[:-1]], ends)):
        dt = max(0., min(float(horizon), upper)-lower)
        total = rates[:,j,:].sum(axis=1)
        if not np.isfinite(total).all(): raise ValueError('Summed hazards overflowed')
        event_mass = -np.expm1(-total*dt)
        shares = np.divide(rates[:,j,:], total[:,None], out=np.zeros_like(cif), where=total[:,None]>0)
        cif += survival[:,None]*event_mass[:,None]*shares
        survival *= np.exp(-total*dt)
    probabilities = np.column_stack([survival, cif])
    if np.max(np.abs(probabilities.sum(axis=1)-1), initial=0) > 1e-10:
        raise ArithmeticError('Competing probabilities do not conserve mass')
    return probabilities


def predict(frame, panel, horizon):
    z = transform(frame, panel['preprocessing'])
    rates = np.stack([np.exp(z@np.asarray(panel['causes'][c]['coefficients'])[:,None] +
                            np.asarray(panel['causes'][c]['baseline_log_hazards'])[None,:]) for c in GROUPS],axis=2)
    return integrate_hazards(rates, horizon, panel['interval_ends_years'])


def labels_at(frame, horizon):
    if ((frame.dead == 0) & (frame.time_years < horizon)).any():
        raise ValueError('Early censoring: fixed-horizon outcome unavailable')
    groups = classify(frame.dead, frame.cause)
    return np.where((groups >= 0) & (frame.time_years <= horizon), groups+1, 0)


def multiclass_brier(labels, probabilities, weights):
    labels = np.asarray(labels); p = np.asarray(probabilities, float); w = np.asarray(weights, float)
    if (p.ndim != 2 or labels.shape != (len(p),) or w.shape != labels.shape
            or not np.issubdtype(labels.dtype, np.integer) or np.any(labels < 0) or np.any(labels >= p.shape[1])
            or not np.isfinite(p).all() or not np.isfinite(w).all() or np.any(w < 0) or w.sum() <= 0
            or np.any(p < 0) or np.any(p > 1) or not np.allclose(p.sum(axis=1), 1, atol=1e-10, rtol=0)):
        raise ValueError('Invalid multistate probability/label/weight')
    target = np.eye(p.shape[1])[labels]
    return float(np.average(np.sum((p-target)**2,axis=1),weights=w))


def panel_select(frame, features):
    rows = []
    for penalty in PENALTIES:
        for cycle in (2005, 2007):
            development = frame[frame.cycle != cycle]; validation = frame[frame.cycle == cycle]
            model = fit_panel(development, features, penalty)
            p = predict(validation, model, 5)
            rows.append({'penalty':penalty, 'held_out_development_cycle':cycle,
                         'multiclass_brier_5y':multiclass_brier(labels_at(validation, 5), p, validation.weight)})
    scores = {p:np.mean([r['multiclass_brier_5y'] for r in rows if r['penalty']==p]) for p in PENALTIES}
    best = min(PENALTIES, key=lambda p:(scores[p],-p))
    return fit_panel(frame, features, best), rows


def load_test(source, features):
    manifest = json.loads((source/'manifest.json').read_text(encoding='utf-8'))
    if manifest['analysis_lock_commit'] != LOCK: raise ValueError('Wrong lock provenance')
    needed = {f'2009_{c}_F.xpt' for c in COMPONENTS}|{'2009_mortality.dat'}
    seen = set()
    for row in manifest['files']:
        name = row['file']
        if Path(name).name != name or name in seen: raise ValueError('Invalid manifest entry')
        seen.add(name)
        if row['status'] != 'retrieved':
            if name in needed: raise ValueError('Required data missing: '+name)
            continue
        if sha(source/name) != row['sha256']: raise ValueError('Changed source: '+name)
    if needed-seen: raise ValueError('Incomplete source manifest')
    demo = read_xpt_checked(source/'2009_DEMO_F.xpt'); ids = set(demo.SEQN)
    frame = demo.loc[demo.RIDAGEYR.between(40,79), ['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']].copy()
    n_age = len(frame)
    for component in COMPONENTS[1:]:
        data = read_xpt_checked(source/f'2009_{component}_F.xpt')
        if not set(data.SEQN).issubset(ids): raise ValueError('Foreign ID in covariates')
        cols = ['SEQN']+[c for c in data if c not in frame]
        frame = frame.merge(data[cols],on='SEQN',how='left',validate='1:1')
    reader = audit_mortality_reader(source/'2009_mortality.dat')
    mortality = independent_mortality_reader(source/'2009_mortality.dat')
    if set(mortality.SEQN) != ids: raise ValueError('Demographic/linkage IDs differ')
    frame = frame.merge(mortality,on='SEQN',how='left',validate='1:1')
    eligible = (frame.elig.eq(1)&frame.dead.isin([0,1])&frame.t_exam.ge(0)&frame.WTMEC2YR.gt(0)
                &np.isfinite(frame.WTMEC2YR)&frame.RIAGENDR.isin([1,2]))
    frame = frame[eligible].copy(); n_eligible=len(frame)
    # BIOPRO_F is IDMS-standardized, no 2005 correction and no withdrawn HbA1c adjustment.
    frame = engineer_new(frame, 2009)
    frame['cycle']=2009;frame['weight']=frame.WTMEC2YR
    frame['time_years']=frame.t_exam.clip(lower=.5)/12
    frame['stratum']='2009:'+frame.SDMVSTRA.astype(str)
    complete = np.isfinite(frame[features].to_numpy(float)).all(axis=1)
    missing = {f:int((~np.isfinite(frame[f])).sum()) for f in features}
    frame = frame[complete].copy().reset_index(drop=True)
    if frame[['SDMVSTRA','SDMVPSU']].isna().any().any(): raise ValueError('Missing survey design')
    if ((frame.dead==0)&(frame.time_years<8)).any(): raise ValueError('Incomplete8year followup')
    return frame, {'survey_n':len(demo),'age40_79_n':n_age,'eligible_n':n_eligible,
                   'complete_profile_n':len(frame),'excluded_incomplete_n':n_eligible-len(frame),
                   'missing_by_feature':missing,'reader_audit':reader,
                   'minimum_survivor_followup':float(frame.loc[frame.dead==0,'time_years'].min()),
                   'zero_month_events':int(((frame.dead==1)&(frame.t_exam==0)).sum())}


def binary_metrics(y, p, weights):
    y, p, w = np.asarray(y,int), np.asarray(p,float),np.asarray(weights,float)
    defined = len(np.unique(y[w>0])) == 2
    auc = weighted_auc(y,p,w) if defined else None
    observed=float(np.average(y,weights=w));predicted=float(np.average(p,weights=w))
    return {'n':len(y),'events':int(y.sum()),'auc':auc,'brier':float(np.average((y-p)**2,weights=w)),
            'observed':observed,'predicted':predicted,'observed_expected':observed/predicted if predicted else None,
            'bias_pp':100*(predicted-observed)}


def main(development_source, test_source, frozen_path, out):
    if out.exists(): raise ValueError('New output directory required')
    if sha(frozen_path)!=MODEL_SHA: raise ValueError('Reference feature list changed')
    verify_new_sources(development_source)
    development, flow, readers = new_cohort(development_source)
    reference = json.loads(frozen_path.read_text(encoding='utf-8'))
    feature_lists = {label:reference['models'][old]['features'] for label,old in
                     [('clinical','M0c_clinical_only_posthoc'),('routine','M1_routine')]}
    complete = np.isfinite(development[feature_lists['routine']].to_numpy(float)).all(axis=1)
    development = development[complete].copy().reset_index(drop=True)
    if len(development)!=5287: raise ValueError('Development differs from frozen complete-profile domain')
    models={};cv=[]
    for label,features in feature_lists.items():
        models[label], rows = panel_select(development, features)
        cv.extend([dict(panel=label,**r) for r in rows])
        print('TRAINING_ONLY_SELECTED',label,models[label]['penalty'],flush=True)
    # Write and hash development outputs BEFORE test outcomes are read by this program.
    out.mkdir(parents=True)
    bundle={'version':'0.9_competing_risks','analysis_lock_commit':LOCK,'groups':list(GROUPS),
            'models':models,'country':'historical_USA','age_range':[40,79],
            'supported_horizons_years':list(HORIZONS),'target':'three_coarse_death_groups',
            'clinical_use_ready':False,'noninfectious_risk_identified':False,
            'development_cycles':[2005,2007], 'development_n':len(development),
            'reference_feature_definition_sha256':MODEL_SHA}
    dump(out/'model_bundle09.json',bundle);model_sha=sha(out/'model_bundle09.json')
    pd.DataFrame(cv).to_csv(out/'training_cv.csv',index=False)
    test,test_flow=load_test(test_source,feature_lists['routine'])
    if set(development.SEQN)&set(test.SEQN): raise ValueError('Training/test overlap')
    records=[];probabilities={};quality=[];cause_counts=[]
    for cohort,frame in [('development',development),('new_test',test)]:
        for horizon in HORIZONS:
            deceased=frame[(frame.dead==1)&(frame.time_years<=horizon)]
            for code in list(range(1,11))+[None]:
                n=int(deceased.cause.isna().sum()) if code is None else int(deceased.cause.eq(code).sum())
                cause_counts.append({'cohort':cohort,'horizon_years':horizon,'ucod_leading':code,'deaths':n})
    masks={'all':np.ones(len(test),bool),'female':test.RIAGENDR.eq(2).to_numpy(),
           'male':test.RIAGENDR.eq(1).to_numpy(),'age40_59':test.RIDAGEYR.lt(60).to_numpy(),
           'age60_79':test.RIDAGEYR.ge(60).to_numpy()}
    for label,model in models.items():
        for i,f in enumerate(model['preprocessing']['features']):
            st=model['preprocessing'];x=test[f].to_numpy(float)
            quality.append({'panel':label,'feature':f,'clipped_n':int(((x<st['lower'][i])|(x>st['upper'][i])).sum())})
        previous=None
        for horizon in HORIZONS:
            p=predict(test,model,horizon);probabilities[label,horizon]=p
            if previous is not None and (np.any(p[:,1:]<previous[:,1:]-1e-12) or np.any(p[:,0]>previous[:,0]+1e-12)):
                raise ArithmeticError('CIF/survival monotonicity failed')
            previous=p;target=labels_at(test,horizon)
            for subgroup,mask in masks.items():
                w=test.weight.to_numpy()[mask]
                records.append({'panel':label,'horizon_years':horizon,'subgroup':subgroup,'outcome':'four_state',
                                'n':int(mask.sum()),'events':int((target[mask]>0).sum()),
                                'multiclass_brier':multiclass_brier(target[mask],p[mask],w)})
                for k,name in [(1,GROUPS[0]),(2,GROUPS[1]),(3,GROUPS[2]),(0,'all_cause')]:
                    y=(target==k).astype(int) if k else (target>0).astype(int)
                    risk=p[:,k] if k else 1-p[:,0]
                    metrics=binary_metrics(y[mask],risk[mask],w)
                    if subgroup=='all' and metrics['auc'] is not None:
                        if abs(metrics['auc']-roc_auc_score(y,risk,sample_weight=test.weight))>1e-12:
                            raise ArithmeticError('Independent AUC mismatch')
                    records.append({'panel':label,'horizon_years':horizon,'subgroup':subgroup,'outcome':name,**metrics})
    metrics=pd.DataFrame(records)
    print(metrics.query("horizon_years==5 and subgroup=='all'").to_string(index=False),flush=True)
    # Paired bootstrap of the frozen test probabilities; undefined AUC is retained as missing.
    rng=np.random.default_rng(SEED);draws=[]
    for rep in range(1000):
        wr=design_bootstrap_weights(test,rng);row={}
        for horizon in (5.,8.):
            y=labels_at(test,horizon)
            for label in models:
                p=probabilities[label,horizon]
                row[f'{horizon}|{label}|four_state|brier']=multiclass_brier(y,p,wr)
                for k,name in [(1,GROUPS[0]),(2,GROUPS[1]),(3,GROUPS[2]),(0,'all_cause')]:
                    binary=(y==k) if k else y>0;risk=p[:,k] if k else 1-p[:,0]
                    row[f'{horizon}|{label}|{name}|brier']=float(np.average((binary-risk)**2,weights=wr))
                    row[f'{horizon}|{label}|{name}|auc']=weighted_auc(binary,risk,wr) if len(np.unique(binary[wr>0]))==2 else np.nan
            for name in ['four_state',*GROUPS,'all_cause']:
                for metric in (['brier'] if name=='four_state' else ['brier','auc']):
                    row[f'{horizon}|delta_routine_clinical|{name}|{metric}']=row[f'{horizon}|routine|{name}|{metric}']-row[f'{horizon}|clinical|{name}|{metric}']
        draws.append(row)
    bootstrap=pd.DataFrame(draws);intervals=[]
    for key in bootstrap:
        horizon,label,outcome,metric=key.split('|');x=bootstrap[key];defined=np.isfinite(x)
        lo,hi=np.quantile(x[defined],[.025,.975]) if defined.any() else (None,None)
        intervals.append({'horizon_years':float(horizon),'panel_or_contrast':label,'outcome':outcome,'metric':metric,
                          'lower':lo,'upper':hi,'valid_replicates':int(defined.sum()),'undefined_replicates':int((~defined).sum())})
    metrics.to_csv(out/'metrics.csv',index=False)
    pd.DataFrame(intervals).to_csv(out/'intervals.csv',index=False)
    pd.DataFrame(cause_counts).to_csv(out/'source_cause_counts.csv',index=False)
    pd.DataFrame(quality).to_csv(out/'feature_transport.csv',index=False)
    summary={'date':'2026-09-10','model_sha256':model_sha,'analysis_lock_commit':LOCK,
             'development_n':len(development),'test_n':len(test),'test_flow':test_flow,
             'development_flow':flow,'development_readers':readers,'participant_overlap':0,
             'selected_penalties':{k:v['penalty'] for k,v in models.items()},
             'test_deaths_by_horizon':{str(t):int((labels_at(test,t)>0).sum()) for t in HORIZONS},
             'bootstrap_replicates':1000,'bootstrap_seed':SEED,'bootstrap_conditional_on_fitted_models':True,
             'training_uncertainty_included':False,'clinical_use_ready':False,
             'default_workbench_changed':False,'participant_rows_exported':False,
             'maximum_probability_mass_error':max(float(np.abs(p.sum(axis=1)-1).max()) for p in probabilities.values()),
             'test_source_manifest_sha256':sha(test_source/'manifest.json'),
             'development_source_manifest_sha256':sha(development_source/'manifest.json'),
             'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
    dump(out/'results.json',summary)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('development-source','test-source','frozen-model','out'):parser.add_argument('--'+key,type=Path,required=True)
    a=parser.parse_args();print(json.dumps(main(a.development_source,a.test_source,a.frozen_model,a.out),indent=2))
