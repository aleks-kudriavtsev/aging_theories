"""Executed detailed-cause study, fixed protocol and availability amendment.

Individual rows and predictions are never exported. The test split holds out
entire design strata; this is internal validation, not RU/DE or temporal external.
"""
from pathlib import Path
import argparse,json,math,hashlib,datetime
import numpy as np
import pandas as pd
from .cohort import load,TARGETS,IDS,CLINICAL,CORE,HEPATIC,sha
from .model import fit,predict
from .scoring import outcomes,metrics,aalen_johansen,resample,calibration

PENALTIES=[1.,10.,100.]

def dump(path,x):
    Path(path).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def arr(frame,causes):
    mapping={c:i+1 for i,c in enumerate(causes)}
    c=np.where(frame.dead.eq(1),frame.model_cause.map(mapping),0).astype(int)
    return frame.time.to_numpy(float),frame.dead.to_numpy(int),c

def evaluate(frame,model,h,w=None):
    if w is None:w=frame.weight.to_numpy(float)
    t,d,c=arr(frame,model['cause_order']);p=predict(frame,model,h);y,ew,info=outcomes(t,d,c,w,h,len(model['cause_order']))
    rows=metrics(y,p,w,ew)
    for r in rows:
        if isinstance(r['outcome'],int):r['outcome']=model['cause_order'][r['outcome']-1]
    return rows,p,(y,ew,info)

def study(eligible,features,split,selected,out,name,bootstrap=True):
    mask=np.isfinite(eligible[features['laboratory']]).all(axis=1).to_numpy()
    d=eligible.loc[mask].copy();train=d[~d.stratum.isin(split)].copy();test=d[d.stratum.isin(split)].copy()
    test_eligible=eligible[eligible.stratum.isin(split)].copy();test_mask=np.isfinite(test_eligible[features['laboratory']]).all(axis=1).to_numpy()
    causes=selected+['infection','external','other_or_unresolved']
    for f in (train,test,test_eligible):f['model_cause']=f.cause.where(f.cause.isin(causes),'other_or_unresolved')
    if len(train)<100 or len(test)<100:raise ValueError('Insufficient prespecified comparison domain')
    if set(train.id)&set(test.id) or set(train.stratum)&set(test.stratum):raise ValueError('Partition leakage')
    folds=np.array_split(np.random.default_rng(20260916).permutation(sorted(train.stratum.unique())),4)
    cv=[];models={}
    for panel,variables in features.items():
        for penalty in PENALTIES:
            for i,fold in enumerate(folds):
                training=train[~train.stratum.isin(fold)];validation=train[train.stratum.isin(fold)]
                model=fit(training,variables,causes,penalty)
                rows,_,_=evaluate(validation,model,15)
                loss=rows[0]['brier'];cv.append(dict(study=name,panel=panel,penalty=penalty,fold=i,n=len(validation),multistate_brier=loss))
            print('CV',name,panel,penalty,float(np.mean([r['multistate_brier'] for r in cv if r['panel']==panel and r['penalty']==penalty])),flush=True)
        best=min(PENALTIES,key=lambda v:(np.mean([r['multistate_brier'] for r in cv if r['panel']==panel and r['penalty']==v]),-v))
        models[panel]=fit(train,variables,causes,best)
    # Both full-development models now frozen; no held-out results used above.
    model_artifact={'study':name,'version':'0.20.0','model_kind':'sex_stratified_penalized_piecewise_cause_specific_hazards',
        'data_source':'NHANES I baseline 1971-1975 / NHEFS1992 public followup',
        'country':'historical_US','clinical_use_ready':False,'horizons':[10,15],
        'models':models,'target_taxonomy_size':20,'estimated_noninfectious_families':selected,
        'additional_competitors':['infection','external','other_or_unresolved'],
        'validation_design':'heldout_survey_strata_internal','heldout_strata':list(map(int,split)),
        'training_n':len(train),'test_n':len(test),'source_plan_commit':'858c7e5a7a13799a093d64eb1295b16af2a65a86',
        'availability_amendment_commit':'95cc634c8492e26079f43766e71c35648a1e602d'}
    out.mkdir(parents=True);dump(out/'model.json',model_artifact)
    table=[];subgroup=[];prediction={};audit=[]
    for h in [10,15]:
        for panel,model in models.items():
            rows,p,o=evaluate(test,model,h);prediction[panel,h]=p
            y,ew,info=o;t,dead,c=arr(test,causes);w=test.weight.to_numpy(float)
            aj=aalen_johansen(t,dead,c,w,h,len(causes));ip=np.bincount(y,weights=ew,minlength=len(causes)+1)/w.sum()
            error=float(np.max(abs(aj-ip)))
            if error>1e-9:raise ArithmeticError('Aalen-Johansen/IPC weighting mismatch')
            audit.append(dict(horizon=h,panel=panel,AJ_IPCW_difference=error,**info))
            for r in rows:table.append(dict(study=name,panel=panel,horizon=h,n=len(test),**r))
            masks={'female':test.sex.eq(2),'male':test.sex.eq(1),'age25_49':test.age.lt(50),'age50_64':test.age.between(50,64),'age65_74':test.age.ge(65)}
            if h==15:
                for label,ms in masks.items():
                    sub=test[ms];rr,_,_=evaluate(sub,model,h)
                    for r in rr:subgroup.append(dict(study=name,subgroup=label,panel=panel,horizon=h,n=len(sub),**r))
                calc=calibration(y>0,1-p[:,0],ew)
                audit.append(dict(panel=panel,horizon=h,calibration=calc))
    print('HELDOUT_RESULT',name,flush=True)
    print(pd.DataFrame(table).query("horizon==15 and outcome in ['all_cause','multistate']").to_string(index=False),flush=True)
    draws=[];failed=0;rng=np.random.default_rng(20260920)
    t,dead,c=arr(test,causes)
    for rep in range(1000 if bootstrap else 0):
        w=resample(test_eligible,rng)[test_mask]
        try:y,ew,_=outcomes(t,dead,c,w,15,len(causes))
        except ValueError:draws.append({});failed+=1;continue
        row={}
        for panel in features:
            for r in metrics(y,prediction[panel,15],w,ew):
                outcome=causes[r['outcome']-1] if isinstance(r['outcome'],int) else r['outcome']
                for key in ['brier','auc','mean_bias_pp']:
                    if key in r:row[f'{panel}|{outcome}|{key}']=r[key]
        for key in list(row):
            if key.startswith('laboratory|'):
                control=key.replace('laboratory|','clinical|',1)
                a,b=row[key],row[control]
                row[key.replace('laboratory|','delta_lab_minus_clinical|')]=a-b if a is not None and b is not None else None
        draws.append(row)
    intervals=[]
    if draws:
        boot=pd.DataFrame(draws)
        for col in boot:
            valid=boot[col].dropna();lo,hi=np.quantile(valid,[.025,.975]) if len(valid) else (None,None)
            panel,cause,metric=col.split('|')
            intervals.append(dict(panel_or_contrast=panel,outcome=cause,metric=metric,horizon=15,lower=lo,upper=hi,
                valid_replicates=len(valid),undefined_replicates=1000-len(valid),
                interpretation='conditional_fixed_models' if len(valid)==1000 else 'defined_draws_only_not_nominal_95pct'))
    counts=[]
    for key,label,icd9,icd10 in TARGETS:
        counts.append(dict(cause=key,label=label,ICD9=icd9,ICD10_target=icd10,
            development_events20=int(((train.cause==key)&(train.time<=20)&train.dead.eq(1)).sum()),
            test_events15=int(((test.cause==key)&(test.time<=15)&test.dead.eq(1)).sum()),
            selected=key in selected))
    for fn,data in [('metrics',table),('intervals',intervals),('subgroups',subgroup),('training_cv',cv),('cause_counts',counts)]:
        pd.DataFrame(data).to_csv(out/(fn+'.csv'),index=False)
    result=dict(study=name,training_n=len(train),test_n=len(test),test_events15=int(((test.time<=15)&test.dead.eq(1)).sum()),
        estimated_noninfectious_families=len(selected),total_estimated_hazards=len(causes),
        outcome_states=len(causes)+1,heldout_strata=list(map(int,split)),training_strata=int(train.stratum.nunique()),
        test_strata=int(test.stratum.nunique()),raw_participant_data_exported=False,
        feature_coefficients_fitted_on_test=False,model_sha256=sha(out/'model.json'),
        bootstrap_replicates=1000 if bootstrap else 0,bootstrap_censoring_failures=failed,
        clinical_use_ready=False,geographic_external_validation=False,contemporary_risk=False,audit=audit)
    dump(out/'results.json',result)
    return result

def main(source,out):
    out=Path(out)
    if out.exists():raise ValueError('Choose a new immutable output directory')
    eligible,audit=load(source)
    strata=np.array(sorted(eligible.stratum.unique()))
    split=np.random.default_rng(20260915).permutation(strata)[:math.ceil(.3*len(strata))]
    complete=np.isfinite(eligible[CORE]).all(axis=1)
    development=eligible[complete&~eligible.stratum.isin(split)]
    count=development[development.dead.eq(1)&development.time.le(20)].cause.value_counts()
    selected=[k for k in IDS if count.get(k,0)>=10]
    if len(selected)<10:raise ValueError('Fewer than10 supported noninfectious heads; do not relabel competitors')
    print('SELECTED_DEVELOPMENT_ONLY',len(selected),selected,flush=True)
    out.mkdir(parents=True);dump(out/'cohort_audit.json',audit)
    results=[study(eligible,{'clinical':CLINICAL,'laboratory':CORE},split,selected,out/'primary','primary',True)]
    detail=eligible[eligible.subsample.eq(2)&eligible.weight_detail.gt(0)].copy();detail['weight']=detail.weight_detail
    # Same stratum split and same cause partition, not chosen after reading primary results.
    features={'clinical':CLINICAL+['smoking_current','smoking_former'],
              'laboratory':HEPATIC+['smoking_current','smoking_former']}
    results.append(study(detail,features,split,selected,out/'hepatic','secondary_hepatic_smoking',True))
    dump(out/'run_summary.json',{'results':results,'analysis_plan_commit':'858c7e5a7a13799a093d64eb1295b16af2a65a86',
       'amendment_commit':'95cc634c8492e26079f43766e71c35648a1e602d','utc_completed':datetime.datetime.now(datetime.timezone.utc).isoformat(),
       'no_raw_data_exported':True,'new_models':True,'external_RU_DE_validation':False,
       'code_sha256':{p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}})
    (out/'SOURCE_MANIFEST.json').write_bytes((Path(source)/'manifest.json').read_bytes())
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.source,a.out)
