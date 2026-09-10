"""Locked compact-panel study. Frozen comparators; no evaluation-time fitting.

Training code depends on existing reviewed iteration12 methods. The late public
cause partition has THREE, not five, identifiable death groups. Only aggregates
are exported; complete-profile comparisons use the same participant domain.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, json
import numpy as np
import pandas as pd
from .transport_no_crp12 import (cohort, CLINICAL, ROUTINE, FOLDS, PENALTIES, GROUPS,
    fit_panel, predict, labels, multiclass_brier, binary_metrics, sha, dump,
    engineer_no_crp, COMPONENTS, read_xpt_checked, independent_mortality_reader,
    audit_mortality_reader, design_bootstrap_weights, diagnostics)

LOCK='2389c6edeb806a3f6c8d3b19faff1ba0b446e0c5'
FULL_SHA='40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094'
SIX_SHA='7c7dc918f042ff44f141a9ad3b49e6573db7f9732c5eb9d93e29b09b0686fa24'
COMPACT=CLINICAL+['hba1c','egfr_low','egfr_high','log_uacr','albumin']
NAMES=('clinical','compact4','six6','full8')
OUTCOMES=('heart_diseases','malignant_neoplasms','other_or_unknown')
HORIZONS=(1,3)
SEED=20260914
ROOT=Path(__file__).parent


def aggregate_states(probabilities):
    p=np.asarray(probabilities,float)
    if (p.ndim!=2 or p.shape[1]!=6 or not np.isfinite(p).all()
        or np.any((p<0)|(p>1)) or not np.allclose(p.sum(axis=1),1,atol=1e-10,rtol=0)):
        raise ValueError('Six coherent states required')
    return np.column_stack([p[:,:3],p[:,3:].sum(axis=1)])


def late_labels(dead,cause,time,horizon):
    d,c,t=map(lambda v:np.asarray(v,float),(dead,cause,time))
    if (d.ndim!=1 or c.shape!=d.shape or t.shape!=d.shape or not len(d)
        or not np.isin(d,[0,1]).all() or not np.isfinite(t).all() or np.any(t<0)
        or np.isinf(c).any() or type(horizon) not in (int,float) or horizon not in HORIZONS):
        raise ValueError('Invalid late-cycle observations or horizon')
    known=np.isfinite(c)
    if np.any((d==0)&known) or not np.isin(c[known],[1,2,10]).all():
        raise ValueError('Late public release does not identify ten detailed causes')
    if np.any((d==0)&(t<horizon)):
        raise ValueError('Early censoring: fixed-horizon analysis stopped')
    group=np.where(c==1,1,np.where(c==2,2,3))
    return np.where((d==1)&(t<=horizon),group,0).astype(int)


def comparator_models():
    p=ROOT/'transport12/model_bundle12.json';s=ROOT/'transport12/model_sensitivity12.json'
    if sha(p)!=FULL_SHA or sha(s)!=SIX_SHA:raise ValueError('Frozen comparator changed')
    b=json.loads(p.read_text());six=json.loads(s.read_text())
    return {'clinical':b['models']['clinical'],'six6':six['model'],'full8':b['models']['routine_no_crp']}


def train(old_source,source08,source09,out):
    if out.exists():raise ValueError('Use a new training output directory')
    roots={y:old_source if y<2005 else source08 if y<2009 else source09 for y in (1999,2001,2003,2005,2007,2009)}
    d,flow,used,reader=cohort(roots)
    expected=json.loads((ROOT/'transport12/expected_development_sources.json').read_text())
    for row in expected:
        if used.get(row['file'],{}).get('sha256')!=row['sha256']:raise ValueError('Development source changed')
    d=d[np.isfinite(d[ROUTINE]).all(axis=1)].copy().reset_index(drop=True)
    if len(d)!=15280:raise ValueError('Development comparison domain changed')
    d['time_years']=d.t_exam.clip(lower=.5)/12;d['weight']=d.pooled_weight
    scores={};cv=[]
    for penalty in PENALTIES:
        vals=[]
        for block in FOLDS:
            mask=d.cycle.isin(block);m=fit_panel(d[~mask],COMPACT,penalty)
            v=multiclass_brier(labels(d[mask],5),predict(d[mask],m,5),d.loc[mask,'weight'])
            vals.append(v);cv.append({'penalty':penalty,'heldout_cycles':str(block),'n':int(mask.sum()),'six_state_brier5y':v})
        scores[penalty]=float(np.mean(vals));print('TRAINING_CV',penalty,scores[penalty],flush=True)
    best=min(PENALTIES,key=lambda p:(scores[p],-p));m=fit_panel(d,COMPACT,best)
    comparators=comparator_models()
    spec={'version':'0.14_compact_candidate','analysis_lock_commit':LOCK,
      'base_commit':'a50a9062a5527fcd5cc146883a8974326ce718b6','model':m,'groups':list(GROUPS),
      'reported_measurements':['HbA1c','creatinine/eGFR','UACR','serum_albumin'],
      'laboratory_reported_measure_count':4,'minimum_underlying_analyte_determinations':5,
      'development_n':len(d),'development_cycles':[1999,2001,2003,2005,2007,2009],
      'comparator_hashes':{'full12':FULL_SHA,'six12':SIX_SHA},'training_source_count':len(used),
      'clinical_use_ready':False,'test_cycle':2015,'horizons':[1,3],
      'evaluation_data_seen':False,'fitted_at_utc':datetime.now(timezone.utc).isoformat(),
      'training_code_sha256':sha(__file__)}
    out.mkdir(parents=True);dump(out/'model_compact14.json',spec)
    pd.DataFrame(cv).to_csv(out/'training_cv.csv',index=False)
    dump(out/'development_sources.json',list(used.values()))
    dump(out/'development_audit.json',{'n':len(d),'events8':int(((d.dead==1)&(d.time_years<=8)).sum()),
        'reader_cells':sum(r['cell_comparisons'] for r in reader),'source_count':len(used),'flow':flow})
    np.save(out/'LOCAL_ONLY_DEVELOPMENT_IDS.npy',d.SEQN.to_numpy())
    print('FROZEN',sha(out/'model_compact14.json'),flush=True)
    return spec


def load_new(source,expected_model_sha):
    doc=json.loads((source/'manifest.json').read_text())
    if doc.get('analysis_lock_commit')!=LOCK or doc.get('model_sha256_before_retrieval')!=expected_model_sha:
        raise ValueError('Unfrozen retrieval or wrong study')
    manifests={};verified=[]
    for r in doc['files']:
        name=r['file']
        if Path(name).name!=name or name in manifests:raise ValueError('Unsafe/duplicate source name')
        manifests[name]=r
        if r['status']=='retrieved':
            if sha(source/name)!=r['sha256']:raise ValueError('Source changed: '+name)
            verified.append(name)
    def read(name):
        if name not in verified:raise ValueError('Required data absent: '+name)
        return read_xpt_checked(source/name) if name.endswith('.xpt') else independent_mortality_reader(source/name)
    demo=read('2015_DEMO_I.xpt');ids=set(demo.SEQN)
    if demo.SEQN.duplicated().any():raise ValueError('Duplicate demographic ID')
    d=demo.loc[demo.RIDAGEYR.between(40,79),['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']].copy();age_n=len(d)
    for comp in COMPONENTS[1:]:
        raw=read(f'2015_{comp}_I.xpt')
        if not set(raw.SEQN)<=ids:raise ValueError('Foreign cycle IDs')
        d=d.merge(raw[['SEQN']+[c for c in raw if c not in d]],on='SEQN',how='left',validate='1:1')
    mort=read('2015_mortality.dat')
    if set(mort.SEQN)!=ids:raise ValueError('Mortality/demo ID mismatch')
    d=d.merge(mort,on='SEQN',how='left',validate='1:1')
    valid=d.elig.eq(1)&d.dead.isin([0,1])&d.t_exam.notna()&d.t_exam.ge(0)&d.WTMEC2YR.gt(0)&np.isfinite(d.WTMEC2YR)&d.RIAGENDR.isin([1,2])
    d=engineer_no_crp(d[valid].copy()).reset_index(drop=True)
    if d[['SDMVSTRA','SDMVPSU']].isna().any().any():raise ValueError('Missing survey design')
    d['weight']=d.WTMEC2YR;d['stratum']='2015:'+d.SDMVSTRA.astype(str)
    d['time_years']=d.t_exam.clip(lower=.5)/12;d['cycle']=2015
    late_labels(d.dead,d.cause,d.time_years,3) # Coverage on entire eligible sample, not only convenient complete cases.
    flow={'survey_n':len(demo),'age40_79_n':age_n,'eligible_n':len(d),
      'minimum_survivor_followup_years':float(d.loc[d.dead==0,'time_years'].min()),
      'unknown_cause_deaths':int((d.dead.eq(1)&d.cause.isna()).sum()),
      'zero_month_deaths':int((d.dead.eq(1)&d.t_exam.eq(0)).sum()),
      'verified_input_files':len(verified),'scientific_data_files':sum(n.endswith(('.dat','.xpt')) for n in verified),
      'reader_audit':audit_mortality_reader(source/'2015_mortality.dat')}
    return d,flow


def evaluate(source,trained,out,expected_model_sha):
    path=trained/'model_compact14.json'
    if out.exists() or sha(path)!=expected_model_sha:raise ValueError('Output exists or frozen model changed')
    spec=json.loads(path.read_text());models=comparator_models();models['compact4']=spec['model']
    eligible,flow=load_new(source,expected_model_sha)
    if set(eligible.SEQN)&set(np.load(trained/'LOCAL_ONLY_DEVELOPMENT_IDS.npy')):raise ValueError('Participant overlap')
    common=np.isfinite(eligible[ROUTINE]).all(axis=1).to_numpy()
    d=eligible.loc[common].copy();w=d.weight.to_numpy(float)
    if not len(d):raise ValueError('Empty common evaluation domain')
    coverage=[]
    for name in NAMES:
        good=np.isfinite(eligible[models[name]['preprocessing']['features']]).all(axis=1).to_numpy()
        coverage.append({'panel':name,'eligible_n':len(eligible),'complete_n':int(good.sum()),
            'additional_complete_vs_common':int(good.sum()-common.sum()),
            'weighted_coverage':float(np.average(good,weights=eligible.weight)),
            'performance_compared_on_common_domain_only':True})
    metrics=[];bins=[];predictions={};counts=[];quality=[]
    masks={'all':np.ones(len(d),bool),'female':d.RIAGENDR.eq(2).to_numpy(),'male':d.RIAGENDR.eq(1).to_numpy(),
      'age40_59':d.RIDAGEYR.lt(60).to_numpy(),'age60_79':d.RIDAGEYR.ge(60).to_numpy()}
    for h in HORIZONS:
        target=late_labels(d.dead,d.cause,d.time_years,h)
        counts.extend({'horizon':h,'outcome':name,'n':len(d),'events':int((target==k).sum())}
            for k,name in enumerate(('alive',)+OUTCOMES))
        for name in NAMES:
            p=aggregate_states(predict(d,models[name],h));predictions[name,h]=p
            for group,mask in masks.items():
                metrics.append({'panel':name,'horizon':h,'subgroup':group,'outcome':'four_state','n':int(mask.sum()),
                    'events':int((target[mask]>0).sum()),'brier':multiclass_brier(target[mask],p[mask],w[mask])})
                for k,label in list(enumerate(OUTCOMES,1))+[(0,'all_cause')]:
                    y=target==k if k else target>0;risk=p[:,k] if k else 1-p[:,0]
                    met=binary_metrics(y[mask],risk[mask],w[mask])
                    if group=='all' and h==3:
                        try:met.update(diagnostics(y,risk,w));met['calibration_status']='computed'
                        except ValueError:met['calibration_status']='not_converged_no_refit'
                    metrics.append({'panel':name,'horizon':h,'subgroup':group,'outcome':label,**met})
                    if group=='all' and h==3:
                        cuts=np.quantile(risk,np.linspace(0,1,6));assign=np.digitize(risk,cuts[1:-1],right=True)
                        for b in range(5):
                            m=assign==b
                            if m.any():bins.append({'panel':name,'outcome':label,'bin':b+1,**binary_metrics(y[m],risk[m],w[m])})
            for j,f in enumerate(models[name]['preprocessing']['features']):
                pp=models[name]['preprocessing'];quality.append({'panel':name,'horizon':h,'feature':f,
                  'clipped_n':int(((d[f]<pp['lower'][j])|(d[f]>pp['upper'][j])).sum())})
    print(pd.DataFrame(metrics).query('horizon==3 and subgroup=="all"')[['panel','outcome','n','events','brier','auc','weighted_predicted']].to_string(index=False),flush=True)
    target=late_labels(d.dead,d.cause,d.time_years,3);rng=np.random.default_rng(SEED);draws=[]
    comparisons=(('compact4','clinical'),('full8','clinical'),('compact4','full8'),('compact4','six6'))
    for _ in range(1000):
        # Survey DOMAIN variance: first resample eligible PSUs, then restrict.
        wr=design_bootstrap_weights(eligible,rng)[common];row={}
        for name in NAMES:
            p=predictions[name,3];row[name+'|four_state|brier']=multiclass_brier(target,p,wr)
            for k,label in list(enumerate(OUTCOMES,1))+[(0,'all_cause')]:
                y=target==k if k else target>0;risk=p[:,k] if k else 1-p[:,0]
                r=binary_metrics(y,risk,wr)
                for metric in ('auc','brier','bias_pp'):row[name+'|'+label+'|'+metric]=r[metric]
        for left,right in comparisons:
            for outcome in ('four_state',)+OUTCOMES+('all_cause',):
                for metric in (['brier'] if outcome=='four_state' else ['auc','brier','bias_pp']):
                    a,b=row[left+'|'+outcome+'|'+metric],row[right+'|'+outcome+'|'+metric]
                    row[left+'_minus_'+right+'|'+outcome+'|'+metric]=a-b if a is not None and b is not None else None
        draws.append(row)
    boot=pd.DataFrame(draws);intervals=[]
    for key in boot:
        panel,outcome,metric=key.split('|');x=boot[key].dropna();lo,hi=np.quantile(x,[.025,.975]) if len(x) else (None,None)
        intervals.append({'panel_or_contrast':panel,'outcome':outcome,'metric':metric,'horizon':3,'lower':lo,'upper':hi,
          'valid_replicates':len(x),'undefined_replicates':1000-len(x),
          'interpretation':'conditional_on_fixed_models' if len(x)==1000 else 'defined_draws_only_not_nominal_coverage'})
    out.mkdir(parents=True)
    for name,data in [('metrics',metrics),('intervals',intervals),('coverage',coverage),('calibration_bins',bins),('counts',counts),('feature_quality',quality)]:
        pd.DataFrame(data).to_csv(out/(name+'.csv'),index=False)
    result={'analysis_lock_commit':LOCK,'model_sha256':expected_model_sha,'status':'new_temporal_validation',
      'primary_comparison':'compact4_minus_clinical','primary_endpoint':'four_state_Brier_at3years',
      'flow':flow,'common_domain_n':len(d),'deaths3':int((target>0).sum()),'development_n':spec['development_n'],
      'participant_overlap':0,'model_coefficients_reestimated_on_evaluation':False,'fitting_uncertainty_included':False,
      'bootstrap_replicates':1000,'seed':SEED,'distinct_CAUSES_validated':3,'detailed_respiratory_stroke_validated':False,
      'clinical_use_ready':False,'noninferiority_claim':False,'laboratory_costs_measured':False,
      'method_equivalence_established':False,'new_source_manifest_sha256':sha(source/'manifest.json')}
    dump(out/'RESULTS.json',result);dump(out/'SOURCE_MANIFEST.json',json.loads((source/'manifest.json').read_text()))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    t=s.add_parser('train')
    for f in ('old-source','source08','source09','out'):t.add_argument('--'+f,type=Path,required=True)
    t=s.add_parser('evaluate')
    for f in ('source','trained','out'):t.add_argument('--'+f,type=Path,required=True)
    t.add_argument('--expected-model-sha',required=True)
    a=p.parse_args();kw=vars(a);cmd=kw.pop('command');globals()[cmd](**kw)
