"""Iteration14 amended IPCW evaluation; original exact training script unchanged.

Amendment b105c62 preceded any comparative performance calculation. No model
selection, coefficients or calibration are fitted on this evaluation dataset.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from .compact_panel14 import (LOCK,NAMES,OUTCOMES,HORIZONS,SEED,ROUTINE,COMPONENTS,
    ROOT,comparator_models,sha,dump,predict,aggregate_states,engineer_no_crp,
    read_xpt_checked,independent_mortality_reader,audit_mortality_reader,
    design_bootstrap_weights,diagnostics)
from .ipcw14 import outcome_weights,brier_states,binary_scores,aalen_johansen_reference
AMENDMENT='b105c62e23337be742c28bfb6c8081e71ad9a4ff'

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
    outcome_weights(d.dead,d.cause,d.time_years,d.weight,3) # Amendment01: retain early censors.
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
    common=np.isfinite(eligible[ROUTINE]).all(axis=1).to_numpy();d=eligible.loc[common].copy();w=d.weight.to_numpy(float)
    if not len(d):raise ValueError('Empty common domain')
    coverage=[]
    for name in NAMES:
        good=np.isfinite(eligible[models[name]['preprocessing']['features']]).all(axis=1).to_numpy()
        coverage.append({'panel':name,'eligible_n':len(eligible),'complete_n':int(good.sum()),
            'additional_complete_vs_common':int(good.sum()-common.sum()),
            'weighted_coverage':float(np.average(good,weights=eligible.weight)),
            'performance_compared_on_common_domain_only':True})
    metrics=[];bins=[];predictions={};counts=[];quality=[];aj_checks=[]
    masks={'all':np.ones(len(d),bool),'female':d.RIAGENDR.eq(2).to_numpy(),'male':d.RIAGENDR.eq(1).to_numpy(),
           'age40_59':d.RIDAGEYR.lt(60).to_numpy(),'age60_79':d.RIDAGEYR.ge(60).to_numpy()}
    def obs(mask,ww,h):return outcome_weights(d.dead.to_numpy()[mask],d.cause.to_numpy()[mask],d.time_years.to_numpy()[mask],ww,h)
    for h in HORIZONS:
        o=obs(masks['all'],w,h);target=o['target']
        reference=aalen_johansen_reference(d.dead,d.cause,d.time_years,w,h)
        ipcw_prob=np.bincount(target,weights=o['effective_weights'],minlength=4)/w.sum()
        difference=float(np.max(abs(reference-ipcw_prob)))
        if difference>1e-10:raise ArithmeticError('Independent Aalen-Johansen check failed')
        aj_checks.append({'horizon':h,'max_abs_difference':difference,'sum_probability':float(reference.sum())})
        for k,name in enumerate(('survival',)+OUTCOMES):
            counts.append({'horizon':h,'outcome':name,'n':len(d),'known_status_count':int(((target==k)&o['known']).sum()),
                 'early_censored_n':o['early_censored_n'],'weighted_cumulative_incidence':float(reference[k])})
        for name in NAMES:
            p=aggregate_states(predict(d,models[name],h));predictions[name,h]=p
            for group,mask in masks.items():
                og=obs(mask,w[mask],h);yg=og['target']
                metrics.append({'panel':name,'horizon':h,'subgroup':group,'outcome':'four_state','n':int(mask.sum()),
                   'events':int((yg>0).sum()),'early_censored_n':og['early_censored_n'],
                   'brier':brier_states(yg,p[mask],og['effective_weights'],w[mask].sum())})
                for k,label in list(enumerate(OUTCOMES,1))+[(0,'all_cause')]:
                    risk=p[:,k] if k else 1-p[:,0];met=binary_scores(yg,risk[mask],w[mask],og,k)
                    if group=='all' and h==3:
                        try:
                            valid=og['effective_weights']>0;y=yg==k if k else yg>0
                            met.update(diagnostics(y[valid],risk[valid],og['effective_weights'][valid]));met['calibration_status']='computed'
                        except ValueError:met['calibration_status']='not_converged_no_refit'
                    metrics.append({'panel':name,'horizon':h,'subgroup':group,'outcome':label,**met})
                    if group=='all' and h==3:
                        assign=np.digitize(risk,np.quantile(risk,np.linspace(0,1,6))[1:-1],right=True)
                        for b in range(5):
                            m=assign==b
                            if m.any():
                                ob=obs(m,w[m],h)
                                bins.append({'panel':name,'outcome':label,'bin':b+1,**binary_scores(ob['target'],risk[m],w[m],ob,k)})
            for j,f in enumerate(models[name]['preprocessing']['features']):
                pp=models[name]['preprocessing'];quality.append({'panel':name,'horizon':h,'feature':f,
                   'clipped_n':int(((d[f]<pp['lower'][j])|(d[f]>pp['upper'][j])).sum())})
    print('PRIMARY_COMPARISON_NOW_CALCULATED',flush=True)
    print(pd.DataFrame(metrics).query('horizon==3 and subgroup=="all"')[['panel','outcome','n','events','brier','auc','predicted','observed']].to_string(index=False),flush=True)
    rng=np.random.default_rng(SEED);draws=[];failed_draws=0
    comparisons=(('compact4','clinical'),('full8','clinical'),('compact4','full8'),('compact4','six6'))
    for _ in range(1000):
        wr=design_bootstrap_weights(eligible,rng)[common];row={}
        try:ob=obs(masks['all'],wr,3)
        except ValueError:failed_draws+=1;draws.append(row);continue
        target=ob['target']
        for name in NAMES:
            p=predictions[name,3];row[name+'|four_state|brier']=brier_states(target,p,ob['effective_weights'],wr.sum())
            for k,label in list(enumerate(OUTCOMES,1))+[(0,'all_cause')]:
                r=binary_scores(target,p[:,k] if k else 1-p[:,0],wr,ob,k)
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
          'interpretation':'conditional_on_fixed_models_censoring_refitted' if len(x)==1000 else 'defined_draws_only_not_nominal_coverage'})
    out.mkdir(parents=True)
    for name,data in [('metrics',metrics),('intervals',intervals),('coverage',coverage),('calibration_bins',bins),('counts',counts),('feature_quality',quality)]:
        pd.DataFrame(data).to_csv(out/(name+'.csv'),index=False)
    urine=read_xpt_checked(source/'2015_ALB_CR_I.xpt')
    flags=urine.set_index('SEQN').reindex(d.SEQN)
    lod={col:int(flags[col].eq(1).sum()) for col in ('URDUMALC','URDUCRLC') if col in flags}
    result={'analysis_lock_commit':LOCK,'amendment_commit':AMENDMENT,'model_sha256':expected_model_sha,'status':'new_temporal_validation',
      'primary_comparison':'compact4_minus_clinical','primary_endpoint':'four_state_IPCW_Brier_at3years',
      'flow':flow,'common_domain_n':len(d),'deaths3':int((obs(masks['all'],w,3)['target']>0).sum()),
      'early_censored_common3':obs(masks['all'],w,3)['early_censored_n'],'development_n':spec['development_n'],
      'G_common_at3_left':obs(masks['all'],w,3)['G_horizon_left'],
      'participant_overlap':0,'model_coefficients_reestimated_on_evaluation':False,'fitting_uncertainty_included':False,
      'bootstrap_replicates':1000,'bootstrap_censoring_failures':failed_draws,'seed':SEED,'distinct_causes_validated':3,
      'independent_Aalen_Johansen_checks':aj_checks,'source_LOD_substitutions_in_common_domain':lod,
      'detailed_respiratory_stroke_validated':False,'clinical_use_ready':False,'noninferiority_claim':False,
      'laboratory_costs_measured':False,'method_equivalence_established':False,
      'new_source_manifest_sha256':sha(source/'manifest.json')}
    dump(out/'RESULTS.json',result);dump(out/'SOURCE_MANIFEST.json',json.loads((source/'manifest.json').read_text()))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('source','trained','out'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--expected-model-sha',required=True);a=p.parse_args();evaluate(**vars(a))
