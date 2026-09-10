"""Secondary CBC-omission sensitivity, metadata-triggered before evaluation outcomes.

The primary no-CRP model stays frozen. Same development/evaluation participants;
no population-mean assay correction, no evaluation refit, no winner selection.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json
import numpy as np
import pandas as pd
from .transport_no_crp12 import (ROUTINE,FOLDS,PENALTIES,GROUPS,LOCK,SEED,sha,dump,
 cohort,fit_panel,predict,labels,load_evaluation,design_bootstrap_weights,
 multiclass_brier,binary_metrics,diagnostics)
PRIMARY_SHA='40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094'
AMENDMENT='d52e6b472a6c069b52a1cee23185641fae897afd'
FEATURES=[f for f in ROUTINE if f not in {'rdw','log_wbc'}]


def train(old_source,source08,source09,primary,out):
    if out.exists():raise ValueError('Output already exists')
    if sha(primary/'model_bundle12.json')!=PRIMARY_SHA:raise ValueError('Primary model changed')
    d,_,used,_=cohort({y:old_source if y<2005 else source08 if y<2009 else source09 for y in (1999,2001,2003,2005,2007,2009)})
    expected=json.loads((primary/'development_sources.json').read_text())
    if any(used[r['file']]['sha256']!=r['sha256'] for r in expected):raise ValueError('Development source changed')
    d=d[np.isfinite(d[ROUTINE]).all(axis=1)].copy();d['time_years']=d.t_exam.clip(lower=.5)/12;d['weight']=d.pooled_weight
    if set(d.SEQN)!=set(np.load(primary/'LOCAL_ONLY_DEVELOPMENT_IDS.npy')):raise ValueError('Different development participants')
    scores={};rows=[]
    for penalty in PENALTIES:
        vals=[]
        for block in FOLDS:
            mask=d.cycle.isin(block);m=fit_panel(d[~mask],FEATURES,penalty)
            v=multiclass_brier(labels(d[mask],5),predict(d[mask],m,5),d.loc[mask,'weight']);vals.append(v)
            rows.append({'penalty':penalty,'heldout':str(block),'brier':v})
        scores[penalty]=float(np.mean(vals));print('SENSITIVITY_CV',penalty,scores[penalty],flush=True)
    selected=min(PENALTIES,key=lambda p:(scores[p],-p));m=fit_panel(d,FEATURES,selected)
    out.mkdir(parents=True)
    dump(out/'model_sensitivity12.json',{'version':'0.12_metadata_sensitivity','amendment_commit':AMENDMENT,
      'primary_model_sha256':PRIMARY_SHA,'groups':list(GROUPS),'model':m,'development_n':len(d),
      'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'evaluation_outcomes_seen':False,
      'clinical_use_ready':False,'source_code_sha256':sha(__file__),'same_no_crp_complete_profile':True})
    pd.DataFrame(rows).to_csv(out/'sensitivity_cv.csv',index=False)
    print('SENSITIVITY_MODEL_SHA',sha(out/'model_sensitivity12.json'),flush=True)


def evaluate(source,trained,primary,out,expected_sha):
    if out.exists():raise ValueError('Output already exists')
    path=trained/'model_sensitivity12.json'
    if sha(path)!=expected_sha or sha(primary/'model_bundle12.json')!=PRIMARY_SHA:raise ValueError('Frozen model changed')
    spec=json.loads(path.read_text());pb=json.loads((primary/'model_bundle12.json').read_text())
    if spec['amendment_commit']!=AMENDMENT or spec['primary_model_sha256']!=PRIMARY_SHA:raise ValueError('Inconsistent amendment')
    models={**pb['models'],'no_crp_no_cbc_sensitivity':spec['model']};metrics=[];intervals=[]
    oldids=set(np.load(primary/'LOCAL_ONLY_DEVELOPMENT_IDS.npy'));seen=set()
    for year,h in [(2011,5),(2013,4)]:
        d,_=load_evaluation(source,year,PRIMARY_SHA)
        if set(d.SEQN)&(oldids|seen):raise ValueError('Participant overlap')
        seen|=set(d.SEQN);target=labels(d,h);w=d.weight.to_numpy();preds={n:predict(d,m,h) for n,m in models.items()}
        def get_metrics(weights):
            vals={}
            for name,p in preds.items():
                vals[name+'|six_state|brier']=multiclass_brier(target,p,weights)
                for k,outcome in list(enumerate(GROUPS,1))+[(0,'all_cause')]:
                    y=target==k if k else target>0;risk=p[:,k] if k else 1-p[:,0]
                    met=binary_metrics(y,risk,weights)
                    for metric in ('auc','brier','bias_pp'):vals[name+'|'+outcome+'|'+metric]=met[metric]
            return vals
        for name,p in preds.items():
            metrics.append({'year':year,'horizon':h,'panel':name,'outcome':'six_state','n':len(d),'events':int((target>0).sum()),'brier':multiclass_brier(target,p,w)})
            for k,outcome in list(enumerate(GROUPS,1))+[(0,'all_cause')]:
                y=target==k if k else target>0;risk=p[:,k] if k else 1-p[:,0];met=binary_metrics(y,risk,w)
                try:met.update(diagnostics(y,risk,w))
                except ValueError:met['calibration_diagnostic_status']='not_converged'
                metrics.append({'year':year,'horizon':h,'panel':name,'outcome':outcome,**met})
        draws=[];rng=np.random.default_rng(SEED+year-2011)
        for rep in range(1000):
            vals=get_metrics(design_bootstrap_weights(d,rng))
            for base in ('clinical','routine_no_crp'):
                for outcome in ('six_state',*GROUPS,'all_cause'):
                    for metric in (['brier'] if outcome=='six_state' else ['brier','auc','bias_pp']):
                        a,b=vals['no_crp_no_cbc_sensitivity|'+outcome+'|'+metric],vals[base+'|'+outcome+'|'+metric]
                        vals['sensitivity_minus_'+base+'|'+outcome+'|'+metric]=a-b if a is not None and b is not None else None
            draws.append(vals)
        frame=pd.DataFrame(draws)
        for key in frame:
            name,outcome,metric=key.split('|');v=frame[key].dropna();lo,hi=np.quantile(v,[.025,.975]) if len(v) else (None,None)
            intervals.append({'year':year,'horizon':h,'panel_or_contrast':name,'outcome':outcome,'metric':metric,
             'lower':lo,'upper':hi,'valid_replicates':len(v),'undefined_replicates':1000-len(v),'interval_kind':'conditional_PSU_percentiles'})
    out.mkdir(parents=True);pd.DataFrame(metrics).to_csv(out/'sensitivity_metrics.csv',index=False);pd.DataFrame(intervals).to_csv(out/'sensitivity_intervals.csv',index=False)
    dump(out/'sensitivity_status.json',{'amendment_commit':AMENDMENT,'model_sha256':expected_sha,'primary_unchanged':True,
      'selection_based_on_evaluation':False,'scope':'metadata_triggered_pre_outcome_secondary_analysis','clinical_use_ready':False})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('train')
    for f in ('old-source','source08','source09','primary','out'):a.add_argument('--'+f,type=Path,required=True)
    b=s.add_parser('evaluate')
    for f in ('source','trained','primary','out'):b.add_argument('--'+f,type=Path,required=True)
    b.add_argument('--expected-sha',required=True);x=p.parse_args();kw=vars(x);cmd=kw.pop('command');globals()[cmd](**kw)
