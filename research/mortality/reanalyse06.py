"""Corrective reanalysis on the ALREADY EXAMINED holdout, not new validation.

No source download; no individual records exported. This corrects XPORT zero
handling without modifying the frozen version-04 model bundle or its reader.
Dependencies are only needed to REBUILD, never to run the browser/CLI workbench.
"""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
from .nhanes_benchmark import (BASE,CLINICAL,PANELS,Preprocessor,fit_piecewise,
                              predict,evaluate,calibration,verify_sources)
from .nhanes_cohort06 import make_cohort_checked
from .xpt_checked import read_xpt_checked


def reanalyse(source: Path,out: Path):
    if out.exists():raise ValueError('Use a NEW reanalysis output directory')
    source_audit=verify_sources(source)
    data,flow=make_cohort_checked(source);train=data[data.cycle<2003];test=data[data.cycle==2003]
    if set(train.SEQN)&set(test.SEQN):raise ValueError('Participant leakage')
    panels={'M0_age_sex':BASE,'M0c_clinical_only_posthoc':BASE+CLINICAL,
            **{k:v for k,v in PANELS.items() if k!='M0_age_sex'}}
    models={};metrics=[];predictions={};subgroups=[];common_stats={}
    for name,features in panels.items():
        pp=Preprocessor().fit(train,features,train.weight.to_numpy())
        model=fit_piecewise(pp.transform(train),train.time_years.to_numpy(),train.dead.to_numpy(),train.weight.to_numpy())
        model.update(columns=pp.columns,features=pp.features,endpoint='all_cause_death',country='US',
            baseline_age_range=[40,79],interval_ends_years=[1,5,10],clinical_use=False)
        models[name]=model;common_stats.update(pp.stats);predictions[name]={}
        for horizon in [1,5,10]:
            p=predict(pp.transform(test),model,horizon);predictions[name][horizon]=p
            y=((test.dead==1)&(test.time_years<=horizon)).to_numpy(int)
            metrics.append({'model':name,'horizon_years':horizon,**evaluate(y,p,test.weight),**calibration(y,p,test.weight)})
            if horizon==10:
                for group,mask in [('female',test.male.eq(0)),('male',test.male.eq(1)),('age40_59',test.RIDAGEYR.lt(60)),('age60_79',test.RIDAGEYR.ge(60))]:
                    subgroups.append({'model':name,'subgroup':group,**evaluate(y[mask],p[mask],test.weight[mask])})
        print(name,metrics[-1],flush=True)
    # Conditional validation uncertainty; the training and model selection are fixed.
    rng=np.random.default_rng(20260910);yw=((test.dead==1)&(test.time_years<=10)).to_numpy(int)
    testw=test.weight.to_numpy();strata=[(sid,sorted(g.SDMVPSU.unique())) for sid,g in test.groupby('SDMVSTRA')]
    rep=[]
    for b in range(300):
        multiplier=np.zeros(len(test))
        for sid,psus in strata:
            m=len(psus)
            if m<2:raise ValueError('Single PSU stratum')
            for psu in rng.choice(psus,m-1,replace=True):
                multiplier[(test.SDMVSTRA.to_numpy()==sid)&(test.SDMVPSU.to_numpy()==psu)]+=m/(m-1)
        es={k:evaluate(yw,v[10],testw*multiplier) for k,v in predictions.items()}
        for extended,base in [('M1_routine','M0c_clinical_only_posthoc'),('M4_cystatinC','M1_routine')]:
            for metric in ['auc','brier']:
                rep.append({'extended':extended,'base':base,'metric':'delta_'+metric,'value':es[extended][metric]-es[base][metric]})
    intervals=[]
    for keys,g in pd.DataFrame(rep).groupby(['extended','base','metric']):
        intervals.append({'extended':keys[0],'base':keys[1],'metric':keys[2],
            'p025':float(g.value.quantile(.025)),'p975':float(g.value.quantile(.975)),
            'replicates':300,'type':'conditional_validation_PSU_percentiles_excludes_training_and_selection'})
    out.mkdir(parents=True)
    origin={'kind':'local_corrective_reanalysis_not_independent_validation','code_base_commit':'fe414dc1957d8c6041dce48e77cf58ee65b30d28',
            'new_remote_commit':None,'reanalysis_date':'2026-09-10','pandas_version':pd.__version__,
            'training_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'reader_sha256':hashlib.sha256(Path(__file__).with_name('xpt_checked.py').read_bytes()).hexdigest()}
    bundle={'version':'0.6_reader_corrected','country':'US','baseline_age_range':[40,79],
        'endpoint':'all_cause_death','clinical_use':False,'interval_ends_years':[1,5,10],
        'preprocessing':common_stats,'models':models,'origin':origin}
    (out/'model_bundle.json').write_text(json.dumps(bundle,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    audit=[]
    for path in sorted(source.glob('*.xpt')):audit.append(read_xpt_checked(path).attrs['xpt_zero_audit'])
    status={'status':'corrective_reanalysis_of_previously_examined_holdout_NOT_independent_validation',
      'source_integrity':source_audit,'origin':origin,'training_n':len(train),'holdout_n':len(test),
      'training_deaths_10y':int(((train.dead==1)&(train.time_years<=10)).sum()),
      'holdout_deaths_10y':int(((test.dead==1)&(test.time_years<=10)).sum()),
      'legacy_training_n':4552,'legacy_holdout_n':2384,
      'legacy_training_deaths_10y':756,'legacy_holdout_deaths_10y':395,
      'trop_all_rows_raw_zero_count':next(a['numeric_zero_counts']['SSTNT'] for a in audit if a['filename']=='1999_SSTROP_A.xpt'),
      'total_repaired_numeric_cells':sum(sum(a['corrected_counts'].values()) for a in audit),
      'files_with_repairs':sum(bool(a['corrected_counts']) for a in audit),
      'frozen_legacy_bundle_unchanged':True,'individual_records_exported':False,'clinical_use_ready':False}
    for fn,content in [('status.json',status),('source_zero_audit.json',audit)]:
        (out/fn).write_text(json.dumps(content,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    for fn,rows in [('metrics.csv',metrics),('subgroups.csv',subgroups),('cohort_flow.csv',flow),('conditional_intervals.csv',intervals)]:
        pd.DataFrame(rows).to_csv(out/fn,index=False)
    return status

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();print(json.dumps(reanalyse(a.sources,a.out),indent=2))
