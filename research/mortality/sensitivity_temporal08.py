"""Prespecified complete-case sensitivity plus explicitly post-hoc missingness audit.

No model or calibration factor changes. PSU draws are taken in the full survey
analytic sample before restriction to domains. Conditional intervals exclude
training and calibration estimation uncertainty. No participant outputs.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from .validate_temporal08 import (MODEL_SHA,MODELS,SEED,sha,verify_new_sources,
 new_cohort,frozen_design,cumulative_hazard,design_bootstrap_weights,weighted_auc,write_json)
from .nhanes_benchmark import evaluate


def main(source: Path,model_path: Path,results: Path,out: Path):
    if out.exists():raise ValueError('New sensitivity output directory required')
    if sha(model_path)!=MODEL_SHA:raise ValueError('Changed model')
    verify_new_sources(source)
    d,_,_=new_cohort(source)
    bundle=json.loads(model_path.read_text());spec=json.loads((results/'results.json').read_text())
    if spec['evaluation_n']!=len(d) or spec['frozen_model_sha256']!=MODEL_SHA:raise ValueError('Analysis versions differ')
    full_features=bundle['models']['M1_routine']['features']
    complete=np.isfinite(d[full_features]).all(axis=1).to_numpy()
    y=((d.dead==1)&(d.time_years<=10)).to_numpy(int);w=d.weight.to_numpy(float)
    predictions={};metrics=[];missing_contributions=[]
    for name in MODELS:
        z=frozen_design(d,bundle,name);model=bundle['models'][name]
        h=cumulative_hazard(z,model,10)
        factors=spec['calibration_factors'][name]
        scale=np.where(d.RIAGENDR.eq(2),factors['female'],factors['male'])
        for version,mult in [('frozen',1),('sex_recalibrated',scale)]:
            p=-np.expm1(-h*mult);predictions[name,version]=p
            for label,mask in [('complete_case_prespecified',complete),('missing_profile_posthoc',~complete)]:
                metrics.append({'model':name,'version':version,'domain':label,**evaluate(y[mask],p[mask],w[mask])})
        # Mechanical decomposition, not removal/reweighting or a causal effect.
        for j,column in enumerate(model['columns']):
            if column.endswith('__missing'):
                contribution=z[:,j]*model['coefficients'][j]
                missing_contributions.append({'model':name,'indicator':column,'coefficient':model['coefficients'][j],
                    'activated_n':int(np.count_nonzero(z[:,j])),
                    'analysis_type':'posthoc_mechanical_decomposition_not_new_prediction'})
    draws=[];rng=np.random.default_rng(SEED)
    for rep in range(1000):
        wr=design_bootstrap_weights(d,rng)[complete];yr=y[complete]
        row={}
        for version in ['frozen','sex_recalibrated']:
            for name in MODELS:
                p=predictions[name,version][complete];prefix=f'{name}|{version}'
                row[prefix+'|auc']=weighted_auc(yr,p,wr)
                row[prefix+'|brier']=float(np.average((yr-p)**2,weights=wr))
                row[prefix+'|calibration_bias']=float(np.average(p-yr,weights=wr))
            for metric in ['auc','brier']:
                row[f'delta_M1_M0c|{version}|{metric}']=row[f'M1_routine|{version}|{metric}']-row[f'M0c_clinical_only_posthoc|{version}|{metric}']
        draws.append(row)
    draws=pd.DataFrame(draws);intervals=[]
    for key in draws:
        name,version,metric=key.split('|');lo,hi=np.quantile(draws[key],[.025,.975])
        intervals.append({'model_or_contrast':name,'version':version,'metric':metric,
            'domain':'complete_case_prespecified','lower':float(lo),'upper':float(hi),
            'replicates':1000,'interval_kind':'conditional_PSU_percentiles'})
    out.mkdir(parents=True)
    pd.DataFrame(metrics).to_csv(out/'domain_metrics.csv',index=False)
    pd.DataFrame(intervals).to_csv(out/'complete_case_intervals.csv',index=False)
    pd.DataFrame(missing_contributions).to_csv(out/'missing_indicator_audit.csv',index=False)
    write_json(out/'sensitivity_status.json',{'complete_case_analysis':'prespecified',
        'missing_profile_metrics_and_coefficient_audit':'post_hoc',
        'frozen_model_sha256':MODEL_SHA,'primary_results_sha256':sha(results/'results.json'),
        'no_model_changes':True,'no_new_validation_sample':True,'replicates':1000,
        'conditional_on_fixed_training_and_calibration':True,'raw_records_exported':False})
    print(pd.DataFrame(metrics).to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','model','results','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.source,a.model,a.results,a.out)
