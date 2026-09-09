"""Post-hoc clinical-only comparator. Does not alter the locked primary models.

Added after inspecting the primary temporal results to distinguish added clinical
history from added routine lab tests. All original and added models are reported;
no test-driven selection or refitting of primary predictions is performed.
"""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from research.mortality.nhanes_benchmark import (BASE, CLINICAL, INTERVAL_ENDS,
    Preprocessor,make_cohort,fit_piecewise,predict,evaluate,calibration,research_predict,verify_sources)

def run(source,out,replicates=300):
    verify_sources(source); data,flow=make_cohort(source)
    tr=data[data.cycle<2003];te=data[data.cycle==2003]
    pp=Preprocessor().fit(tr,BASE+CLINICAL,tr.weight.to_numpy())
    m=fit_piecewise(pp.transform(tr),tr.time_years.to_numpy(),tr.dead.to_numpy(),tr.weight.to_numpy())
    m.update({'features':pp.features,'preprocessing':pp.stats,'columns':pp.columns,
              'endpoint':'all_cause_death','country':'US','baseline_age_range':[40,79],
              'interval_ends_years':INTERVAL_ENDS.tolist(),'clinical_use':False,
              'analysis_status':'post_hoc_sensitivity_not_prespecified_primary'})
    base=json.loads((out/'fitted_models.json').read_text()); y=((te.dead==1)&(te.time_years<=10)).to_numpy(int);w=te.weight.to_numpy()
    p=predict(pp.transform(te),m,10)
    p1=research_predict(te,base['M1_routine'],10,country='US',acknowledge_research_only=True)
    p4=research_predict(te,base['M4_cystatinC'],10,country='US',acknowledge_research_only=True)
    rows=[]
    for name,pr in [('M0c_clinical_only_posthoc',p),('M1_routine',p1),('M4_cystatinC',p4)]:
        rows.append({'model':name,'horizon_years':10,**evaluate(y,pr,w),**calibration(y,pr,w)})
    rng=np.random.default_rng(20260909);res=[]
    strata=[(sid,sorted(g.SDMVPSU.unique()))for sid,g in te.groupby('SDMVSTRA')]
    for b in range(replicates):
        mult=np.zeros(len(te))
        for sid,psus in strata:
            chosen=rng.choice(psus,len(psus)-1,replace=True)
            for psu in chosen:
                mult[(te.SDMVSTRA.to_numpy()==sid)&(te.SDMVPSU.to_numpy()==psu)]+=len(psus)/(len(psus)-1)
        val=[evaluate(y,pr,w*mult)for pr in [p,p1,p4]]
        res.append({'auc_clinical':val[0]['auc'],'delta_auc_routine_minus_clinical':val[1]['auc']-val[0]['auc'],
                    'delta_brier_routine_minus_clinical':val[1]['brier']-val[0]['brier'],
                    'delta_auc_expanded_minus_routine':val[2]['auc']-val[1]['auc']})
    frame=pd.DataFrame(res)
    intervals=[{'metric':k,'p025':float(frame[k].quantile(.025)),'p975':float(frame[k].quantile(.975)),
                'replicates':replicates,'status':'post_hoc_conditional_validation_PSU_resampling'}for k in frame]
    pd.DataFrame(rows).to_csv(out/'clinical_sensitivity.csv',index=False)
    pd.DataFrame(intervals).to_csv(out/'clinical_sensitivity_intervals.csv',index=False)
    (out/'clinical_sensitivity_model.json').write_text(json.dumps(m,indent=2,allow_nan=False))
    print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(intervals).to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();run(args.source,args.out)
