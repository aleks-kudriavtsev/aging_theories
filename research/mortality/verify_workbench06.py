"""Numerical parity against existing historical holdout, NOT new validation."""
from pathlib import Path
import json, argparse
import numpy as np
from research.mortality.nhanes_benchmark import Preprocessor,predict
from research.mortality.nhanes_cohort06 import make_cohort_checked as make_cohort
from research.mortality.workbench.runtime import load_bundle,example,run,MODELS,BUNDLE_SHA256

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    source=args.sources;data,_=make_cohort(source);te=data.loc[data.cycle.eq(2003)]
    bundle=load_bundle();results=[]
    for name in MODELS:
     m=bundle['models'][name];pp=Preprocessor();pp.features=m['features'];pp.stats=bundle['preprocessing'];pp.columns=m['columns']
     complete=te.loc[np.isfinite(te[m['features']]).all(axis=1)].copy()
     old=np.column_stack([predict(pp.transform(complete),m,h) for h in [1,5,10]])
     actual=[]
     for row in complete.to_dict('records'):
      p=example(name);p.update(age_years=float(row['RIDAGEYR']),sex='male' if row['male']==1 else 'female',dataset_context='historical_US_research')
      if name!='M0_age_sex':
       p['clinical']={'smoking':'current' if row['smoker_current']==1 else 'former' if row['smoker_former']==1 else 'never',
         **{k:bool(row[k]) for k in ['bp_treatment','diabetes_history','cvd_history','cancer_history']}}
      vals={'bmi':row['BMXBMI'],'sbp':(row['sbp']*20+130),
        'total_cholesterol':row['LBXTC'],'hdl':row['LBXHDD'],'hba1c':row['LBXGH'],'creatinine':row['creatinine_standardized'],
        'uacr':100*row['URXUMA']/row['URXUCR'],'albumin':row['LBXSAL']*10,'rdw':row['LBXRDW'],
        'wbc':row['LBXWBCSI'],'crp':row['LBXCRP']*10,'ntprobnp':row['SSBNP'],'hs_ctnt':row['SSTNT'],'cystatin_c':row['SSCYST']}
      for key in p['measurements']:p['measurements'][key]['value']=float(vals[key])
      r=run(p)
      if r['status']=='blocked':raise AssertionError((name,r['errors']))
      actual.append([x['probability'] for x in r['probabilities']])
     actual=np.array(actual);delta=np.abs(old-actual);np.testing.assert_allclose(old,actual,rtol=1e-12,atol=1e-13)
     results.append({'model':name,'complete_holdout_profiles':len(complete),'horizons':[1,5,10],
      'probabilities_compared':int(delta.size),'max_absolute_difference':float(delta.max()),'pass':True})
     print(results[-1],flush=True)
    summary={'purpose':'raw_input_runtime_parity_NOT_new_clinical_validation','model_sha256':BUNDLE_SHA256,
     'previously_examined_holdout':'NHANES_2003_2004','holdout_n':len(te),
     'models':results,'total_probabilities_compared':sum(r['probabilities_compared'] for r in results),
     'new_training':False,'recalibration':False,'individual_records_exported':False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(summary,indent=2)+'\n')
