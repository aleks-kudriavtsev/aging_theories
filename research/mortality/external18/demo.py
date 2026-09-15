"""Generate unmistakably synthetic local fixtures, not a foreign cohort."""
from copy import deepcopy
import json
from pathlib import Path
import random
from .runner import ARTIFACTS, MEASURES, FIELDS, digest, ATTESTATIONS


def dictionary():
    measurements={}
    for k in MEASURES:
        measurements[k]={'unit':FIELDS[k]['unit'],'matrix':FIELDS[k]['matrix'],
            'method_id':'SYNTHETIC_METHOD', 'source_variable':'SYNTHETIC_'+k,
            'standardization':'NGSP_traceable' if k=='hba1c' else 'IDMS_standardized' if k=='creatinine' else 'not_applicable'}
    return {'measurements':measurements,'clinical_mapping':'Synthetic boolean history only; not a real cohort mapping'}


def rows(n=240):
    rng=random.Random(1801);out=[]
    for i in range(n):
        age=rng.randint(40,79)
        out.append({'record_id':f'SYNTHETIC_{i:04d}','age_years':age,'sex':'male' if i%2 else 'female',
          'clinical':{'smoking':rng.choice(['never','former','current']),'bp_treatment':i%3==0,
            'diabetes_history':i%6==0,'cvd_history':i%8==0,'cancer_history':i%11==0},
          'measurements':{'bmi':rng.uniform(20,34),'sbp':rng.uniform(100,165),'hba1c':rng.uniform(4.7,8.5),
            'creatinine':rng.uniform(.6,1.8),'uacr':rng.uniform(2,70),'albumin':rng.uniform(34,47)},
          'death':i%4==0,'followup_years':rng.uniform(.1,2.9) if i%4==0 else 4.0,
          'weight':1,'stratum':None,'cluster':None})
    out[-1]['death']=False;out[-1]['followup_years']=1.5
    out[-2]['measurements']['uacr']=None
    return out


def generate(out):
    out=Path(out)
    if out.exists():raise ValueError('New directory required')
    out.mkdir(parents=True)
    (out/'protocol.txt').write_text('SYNTHETIC SOFTWARE TEST ONLY. No real participants.\n',encoding='utf-8')
    (out/'dictionary.json').write_text(json.dumps(dictionary(),indent=2),encoding='utf-8')
    (out/'synthetic.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows()),encoding='utf-8')
    config={'schema':'external18.1','study_id':'SYNTHETIC_NOT_A_COHORT','country':'DE','data_origin':'synthetic',
        'endpoint':'all_cause_mortality','horizon_years':3,'time_origin':'baseline_measurements','delayed_entry':False,
        'design':'iid','censoring_assumption':'independent_within_domain','bootstrap_replicates':100,
        'seed':20260913,'minimum_cell_count':10,'protocol_sha256':digest(out/'protocol.txt'),
        'dictionary_sha256':digest(out/'dictionary.json'),'data_sha256':digest(out/'synthetic.jsonl'),
        'model_sha256':ARTIFACTS['compact'][1],'clinical_sha256':ARTIFACTS['full'][1],
        'approval_sha256':None,'attestations':dict.fromkeys(ATTESTATIONS,False)}
    (out/'config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    return config

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();generate(a.out)
