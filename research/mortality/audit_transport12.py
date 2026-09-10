"""Post-computation integrity audit; never changes fitted predictions or endpoints.

Intervals with undefined bootstrap draws are flagged, not upgraded to nominal
coverage by silently deleting undefined replicates. No participant data export.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from .transport_no_crp12 import sha,dump,load_evaluation,predict,GROUPS
from .replay_transport12 import MODEL_SHA,score_engineered


def interval_status(valid,undefined):
    if type(valid) is not int or type(undefined) is not int or valid<0 or undefined<0 or valid+undefined!=1000:
        raise ValueError('Expected complete accounting of1000 bootstrap replicates')
    return 'defined_draws_only_not_nominal_interval' if undefined else 'conditional_on_fixed_models'


def audit(source,trained,evaluated,out):
    if out.exists():raise ValueError('Audit output exists')
    bundle=json.loads((trained/'model_bundle12.json').read_text())
    if sha(trained/'model_bundle12.json')!=MODEL_SHA:raise ValueError('Model changed')
    files=json.loads((source/'manifest.json').read_text())['files'];names=set()
    for row in files:
        n=row['file']
        if Path(n).name!=n or n in names:raise ValueError('Unsafe/duplicate manifest filename')
        names.add(n)
        if row['status']!='retrieved' or sha(source/n)!=row['sha256']:raise ValueError('Source hash mismatch')
    comparisons=0;maxerror=0.;total=0;counts_checks=0;cohort_ids=set(np.load(trained/'LOCAL_ONLY_DEVELOPMENT_IDS.npy'));seen=set()
    metrics=pd.read_csv(evaluated/'metrics.csv');counts=pd.read_csv(evaluated/'cause_counts.csv')
    for year,hs in [(2011,[1,5]),(2013,[1,4])]:
        d,flow=load_evaluation(source,year,MODEL_SHA);total+=len(d);ids=set(d.SEQN)
        if ids&cohort_ids or ids&seen:raise ValueError('Participant overlap')
        seen|=ids
        for panel,model in bundle['models'].items():
            vectors={h:predict(d,model,h) for h in hs}
            for i,record in enumerate(d[model['preprocessing']['features']].to_dict('records')):
                result,_=score_engineered(record,model,hs)
                for h,row in zip(hs,result):
                    scalars=np.array([row['survival']]+[row['cause_specific_cif'][g] for g in GROUPS])
                    maxerror=max(maxerror,float(abs(scalars-vectors[h][i]).max()));comparisons+=6
            for h in hs:
                block=metrics.query('year==@year and horizon==@h and panel==@panel and subgroup=="all"')
                deaths=int(block[block.outcome=='all_cause'].events.iloc[0])
                if counts.query('year==@year and horizon==@h').events.sum()!=deaths:raise ValueError('Cause sum mismatch')
                for prefix in ('','female_','male_'):
                    if prefix=='':groups=['female','male']
                    else:groups=[prefix+'age40_59',prefix+'age60_79']
                    sub=metrics.query('year==@year and horizon==@h and panel==@panel and outcome=="all_cause"')
                    target='all' if not prefix else prefix.rstrip('_')
                    ref=sub[sub.subgroup==target].iloc[0];parts=sub[sub.subgroup.isin(groups)]
                    if int(parts.n.sum())!=int(ref.n) or int(parts.events.sum())!=int(ref.events):raise ValueError('Subgroup sum mismatch')
                    counts_checks+=2
    if maxerror>1e-12:raise ArithmeticError('Scalar/vector replay disagreement')
    ci=pd.read_csv(evaluated/'intervals.csv')
    ci['interpretation_status']=[interval_status(int(v),int(u)) for v,u in zip(ci.valid_replicates,ci.undefined_replicates)]
    out.mkdir(parents=True);ci.to_csv(out/'interval_interpretation.csv',index=False)
    result={'files_verified_including_documentation':len(names),'scientific_data_files':sum(n.endswith(('.dat','.xpt')) for n in names),
      'evaluation_n':total,'probabilities_compared':comparisons,'max_absolute_replay_difference':maxerror,
      'subgroup_arithmetic_checks':counts_checks,'intervals_with_undefined_draws':int(ci.undefined_replicates.gt(0).sum()),
      'guard_added_after_evaluation':'interval reporting only; no coefficient, point metric or primary interval modification',
      'clinical_use_ready':False,'raw_participant_records_exported':False,'model_sha256':MODEL_SHA}
    dump(out/'integrity_audit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for f in ('source','trained','evaluated','out'):p.add_argument('--'+f,type=Path,required=True)
    a=p.parse_args();audit(a.source,a.trained,a.evaluated,a.out)
