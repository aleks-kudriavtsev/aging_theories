"""Locked calibration diagnostics of saved models; never modifies predictions.

Eight quantile bins are an explicitly descriptive display choice. Slopes and
intercepts are evaluation statistics, not newly applied calibration parameters.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from .competing_risks09 import load_test,predict,labels_at,GROUPS,binary_metrics,sha,dump
from .validate_temporal08 import diagnostics
from .replay_competing09 import score_engineered,BUNDLE_SHA256


def main(source:Path,bundle_path:Path,out:Path):
    if out.exists():raise ValueError('Use a new diagnostics directory')
    if sha(bundle_path)!=BUNDLE_SHA256:raise ValueError('Changed model')
    bundle=json.loads(bundle_path.read_text(encoding='utf-8'))
    data,flow=load_test(source,bundle['models']['routine']['preprocessing']['features'])
    rows=[];bins=[];maximum=0.;comparisons=0
    for panel,model in bundle['models'].items():
        all_p={h:predict(data,model,h) for h in [1,5,8]}
        # Verify all 3182 rows, both panels, every cause/survival and horizon with
        # independently implemented stdlib integration of saved coefficients.
        for i,(_,record) in enumerate(data.iterrows()):
            results,_=score_engineered(record.to_dict(),model,[1,5,8])
            for h,r in zip([1,5,8],results):
                actual=np.array([r['survival']]+[r['cause_specific_cif'][g] for g in GROUPS])
                maximum=max(maximum,float(np.max(np.abs(actual-all_p[h][i]))));comparisons+=4
        for horizon in [5,8]:
            p=all_p[horizon];target=labels_at(data,horizon)
            for k,name in [(1,GROUPS[0]),(2,GROUPS[1]),(3,GROUPS[2]),(0,'all_cause')]:
                y=(target==k) if k else target>0;risk=p[:,k] if k else 1-p[:,0]
                rows.append({'panel':panel,'horizon_years':horizon,'outcome':name,
                    **binary_metrics(y,risk,data.weight),**diagnostics(y,risk,data.weight)})
                cuts=np.quantile(risk,np.linspace(0,1,9))
                assignment=np.digitize(risk,cuts[1:-1],right=True)
                for b in range(8):
                    idx=assignment==b
                    if idx.any():bins.append({'panel':panel,'horizon_years':horizon,'outcome':name,
                        'quantile_bin':b+1,**binary_metrics(y[idx],risk[idx],data.weight.to_numpy()[idx])})
    if maximum>1e-12:raise ArithmeticError('Independent runtime replay mismatch')
    out.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(out/'calibration_diagnostics.csv',index=False)
    pd.DataFrame(bins).to_csv(out/'calibration_bins.csv',index=False)
    dump(out/'replay_validation.json',{'number_of_probabilities_compared':comparisons,
        'maximum_absolute_difference':maximum,'model_sha256':BUNDLE_SHA256,
        'all_test_rows':len(data),'feature_coefficients_changed':False,
        'calibration_parameters_applied':False,'quantile_bins':'8; descriptive display choice',
        'clinical_use_ready':False})
    print(pd.DataFrame(rows)[['panel','horizon_years','outcome','calibration_slope','calibration_in_the_large_logit']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','model','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.source,a.model,a.out)
