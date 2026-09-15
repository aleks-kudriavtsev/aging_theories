"""Balanced two-run repeated-measure precision, NOT automatic method correction.

Nested ANOVA of one material/concentration level, day > run > replicate.
Variance estimates truncated at zero are reported explicitly. Constant Deming
variance ratio across concentration levels is not assumed or selected here.
"""
from collections import defaultdict
import math
import numpy as np


def components(a):
    a=np.asarray(a,float)
    if a.ndim!=3 or any(n<2 for n in a.shape) or not np.isfinite(a).all() or np.any(a<=0):
        raise ValueError('Balanced positive day/run/replicate array; >=2 per axis required')
    d,r,k=a.shape;grand=float(a.mean());dm=a.mean((1,2));rm=a.mean(2)
    ms_repeat=float(np.sum((a-rm[:,:,None])**2)/(d*r*(k-1)))
    ms_run=float(k*np.sum((rm-dm[:,None])**2)/(d*(r-1)))
    ms_day=float(r*k*np.sum((dm-grand)**2)/(d-1))
    untruncated={'repeatability':ms_repeat,'between_run':(ms_run-ms_repeat)/k,'between_day':(ms_day-ms_run)/(r*k)}
    v={key:max(0,val) for key,val in untruncated.items()};total=sum(v.values())
    return {'days':d,'runs_per_day':r,'replicates_per_run':k,'mean':grand,
      'variance_components':v,'untruncated_variance_components':untruncated,
      'truncated_components':[key for key,val in untruncated.items() if val<0],
      'within_lab_variance':total,'within_lab_sd':math.sqrt(total),
      'within_lab_cv_percent':100*math.sqrt(total)/grand,
      'mathematical_minimum_is_not_a_sufficient_validation_design':True}


def analyse(doc):
    if not isinstance(doc,dict) or doc.get('schema')!='precision18.1':raise ValueError('Schema')
    if doc.get('data_origin') not in ('synthetic','paired_control_measurements'):raise ValueError('Origin')
    for key in ('analyte','unit','matrix','reference_method','candidate_method','reference_lot','candidate_lot'):
        if not isinstance(doc.get(key),str) or not doc[key]:raise ValueError('Metadata missing')
    grouped=defaultdict(dict)
    for row in doc.get('observations',[]):
        if set(row)!={'level','day','run','replicate','reference','candidate'}:raise ValueError('Observation fields')
        key=(row['day'],row['run'],row['replicate'])
        if any(type(x) is not int or x<1 for x in key):raise ValueError('Index')
        if key in grouped[row['level']]:raise ValueError('Duplicate replicate')
        if any(type(row[k]) not in (int,float) or not math.isfinite(row[k]) or row[k]<=0 for k in ('reference','candidate')):raise ValueError('Numeric value')
        grouped[row['level']][key]=row
    if not grouped:raise ValueError('Empty study')
    output=[]
    for level,observations in grouped.items():
        axes=[sorted({key[j] for key in observations}) for j in range(3)]
        if len(observations)!=math.prod(len(x) for x in axes):raise ValueError('Unbalanced nested design')
        data={method:np.array([[[observations[d,r,k][method] for k in axes[2]] for r in axes[1]] for d in axes[0]]) for method in ('reference','candidate')}
        a,b=components(data['reference']),components(data['candidate'])
        ratio=b['within_lab_variance']/a['within_lab_variance'] if a['within_lab_variance']>0 else None
        output.append({'level':str(level),'reference':a,'candidate':b,'candidate_to_reference_variance_ratio':ratio})
    return {'status':'synthetic_precision_check' if doc['data_origin']=='synthetic' else 'precision_for_analytical_review',
        'levels':output,'data_origin':doc['data_origin'],'analyte':doc['analyte'],'unit':doc['unit'],'matrix':doc['matrix'],
        'clinical_use_ready':False,'automatic_Deming_lambda':None,'automatic_correction_applied':False,
        'limitations':['No between-lot component from a single lot pair','No proof of commutability from control pools',
            'No mortality validation','Deming error independence and variance ratio stability require review',
            'No confidence interval for variance components in this implementation']}

if __name__=='__main__':
    import argparse, json
    from pathlib import Path
    from ..workbench.runtime import strict_json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise SystemExit('Output exists')
    result=analyse(strict_json(args.input.read_bytes()))
    with args.out.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps({'status':result['status'],'clinical_use_ready':False}))
