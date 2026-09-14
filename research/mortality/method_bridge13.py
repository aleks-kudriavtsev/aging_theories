"""Paired quantitative method comparison, research only; standard library.

Candidate = intercept + slope * reference, constant-variance Deming regression.
The caller MUST specify var(candidate error)/var(reference error), in the same
units as paired measurements. Independent specimen IDs cannot cross fit/validation.
Numerical comparability is NOT interchangeability or mortality-model validation.
No actual laboratory pairs are bundled; --demo produces labelled synthetic data.
"""
from __future__ import annotations
import argparse
import csv
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics
from .workbench.runtime import strict_json


def numeric(v):
    if type(v) not in (int, float): return False
    try: return math.isfinite(v)
    except OverflowError: return False


def positive_vector(values):
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError('Nonempty numerical list required')
    if any(not numeric(v) or v <= 0 or v > 1e100 for v in values):
        raise ValueError('Finite positive measurements within numerical support required')
    return [float(v) for v in values]


def paired(x, y, minimum=3):
    x, y = positive_vector(x), positive_vector(y)
    if len(x) != len(y) or len(x) < minimum:
        raise ValueError('Equal paired arrays and adequate arithmetic support required')
    return x, y


def deming(x, y, *, variance_ratio):
    """Positive affine relationship, known homoscedastic error-variance ratio."""
    x, y = paired(x, y)
    if not numeric(variance_ratio) or not 1e-12 <= variance_ratio <= 1e12:
        raise ValueError('Explicit positive variance ratio var(y error)/var(x error) required')
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sx = math.fsum((v-mx)**2 for v in x)/len(x)
    sy = math.fsum((v-my)**2 for v in y)/len(y)
    sxy = math.fsum((a-mx)*(b-my) for a,b in zip(x,y))/len(x)
    if sx <= 0 or sy <= 0 or sxy <= 0:
        raise ValueError('Insufficient variation or nonpositive association for analyte bridge')
    lam = float(variance_ratio); diff = sy-lam*sx
    radical = math.hypot(diff, 2*math.sqrt(lam)*sxy)
    slope = (diff+radical)/(2*sxy) if diff >= 0 else 2*lam*sxy/(radical-diff)
    intercept = my-slope*mx
    if not all(numeric(v) for v in (slope,intercept)) or slope <= 0:
        raise ValueError('Unstable bridge fit')
    return {'intercept':intercept,'slope':slope,'variance_ratio':lam,
       'equation':'candidate = intercept + slope * reference',
       'reference_range':[min(x),max(x)],'candidate_range':[min(y),max(y)],'fit_n':len(x)}


def correct(value, model):
    if not numeric(value) or value <= 0: raise ValueError('Positive measurement required')
    if not isinstance(model,dict): raise ValueError('Model mapping required')
    slope,intercept=model.get('slope'),model.get('intercept');bounds=model.get('candidate_range')
    if (not numeric(slope) or slope<=0 or not numeric(intercept) or not isinstance(bounds,list)
        or len(bounds)!=2 or any(not numeric(v) for v in bounds) or not 0<bounds[0]<=bounds[1]):
        raise ValueError('Invalid bridge specification')
    if not bounds[0]<=value<=bounds[1]: raise ValueError('Outside fitted candidate range; no extrapolation')
    result=(value-intercept)/slope
    if not numeric(result) or result<=0: raise ValueError('Correction outside positive measurement domain')
    return result


def quantile(values,p):
    if not values or not numeric(p) or not 0<=p<=1:raise ValueError('Invalid quantile')
    v=sorted(values);at=(len(v)-1)*p;lo=math.floor(at);hi=math.ceil(at)
    return v[lo]+(v[hi]-v[lo])*(at-lo)


def agreement(reference,candidate):
    x,y=paired(reference,candidate);dif=[b-a for a,b in zip(x,y)]
    mean=statistics.fmean(dif);sd=statistics.stdev(dif)
    mx,my=statistics.fmean(x),statistics.fmean(y)
    # Separate square roots prevent overflow of a product of valid variances.
    denominator=math.sqrt(math.fsum((v-mx)**2 for v in x))*math.sqrt(math.fsum((v-my)**2 for v in y))
    corr=math.fsum((a-mx)*(b-my) for a,b in zip(x,y))/denominator if denominator else None
    return {'n':len(x),'mean_difference':mean,'sd_difference':sd,
        'mean_absolute_difference':statistics.fmean(abs(v) for v in dif),
        'p95_absolute_difference':quantile([abs(v) for v in dif],.95),
        'approximate_limits_of_agreement':[mean-1.96*sd,mean+1.96*sd],
        'limits_definition':'mean +/- 1.96 SD of paired differences; NOT confidence limits',
        'limits_assumptions':'independent specimens; approximately normal constant-variance differences',
        'pearson_correlation':corr,'correlation_establishes_interchangeability':False}


def validate_document(doc):
    if not isinstance(doc,dict):raise ValueError('Document mapping required')
    expected={'schema_version','data_origin','study_id','model_sha256','reference_method','candidate_method',
              'variance_ratio','variance_ratio_basis','pairs'}
    if set(doc)!=expected or doc['schema_version']!='0.13_paired_methods':
        raise ValueError('Unsupported or extra top-level fields')
    if doc['data_origin'] not in ('synthetic','paired_specimens_declared'):
        raise ValueError('Explicit origin declaration required')
    if not isinstance(doc['study_id'],str) or not re.fullmatch(r'[A-Za-z0-9_.:-]{1,80}',doc['study_id']):
        raise ValueError('Non-identifying study code required')
    if not isinstance(doc['model_sha256'],str) or not re.fullmatch(r'[0-9a-f]{64}',doc['model_sha256']):
        raise ValueError('Explicit target risk-model hash required')
    for label in ('reference_method','candidate_method'):
        method=doc[label]
        if not isinstance(method,dict) or set(method)!={'analyte','unit','matrix','method_id','lot_id'}:
            raise ValueError('Exact analyte, unit, matrix, method and lot required')
        if any(not isinstance(v,str) or not v.strip() or v.lower() in ('unknown','not_reported') for v in method.values()):
            raise ValueError('Unresolved method metadata')
    for key in ('analyte','unit','matrix'):
        if doc['reference_method'][key]!=doc['candidate_method'][key]:
            raise ValueError('Paired methods must have identical '+key+'; no automatic unit or matrix conversion')
    if doc['reference_method']['method_id']==doc['candidate_method']['method_id']:
        raise ValueError('Declare distinct method/version IDs for a comparison')
    if not isinstance(doc['variance_ratio_basis'],str) or not doc['variance_ratio_basis'].strip():
        raise ValueError('Precision study supporting error-variance ratio must be identified')
    records=doc['pairs']
    if not numeric(doc['variance_ratio']) or not 1e-12<=doc['variance_ratio']<=1e12:
        raise ValueError('Invalid declared error-variance ratio')
    if not isinstance(records,list):raise ValueError('Pair list required')
    splits={'fit':[],'validation':[]};seen=set()
    for row in records:
        if not isinstance(row,dict) or set(row)!={'specimen_id','split','reference','candidate'}:
            raise ValueError('Exactly one aggregate pair per specimen; raw identifiers/extra fields prohibited')
        sid=row['specimen_id']
        if not isinstance(sid,str) or not re.fullmatch(r'[A-Za-z0-9_.:-]{1,80}',sid) or sid in seen:
            raise ValueError('Unique non-identifying specimen IDs; no split leakage or duplicate weighting')
        seen.add(sid)
        if not isinstance(row['split'],str) or row['split'] not in splits:raise ValueError('Explicit fit/validation split required')
        positive_vector([row['reference'],row['candidate']]);splits[row['split']].append(row)
    if min(map(len,splits.values()))<5:
        raise ValueError('Each split needs >=5 for software operation, NOT a validation sample-size recommendation')
    return splits


def analyse(doc, *, bootstrap_replicates=1000, seed=20260913):
    splits=validate_document(doc)
    if type(bootstrap_replicates)!=int or not 100<=bootstrap_replicates<=10000:
        raise ValueError('Declare 100-10000 bootstrap replicates')
    x=[r['reference'] for r in splits['fit']];y=[r['candidate'] for r in splits['fit']]
    xv=[r['reference'] for r in splits['validation']];yv=[r['candidate'] for r in splits['validation']]
    model=deming(x,y,variance_ratio=doc['variance_ratio']);raw=agreement(xv,yv)
    outside=sum(not(model['reference_range'][0]<=a<=model['reference_range'][1] and
                    model['candidate_range'][0]<=b<=model['candidate_range'][1]) for a,b in zip(xv,yv))
    corrected=None;correction_errors=0
    if not outside:
        try: corrected=agreement(xv,[correct(v,model) for v in yv])
        except ValueError:correction_errors=1
    rng=random.Random(seed);estimates=[];undefined=0
    for _ in range(bootstrap_replicates):
        idx=[rng.randrange(len(x)) for _ in x]
        try:
            m=deming([x[i] for i in idx],[y[i] for i in idx],variance_ratio=doc['variance_ratio'])
            estimates.append((m['intercept'],m['slope']))
        except ValueError:undefined+=1
    interval={}
    for j,key in enumerate(('intercept','slope')):
        v=[e[j] for e in estimates]
        interval[key]={'lower':quantile(v,.025) if v else None,'upper':quantile(v,.975) if v else None,
            'valid_replicates':len(v),'undefined_replicates':undefined,
            'interval_kind':'paired_specimen_percentiles_conditional_on_variance_ratio',
            'nominal_coverage_claimed':False}
    serial=json.dumps(doc,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()
    return {'schema_version':'0.13_paired_methods_result','data_origin':doc['data_origin'],
      'study_id':doc['study_id'],'input_sha256':hashlib.sha256(serial).hexdigest(),
      'target_model_sha256':doc['model_sha256'],'reference_method':deepcopy(doc['reference_method']),
      'candidate_method':deepcopy(doc['candidate_method']),'variance_ratio_basis':doc['variance_ratio_basis'],
      'fit':model,'fit_parameter_intervals':interval,'validation_raw':raw,'validation_corrected':corrected,
      'validation_outside_fitted_range':outside,'validation_correction_errors':correction_errors,
      'calibration_validation_overlap':0,'bootstrap_replicates':bootstrap_replicates,'seed':seed,
      'source_independently_verified':False,'clinical_use_ready':False,
      'interchangeability_established':False,'mortality_transport_validated':False,
      'automatic_risk_model_activation':False,'specimen_rows_exported':False,
      'uncertainty_exclusions':'variance-ratio estimation, between-lot/batch effects and risk-model uncertainty',
      'status':'synthetic_method_demo' if doc['data_origin']=='synthetic' else 'paired_method_analysis_for_review'}


def application_gate(report, *, analyte, unit, matrix, method_id, lot_id, model_sha256):
    """Reject automatic use; report is analytic evidence, never an activation token."""
    reasons=[]
    if not isinstance(report,dict):return {'ready':False,'reasons':['invalid_report']}
    method=report.get('candidate_method',{})
    for key,value in [('analyte',analyte),('unit',unit),('matrix',matrix),('method_id',method_id),('lot_id',lot_id)]:
        if not isinstance(method,dict) or method.get(key)!=value:reasons.append('different_'+key)
    if report.get('target_model_sha256')!=model_sha256:reasons.append('different_risk_model')
    if report.get('data_origin')!='paired_specimens_declared':reasons.append('synthetic_or_undeclared_data')
    reasons+=['independent_assay_review_required','separate_risk_transport_validation_required']
    return {'ready':False,'reasons':reasons,'clinical_use_ready':False}


def demo():
    from .replay_transport12 import MODEL_SHA
    rng=random.Random(1301)
    pairs=[]
    for split,values in [('fit',[10+i*.5 for i in range(61)]),('validation',[11+i*.6 for i in range(46)])]:
        for j,x in enumerate(values):
            pairs.append({'specimen_id':f'SYNTHETIC_{split}_{j:03d}','split':split,
                          'reference':x,'candidate':1.5+1.08*x+rng.gauss(0,.25)})
    common={'analyte':'synthetic_analyte','unit':'arbitrary_unit','matrix':'synthetic_matrix'}
    return {'schema_version':'0.13_paired_methods','data_origin':'synthetic','study_id':'SYNTHETIC_DEMO',
       'model_sha256':MODEL_SHA,'reference_method':{**common,'method_id':'synthetic_reference','lot_id':'DEMO_ONLY'},
       'candidate_method':{**common,'method_id':'synthetic_candidate','lot_id':'DEMO_ONLY'},
       'variance_ratio':1.,'variance_ratio_basis':'synthetic assumption, not an empirical precision study',
       'pairs':pairs}


def load_pairs_csv(path):
    """Strict four-column input; no defaulting of empty, censored or repeated rows."""
    expected=['specimen_id','split','reference','candidate']
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=expected:raise ValueError('CSV headers must be '+','.join(expected))
        rows=[]
        for row in reader:
            if set(row)!=set(expected) or any(v is None for v in row.values()):
                raise ValueError('Malformed CSV row')
            try:
                row['reference']=float(row['reference']);row['candidate']=float(row['candidate'])
            except (TypeError,ValueError):raise ValueError('CSV requires quantified finite numerical measurements') from None
            positive_vector([row['reference'],row['candidate']]);rows.append(row)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path,nargs='?')
    p.add_argument('--demo',action='store_true');p.add_argument('--pairs',type=Path)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if bool(a.input)==a.demo:p.error('Exactly one file or --demo required')
    doc=demo() if a.demo else strict_json(a.input.read_bytes())
    if a.pairs:
        if a.demo or not isinstance(doc,dict) or doc.get('pairs')!=[]:
            p.error('--pairs requires a metadata JSON document with an empty pairs list')
        doc['pairs']=load_pairs_csv(a.pairs)
    result=analyse(doc)
    with a.out.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False)
    print(json.dumps({'status':result['status'],'clinical_use_ready':False}))

if __name__=='__main__':main()
