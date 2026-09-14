"""Opt-in standard-library replay of frozen five-cause models without CRP.

The old runtime is used for clinical input validation and unit definitions, not
its coefficients. No dummy CRP is introduced. Both panels require the same full
no-CRP profile to match the comparative evaluation domain.
"""
from copy import deepcopy
import argparse
import hashlib
import json
import math
from pathlib import Path
from .workbench.runtime import (normalize as normalize_clinical, engineer, example as legacy_example,
                                FIELDS, finite, strict_json)

MODEL_SHA='40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094'
MODEL_PATH=Path(__file__).parent/'transport12/model_bundle12.json'
GROUPS=('heart_diseases','malignant_neoplasms','chronic_lower_respiratory','cerebrovascular','other_or_unknown')
REQUIRED=('bmi','sbp','total_cholesterol','hdl','hba1c','creatinine','uacr','albumin','rdw','wbc')
HORIZONS=(1,4,5)


def load_bundle(path=MODEL_PATH):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=MODEL_SHA:raise ValueError('Five-cause model checksum mismatch')
    b=strict_json(raw)
    if b.get('groups')!=list(GROUPS) or b.get('clinical_use_ready') is not False or b.get('CRP_used') is not False:
        raise ValueError('Invalid model scope')
    return b


def example(panel='routine_no_crp'):
    if panel not in ('clinical','routine_no_crp'):raise ValueError('Unknown panel')
    p=legacy_example();p.pop('model');p['measurements'].pop('crp')
    p.update(schema_version='0.12.0',panel=panel,endpoint='five_competing_death_groups',horizons_years=[1,4,5])
    return p


def normalize(payload):
    p=deepcopy(payload) if isinstance(payload,dict) else {}
    errors=[]
    if not isinstance(payload,dict):errors.append({'field':'root','code':'expected_object'})
    for key,value in [('schema_version','0.12.0'),('endpoint','five_competing_death_groups')]:
        if p.get(key)!=value:errors.append({'field':key,'code':'unsupported_definition'})
    panel=p.get('panel')
    if not isinstance(panel,str) or panel not in ('clinical','routine_no_crp'):errors.append({'field':'panel','code':'explicit_panel_required'})
    if 'model' in p:errors.append({'field':'model','code':'legacy_model_not_allowed'})
    h=p.get('horizons_years')
    if (not isinstance(h,list) or not h or len(h)>3 or any(not finite(v) or v not in HORIZONS for v in h)
          or (all(finite(v) for v in h) and len(h)!=len(set(h)))):
        errors.append({'field':'horizons_years','code':'only_unique_1_4_5_years'})
    original=p.get('measurements');m=original if isinstance(original,dict) else {}
    a=deepcopy(p);a.pop('panel',None)
    a.update(schema_version='0.6.0',endpoint='all_cause_death',model='M0c_clinical_only_posthoc',horizons_years=[5])
    a['measurements']={k:m[k] for k in ('bmi','sbp') if k in m} if isinstance(original,dict) else original
    values,trace,clinical_errors,warnings=normalize_clinical(a);errors+=clinical_errors
    for key in set(m)-set(REQUIRED):
        errors.append({'field':'measurements.'+key,'code':'not_in_no_crp_profile'})
    for key in REQUIRED[2:]:
        spec=FIELDS[key];row=m.get(key);path='measurements.'+key
        if not isinstance(row,dict):errors.append({'field':path,'code':'required_measurement_object'});continue
        if set(row)-{'value','unit','matrix','assay'}:errors.append({'field':path,'code':'unknown_measurement_metadata'})
        value=row.get('value');unit=row.get('unit')
        if not finite(value):errors.append({'field':path+'.value','code':'finite_numeric_value_required'});continue
        if not isinstance(unit,str) or unit not in spec['conversions']:
            errors.append({'field':path+'.unit','code':'incompatible_or_missing_unit'});continue
        value=value*spec['conversions'][unit]
        if not finite(value) or not spec['software_limits'][0]<=value<=spec['software_limits'][1]:
            errors.append({'field':path+'.value','code':'outside_software_support'});continue
        if row.get('matrix')!=spec['matrix']:errors.append({'field':path+'.matrix','code':'matrix_mismatch'})
        if spec['required_assay'] and row.get('assay')!=spec['required_assay']:
            errors.append({'field':path+'.assay','code':'standardization_not_confirmed'})
        elif not spec['required_assay']:warnings.append({'field':path+'.assay','code':'assay_transport_not_independently_validated'})
        values[key]=value;trace.append({'field':key,'input_value':row['value'],'input_unit':unit,
            'canonical_value':value,'canonical_unit':spec['unit'],'matrix':row.get('matrix'),'assay':row.get('assay')})
    if values.get('hdl',0)>values.get('total_cholesterol',float('inf')):
        errors.append({'field':'measurements.hdl','code':'HDL_exceeds_total_cholesterol'})
    return values,trace,errors,warnings


def score_engineered(features, model, horizons):
    if (not isinstance(horizons,list) or not horizons or len(horizons)>3
        or any(not finite(v) or v not in HORIZONS for v in horizons) or len(set(horizons))!=len(horizons)):
        raise ValueError('Unsupported horizon list')
    pp=model['preprocessing'];names=pp['features'];z=[];trace=[]
    if any(len(pp[k])!=len(names) for k in ('lower','upper','mean','sd')):raise ValueError('Preprocessing dimensions differ')
    if 'log_crp' in names:raise ValueError('CRP is not a model12 predictor')
    for j,name in enumerate(names):
        x=features.get(name)
        if not finite(x) or not finite(pp['sd'][j]) or pp['sd'][j]<=0:raise ValueError('Invalid feature/scale')
        v=min(max(x,pp['lower'][j]),pp['upper'][j]);standard=(v-pp['mean'][j])/pp['sd'][j]
        if not finite(standard):raise ValueError('Feature transform overflow')
        z.append(standard);trace.append({'feature':name,'raw_value':x,'used_value':v,'standardized_value':standard,'clipped':v!=x})
    if model['interval_ends_years']!=[5.,8.]:raise ValueError('Unexpected interval boundaries')
    rates=[]
    for name in GROUPS:
        c=model['causes'][name]
        if len(c['coefficients'])!=len(z) or len(c['baseline_log_hazards'])!=2:raise ValueError('Invalid coefficient dimensions')
        lp=math.fsum(a*b for a,b in zip(z,c['coefficients']))
        lograte=[lp+a for a in c['baseline_log_hazards']]
        if any(not finite(v) or v>700 for v in lograte):raise ValueError('Hazard overflow')
        rates.append([math.exp(v) for v in lograte])
    result=[]
    for h in horizons:
        s=1.;cif=[0.]*5
        for j,(lo,hi) in enumerate(((0,5),(5,8))):
            dt=max(0.,min(h,hi)-lo);total=math.fsum(row[j] for row in rates)
            if not finite(total):raise ValueError('Summed hazards overflow')
            mass=-math.expm1(-total*dt)
            if total>0:
                for k in range(5):cif[k]+=s*mass*rates[k][j]/total
            s*=math.exp(-total*dt)
        if abs(s+math.fsum(cif)-1)>1e-10:raise ArithmeticError('Probability mass mismatch')
        result.append({'horizon_years':h,'survival':s,'cause_specific_cif':dict(zip(GROUPS,cif)),
            'all_cause_death_probability':math.fsum(cif)})
    return result,trace


def run(payload,*,bundle_path=MODEL_PATH):
    bundle=load_bundle(bundle_path);values,trace,errors,warnings=normalize(payload)
    r={'version':'0.12.0','status':'blocked' if errors else 'calculated_research_only','clinical_use_ready':False,
        'model_sha256':MODEL_SHA,'endpoint':'five_competing_death_groups','probabilities':None,
        'noninfectious_natural_probability':None,'individual_uncertainty_interval':None,'imputation_applied':False,
        'CRP_used':False,'fit_on_request':False,'errors':errors,'warnings':warnings,'measurement_trace':trace,
        'disclaimer':'Исследовательская модель исторической смертности США, 40–79 лет. Пять широких групп первоначальной причины, не точные диагнозы. Остаток сохраняет инфекции и внешние причины. Не для выбора лечения.'}
    if errors:return r
    f=engineer(payload['age_years'],payload['sex'],payload['clinical'],values)
    r['probabilities'],r['feature_trace']=score_engineered(f,bundle['models'][payload['panel']],payload['horizons_years'])
    r['panel']=payload['panel'];r['declaration']='synthetic_not_a_person' if payload['dataset_context']=='synthetic_example' else 'historical_research'
    if payload['panel']=='clinical':r['warnings'].append({'code':'same_complete_laboratory_domain_required_but_not_used_as_predictors'})
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path,nargs='?')
    p.add_argument('--demo',action='store_true');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if bool(a.input)==a.demo:p.error('Provide exactly one input file or --demo')
    if a.out.exists():p.error('Output already exists')
    r=run(example() if a.demo else strict_json(a.input.read_bytes()))
    with a.out.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps({'status':r['status'],'clinical_use_ready':False,'output':str(a.out)}))
    return 2 if r['status']=='blocked' else 0

if __name__=='__main__':raise SystemExit(main())
