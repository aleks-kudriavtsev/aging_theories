"""Unified, fail-closed adapters. Frozen coefficients; no fitting on requests."""
from copy import deepcopy
from hashlib import sha256
import json
import math
from .registry import REGISTRY, FULL, COARSE, COMMON_WARNING, ARTIFACTS, model_spec, evidence
from ..workbench.runtime import FIELDS, finite, normalize as old_normalize, engineer, example as old_example
from ..replay_compact14 import score_engineered as score_coarse
from ..replay_transport12 import score_engineered as score_detailed

VERSION='0.16.0'
KEYS={'schema_version','model_id','research_only','acknowledge_limitations','dataset_context',
      'country','endpoint','horizons_years','age_years','sex','clinical','measurements'}


def example(model_id):
    if model_id not in REGISTRY: raise ValueError('Unknown model_id')
    e=REGISTRY[model_id];p=old_example()
    p.pop('model');p.update(schema_version=VERSION,model_id=model_id,
        endpoint=e['endpoint'],horizons_years=deepcopy(e['horizons']))
    p['measurements']={k:v for k,v in p['measurements'].items() if k in e['required']}
    return p


def normalize(payload):
    errors=[];warnings=[];trace=[];values={}
    def err(field,code):errors.append({'field':field,'code':code})
    if not isinstance(payload,dict):return None,values,trace,[{'field':'root','code':'expected_object'}],warnings
    for k in set(payload)-KEYS:err(str(k),'unknown_field_identifiers_not_allowed')
    if payload.get('schema_version')!=VERSION:err('schema_version','expected_0.16.0')
    mid=payload.get('model_id');entry=REGISTRY.get(mid) if isinstance(mid,str) else None
    if entry is None:err('model_id','explicit_registered_model_required');return None,values,trace,errors,warnings
    if payload.get('endpoint')!=entry['endpoint']:err('endpoint','endpoint_mismatch')
    horizons=payload.get('horizons_years')
    if (not isinstance(horizons,list) or not horizons or len(horizons)>3
        or any(not finite(h) or h not in entry['horizons'] for h in horizons)
        or (all(finite(h) for h in horizons) and len(set(horizons))!=len(horizons))):
        err('horizons_years','unsupported_or_duplicate_horizon')
    measurements=payload.get('measurements')
    m=measurements if isinstance(measurements,dict) else {}
    # Reuse the old clinical/anthropometry validator ONLY, never its predictions.
    adapted={k:deepcopy(v) for k,v in payload.items() if k in KEYS-{'model_id'}}
    adapted.update(schema_version='0.6.0',model='M0c_clinical_only_posthoc',
                   endpoint='all_cause_death',horizons_years=[1])
    adapted['measurements']={k:m[k] for k in ['bmi','sbp'] if k in m} if isinstance(measurements,dict) else measurements
    values,trace,ce,cw=old_normalize(adapted);errors+=ce;warnings+=cw
    for k in set(m)-set(entry['required']):err('measurements.'+str(k),'unused_or_unknown_measurement')
    for k in entry['required']:
        if k in ('bmi','sbp'):continue
        row=m.get(k);field='measurements.'+k;sp=FIELDS[k]
        if not isinstance(row,dict):err(field,'required_measurement_missing');continue
        if set(row)-{'value','unit','matrix','assay'}:err(field,'unknown_measurement_metadata')
        v=row.get('value');u=row.get('unit')
        if not finite(v):err(field+'.value','finite_number_required_no_LOD_imputation');continue
        if not isinstance(u,str) or u not in sp['conversions']:err(field+'.unit','incompatible_unit');continue
        v=v*sp['conversions'][u]
        if not finite(v) or not sp['software_limits'][0]<=v<=sp['software_limits'][1]:err(field+'.value','outside_software_support');continue
        if row.get('matrix')!=sp['matrix']:err(field+'.matrix','matrix_mismatch')
        if sp['required_assay'] and row.get('assay')!=sp['required_assay']:err(field+'.assay','standardization_not_confirmed')
        elif not sp['required_assay']:warnings.append({'field':field,'code':'method_transport_not_validated'})
        values[k]=v;trace.append({'field':k,'input_value':row['value'],'input_unit':u,
            'canonical_value':v,'canonical_unit':sp['unit'],'matrix':row.get('matrix'),'assay':row.get('assay')})
    if values.get('hdl',0)>values.get('total_cholesterol',float('inf')):err('measurements.hdl','HDL_exceeds_total_cholesterol')
    return mid,values,trace,errors,warnings


def blocked(errors):
    return {'version':VERSION,'status':'blocked','clinical_use_ready':False,'probabilities':None,
            'errors':errors,'disclaimer':COMMON_WARNING}


def verify_probabilities(rows):
    for r in rows:
        p=[r['survival'],*r['cause_specific_cif'].values()]
        if any(not finite(v) or not 0<=v<=1 for v in p):raise ArithmeticError('Invalid probability')
        if abs(math.fsum(p)-1)>1e-10:raise ArithmeticError('Probability mass failure')
        if abs(r['all_cause_death_probability']-math.fsum(p[1:]))>1e-10:raise ArithmeticError('Inconsistent all-cause risk')


def run(payload):
    mid,values,trace,errors,warnings=normalize(payload)
    if errors:return blocked(errors)
    e=REGISTRY[mid]
    # Integrity failures abort the request; never return stale or cached predictions.
    m=model_spec(mid)
    f=engineer(payload['age_years'],payload['sex'],payload['clinical'],values)
    score=score_coarse if e['endpoint']==COARSE else score_detailed
    rows,features=score(f,m,payload['horizons_years']);verify_probabilities(rows)
    request_hash=sha256(json.dumps(payload,sort_keys=True,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode()).hexdigest()
    return {'version':VERSION,'status':'calculated_research_only','clinical_use_ready':False,
        'model_id':mid,'model_title':e['title'],'model_sha256':ARTIFACTS[e['artifact']][1],
        'endpoint':e['endpoint'],'request_sha256':request_hash,'probabilities':rows,
        'validation':evidence(mid,payload['horizons_years']), 'errors':[], 'warnings':warnings,
        'measurement_trace':trace,'feature_trace':features,
        'clipped_features':[x['feature'] for x in features if x['clipped']],
        'country':payload['country'],'age_years':payload['age_years'],'sex':payload['sex'],
        'dataset_context':payload['dataset_context'], 'laboratory_predictor_count':e['predictors'],
        'minimum_input_determinations':e['determinations'],
        'imputation_applied':False,'refitting_on_request':False,'biological_age':None,
        'noninfectious_natural_probability':None,'individual_uncertainty_interval':None,
        'disclaimer':COMMON_WARNING,
        'next_measurement_needed':False,
        'not_estimated':['treatment_effect','years_of_life_regained','contemporary_country_risk']}


def compare(payload, model_ids):
    # Display only, not a new validation or a model selection algorithm.
    if (not isinstance(model_ids,list) or not 2<=len(model_ids)<=4
        or any(not isinstance(x,str) or x not in REGISTRY for x in model_ids)
        or len(set(model_ids))!=len(model_ids)):
        return blocked([{'field':'models','code':'two_to_four_unique_models_required'}])
    mid,_,_,errors,_=normalize(payload)
    if errors:return blocked(errors)
    if REGISTRY[mid]['required']!=FULL or REGISTRY[mid]['endpoint']!=COARSE:
        return blocked([{'field':'input','code':'complete_full8_coarse_input_required_for_comparison'}])
    if any(REGISTRY[x]['endpoint']!=COARSE or not set(payload['horizons_years'])<=set(REGISTRY[x]['horizons']) for x in model_ids):
        return blocked([{'field':'models','code':'different_endpoint_or_horizon'}])
    results=[]
    for model_id in model_ids:
        p=deepcopy(payload);p['model_id']=model_id
        fields=REGISTRY[model_id]['required']
        excluded=sorted(set(p['measurements'])-set(fields))
        p['measurements']={k:v for k,v in p['measurements'].items() if k in fields}
        r=run(p)
        if r['status']=='blocked':return blocked([{'field':'comparison','code':'component_blocked'}])
        r['explicitly_unused_input_measurements']=excluded;results.append(r)
    return {'version':VERSION,'status':'compared_research_only','clinical_use_ready':False,
        'endpoint':COARSE,'horizons_years':deepcopy(payload['horizons_years']),
        'results':results,'probabilities':None,'errors':[],
        'input_domain':'same_complete_full8_profile','automatic_winner':None,'averaged_risk':None,
        'disclaimer':COMMON_WARNING+' Различие прогнозов для одного профиля не является доказательством качества или эквивалентности моделей.'}
