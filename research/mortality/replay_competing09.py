"""Opt-in, stdlib-only research replay of the three-cause model 09.

This is a new fitted model, not an allocation of the old all-cause prediction.
Both panels require the complete clinical/routine profile used in validation.
The old runtime is imported for unit validation and feature algebra ONLY.
"""
from __future__ import annotations
from copy import deepcopy
from hashlib import sha256
import argparse
import json
import math
from pathlib import Path
from .workbench.runtime import normalize, engineer, example as old_example, strict_json, finite

BUNDLE_PATH=Path(__file__).parent/'competing09/model_bundle09.json'
BUNDLE_SHA256='e6c0a9b3ca57ef9556c1556a28da03271a6caf692f5649bfc804e9fba8747f31'
GROUPS=('heart_diseases','malignant_neoplasms','other_or_unknown')
HORIZONS=(1,5,8)
DISCLAIMER=('Исследовательская модель по историческим данным США; возраст 40–79 лет. '
    'Первоначальные причины объединены в болезни сердца, злокачественные новообразования '
    'и остальные/неуточнённые причины. Инфекции и внешние причины не исключены. '
    'Болезни сердца не означают только ИБС. Требуется полный клинико-лабораторный профиль. '
    'Расчёт не предназначен для диагноза или выбора лечения.')


def load_bundle(path=BUNDLE_PATH):
    raw=Path(path).read_bytes()
    if sha256(raw).hexdigest()!=BUNDLE_SHA256:
        raise ValueError('Контрольная сумма модели 09 не совпадает')
    bundle=strict_json(raw)
    if (bundle.get('clinical_use_ready') is not False or bundle.get('groups')!=list(GROUPS)
            or bundle.get('noninfectious_risk_identified') is not False):
        raise ValueError('Несовместимая область модели')
    return bundle


def example(panel='routine'):
    if panel not in ('clinical','routine'):raise ValueError('Unknown panel')
    payload=old_example('M1_routine')
    payload.update(schema_version='0.9.0',endpoint='competing_death_groups',
                   panel=panel,horizons_years=list(HORIZONS))
    payload.pop('model')
    return payload


def score_engineered(features, model, horizons):
    """No fitting and no imputation. Exact coherent piecewise integration."""
    if (not isinstance(horizons,list) or not horizons or len(horizons)>3
            or any(not finite(h) or h not in HORIZONS for h in horizons)
            or len(horizons)!=len(set(horizons))):
        raise ValueError('Only unique 1/5/8-year horizons are supported')
    pp=model['preprocessing'];z=[];trace=[]
    for j,name in enumerate(pp['features']):
        raw=features.get(name)
        if not finite(raw):raise ValueError('Missing/nonfinite engineered predictor: '+name)
        used=min(max(raw,pp['lower'][j]),pp['upper'][j])
        standardized=(used-pp['mean'][j])/pp['sd'][j]
        z.append(standardized)
        trace.append({'feature':name,'raw_value':raw,'used_value':used,
                      'standardized_value':standardized,'clipped':used!=raw})
    ends=model['interval_ends_years']
    if ends!=[1.,5.,8.]:raise ValueError('Unsupported interval specification')
    rates=[]
    for group in GROUPS:
        spec=model['causes'][group]
        if len(spec['coefficients'])!=len(z):raise ValueError('Coefficient dimensions differ')
        lp=math.fsum(x*b for x,b in zip(z,spec['coefficients']))
        log_rates=[lp+a for a in spec['baseline_log_hazards']]
        if any(not finite(a) or a>700 for a in log_rates):raise ValueError('Hazard outside numeric support')
        rates.append([math.exp(a) for a in log_rates])
    rows=[]
    for horizon in horizons:
        survival=1.;cif=[0.,0.,0.]
        for j,(start,end) in enumerate(zip([0.,1.,5.],ends)):
            dt=max(0.,min(horizon,end)-start)
            total=math.fsum(r[j] for r in rates)
            if not finite(total):raise ValueError('Summed hazard overflow')
            mass=-math.expm1(-total*dt)
            if total>0:
                for k in range(3):cif[k]+=survival*mass*(rates[k][j]/total)
            survival*=math.exp(-total*dt)
        if abs(survival+math.fsum(cif)-1)>1e-10:raise ArithmeticError('Probability mass mismatch')
        rows.append({'horizon_years':horizon,'survival':survival,
            'cause_specific_cif':dict(zip(GROUPS,cif)),'all_cause_death_probability':math.fsum(cif)})
    return rows,trace


def run(payload,*,bundle_path=BUNDLE_PATH):
    bundle=load_bundle(bundle_path)
    errors=[]
    if not isinstance(payload,dict):payload={};errors.append({'field':'root','code':'expected_object'})
    if payload.get('schema_version')!='0.9.0':errors.append({'field':'schema_version','code':'expected_0.9.0'})
    if payload.get('endpoint')!='competing_death_groups':errors.append({'field':'endpoint','code':'coarse_competing_groups_only'})
    panel=payload.get('panel')
    if not isinstance(panel,str) or panel not in ('clinical','routine'):
        errors.append({'field':'panel','code':'explicit_clinical_or_routine_required'})
    if 'model' in payload:errors.append({'field':'model','code':'use_explicit_panel_not_legacy_model'})
    horizons=payload.get('horizons_years')
    if (not isinstance(horizons,list) or not horizons or len(horizons)>3
            or any(not finite(h) or h not in HORIZONS for h in horizons)
            or (all(finite(h) for h in horizons) and len(horizons)!=len(set(horizons)))):
        errors.append({'field':'horizons_years','code':'only_unique_1_5_8_year_horizons'})
    # Adapter routing does not change the requested scientific endpoint: no old
    # score is ever computed. The M1 contract checks full profile in BOTH panels.
    adapter=deepcopy(payload);adapter.pop('panel',None)
    adapter.update(schema_version='0.6.0',endpoint='all_cause_death',model='M1_routine',horizons_years=[5])
    values,measurement_trace,input_errors,warnings=normalize(adapter)
    errors.extend(input_errors)
    result={'version':'0.9.0','status':'blocked' if errors else 'calculated_research_only',
        'clinical_use_ready':False,'endpoint':'competing_death_groups','panel':panel if isinstance(panel,str) else None,
        'model_sha256':BUNDLE_SHA256,'probabilities':None,'individual_uncertainty_interval':None,
        'noninfectious_natural_probability':None,'imputation_applied':False,'fit_on_request':False,
        'declaration':'synthetic_not_a_person' if payload.get('dataset_context')=='synthetic_example' else 'historical_research',
        'warnings':warnings,'errors':errors,'disclaimer':DISCLAIMER,
        'validation':{'cycles':[2009,2010],'n':3182,'baseline_age_range':[40,79],
                      'profile':'complete_clinical_and_routine','deaths_5y':150,'deaths_8y':302},
        'measurement_trace':measurement_trace}
    if errors:return result
    features=engineer(payload['age_years'],payload['sex'],payload['clinical'],values)
    probabilities,trace=score_engineered(features,bundle['models'][panel],horizons)
    result.update(probabilities=probabilities,feature_trace=trace,
        clipped_features=[r['feature'] for r in trace if r['clipped']])
    if panel=='clinical':result['warnings'].append({'code':'complete_laboratory_profile_required_for_validation_domain_not_used_in_clinical_score'})
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,nargs='?')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--demo',action='store_true')
    args=parser.parse_args()
    if bool(args.input)==args.demo:parser.error('Provide exactly an input file or --demo')
    if args.out.exists():parser.error('Output exists; use a new file')
    payload=example() if args.demo else strict_json(args.input.read_bytes())
    result=run(payload)
    with args.out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps({'status':result['status'],'output':str(args.out),'clinical_use_ready':False}))
    return 2 if result['status']=='blocked' else 0

if __name__=='__main__':raise SystemExit(main())
