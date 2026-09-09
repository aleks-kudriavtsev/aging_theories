"""Replay the reader-corrected, frozen NHANES model with typed raw inputs, using stdlib only.

This request path does not fit or calibrate. The bundled artifact was refitted
after a source-reader correction on the same previously examined data.
Country and endpoint are declarations of intended use, not inferred identities.
Coefficient artifact: reanalysis06/model_bundle.json. Legacy algorithm base:
aging_theories@fe414dc; this local corrective artifact has no remote commit yet.
"""
from __future__ import annotations
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
from ..engine import egfr_ckd_epi_2021

BUNDLE_SHA256 = '7b287f3f4778b2127f2226f79d2e503e1e39087306a172e5f524a561964b07d0'
BUNDLE_PATH = Path(__file__).resolve().parents[1] / 'reanalysis06/model_bundle.json'
SOURCE_COMMIT = None  # Content-addressed by BUNDLE_SHA256; release commit is separate provenance.
BASE_CODE_COMMIT = 'fe414dc1957d8c6041dce48e77cf58ee65b30d28'
MAX_BODY = 1_048_576
MODELS = {
    'M0_age_sex': 'Возраст и пол — контроль',
    'M0c_clinical_only_posthoc': 'Клинический контроль — post hoc',
    'M1_routine': 'Клинические данные + рутинная лаборатория',
    'M2_NTproBNP': 'Рутинная панель + NT-proBNP',
    'M3_troponin': 'Предыдущая панель + hs-cTnT',
    'M4_cystatinC': 'Предыдущая панель + цистатин C',
}
# Input sanity limits below are software limits, NOT clinical reference intervals.
# Known out-of-training feature values are winsorized exactly as in the frozen model.
def field(label, unit, matrix, limits, conversions=None, assay=None):
    return {'label': label, 'unit': unit, 'matrix': matrix, 'software_limits': limits,
            'conversions': conversions or {unit: 1.0}, 'required_assay': assay}
FIELDS = {
    'bmi': field('Индекс массы тела', 'kg/m2', 'clinical', [5,150]),
    'sbp': field('Среднее систолическое АД', 'mmHg', 'clinical', [40,300]),
    'total_cholesterol': field('Общий холестерин', 'mg/dL', 'serum', [1,2000], {'mg/dL':1., 'mmol/L':38.67}),
    'hdl': field('Холестерин HDL', 'mg/dL', 'serum', [1,500], {'mg/dL':1., 'mmol/L':38.67}),
    'hba1c': field('HbA1c (NGSP)', '%', 'whole_blood', [1,30]),
    'creatinine': field('Креатинин, стандартизованный по IDMS', 'mg/dL', 'serum', [.01,30], {'mg/dL':1., 'umol/L':1/88.4}, 'IDMS_standardized'),
    'uacr': field('Альбумин/креатинин мочи (UACR)', 'mg/g', 'urine', [.0001,100000], {'mg/g':1.,'mg/mmol':8.84}),
    'albumin': field('Альбумин', 'g/L', 'serum', [1,100], {'g/L':1.,'g/dL':10.}),
    'rdw': field('RDW-CV (не RDW-SD)', '%', 'whole_blood', [1,50]),
    'wbc': field('Лейкоциты', '10^9/L', 'whole_blood', [.001,500], {'10^9/L':1.,'10^3/uL':1.,'cells/uL':.001}),
    'crp': field('C-реактивный белок', 'mg/L', 'serum', [.0001,1000], {'mg/L':1.,'mg/dL':10.}),
    'ntprobnp': field('NT-proBNP (не BNP)', 'pg/mL', 'serum', [.001,1000000], {'pg/mL':1.,'ng/L':1.}, 'Roche_Cobas_e601'),
    'hs_ctnt': field('Высокочувствительный тропонин T (не I)', 'ng/L', 'serum', [.001,100000], {'ng/L':1.,'pg/mL':1.}, 'Roche_Elecsys_TnT_Gen5_e601'),
    'cystatin_c': field('Цистатин C', 'mg/L', 'serum', [.001,100], assay='Siemens_Dimension_Vista_1500'),
}
CLINICAL = ['smoking','bp_treatment','diabetes_history','cvd_history','cancer_history']
ROUTINE = ['total_cholesterol','hdl','hba1c','creatinine','uacr','albumin','rdw','wbc','crp']
REQUIRED = {
    'M0_age_sex': [], 'M0c_clinical_only_posthoc': ['bmi','sbp'],
    'M1_routine': ['bmi','sbp']+ROUTINE,
    'M2_NTproBNP': ['bmi','sbp']+ROUTINE+['ntprobnp'],
    'M3_troponin': ['bmi','sbp']+ROUTINE+['ntprobnp','hs_ctnt'],
    'M4_cystatinC': ['bmi','sbp']+ROUTINE+['ntprobnp','hs_ctnt','cystatin_c'],
}
ROOT_FIELDS = {'schema_version','research_only','acknowledge_limitations','dataset_context',
               'country','endpoint','model','horizons_years','age_years','sex','clinical','measurements'}
DISCLAIMER = (
    'Исследовательский расчёт общей смертности по исторической модели США. '
    'Это не диагноз, не рекомендация по лечению и не клинически валидированный индивидуальный риск. '
    'Инфекции и внешние причины включены; разделение по причинам не выполнено. '
    'Временная проверка: NHANES 2003–2004, возраст 40–79 лет; исходы до 2019 года. '
    'Выполнен исправленный повторный анализ уже просмотренной выборки, не новая независимая валидация. '
    'Калибровка для России, Германии и современного населения США не установлена.'
)


def finite(value):
    """JSON integers are unbounded; converting a huge integer must not crash."""
    if type(value) not in (float, int):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def strict_json(raw):
    """Refuse duplicates, nonfinite JSON extensions and oversized documents."""
    if not isinstance(raw,(str,bytes)) or len(raw if isinstance(raw,bytes) else raw.encode()) > MAX_BODY:
        raise ValueError('Размер JSON превышает допустимый предел')
    def pairs(items):
        obj={}
        for k,v in items:
            if k in obj: raise ValueError('Повторяющееся поле JSON: '+k)
            obj[k]=v
        return obj
    def constant(value):
        raise ValueError('Недопустимое число JSON: '+value)
    def numeric(text, convert):
        value = convert(text)
        if not finite(value):
            raise ValueError('Число вне конечного диапазона JSON')
        return value
    try:
        return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant,
                          parse_int=lambda text: numeric(text, int),
                          parse_float=lambda text: numeric(text, float))
    except (UnicodeDecodeError, RecursionError) as exc:
        raise ValueError('Некорректная структура или кодировка JSON') from exc


def load_bundle(path=BUNDLE_PATH):
    content=Path(path).read_bytes()
    if sha256(content).hexdigest()!=BUNDLE_SHA256:
        raise ValueError('Контрольная сумма модели не совпадает: расчёт остановлен')
    bundle=strict_json(content)
    if bundle['clinical_use'] is not False or bundle['endpoint']!='all_cause_death':
        raise ValueError('Несовместимое назначение модели')
    return bundle


def schema():
    return {'version':'0.7.0','input_schema_version':'0.6.0','models':MODELS,'required_measurements':REQUIRED,
            'measurements':FIELDS,'clinical_fields':CLINICAL,'country':'US',
            'endpoint':'all_cause_death','horizons_years':[1,5,10],
            'age_range':[40,79], 'sex_values':['female','male'],
            'imputation':'disabled: missing required inputs block the calculation',
            'clinical_use_ready':False,'disclaimer':DISCLAIMER}


def example(model='M1_routine'):
    if model not in MODELS: raise ValueError('Неизвестная модель')
    values={'bmi':26.,'sbp':125.,'total_cholesterol':195.,'hdl':55.,'hba1c':5.4,
            'creatinine':.9,'uacr':8.,'albumin':43.,'rdw':12.8,'wbc':6.5,
            'crp':1.5,'ntprobnp':60.,'hs_ctnt':6.,'cystatin_c':.8}
    return {'schema_version':'0.6.0','research_only':True,'acknowledge_limitations':True,
            'dataset_context':'synthetic_example','country':'US','endpoint':'all_cause_death',
            'model':model,'horizons_years':[1,5,10], 'age_years':55, 'sex':'female',
            'clinical':({} if model=='M0_age_sex' else {'smoking':'never','bp_treatment':False,
                'diabetes_history':False,'cvd_history':False,'cancer_history':False}),
            'measurements':{k:{'value':values[k], 'unit':FIELDS[k]['unit'],
                'matrix':FIELDS[k]['matrix'],'assay':FIELDS[k]['required_assay'] or 'not_reported'} for k in REQUIRED[model]}}


def normalize(payload):
    errors=[]; warnings=[]; canonical={}; trace=[]
    def err(path,code): errors.append({'field':path,'code':code})
    if not isinstance(payload,dict):
        return {},[],[{'field':'root','code':'expected_object'}],[]
    for key in set(payload)-ROOT_FIELDS: err(key,'unknown_field_or_identifier_not_allowed')
    if payload.get('schema_version')!='0.6.0': err('schema_version','expected_0.6.0')
    for flag in ['research_only','acknowledge_limitations']:
        if payload.get(flag) is not True: err(flag,'explicit_true_required')
    if payload.get('dataset_context') not in ('synthetic_example','historical_US_research'):
        err('dataset_context','historical_research_or_synthetic_only')
    if payload.get('country')!='US': err('country','country_transport_not_validated')
    if payload.get('endpoint')!='all_cause_death': err('endpoint','cause_specific_model_not_available')
    model=payload.get('model')
    if not isinstance(model,str) or model not in MODELS:
        err('model','unknown_model'); model=None
    a=payload.get('age_years')
    if not finite(a) or not 40<=a<=79: err('age_years','baseline_age_must_be_40_to_79')
    if payload.get('sex') not in ('female','male'): err('sex','published_sex_variable_required')
    h=payload.get('horizons_years')
    if (not isinstance(h,list) or not h or len(h)>3 or
        any(not finite(x) or x not in (1,5,10) for x in h) or (isinstance(h,list) and all(finite(x) for x in h) and len(h)!=len(set(h)))):
        err('horizons_years','only_unique_1_5_10_year_horizons')
    clinical=payload.get('clinical')
    if not isinstance(clinical,dict):
        err('clinical','expected_object');clinical={}
    for key in set(clinical)-set(CLINICAL): err('clinical.'+key,'unknown_field')
    if model and model!='M0_age_sex':
        if clinical.get('smoking') not in ('never','former','current'): err('clinical.smoking','required_never_former_current')
        for key in CLINICAL[1:]:
            if type(clinical.get(key)) is not bool: err('clinical.'+key,'explicit_boolean_required')
    measurements=payload.get('measurements')
    if not isinstance(measurements,dict): err('measurements','expected_object'); measurements={}
    for k in set(measurements)-set(FIELDS): err('measurements.'+k,'unknown_marker_no_substitution')
    for k in REQUIRED.get(model,[]):
        spec=FIELDS[k]; row=measurements.get(k); path='measurements.'+k
        if not isinstance(row,dict): err(path,'required_measurement_missing');continue
        for key in set(row)-{'value','unit','matrix','assay'}: err(path+'.'+key,'unknown_field')
        value=row.get('value'); unit=row.get('unit')
        if not finite(value): err(path+'.value','finite_number_required_no_LOD_imputation');continue
        if not isinstance(unit,str) or unit not in spec['conversions']:
            err(path+'.unit','incompatible_or_missing_units');continue
        value=float(value)*spec['conversions'][unit]
        if not finite(value) or not spec['software_limits'][0]<=value<=spec['software_limits'][1]:
            err(path+'.value','outside_software_support_not_a_clinical_reference_range');continue
        if row.get('matrix')!=spec['matrix']: err(path+'.matrix','matrix_mismatch_or_missing')
        required=spec['required_assay']
        if required and row.get('assay')!=required: err(path+'.assay','reference_assay_or_standardization_not_confirmed')
        if not required and spec['matrix']!='clinical':
            warnings.append({'field':path+'.assay','code':'routine_assay_transport_not_validated'})
        canonical[k]=value
        trace.append({'field':k,'input_value':row['value'],'input_unit':unit,
                      'canonical_value':value,'canonical_unit':spec['unit'],
                      'matrix':row.get('matrix'),'assay':row.get('assay')})
    ignored=sorted(set(measurements)&set(FIELDS)-set(REQUIRED.get(model,[])))
    if ignored: warnings.append({'code':'measurements_not_used_by_selected_model','fields':ignored})
    if canonical.get('hdl',0)>canonical.get('total_cholesterol',float('inf')):
        err('measurements.hdl','HDL_exceeds_total_cholesterol_check_units')
    return canonical,trace,errors,warnings


def engineer(age, sex, clinical, values):
    """Same raw-feature algebra as frozen make_cohort; no 1999 lab correction.

    The caller provides already IDMS-standardized creatinine. The cycle-specific
    1999 raw calibration belongs ONLY to reading that historical source dataset.
    """
    male=1. if sex=='male' else 0.
    age_scaled=(age-60)/10
    f={'age':age_scaled,'male':male,'age_male':age_scaled*male,
       'age_rcs':(max(age-40,0)**3-(39/19)*max(age-60,0)**3+(20/19)*max(age-79,0)**3)/(39**2)}
    if 'smoking' in clinical:
        f['smoker_current']=float(clinical['smoking']=='current')
        f['smoker_former']=float(clinical['smoking']=='former')
        f.update({k:float(clinical[k]) for k in CLINICAL[1:]})
    if 'bmi' in values:
        f['bmi']=(values['bmi']-27)/5; f['bmi_sq']=f['bmi']**2
    if 'sbp' in values:f['sbp']=(values['sbp']-130)/20
    for k in ('total_cholesterol','hdl','hba1c','albumin','rdw'):
        if k in values:f[k]=values[k]
    if 'creatinine' in values:
        egfr=egfr_ckd_epi_2021(age,sex,values['creatinine'])
        f['egfr_low']=min(egfr,60)/15;f['egfr_high']=max(egfr-60,0)/30
    for source,target in [('uacr','log_uacr'),('wbc','log_wbc'),('crp','log_crp'),('ntprobnp','log_ntprobnp'),('hs_ctnt','log_troponin'),('cystatin_c','log_cystatin')]:
        if source in values:f[target]=math.log2(values[source])
    return f


def score_features(features, model, preprocessing, horizons):
    """Exact finite-data replay. Missingness is blocked at the raw-input layer."""
    z=[]; columns=[]; trace=[]
    for name in model['features']:
        raw=features[name]
        if not finite(raw): raise ValueError('Missing/nonfinite engineered feature')
        p=preprocessing[name];used=min(max(raw,p['lower']),p['upper'])
        z.append((used-p['mean'])/p['sd']);columns.append(name)
        trace.append({'feature':name,'raw_value':raw,'used_value':used,
                      'standardized_value':z[-1], 'clipped':used!=raw})
        if p['missing_indicator']:z.append(0.);columns.append(name+'__missing')
    if columns!=model['columns'] or len(z)!=len(model['coefficients']):
        raise ValueError('Несовместимый порядок признаков')
    lp=math.fsum(x*b for x,b in zip(z,model['coefficients']))
    if not math.isfinite(lp) or lp>700:raise ValueError('Расчёт за пределами численной поддержки')
    risks=[]
    for h in horizons:
        dt=[max(0.,min(h,b)-a) for a,b in [(0,1),(1,5),(5,10)]]
        hazard=math.exp(lp)*math.fsum(d*math.exp(b) for d,b in zip(dt,model['baseline_log_hazards']))
        risks.append({'horizon_years':h,'probability':-math.expm1(-hazard)})
    return risks,trace


def run(payload, *, bundle_path=BUNDLE_PATH):
    bundle=load_bundle(bundle_path)  # Hash check before every request; no arbitrary model upload.
    values,trace,errors,warnings=normalize(payload)
    base={'application_version':'0.7.0','status':'blocked' if errors else 'calculated_research_only',
          'clinical_use_ready':False,'disclaimer':DISCLAIMER,
          'model_sha256':BUNDLE_SHA256,'model_source_commit':SOURCE_COMMIT,
          'endpoint':'all_cause_death','cause_specific_probabilities':None,
          'individual_uncertainty_interval':None, 'fit_on_request':False,
          'model_origin':deepcopy(bundle['origin']),'model_code_base_commit':BASE_CODE_COMMIT,
          'model_artifact_version':bundle['version'],
          'errors':errors,'warnings':warnings,'probabilities':None}
    if errors:return base
    name=payload['model'];clinical=payload['clinical'] if name!='M0_age_sex' else {}
    features=engineer(payload['age_years'],payload['sex'],clinical,values)
    risks,ft=score_features(features,bundle['models'][name],bundle['preprocessing'],sorted(payload['horizons_years']))
    clipped=[r['feature'] for r in ft if r['clipped']]
    if clipped:warnings.append({'code':'training_winsorization_applied_not_reference_limits','features':clipped})
    if payload['sex']=='female':warnings.append({'code':'historical_validation_overprediction_in_women'})
    if payload['age_years']>=60:warnings.append({'code':'lower_discrimination_in_older_validation_stratum'})
    warnings.append({'code':'no_clinical_calibration_or_decision_thresholds'})
    if 1 in payload['horizons_years']:
        warnings.append({'code':'one_year_validation_only_33_deaths'})
    # Never emit an age-equivalent, diagnoses, treatment changes or causal contributions.
    base.update(model=name,model_label=MODELS[name],dataset_context=payload['dataset_context'],
                age_years=payload['age_years'],sex=payload['sex'],country='US',
                probabilities=risks,measurement_trace=trace,feature_trace=ft,
                imputation_applied=False,clipped_features=clipped,
                declaration='synthetic_not_a_person' if payload['dataset_context']=='synthetic_example' else 'user_declared_historical_research',
                post_hoc_model=name=='M0c_clinical_only_posthoc')
    if 'creatinine' in values:base['egfr_ckd_epi_2021']=egfr_ckd_epi_2021(payload['age_years'],payload['sex'],values['creatinine'])
    return base
