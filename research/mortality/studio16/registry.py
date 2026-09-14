"""Model identity, endpoint and validation scope are inseparable.

Only frozen artifacts from iterations12/14 are served. This module makes no
clinical-readiness decision and imports no participant-level records.
"""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from ..workbench.runtime import FIELDS, strict_json

ROOT = Path(__file__).resolve().parents[1]
BASE_COMMIT = '739eb5cd4dcb7dff24ca72b42b8571500ceecec5'
ARTIFACTS = {
    'compact': ('compact14/model_compact14.json', '274753ea36728652d10b10fa7d02b6f77b37304084bb2b0b056128b734e394d4'),
    'full': ('transport12/model_bundle12.json', '40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094'),
    'six': ('transport12/model_sensitivity12.json', '7c7dc918f042ff44f141a9ad3b49e6573db7f9732c5eb9d93e29b09b0686fa24'),
}
COMPACT = ['bmi','sbp','hba1c','creatinine','uacr','albumin']
SIX = COMPACT + ['total_cholesterol','hdl']
FULL = SIX + ['rdw','wbc']
COARSE = 'three_coarse_death_groups'
DETAILED = 'five_competing_death_groups'
LABELS = {'heart_diseases':'Болезни сердца (не только ИБС)',
    'malignant_neoplasms':'Злокачественные новообразования',
    'chronic_lower_respiratory':'Хронические болезни нижних дыхательных путей',
    'cerebrovascular':'Цереброваскулярные болезни',
    'other_or_unknown':'Другие и неуточнённые причины, включая инфекции и внешние причины'}
# Metrics are historical, immutable references; never computed on the submitted input.
METRICS14 = {
 'clinical3': (.0432887665481971,.7613219393645184,.0269999066171649),
 'compact4': (.0423250928134296,.787621733785503,.0261815719757362),
 'six6': (.0423212873067137,.7877102467257012,.0258966634301786),
 'full8': (.0426795118682321,.8088243045714316,.032794887793869),
}
REGISTRY = {
 'compact4': dict(title='Компактная панель: 4 показателя', artifact='compact', key=None, endpoint=COARSE, horizons=[1,3], required=COMPACT, predictors=4, determinations=5),
 'six6': dict(title='Панель: 6 показателей, без RDW и лейкоцитов', artifact='six', key=None, endpoint=COARSE, horizons=[1,3], required=SIX, predictors=6, determinations=7),
 'full8': dict(title='Расширенная панель: 8 показателей', artifact='full', key='routine_no_crp', endpoint=COARSE, horizons=[1,3], required=FULL, predictors=8, determinations=9),
 'clinical3': dict(title='Клинический контроль: общий полный профиль', artifact='full', key='clinical', endpoint=COARSE, horizons=[1,3], required=FULL, predictors=0, determinations=9),
 'full8_detailed': dict(title='8 показателей: пять групп причин', artifact='full', key='routine_no_crp', endpoint=DETAILED, horizons=[1,4,5], required=FULL, predictors=8, determinations=9),
}
COMMON_WARNING = ('Исторические данные США, исходный возраст 40–79 лет. Это прогноз общей смертности, '
    'а не диагноз, биологический возраст или эффект лечения. Клиническое применение, перенос на Россию '
    'и Германию и индивидуальные интервалы неопределённости не подтверждены.')


def load_artifact(name):
    path, expected = ARTIFACTS[name]
    raw = (ROOT/path).read_bytes()
    if sha256(raw).hexdigest() != expected:
        raise ValueError('Frozen artifact integrity failure: '+name)
    data = strict_json(raw)
    if data.get('clinical_use_ready') is not False:
        raise ValueError('Unexpected clinical-readiness declaration')
    return data


def model_spec(model_id):
    entry = REGISTRY.get(model_id) if isinstance(model_id,str) else None
    if entry is None: raise ValueError('Unknown model_id')
    data = load_artifact(entry['artifact'])
    return data['models'][entry['key']] if entry['key'] else data['model']


def evidence(model_id, horizons):
    entry=REGISTRY[model_id]
    rows=[]
    for h in horizons:
        row={'horizon_years':h, 'endpoint':entry['endpoint'], 'source_commit':BASE_COMMIT,
             'validation_kind':'historical_US_temporal', 'new_validation_in_this_release':False,
             'individual_uncertainty_interval':None, 'clinical_use_ready':False}
        if entry['endpoint']==COARSE:
            row.update(cohort='NHANES 2015–2016', n=2922, deaths=93 if h==3 else 29,
                population='40–79 лет; общий полный восьмипоказательный клинико-лабораторный профиль',
                source_path='research/mortality/compact14/REPORT_RU.md',
                evaluated_death_groups=3,
                caveat='Сравнение выполнено на общем полном профиле. Это не доказательство эквивалентности панелей.')
            if h==3:
                brier,auc,pred=METRICS14[model_id]
                row.update(multistate_brier=brier, all_cause_auc=auc, mean_predicted=pred,
                           observed_ipcw=.0222660851071608)
            else:
                row['caveat']+=' Годовой анализ вторичный; метрики за три года не переносятся на один год.'
        else:
            row.update(cohort='NHANES 2011–2012 / 2013–2014',
                source_path='research/mortality/ITERATION_12_RU.md', evaluated_death_groups=5,
                caveat='Для редких причин число событий мало. В пятилетнем респираторном исходе — две смерти; устойчивость точности не установлена.')
            if h==5: row.update(cohort='NHANES 2011–2012', n=2666, deaths=130)
            elif h==4: row.update(cohort='NHANES 2013–2014', n=3013, deaths=113)
            else: row['caveat']+=' Один год — вторичный горизонт, отдельные метрики здесь не показаны.'
        rows.append(row)
    return rows


def catalog():
    # All artifacts are checked; no stale cached evidence after file replacement.
    for name in ARTIFACTS: load_artifact(name)
    items=[]
    for model_id,e in REGISTRY.items():
        r=deepcopy(e);r.pop('key');r['model_id']=model_id
        r['model_sha256']=ARTIFACTS[e['artifact']][1]
        r['clinical_use_ready']=False;r['age_range']=[40,79];r['countries']=['US']
        r['measurements']={f:deepcopy(FIELDS[f]) for f in e['required']}
        r['evidence']=evidence(model_id,e['horizons'])
        r['input_domain_note']=('Контроль не использует анализы как предикторы, но требует общий профиль для сопоставимости.'
             if model_id=='clinical3' else 'Кратность лабораторных определений не равна измеренной стоимости панели.')
        items.append(r)
    return {'version':'0.16.0','models':items,'warning':COMMON_WARNING,'cause_labels':LABELS,
            'legacy_available_separately':True, 'default_model_selected':False,
            'model_selection_optimized_on_user_input':False}
