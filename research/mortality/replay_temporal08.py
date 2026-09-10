"""Opt-in research replay of the separately evaluated sex recalibration.

No change to the Workbench default or frozen feature coefficients. The complete
routine profile is mandatory, consistent with the informative validation domain.
The new specification is a historical-US research result, not clinical approval.
"""
from __future__ import annotations
from copy import deepcopy
from hashlib import sha256
import argparse
import json
import math
from pathlib import Path
from .workbench.runtime import run, strict_json, BUNDLE_SHA256

SPEC_PATH=Path(__file__).parent/'validation08/recalibration_spec.json'
SPEC_SHA256='ed8d8b8d6aa464bc42b5d47b70176677895b5823fe899ee81494e9c6804a353c'


def recalibrate_probability(probability: float,factor: float) -> float:
    if type(probability) not in (int,float) or type(factor) not in (int,float):
        raise ValueError('Finite numerical probability and factor required')
    if not math.isfinite(probability) or not math.isfinite(factor) or not 0<=probability<=1 or factor<=0:
        raise ValueError('Invalid probability or factor')
    return 1. if probability==1 else -math.expm1(factor*math.log1p(-probability))


def replay(payload: dict,*,acknowledge_temporal_research: bool=False,spec_path: Path=SPEC_PATH) -> dict:
    if acknowledge_temporal_research is not True:
        raise ValueError('Explicit temporal-research acknowledgement required')
    if not isinstance(payload,dict) or payload.get('model')!='M1_routine':
        raise ValueError('This opt-in path requires M1 and the complete routine profile')
    raw=spec_path.read_bytes()
    if sha256(raw).hexdigest()!=SPEC_SHA256:raise ValueError('Changed recalibration specification')
    spec=strict_json(raw)
    if spec['parent_model_sha256']!=BUNDLE_SHA256 or spec['clinical_use_ready'] is not False:
        raise ValueError('Incompatible model provenance or scope')
    result=run(payload)  # Retains strict original units, country, age and completeness gates.
    if result['status']=='blocked':return result
    factor=spec['factors']['M1_routine'][payload['sex']]
    result=deepcopy(result)
    result['frozen_probabilities']=result['probabilities']
    result['probabilities']=[{**row,'probability':recalibrate_probability(row['probability'],factor)}
                             for row in result['frozen_probabilities']]
    result['recalibration']={'specification_sha256':SPEC_SHA256,'factor':factor,
        'estimation_cycle':'NHANES2003_2004','evaluation_cycles':['NHANES2005_2006','NHANES2007_2008'],
        'validation_domain':'complete_clinical_and_routine_measurements',
        'analysis_type':'prespecified_complete_case_sensitivity','evaluation_n':5287,
        'evaluation_events_10y':760,'individual_prediction_interval':None,
        'clinical_use_ready':False,'changes_to_feature_coefficients':False}
    # Preserve earlier warnings as history, not silently erase known model limitations.
    result['historical_validation_warnings']=result['warnings']
    result['warnings']=[{'code':'historical_US_temporal_research_only'},
        {'code':'complete_profile_validation_does_not_generalize_to_missing_data'},
        {'code':'no_individual_uncertainty_interval_or_clinical_decision_threshold'}]
    result['disclaimer']=('Исследовательский прогноз общей смертности по историческим данным США, '
        'возраст 40–79 лет. Половая рекалибровка оценена на NHANES 2003–2004, '
        'проверена на циклах 2005–2008. Требуется полный клинико-лабораторный профиль. '
        'Все причины смерти включены. Клинические решения, современные национальные риски '
        'и причинные эффекты вмешательств этим расчётом не устанавливаются.')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--acknowledge-temporal-research',action='store_true')
    a=p.parse_args()
    if a.out.exists():raise SystemExit('Output already exists')
    r=replay(strict_json(a.input.read_bytes()),acknowledge_temporal_research=a.acknowledge_temporal_research)
    with a.out.open('x',encoding='utf-8') as handle:json.dump(r,handle,ensure_ascii=False,indent=2,allow_nan=False)
