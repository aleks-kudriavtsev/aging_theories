"""Freeze two final routine29 fits; export parameters and aggregate audit only.

Inner heldout Brier is a training-selection diagnostic, never a new validation
estimate. Participant outcomes/measurements and prediction arrays are transient
local inputs and are not written or printed by this module.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from ..routine29 import model
from ..routine29.cohort import load
from ..routine29.features import CORE
from ..binary28.cohort import PROTEINS
from . import adapter

ROOT=Path(__file__).parent
ARMS=('routine','routine_all3')
INNER_SALT='freeze30:final-inner:v1'
SCHEMA='freeze30-model/1'


def _canonical_bytes(value):
    return (json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode('utf-8')


def fit_final(core,protein,y,w,clusters):
    """Select a penalty in3whole-PSU folds and refit on all supplied training.

    The returned payload contains only learned parameters and aggregate support
    and optimization diagnostics; no fitted/heldout probability arrays or IDs.
    """
    core,protein=model._features(core,protein)
    y,w=model._rows(y,w,len(core));clusters=np.asarray(clusters)
    if clusters.shape!=(len(core),) or any(not isinstance(c,str) or not c for c in clusters.tolist()):
        raise ValueError('Aligned nonempty string PSU identities required')
    if np.unique(y[w>0]).size!=2:raise ValueError('Training needs positive-weight events and nonevents')
    folds=model.psu_partition(clusters,INNER_SALT,folds=3)
    result={}
    for arm in ARMS:
        predictions={penalty:np.full(len(y),np.nan) for penalty in model.PENALTIES}
        aggregate=[];diagnostics=[]
        for fold in range(3):
            training,scoring=folds!=fold,folds==fold
            train_psus=set(clusters[training].tolist());score_psus=set(clusters[scoring].tolist())
            if train_psus&score_psus:raise ValueError('PSU crosses final-inner boundary')
            prep=model.fit_preprocessing(core[training],protein[training],w[training],arm)
            ztrain=model.transform_features(core[training],protein[training],prep)
            zscore=model.transform_features(core[scoring],protein[scoring],prep)
            aggregate.append(dict(fold=fold,training_n=int(training.sum()),scoring_n=int(scoring.sum()),
                training_events=int(y[training].sum()),scoring_events=int(y[scoring].sum()),
                training_positive_weight_n=int((training&(w>0)).sum()),
                training_positive_weight_events=int((training&(w>0)&(y==1)).sum()),
                scoring_positive_weight_n=int((scoring&(w>0)).sum()),
                training_psu_n=len(train_psus),scoring_psu_n=len(score_psus),PSU_overlap_n=0))
            for penalty in model.PENALTIES:
                fitted=model.fit_ridge(ztrain,y[training],w[training],penalty)
                predictions[penalty][scoring]=model.predict_ridge(zscore,fitted)
                diagnostics.append(dict(fold=fold,penalty=float(penalty),
                    gradient_max=fitted['gradient_max'],iterations=fitted['iterations'],
                    convergence_status=fitted['convergence_status'],
                    normalized_weight_sum=fitted['normalized_weight_sum']))
        scores={}
        for penalty,p in predictions.items():
            if not np.isfinite(p).all():raise ArithmeticError('Missing final-inner probability')
            scores[penalty]=float(np.average((y-p)**2,weights=w))
        minimum=min(scores.values())
        chosen=max(p for p in model.PENALTIES if scores[p]<=minimum+1e-12)
        prep=model.fit_preprocessing(core,protein,w,arm)
        fitted=model.fit_ridge(model.transform_features(core,protein,prep),y,w,chosen)
        names=list(model.CORE_NAMES)+([] if arm=='routine' else list(model.PROTEIN_NAMES))
        result[arm]=dict(preprocessing=prep,fit=fitted,feature_names=names,
            inner_brier_by_penalty={str(p):score for p,score in scores.items()},selected_penalty=float(chosen),
            inner_folds_aggregate=aggregate,inner_fit_diagnostics=diagnostics)
    return dict(models=result,training=dict(n=len(y),events8=int(y.sum()),
        positive_weight_n=int((w>0).sum()),weight_sum=float(w.sum()),psu_n=len(set(clusters.tolist())),
        inner_folds=3,inner_salt=INNER_SALT,weight_normalization='sum_to_raw_training_n'),
        selection_scope='training_only_not_validation')


def _check_training_sources(source):
    """Pin all historical29code and non-index scientific/doc source bytes."""
    adapter.check_study()
    prior=adapter.strict_json(ROOT.parent/'routine29/RESULTS_SUMMARY.json')
    for name,digest in prior['scientific_python_sha256'].items():
        if adapter.digest(ROOT.parent/name)!=digest:raise ValueError('Historical scientific code changed')
    d,common,strict,audit,manifest=load(source)
    actual={row['file']:row['sha256'] for row in manifest['files']}
    expected=prior['source_audit']['routine_source_fingerprints']
    if set(actual)!=set(expected):raise ValueError('Historical source file set changed')
    if any(actual[name]!=digest for name,digest in expected.items() if name!='index.htm'):
        raise ValueError('Historical scientific/documentation source bytes changed')
    old_manifest=adapter.strict_json(ROOT.parent/'routine29/SOURCE_MANIFEST.json')
    urls={row['file']:(row['url'],row['resolved_url']) for row in manifest['files']}
    old_urls={row['file']:(row['url'],row['resolved_url']) for row in old_manifest['files']}
    if urls!=old_urls:raise ValueError('Historical resolved release URLs changed')
    source_policy=dict(pinned_non_index_files=len(expected)-1,
        index_is_dynamic_catalog_not_scientific_table=True,
        prior_index_sha256=expected['index.htm'],actual_index_sha256=actual['index.htm'],
        index_changed=actual['index.htm']!=expected['index.htm'],
        every_non_index_source_sha_unchanged=True,every_resolved_release_url_unchanged=True)
    return d,common,audit,manifest,source_policy


def _integer(value):
    if type(value) is bool or not np.isfinite(value) or value!=int(value):
        raise ValueError('Training metadata category is not an integer')
    return int(value)


def _training_document(d):
    """Temporary canonical batch for a numerical adapter replay; never export."""
    contract=adapter.check_study();records=[]
    for row in d.to_dict('records'):
        smoking='never' if row['HAR1']==2 else 'current' if row['HAR3']==1 else 'former' if row['HAR3']==2 else None
        records.append(dict(age_years=_integer(row['HSAGEIR']),sex='male' if row['HSSEX']==1 else 'female',
            smoking=smoking,sbp_mmhg=float(row['HAZMNK1R']),sbp_measurement_count=_integer(row['HAZNOK1R']),
            bmi_kg_m2=float(row['BMPBMI']),bmi_weight_flag=_integer(row['BMPWTFLG']),
            bmi_height_flag=_integer(row['BMPHTFLG']),hba1c_percent=float(row['GHP']),
            hba1c_method=_integer(row['GHPMETH']),creatinine_historical_mg_dl=float(row['CEP']),
            urine_creatinine_mg_dl=float(row['URP']),urine_albumin_ug_ml=float(row['UBP']),
            cystatin_c_reported_mg_l=float(row['SSCYSTAT']),beta2_microglobulin_reported_mg_l=float(row['SSB2M']),
            beta_trace_protein_reported_mg_l=float(row['SSBTP'])))
    return dict(schema=contract['input_schema'],contract_sha256=adapter.digest(ROOT/'CONTRACT.json'),
        mapping_id=contract['mapping_id'],units=contract['units'],records=records)


def run(source,out):
    source,out=Path(source),Path(out)
    if out.exists():raise ValueError('New freeze output directory required')
    lock=adapter.strict_json(ROOT/'LOCK.json');contract=adapter.check_study()
    eligible,common,audit,source_manifest,source_policy=_check_training_sources(source)
    d=eligible.loc[common].copy();y=(d.dead.eq(1)&d.t_exam.le(96)).to_numpy(int)
    if len(d)!=lock['expected_training_n'] or int(y.sum())!=lock['expected_events8']:
        raise ValueError('Frozen training domain changed')
    if (d.dead.eq(0)&d.t_exam.lt(96)).any():raise ValueError('Unexpected early censoring in final training')
    core=d[list(CORE)].to_numpy(float);protein=np.log2(d[list(PROTEINS)].to_numpy(float))
    frozen=fit_final(core,protein,y,d.weight.to_numpy(float),d.cluster.to_numpy())
    frozen['training'].update(weight='WTCYPEX6',strata_n=int(d.stratum.nunique()),
        effective_n=float(d.weight.sum()**2/(d.weight.to_numpy()@d.weight.to_numpy())))
    scientific={str(path.relative_to(ROOT.parent)):adapter.digest(path) for path in sorted(ROOT.parent.rglob('*.py'))}
    artifact=dict(schema=SCHEMA,target=contract['target'],mapping_id=contract['mapping_id'],
        contract_sha256=adapter.digest(ROOT/'CONTRACT.json'),plan_sha256=adapter.digest(ROOT/'FREEZE_PLAN.md'),
        lock_sha256=adapter.digest(ROOT/'LOCK.json'),plan_commit='5cbf90e',
        model_family='eight_year_survey_weighted_ridge_logistic',
        fit_scope='full_training_final_fit_not_new_validation',selection_scope=frozen['selection_scope'],
        models=frozen['models'],training=frozen['training'],
        source_manifest_sha256=adapter.digest(source/'routine29_manifest.json'),
        source_content_sha256=audit['routine_source_content_sha256'],
        source_fingerprints={row['file']:row['sha256'] for row in source_manifest['files']},
        source_policy=source_policy,scientific_python_sha256=scientific,
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
        clinical_use_ready=False,independent_external_validation=False,participant_records_exported=False)
    adapter.validate_artifact(artifact)
    actual=adapter.score(_training_document(d),artifact)
    replay={}
    for arm in ARMS:
        spec=frozen['models'][arm]
        reference=model.predict_ridge(model.transform_features(core,protein,spec['preprocessing']),spec['fit'])
        difference=float(np.max(np.abs(reference-actual[arm])))
        if difference>1e-12:raise ArithmeticError('Frozen adapter numerical replay mismatch')
        replay[arm]=dict(compared_n=len(d),maximum_absolute_probability_difference=difference)
    artifact_bytes=_canonical_bytes(artifact);artifact_sha=hashlib.sha256(artifact_bytes).hexdigest()
    manifest=dict(schema='freeze30-manifest/1',artifact_file='FROZEN_MODELS.json',artifact_sha256=artifact_sha,
        contract_sha256=artifact['contract_sha256'],plan_sha256=artifact['plan_sha256'],lock_sha256=artifact['lock_sha256'],
        source_manifest_sha256=artifact['source_manifest_sha256'],source_content_sha256=artifact['source_content_sha256'],
        source_fingerprints=artifact['source_fingerprints'],source_policy=source_policy,
        scientific_python_sha256=scientific,source_audit=audit,training=frozen['training'],
        numerical_replay=replay,training_performance_exported=False,participant_records_exported=False,
        clinical_use_ready=False)
    out.mkdir(parents=True)
    (out/'FROZEN_MODELS.json').write_bytes(artifact_bytes)
    (out/'FREEZE_MANIFEST.json').write_bytes(_canonical_bytes(manifest))
    for label,value in [('ARTIFACT',artifact),('MANIFEST',manifest)]:
        print('FREEZE30_'+label+'_BEGIN');print(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False));print('FREEZE30_'+label+'_END')
    return artifact,manifest

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run(args.source,args.out)
