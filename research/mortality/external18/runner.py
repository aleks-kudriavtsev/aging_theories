"""Owner-site external validation of frozen models, aggregate output only.

The scope attestation is NOT legal verification. No country is changed to US;
only the frozen numerical core is used. No request to a website is performed.
Independent censoring within the analysed domain is a declared assumption.
"""
from __future__ import annotations
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import numpy as np
from ..studio16.registry import ARTIFACTS, model_spec
from ..replay_compact14 import score_engineered
from ..workbench.runtime import engineer, FIELDS, finite, strict_json
from ..ipcw14 import outcome_weights, binary_scores
from ..validate_temporal08 import diagnostics

MODELS = ('clinical3', 'compact4')
MEASURES = ('bmi','sbp','hba1c','creatinine','uacr','albumin')
CLINICAL = ('smoking','bp_treatment','diabetes_history','cvd_history','cancer_history')
ROW_FIELDS = {'record_id','age_years','sex','clinical','measurements',
              'death','followup_years','weight','stratum','cluster'}
ATTESTATIONS = ('owner_permission_for_model_validation','ethics_scope_reviewed',
 'consent_or_other_basis_reviewed','local_execution_permitted','mortality_linkage_permitted',
 'baseline_mapping_independently_reviewed','no_development_overlap_reviewed',
 'analysis_plan_frozen_before_outcomes','output_review_required')


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def preflight(config: dict, protocol: Path, dictionary: Path, approval: Path | None = None) -> dict:
    """Metadata-only checks; never opens participant data. Pending scope blocks."""
    errors=[]
    if not isinstance(config,dict):return {'ready':False,'errors':['config_not_object']}
    if config.get('schema')!='external18.1':errors.append('schema')
    if config.get('country') not in ('RU','DE'):errors.append('target_country')
    origin=config.get('data_origin')
    if origin not in ('synthetic','authorized_cohort'):errors.append('data_origin')
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,64}',str(config.get('study_id',''))):errors.append('study_id')
    if config.get('endpoint')!='all_cause_mortality':errors.append('endpoint')
    if config.get('horizon_years')!=3 or type(config.get('horizon_years')) is not int:errors.append('primary_horizon_3y')
    if config.get('time_origin')!='baseline_measurements' or config.get('delayed_entry') is not False:errors.append('time_origin_or_delayed_entry')
    design=config.get('design')
    if design not in ('iid','stratified_psu'):errors.append('sampling_design')
    if config.get('censoring_assumption')!='independent_within_domain':errors.append('censoring_assumption')
    nboot=config.get('bootstrap_replicates')
    if type(nboot) is not int or nboot<100 or nboot>5000:errors.append('bootstrap_range')
    if origin=='authorized_cohort' and nboot!=1000:errors.append('real_study_requires_1000_replicates')
    if type(config.get('seed')) is not int:errors.append('seed')
    minimum=config.get('minimum_cell_count')
    if type(minimum) is not int or minimum<10:errors.append('disclosure_policy_minimum')
    for name,path in [('protocol',protocol),('dictionary',dictionary)]:
        if not Path(path).is_file() or digest(path)!=config.get(name+'_sha256'):errors.append(name+'_integrity')
    if config.get('model_sha256')!=ARTIFACTS['compact'][1] or config.get('clinical_sha256')!=ARTIFACTS['full'][1]:errors.append('model_version')
    if not re.fullmatch(r'[0-9a-f]{64}',str(config.get('data_sha256',''))):errors.append('data_version_required')
    if origin=='authorized_cohort':
        if approval is None or not Path(approval).is_file() or digest(approval)!=config.get('approval_sha256'):errors.append('approval_document_missing_or_changed')
        a=config.get('attestations',{})
        if not isinstance(a,dict):a={}
        errors += ['pending_'+k for k in ATTESTATIONS if a.get(k) is not True]
        for key in ('owner_review_reference','mapping_review_reference','execution_environment_reference'):
            if not isinstance(config.get(key),str) or not config[key].strip():errors.append('missing_'+key)
    if not errors:
        for name in MODELS:model_spec(name)
    return {'ready':not errors,'errors':errors,'country':config.get('country'),
        'scope_check':'operator_attestation_and_file_integrity_NOT_legal_verification',
        'participant_data_opened':False,'clinical_use_ready':False}


def measurement(value, key, definition):
    if value is None:return None
    if not finite(value):raise ValueError('nonfinite_or_censored_measurement')
    sp=FIELDS[key]
    if definition.get('matrix')!=sp['matrix']:raise ValueError('specimen_matrix_mismatch')
    if not definition.get('method_id'):raise ValueError('method_id_missing')
    u=definition.get('unit')
    if key=='hba1c':
        if u=='mmol/mol' and definition.get('standardization')=='IFCC_traceable':
            value=.09148*value+2.152
        elif u=='%' and definition.get('standardization')=='NGSP_traceable':pass
        else:raise ValueError('HbA1c_standardization_not_confirmed')
    else:
        if u not in sp['conversions']:raise ValueError('unit_mismatch')
        value*=sp['conversions'][u]
    if key=='creatinine' and definition.get('standardization')!='IDMS_standardized':raise ValueError('creatinine_standardization')
    if not finite(value) or not sp['software_limits'][0]<=value<=sp['software_limits'][1]:raise ValueError('outside_numerical_support')
    return float(value)


def read_cohort(path: Path, dictionary: dict, design: str):
    if set(dictionary.get('measurements',{}))!=set(MEASURES):raise ValueError('dictionary_measurements')
    rows=[];seen=set()
    with Path(path).open('rb') as stream:
        for line in stream:
            r=strict_json(line)
            if not isinstance(r,dict) or set(r)!=ROW_FIELDS:raise ValueError('invalid_columns_or_identifier_fields')
            rid=r['record_id']
            if not isinstance(rid,str) or not re.fullmatch('[A-Za-z0-9_-]{1,64}',rid) or rid in seen:raise ValueError('duplicate_or_invalid_local_record_id')
            seen.add(rid)
            if not finite(r['age_years']) or not 0<=r['age_years']<=120 or r['sex'] not in ('female','male'):raise ValueError('invalid_age_or_recorded_sex')
            if type(r['death']) is not bool or not finite(r['followup_years']) or r['followup_years']<=0:raise ValueError('death_or_time_invalid_no_zero_time_substitution')
            if not finite(r['weight']) or r['weight']<=0:raise ValueError('invalid_weight')
            if design=='iid' and (r['weight']!=1 or r['stratum'] is not None or r['cluster'] is not None):raise ValueError('iid_design_mismatch')
            if design=='stratified_psu' and any(not isinstance(r[k],str) or not r[k] for k in ('stratum','cluster')):raise ValueError('missing_design_field')
            if not isinstance(r['clinical'],dict) or set(r['clinical'])!=set(CLINICAL):raise ValueError('clinical_fields')
            if r['clinical']['smoking'] not in ('never','former','current',None):raise ValueError('smoking_code')
            if any(r['clinical'][k] is not None and type(r['clinical'][k]) is not bool for k in CLINICAL[1:]):raise ValueError('clinical_boolean')
            if not isinstance(r['measurements'],dict) or set(r['measurements'])!=set(MEASURES):raise ValueError('measurement_columns')
            values={k:measurement(r['measurements'][k],k,dictionary['measurements'][k]) for k in MEASURES}
            r['measurements']=values
            r['complete']=all(v is not None for v in values.values()) and all(v is not None for v in r['clinical'].values())
            rows.append(r)
    if not rows:raise ValueError('empty_cohort')
    if design=='stratified_psu':
        strata={}
        for r in rows:strata.setdefault(r['stratum'],set()).add(r['cluster'])
        if any(len(x)<2 for x in strata.values()):raise ValueError('singleton_stratum_requires_protocol_amendment')
    return rows


def resample_weights(rows,design,rng):
    n=len(rows);w=np.array([r['weight'] for r in rows],float)
    if design=='iid':return np.bincount(rng.integers(0,n,n),minlength=n).astype(float)
    groups={};mult=np.zeros(n)
    for i,r in enumerate(rows):groups.setdefault(r['stratum'],{}).setdefault(r['cluster'],[]).append(i)
    for units in groups.values():
        ids=list(units);m=len(ids)
        for index,count in Counter(rng.integers(0,m,m-1)).items():mult[units[ids[index]]]=count*m/(m-1)
    return w*mult


def evaluate(config,rows):
    eligible=[r for r in rows if 40<=r['age_years']<=79]
    if not eligible:raise ValueError('empty_eligible_age_domain')
    mask=np.array([40<=r['age_years']<=79 and r['complete'] for r in rows],bool);d=[r for r in eligible if r['complete']]
    minimum=config['minimum_cell_count']
    if not d:raise ValueError('empty_common_complete_domain')
    y=np.array([r['death'] for r in d],float);t=np.array([r['followup_years'] for r in d]);w=np.array([r['weight'] for r in d])
    # Sentinel10 means only "death" to this all-cause scorer, NOT an inferred cause.
    c=np.where(y==1,10,np.nan);h=config['horizon_years'];obs=outcome_weights(y,c,t,w,h)
    events=int(((y==1)&(t<=h)).sum())
    if len(d)<2*minimum or events<minimum or int(((obs['target']==0)&obs['known']).sum())<minimum:
        return {'status':'aggregate_output_suppressed','reason':'minimum_event_or_known_survivor_cell',
                'clinical_use_ready':False,'rows_exported':False}
    models={k:model_spec(k) for k in MODELS};pred={k:[] for k in MODELS};clipping=Counter()
    for r in d:
        f=engineer(r['age_years'],r['sex'],r['clinical'],r['measurements'])
        for k,m in models.items():
            p,trace=score_engineered(f,m,[h]);pred[k].append(p[0]['all_cause_death_probability'])
            for x in trace:
                if x['clipped']:clipping[k+':'+x['feature']]+=1
    pred={k:np.array(v) for k,v in pred.items()}
    metrics={}
    for k,p in pred.items():
        m=binary_scores(obs['target'],p,w,obs)
        try:
            z=obs['effective_weights']>0
            m.update(diagnostics((obs['target']>0)[z],p[z],obs['effective_weights'][z]))
        except ValueError:m['calibration_diagnostic_status']='failed_no_model_refit'
        metrics[k]=m
    draws=[];rng=np.random.default_rng(config['seed'])
    for _ in range(config['bootstrap_replicates']):
        # Resample the full supplied sample BEFORE the age/complete-profile domain.
        wr=resample_weights(rows,config['design'],rng)[mask]
        try:
            o=outcome_weights(y,c,t,wr,h)
            scores={k:binary_scores(o['target'],p,wr,o) for k,p in pred.items()}
            draws.append({q:scores['compact4'][q]-scores['clinical3'][q]
                if scores['compact4'][q] is not None and scores['clinical3'][q] is not None else None for q in ('brier','auc','bias_pp')})
        except (ValueError,ArithmeticError):draws.append(dict.fromkeys(('brier','auc','bias_pp')))
    intervals={}
    for q in ('brier','auc','bias_pp'):
        good=[r[q] for r in draws if r[q] is not None]
        bounds=np.quantile(good,[.025,.975]).tolist() if good else [None,None]
        intervals[q]={'difference':metrics['compact4'][q]-metrics['clinical3'][q] if metrics['compact4'][q] is not None else None,
            'lower':bounds[0],'upper':bounds[1],'valid':len(good),'undefined':len(draws)-len(good),
            'interval_scope':'conditional_on_frozen_models_censor_distribution_reestimated'}
    # No small-group tables are exported. These require the holder's disclosure review.
    return {'status':'synthetic_pipeline_check' if config['data_origin']=='synthetic' else 'external_analysis_pending_owner_review',
        'clinical_use_ready':False,'country':config['country'],'endpoint':'all_cause_mortality','horizon_years':h,
        'input_domain':'same_complete_compact_profile_for_both_models',
        'flow':{'received':len(rows),'age_eligible':len(eligible),'complete_common':len(d),'deaths_by_horizon':events},
        'metrics':metrics,'compact_minus_clinical':intervals,
        'clipping_counts_suppressed':True,'clipping_present':bool(clipping),
        'country_substituted_with_US':False,'risk_coefficients_refitted':False,'predictions_exported':False,
        'cause_specific_predictions_exported':False,'recalibration_applied':False,'output_requires_owner_review':True,'statistical_disclosure_control_complete':False,
        'assumptions':['independent_censoring_within_complete_domain','sampling_design_as_declared'],
        'uncertainty_excludes':['model_training','sample_selection','measurement_method_uncertainty'],
        'synthetic_data_are_not_external_validation':config['data_origin']=='synthetic'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['preflight','evaluate'])
    for k in ('config','protocol','dictionary'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--approval',type=Path);p.add_argument('--data',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();config=strict_json(a.config.read_bytes())
    if a.out.exists():raise SystemExit('Output exists')
    status=preflight(config,a.protocol,a.dictionary,a.approval)
    if a.command=='preflight' or not status['ready']:result=status
    else:
        if a.data is None or digest(a.data)!=config['data_sha256']:raise SystemExit('Data version mismatch')
        dictionary=strict_json(a.dictionary.read_bytes())
        rows=read_cohort(a.data,dictionary,config['design']);result=evaluate(config,rows)
        result['provenance']={k:config[k] for k in ('model_sha256','clinical_sha256','protocol_sha256','dictionary_sha256')}
        result['source_data_hash_retained_locally']=True
    with a.out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps({'status':result.get('status','preflight_ready' if status['ready'] else 'blocked'),
        'clinical_use_ready':False,'participant_rows_exported':False}))

if __name__=='__main__':main()
