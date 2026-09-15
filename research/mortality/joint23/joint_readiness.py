"""Auditable source readiness; metadata presence is never individual-data access.

Use before importing protected data. The input is a holder-completed metadata
contract, not participant records. No coefficient is trained by this module.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

PROTEINS={'NT_PROBNP':'P16860','CST3':'P01034','GDF15':'Q99988',
          'IL6':'P05231','TNF':'P01375','GFAP':'P14136','NEFL':'P07196'}
BASELINE=('total_cholesterol','hdl','hba1c','creatinine','uacr','albumin',
          'rdw_cv','wbc','alp','ggt','bilirubin_total','platelets','nlr')
CLINICAL=('age','sex','smoking','bmi','sbp','bp_treatment',
          'diabetes_history','cvd_history','cancer_history')
REQUIRED=tuple(PROTEINS)+BASELINE+CLINICAL


def present(x):return isinstance(x,str) and bool(x.strip()) and x.strip().lower() not in {'unknown','none','not_reported','pending'}

def audit(contract):
    if not isinstance(contract,dict):raise ValueError('Expected a metadata object')
    records=contract.get('variables',{})
    if not isinstance(records,dict):raise ValueError('Variables must be a mapping')
    missing=[];invalid=[];unresolved=[];methods=[]
    for key in REQUIRED:
        r=records.get(key)
        if not isinstance(r,dict):missing.append(key);continue
        if r.get('released') is not True:unresolved.append(key+':release_unconfirmed')
        for f in ('variable','source_release','baseline_window','source_locator'):
            if not present(r.get(f)):unresolved.append(key+':missing_'+f)
        if key in PROTEINS:
            if r.get('uniprot')!=PROTEINS[key]:invalid.append(key+':wrong_or_unknown_protein_identity')
            if key=='NT_PROBNP' and r.get('measured_form')!='NT-proBNP':invalid.append('NT_PROBNP:fragment_not_confirmed')
        if key in PROTEINS or key in BASELINE:
            for f in ('matrix','assay','scale','quantification_policy'):
                if not present(r.get(f)):unresolved.append(key+':missing_'+f)
            if r.get('scale') in ('NPX','RFU'):
                methods.append(key+':relative_abundance_requires_own_training_or_paired_bridge')
    # Metadata documents do not imply a shared participant sample.
    linkage=contract.get('linkage',{})
    if not isinstance(linkage,dict):raise ValueError('Linkage must be a mapping')
    gates={
        'access_approved':contract.get('access_approved') is True,
        'raw_records_received':contract.get('raw_records_received') is True,
        'same_person_verified':linkage.get('same_person_verified') is True,
        'same_baseline_verified':linkage.get('same_baseline_verified') is True,
        'outcome_after_baseline_verified':linkage.get('outcome_after_baseline_verified') is True,
        'mortality_followup_verified':linkage.get('mortality_followup_verified') is True,
        'clinical_baseline_complete':not any(x in CLINICAL for x in missing),
        'weights_and_sampling_documented':linkage.get('weights_and_sampling_documented') is True,
        'assay_transport_adjudicated':linkage.get('assay_transport_adjudicated') is True,
    }
    n=linkage.get('joint_complete_participants');events=linkage.get('deaths_after_baseline')
    gates['joint_count_verified']=type(n) is int and n>0
    gates['event_count_verified']=type(events) is int and events>0 and type(n) is int and events<=n
    if not present(contract.get('permission_locator')):gates['access_approved']=False
    ready=all(gates.values()) and not missing and not invalid and not unresolved
    # Even a complete candidate requires a separate analysis plan/sample-size review.
    return {'status':'ready_for_design_review' if ready else 'blocked_before_training',
      'missing_components':missing,'identity_errors':invalid,'unresolved_metadata':unresolved,
      'method_notes':methods,'gates':gates,'joint_participant_count':n if gates['joint_count_verified'] else None,
      'coefficient_training_authorized_by_this_audit':False,'clinical_use_ready':False,
      'existing_NHANES_model_applicable':False,'protected_records_processed':False}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():raise ValueError('New output file required')
    if a.input.stat().st_size>1_000_000:raise ValueError('Metadata contract too large')
    data=json.loads(a.input.read_text(encoding='utf-8'));r=audit(data)
    with a.out.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2,allow_nan=False)
    print(r['status'])
