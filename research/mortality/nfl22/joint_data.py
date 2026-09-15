"""Audit same-person/same-visit overlap, not fit a seven-protein risk model.

The data holder must approve use and disclosure. User-declared authorization is
not a legal verification. Never export participant identifiers from this audit.
"""
from collections import Counter, defaultdict
import math

PROTEINS=('NT_PROBNP','CYSTATIN_C','GDF15','IL6','TNF_ALPHA','GFAP','NFL')
FORBIDDEN_SUBSTITUTES={'TNF_ALPHA':{'TNFRSF1A','TNFRSF1B','sTNFR1','sTNFR2'},
    'NT_PROBNP':{'NPPB','BNP'},'IL6':{'IL6R','sIL6R'}}


def validate_plan(plan):
    if not isinstance(plan,dict):raise ValueError('Plan must be an object')
    for key in ('cohort','baseline_visit','approval_reference','data_release','dictionary_sha256'):
        if not isinstance(plan.get(key),str) or not plan[key].strip():raise ValueError('Missing '+key)
    if plan.get('data_use_approved') is not True:raise ValueError('Explicit approved-use declaration required')
    required=plan.get('required_markers')
    if (not isinstance(required,list) or not required or any(not isinstance(x,str) or not x for x in required)
        or len(set(required))!=len(required) or not set(PROTEINS)<=set(required)):
        raise ValueError('All seven exact protein identities and explicit baseline fields required')
    methods=plan.get('measurements')
    if not isinstance(methods,dict) or set(methods)!=set(required):raise ValueError('Incomplete assay dictionary')
    for target,row in methods.items():
        if not isinstance(row,dict):raise ValueError('Invalid assay record')
        for key in ('source_variable','source_analyte','matrix','unit','assay'):
            if not isinstance(row.get(key),str) or row[key] in {'','unknown','not_reported'}:raise ValueError('Unresolved '+target+' '+key)
        if row['source_analyte'] in FORBIDDEN_SUBSTITUTES.get(target,set()):raise ValueError('Different molecular target: '+target)
        if row.get('identity_verified') is not True:raise ValueError('Assay identity review required')
        # A relative proteomics number is not an absolute concentration.
        if row.get('scale') in {'NPX','RFU'} and row['unit']!=row['scale']:raise ValueError('Undeclared relative-to-absolute conversion')
    return required


def audit_overlap(rows,plan):
    required=validate_plan(plan);cells={};visits=defaultdict(set); counts=Counter(); qualifiers=Counter()
    for row in rows:
        if not isinstance(row,dict):raise ValueError('Each measurement must be an object')
        if row.get('cohort')!=plan['cohort']:raise ValueError('Cross-cohort joining forbidden')
        pid=row.get('participant_id');visit=row.get('visit');marker=row.get('marker')
        if not isinstance(pid,str) or not pid or not isinstance(visit,str) or not visit:raise ValueError('Explicit participant/visit keys required')
        if marker not in required:raise ValueError('Marker outside declared plan')
        if visit!=plan['baseline_visit']:raise ValueError('Different visit; no nearest-person/visit substitution')
        key=(pid,visit,marker)
        if key in cells:raise ValueError('Duplicate assay; adjudicate replicates before joining')
        meta=plan['measurements'][marker]
        if any(row.get(k)!=meta[k] for k in ('matrix','unit','assay')):raise ValueError('Assay metadata mismatch')
        qualifier=row.get('qualifier');value=row.get('value')
        if qualifier not in {'quantified','below_LOQ','above_LOQ','missing'}:raise ValueError('Explicit quantification status required')
        if qualifier=='quantified' and (type(value) not in (int,float) or not math.isfinite(value)):
            raise ValueError('Finite measured value required; no string-to-zero conversion')
        cells[key]=qualifier;qualifiers[qualifier]+=1
        # Missing/censored concentrations remain present but do not qualify as quantified overlap.
        if qualifier=='quantified':visits[pid].add(marker);counts[marker]+=1
        else:visits[pid]
    complete=sum(set(required)<=markers for markers in visits.values())
    return {'cohort':plan['cohort'],'baseline_visit':plan['baseline_visit'],'records':len(cells),
        'participants_seen':len(visits),'complete_quantified_profiles':complete,
        'quantified_per_marker':dict(counts),'qualifier_counts':dict(qualifiers),
        'identifiers_exported':False,'data_holder_disclosure_review_required':True,
        'seven_protein_model_fitted':False,'clinical_use_ready':False,
        'authorization_legally_verified':False}
