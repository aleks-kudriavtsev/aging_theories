"""Add source-checked renal/routine measurements without modifying cohort28."""
from __future__ import annotations
import hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from ..binary28 import cohort as old
from ..cohort27.gate import require_qualification
from ..cohort27.importer import join_one_to_one
from ..nfl22.analyze import sha,independent_mortality_reader
from .retrieve import ROOT,STATIC,DYNAMIC,check_study,resolve_additional,prior
from .reader import read_checked
from .features import engineer,CORE,strict_metadata_mask


def verify_sources(source):
    check_study();require_qualification(json.loads((ROOT.parent/'cohort27/CONTRACT.json').read_text()))
    source=Path(source);doc=json.loads((source/'routine29_manifest.json').read_text());base=json.loads((source/'manifest.json').read_text())
    if (doc.get('plan_commit')!='7800458' or doc.get('lock_sha256')!=sha(ROOT/'LOCK.json')
            or doc.get('base_manifest_sha256')!=sha(source/'manifest.json')):raise ValueError('New sources not bound to routine29 study')
    expected={r['file'] for r in base['files']}|set(STATIC)|set(DYNAMIC);rows={}
    for row in doc['files']:
        name=row['file']
        if name in rows or name not in expected or Path(name).name!=name:raise ValueError('Unexpected/duplicate routine source')
        if row.get('status')!='retrieved' or not prior.allowed_url(row.get('url','')) or not prior.allowed_url(row.get('resolved_url','')):raise ValueError('Unverified public source')
        rows[name]=row
    if set(rows)!=expected:raise ValueError('Missing routine input')
    for name,row in rows.items():
        path=source/name
        if not path.is_file() or path.stat().st_size!=row['bytes'] or sha(path)!=row['sha256']:raise ValueError('Routine source changed: '+name)
    urls={**STATIC,**resolve_additional((source/'index.htm').read_text(encoding='utf-8-sig'))}
    if any(rows[name]['url']!=url for name,url in urls.items()):raise ValueError('Wrong LAB/EXAM release')
    for row in base['files']:
        if rows[row['file']]!=row:raise ValueError('Original source row changed')
    return doc


def load(source):
    source=Path(source);manifest=verify_sources(source)
    eligible,oldcommon,audit,_=old.load(source)
    if int(oldcommon.sum())!=3857:raise ValueError('Old28 common domain changed')
    tables={};readers={}
    mortality_ids=set(independent_mortality_reader(source/'mortality.dat').SEQN)
    for kind,specfile in [('lab','LAB_FIELDS.json'),('exam','EXAM_FIELDS.json')]:
        spec=json.loads((ROOT/specfile).read_text());fields=spec['positions_1_based_inclusive']
        width=spec.get('nominal_content_width',spec.get('nominal_record_content_width'))
        count=spec.get('official_source_rows',spec.get('expected_record_count'))
        tables[kind],readers[kind]=read_checked(source/(kind+'.dat'),source/(kind+'.sas'),fields,width,count)
        if not set(tables[kind].SEQN)<=mortality_ids:raise ValueError('Foreign survey IDs in routine source')
        eligible=join_one_to_one(eligible,tables[kind])
    eligible=engineer(eligible)
    protein=np.isfinite(eligible[list(old.PROTEINS)]).all(axis=1).to_numpy()&eligible[list(old.PROTEINS)].gt(0).all(axis=1).to_numpy()
    routine=np.isfinite(eligible[list(CORE)]).all(axis=1).to_numpy()
    common=protein&routine
    strict=common&strict_metadata_mask(eligible)
    if (common&~oldcommon).any():raise ValueError('Routine comparison must be a subset of old28 common domain')
    reasons={'legacy_clinical':~np.isfinite(eligible[list(old.CORE)]).all(axis=1).to_numpy(),
             'protein':~protein,'BMI':~np.isfinite(eligible.bmi).to_numpy(),
             'HbA1c':~np.isfinite(eligible.hba1c).to_numpy(),
             'eGFR':~np.isfinite(eligible.egfr).to_numpy(),'UACR':~np.isfinite(eligible.log_uacr).to_numpy()}
    patterns=['|'.join(k for k,v in reasons.items() if v[i]) or 'complete' for i in range(len(eligible))]
    common_flag_counts={name:{str(k):int(v) for k,v in eligible.loc[common,name].fillna(-1).value_counts().sort_index().items()}
                        for name in ('GHPMETH','BMPWTFLG','BMPHTFLG')}
    audit.update(exclusion_reason_counts={k:int(v.sum()) for k,v in reasons.items()},
                 exclusion_reason_patterns=dict(Counter(patterns)),
                 exclusion_patterns_in_old28_common=dict(Counter(p for p,keep in zip(patterns,oldcommon) if keep)),
                 source_flag_counts_in_common=common_flag_counts,
                 common_urine_albumin_point4_n=int(eligible.loc[common,'UBP'].eq(.4).sum()))
    audit.update(old28_common_n=int(oldcommon.sum()),old28_common_events8=int((eligible.loc[oldcommon,'dead'].eq(1)&eligible.loc[oldcommon,'t_exam'].le(96)).sum()),
        routine_complete_n=int(routine.sum()),routine_common_n=int(common.sum()),strict_metadata_domain_n=int(strict.sum()),
        old28_cases_removed_by_routine_n=int((oldcommon&~common).sum()),
        missing_routine={f:int((~np.isfinite(eligible[f])).sum()) for f in CORE},routine_readers=readers,
        urine_creatinine_suspect_replacement_n=int(eligible.URP.eq(7.9).sum()),urine_creatinine_below10_n=int(eligible.URP.lt(10).sum()),
        urine_albumin_replacement_point4_n=int(eligible.UBP.eq(.4).sum()),
        glyco_method_counts={str(k):int(v) for k,v in eligible.GHPMETH.fillna(-1).value_counts().sort_index().items()},
        bmi_weight_flag_counts={str(k):int(v) for k,v in eligible.BMPWTFLG.fillna(-1).value_counts().sort_index().items()},
        bmi_height_flag_counts={str(k):int(v) for k,v in eligible.BMPHTFLG.fillna(-1).value_counts().sort_index().items()},
        routine_source_manifest_sha256=sha(source/'routine29_manifest.json'),
        routine_source_fingerprints={r['file']:r['sha256'] for r in manifest['files']},
        routine_source_content_sha256=hashlib.sha256(json.dumps({r['file']:r['sha256'] for r in manifest['files']},sort_keys=True,separators=(',',':')).encode()).hexdigest())
    return eligible,common,strict,audit,manifest
