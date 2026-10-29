"""Source-pinned matched domains for the prespecified reused-cohort B2M study."""
from __future__ import annotations
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from ..xpt_checked import read_xpt_checked
from ..validate_temporal08 import independent_mortality_reader, audit_mortality_reader
from .retrieve import ROOT, PLAN_COMMIT, check_study, allowed_url, sha, static_urls, resolve_sscard, unavailable_documentation
from .features import engineer, CORE

PROTEINS = ('SSCYST', 'SSB2M')


def verify_sources(source):
    check_study(); source=Path(source)
    manifest=json.loads((source/'manifest.json').read_text())
    if (manifest.get('schema')!='b2m31-retrieval/1' or manifest.get('plan_commit')!=PLAN_COMMIT
            or manifest.get('lock_sha256')!=sha(ROOT/'LOCK.json')
            or manifest.get('documentation_not_retrieved',[])!=unavailable_documentation()):
        raise ValueError('Source intake not bound to fixed protocol')
    spec=json.loads((ROOT/'SOURCE_SPEC.json').read_text()); rows={}
    for row in manifest['files']:
        name=row['file']
        if (name in rows or Path(name).name!=name or row.get('status')!='retrieved'
                or not allowed_url(row.get('url','')) or not allowed_url(row.get('resolved_url',''))):
            raise ValueError('Duplicate, unsafe or failed source')
        path=source/name
        if not path.is_file() or path.stat().st_size!=row['bytes'] or sha(path)!=row['sha256']:
            raise ValueError('Source bytes changed before parsing')
        rows[name]=row
    if 'index.htm' not in rows:raise ValueError('Official catalogue required')
    expected={**static_urls(),'index.htm':spec['sscard_index_url'],
        '1999_SSCARD_A.xpt':resolve_sscard((source/'index.htm').read_text(encoding='utf-8-sig'),rows['index.htm']['resolved_url'])}
    if set(rows)!=set(expected) or any(rows[name]['url']!=url for name,url in expected.items()):
        raise ValueError('Unexpected official release file or URL')
    historical={name:digest for digest,name in (s.split() for s in (ROOT.parent/'nhanes04/SOURCE_DATA_SHA256.txt').read_text().splitlines() if s.strip())}
    data=[name for name in rows if name.endswith(('.dat','.xpt'))]
    if any(name not in historical or rows[name]['sha256']!=historical[name] for name in data):
        raise ValueError('Historical public data release fingerprint changed')
    return manifest,dict(files_verified=len(rows),historical_data_hashes_verified=len(data),
        documentation_not_retrieved=unavailable_documentation(),
        source_manifest_sha256=sha(source/'manifest.json'),
        source_fingerprints={name:row['sha256'] for name,row in rows.items()},
        source_content_sha256=hashlib.sha256(json.dumps({name:row['sha256'] for name,row in rows.items()},sort_keys=True,separators=(',',':')).encode()).hexdigest())


def six_year_weights(frame):
    if not frame.SDDSRVYR.isin([1,2,3]).all():raise ValueError('Unknown survey cycle for pooled weight')
    w=pd.Series(np.where(frame.SDDSRVYR.isin([1,2]),(2/3)*frame.WTSSCB4Y,(1/3)*frame.WTSSCB2Y),index=frame.index)
    return w.where(np.isfinite(w)&w.gt(0))


def validate_sscard(frame, spec):
    if len(frame)!=spec['expected_rows_full_release'] or not set(spec['required_fields'])<=set(frame):
        raise ValueError('Wrong corrected SSCARD release schema/count')
    if frame.SEQN.isna().any() or frame.SEQN.duplicated().any():raise ValueError('Invalid SSCARD keys')
    for flag,expected in spec['metadata_marginal_flag_counts'].items():
        if not frame[flag].dropna().isin([0,1,2]).all():raise ValueError('Unknown protein detection flag')
        actual={str(k):int(frame[flag].eq(int(k)).sum()) for k in ('0','1','2')}
        actual['missing']=int(frame[flag].isna().sum())
        if actual!=expected:raise ValueError('Published protein flag counts changed')
    for col,expected in spec['metadata_marginal_pristine_counts'].items():
        actual={str(k):int(frame[col].eq(int(k)).sum()) for k in ('0','1')};actual['missing']=int(frame[col].isna().sum())
        if actual!=expected:raise ValueError('Published pristine metadata changed')
    for col,expected in spec['metadata_weight_counts'].items():
        actual=dict(positive=int((np.isfinite(frame[col])&frame[col].gt(0)).sum()),zero=int(frame[col].eq(0).sum()),missing=int(frame[col].isna().sum()))
        if actual!=expected:raise ValueError('Published corrected weight counts changed')
    for col,flag in [('SSCYST','SSCYSTL'),('SSB2M','SSB2ML')]:
        if frame[col].isna().sum()!=frame[flag].isna().sum() or (frame[col].notna()!=frame[flag].notna()).any():
            raise ValueError('Protein result/detection flag missingness mismatch')
        if not (np.isfinite(frame.loc[frame[col].notna(),col])&frame.loc[frame[col].notna(),col].gt(0)).all():
            raise ValueError('Invalid released protein value')


def validate_design_inventory(full, eligible):
    full_pairs=set(zip(full.SDMVSTRA.astype(int),full.SDMVPSU.astype(int)))
    target_pairs=set(zip(eligible.SDMVSTRA.astype(int),eligible.SDMVPSU.astype(int)))
    if full_pairs!=target_pairs:raise ValueError('Target eligibility removed a full-DEMO masked PSU; separate design protocol required')
    counts=eligible.groupby('SDMVSTRA').SDMVPSU.nunique()
    if (counts<2).any():raise ValueError('Singleton eligible masked stratum')
    return dict(full_DEMO_masked_psu_n=len(full_pairs),eligible_masked_psu_n=len(target_pairs),
                full_DEMO_design_inventory_preserved=True,
                full_DEMO_masked_design={str(st):sorted(p for s,p in full_pairs if s==st) for st in sorted({s for s,p in full_pairs})},
                eligible_masked_design={str(st):sorted(p for s,p in target_pairs if s==st) for st in sorted({s for s,p in target_pairs})},
                eligible_psu_count_distribution={str(k):int(v) for k,v in counts.value_counts().sort_index().items()})


def domain_masks(d):
    routine=np.isfinite(d[list(CORE)]).all(axis=1).to_numpy()
    proteins=(np.isfinite(d[list(PROTEINS)]).all(axis=1)&d[list(PROTEINS)].gt(0).all(axis=1)).to_numpy()
    quantified=d.SSCYSTL.eq(0).to_numpy()&d.SSB2ML.eq(0).to_numpy()
    reported=d.SSCYSTL.isin([0,1]).to_numpy()&d.SSB2ML.isin([0,1]).to_numpy()
    return routine&proteins&quantified,routine&proteins&reported


def checked_mortality_audit(path):
    try:
        return audit_mortality_reader(path)
    except (AssertionError, ValueError, ArithmeticError):
        raise ValueError('Independent mortality readers disagree; row diagnostics suppressed') from None


def load(source):
    source=Path(source);manifest,audit=verify_sources(source);spec=json.loads((ROOT/'SOURCE_SPEC.json').read_text())
    zero_audits=[];reader_audits=[]
    def read(name,fields):
        d=read_xpt_checked(source/name)
        if not set(fields)<=set(d) or any(not pd.api.types.is_numeric_dtype(d[c]) for c in fields):
            raise ValueError('Missing or nonnumeric source measurement field')
        if not (d.SEQN.gt(0)&np.isfinite(d.SEQN)&(d.SEQN%1).eq(0)).all():
            raise ValueError('Invalid source participant key representation')
        zero_audits.append(d.attrs['xpt_zero_audit']);return d[fields].copy()
    sscard=read('1999_SSCARD_A.xpt',spec['sscard']['required_fields']);validate_sscard(sscard,spec['sscard'])
    frames=[];full_demos=[];flow=[];all_ids=set()
    for c in spec['cycles']:
        yr=c['year'];components=c['components'];demo_name=next(n for n in components if n.startswith('DEMO'))
        d=read(f'{yr}_{demo_name}.xpt',components[demo_name]['required_fields']);ids=set(d.SEQN)
        if len(d)!=c['demo_rows'] or not d.SDDSRVYR.eq(c['cycle_code']).all():raise ValueError('Wrong demographic release/cycle')
        if all_ids&ids:raise ValueError('Repeated participant IDs across cycles')
        all_ids|=ids
        lo,hi=c['stratum_range']
        if (not d.SDMVSTRA.between(lo,hi).all() or not d.SDMVPSU.isin(c['allowed_psu']).all()
                or not (d.SDMVSTRA%1).eq(0).all()):raise ValueError('Invalid masked design variables')
        full_demos.append(d[['SDMVSTRA','SDMVPSU']].copy())
        for name,row in components.items():
            if name==demo_name:continue
            raw=read(f'{yr}_{name}.xpt',row['required_fields'])
            if not set(raw.SEQN)<=ids:raise ValueError('Foreign covariate source keys')
            d=d.merge(raw,on='SEQN',how='left',validate='1:1')
        mort=independent_mortality_reader(source/f'{yr}_mortality.dat')
        if set(mort.SEQN)!=ids:raise ValueError('Mortality/demographic key mismatch')
        reader_audits.append(dict(cycle=yr,**checked_mortality_audit(source/f'{yr}_mortality.dat')))
        d=d.merge(mort,on='SEQN',how='left',validate='1:1').merge(sscard,on='SEQN',how='left',validate='1:1')
        d['cycle']=yr;d['raw_creatinine']=d[c['creatinine_variable']];d['weight']=six_year_weights(d);d['weight6y']=d.weight
        age=d.RIDAGEYR.between(40,75);exam=age&d.RIDSTATR.eq(2)
        linked=exam&d.elig.eq(1)&d.dead.isin([0,1])&np.isfinite(d.t_exam)&d.t_exam.ge(0)&d.RIAGENDR.isin([1,2])
        eligible=linked&np.isfinite(d.weight)&d.weight.gt(0)
        flow.append(dict(cycle=yr,survey_n=len(d),age40_75_n=int(age.sum()),MEC_target_n=int(exam.sum()),
            linked_target_n=int(linked.sum()),eligible_weighted_subsample_n=int(eligible.sum()),
            excluded_no_corrected_subsample_weight_n=int((linked&~eligible).sum())))
        frames.append(d.loc[eligible].copy())
    if not set(sscard.SEQN)<=all_ids:raise ValueError('SSCARD contains unmatched survey participants')
    d=pd.concat(frames,ignore_index=True)
    if d.SEQN.duplicated().any():raise ValueError('Repeated participants in pooled target')
    if (d.dead.eq(0)&d.t_exam.lt(96)).any():raise ValueError('Eight-year survivor followup gate failed before comparison')
    design=validate_design_inventory(pd.concat(full_demos,ignore_index=True),d)
    d=engineer(d);d['stratum']='NHANES:'+d.SDMVSTRA.astype(int).astype(str)
    d['cluster']=d.SDMVSTRA.astype(int).astype(str)+':'+d.SDMVPSU.astype(int).astype(str)
    common,reported=domain_masks(d)
    reasons={'routine':~np.isfinite(d[list(CORE)]).all(axis=1).to_numpy(),
        'protein_missing_nonpositive':~(np.isfinite(d[list(PROTEINS)]).all(axis=1)&d[list(PROTEINS)].gt(0).all(axis=1)).to_numpy(),
        'Cys_not_measurable_flag0':~d.SSCYSTL.eq(0).to_numpy(),
        'B2M_not_measurable_flag0':~d.SSB2ML.eq(0).to_numpy()}
    pattern=['|'.join(k for k,v in reasons.items() if v[i]) or 'complete' for i in range(len(d))]
    def counts(mask):
        s=d.loc[mask];return dict(n=len(s),events8=int((s.dead.eq(1)&s.t_exam.le(96)).sum()),
            cycles={str(k):int(v) for k,v in s.cycle.value_counts().sort_index().items()},
            weight_sum=float(s.weight.sum()),effective_n=float(s.weight.sum()**2/(s.weight@s.weight)) if len(s) else None,
            pristine_counts={str(k):int(v) for k,v in s.SSPRIS.value_counts(dropna=False).sort_index().items()},
            flag_counts={f:{str(k):int(v) for k,v in s[f].value_counts(dropna=False).sort_index().items()} for f in ('SSCYSTL','SSB2ML')})
    audit.update(eligible_n=len(d),common_domain=counts(common),reported_sensitivity_domain=counts(reported),cohort_flow=flow,
        corrected_sscard_source_rows=len(sscard),source_all_participants_from_previously_loaded_release=True,
        old_fitted_M4_participant_intersection_not_reconstructed=True,prior_Cys_fitted=True,prior_B2M_loaded=True,
        common_weighted_coverage=float(np.average(common,weights=d.weight)),
        exclusion_reason_counts={k:int(v.sum()) for k,v in reasons.items()},exclusion_reason_patterns=dict(Counter(pattern)),
        missing_routine={k:int((~np.isfinite(d[k])).sum()) for k in CORE},
        minimum_survivor_followup_months=float(d.loc[d.dead.eq(0),'t_exam'].min()),
        zero_month_deaths=int((d.dead.eq(1)&d.t_exam.eq(0)).sum()),
        serum_calibration_applied_by_cycle={str(c['year']):'1.013*raw+0.147' if c['year']==1999 else 'no_additional_adjustment' for c in spec['cycles']},
        mortality_reader_audits=reader_audits,xpt_zero_audits=zero_audits,**design)
    return d,common,reported,audit,manifest
