"""Gate and verify every source before parsing public first-exam participant files."""
from __future__ import annotations
import hashlib,json,re
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from ..hrs24.access import strict_json
from ..nfl22 import analyze as base
from .gate import require_qualification
from .retrieve import STATIC,DYNAMIC,INDEX,resolve_links,allowed_url
ROOT=Path(__file__).parent
FIELDS={'SEQN':(1,5),'HSSEX':(15,15),'HSAGEIR':(18,19),'HSAGEU':(20,20),'SDPPSU6':(43,43),'SDPSTRA6':(44,45)}


def preflight(source):
    # These checks precede any participant-file open, not merely model fitting.
    contract=strict_json((ROOT/'CONTRACT.json').read_text());gate=require_qualification(contract)
    lock=strict_json((ROOT/'LOCK.json').read_text())
    for key,name in [('plan_sha256','ANALYSIS_PLAN.md'),('contract_sha256','CONTRACT.json'),('source_map_sha256','SOURCE_MAP.json')]:
        if base.sha(ROOT/name)!=lock[key]:raise ValueError('Locked study metadata changed')
    source=Path(source);manifest=strict_json((source/'manifest.json').read_text())
    if manifest.get('plan_commit')!='a15532f' or manifest.get('lock_sha256')!=base.sha(ROOT/'LOCK.json'):raise ValueError('Source/protocol binding mismatch')
    expected={'index.htm',*STATIC,*DYNAMIC};rows={}
    for row in manifest['files']:
        name=row['file']
        if name in rows or name not in expected or Path(name).name!=name:raise ValueError('Unexpected/duplicate source')
        if row.get('status')!='retrieved' or not allowed_url(row.get('url','')) or not allowed_url(row.get('resolved_url','')):raise ValueError('Unverified official public source')
        rows[name]=row
    if set(rows)!=expected:raise ValueError('Missing required source')
    if rows['index.htm']['url']!=INDEX or any(rows[n]['url']!=u for n,u in STATIC.items()):raise ValueError('Wrong official release')
    # Verify all bytes first. No table parser or SAS interpreter runs here.
    for name,row in rows.items():
        path=source/name
        if not path.is_file() or path.stat().st_size!=row['bytes'] or base.sha(path)!=row['sha256']:raise ValueError('Source fingerprint mismatch: '+name)
    dynamic=resolve_links((source/'index.htm').read_text(encoding='utf-8-sig'))
    if any(rows[n]['url']!=u for n,u in dynamic.items()):raise ValueError('Data URL not bound to release index')
    fingerprints={n:r['sha256'] for n,r in rows.items()}
    stable=hashlib.sha256(json.dumps(fingerprints,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return gate,manifest,stable


def read_adult_checked(data_path,sas_path):
    sas=Path(sas_path).read_text(encoding='ascii')
    blocks=re.findall(r'(?ims)^[ \t]*INPUT\b(.*?);',sas)
    if len(blocks)!=1:raise ValueError('SAS INPUT dictionary block unresolved')
    for name,(start,end) in FIELDS.items():
        spans=re.findall(r'(?im)^[ \t]*'+name+r'[ \t]+(\d+)(?:[ \t]*-[ \t]*(\d+))?[ \t]*$',blocks[0])
        normalized=[(int(a),int(b or a)) for a,b in spans]
        if normalized!=[(start,end)]:raise ValueError('SAS dictionary position unresolved: '+name)
    infile_blocks=re.findall(r'(?ims)^[ \t]*INFILE\b([^;]*);',sas)
    missover=any(re.search(r'\bMISSOVER\b',block,re.I) for block in infile_blocks)
    selected_width=max(e for _,e in FIELDS.values())
    records=[]
    lines=Path(data_path).read_text(encoding='ascii').splitlines()
    widths=dict(Counter(map(len,lines)))
    for line in lines:
        valid_width=(selected_width<=len(line)<=3346) if missover else len(line)==3346
        if not valid_width:raise ValueError('Unexpected adult record width; aggregate_width_counts='+json.dumps(widths,sort_keys=True))
        row={}
        for name,(start,end) in FIELDS.items():
            token=line[start-1:end].strip()
            if token in ('','.'):row[name]=np.nan
            elif token.isdigit():row[name]=float(token)
            else:raise ValueError('Non-numeric adult field')
        records.append(row)
    manual=pd.DataFrame(records)
    independent=pd.read_fwf(data_path,colspecs=[(s-1,e) for s,e in FIELDS.values()],names=list(FIELDS),header=None,na_values=['.'],dtype=float)
    pd.testing.assert_frame_equal(manual,independent,check_exact=True,check_dtype=False)
    if manual.SEQN.isna().any() or manual.SEQN.duplicated().any():raise ValueError('Invalid adult identifiers')
    return manual,dict(rows=len(manual),fields=len(FIELDS),cell_comparisons=int(manual.size),exact_match=True,
                       source_record_width_counts=widths,SAS_MISSOVER=missover,selected_fields_end=selected_width,
                       rows_with_absent_trailing_unselected_characters=sum(n for width,n in widths.items() if width<3346))


def join_one_to_one(left,right):
    for d in (left,right):
        if d.SEQN.isna().any() or d.SEQN.duplicated().any():raise ValueError('Invalid/duplicate linkage identifiers')
    return left.merge(right,on='SEQN',how='left',validate='1:1')


def load(source):
    source=Path(source);gate,manifest,stable=preflight(source)
    adult,adult_audit=read_adult_checked(source/'adult.dat',source/'adult.sas')
    if len(adult)!=20050:raise ValueError('Adult release row count differs from official20050')
    cys=base.read_xpt_checked(source/'cystatin.xpt');btp=base.read_xpt_checked(source/'btp.xpt')
    mort=base.independent_mortality_reader(source/'mortality.dat')
    if not set(adult.SEQN)<=set(mort.SEQN) or not set(cys.SEQN)<=set(mort.SEQN) or not set(btp.SEQN)<=set(mort.SEQN):raise ValueError('Survey/mortality identifier universe mismatch')
    shared=cys[['SEQN','WTCYPEX6']].merge(btp[['SEQN','WTCYPEX6']],on='SEQN',validate='1:1',suffixes=('_cys','_btp'))
    mask=shared[['WTCYPEX6_cys','WTCYPEX6_btp']].notna().all(axis=1)
    differences=np.abs(shared.loc[mask,'WTCYPEX6_cys']-shared.loc[mask,'WTCYPEX6_btp'])
    if not np.allclose(shared.loc[mask,'WTCYPEX6_cys'],shared.loc[mask,'WTCYPEX6_btp'],rtol=1e-12,atol=1e-8):raise ValueError('First-exam protein weights disagree')
    d=adult.loc[adult.HSAGEIR.between(40,75)&adult.HSAGEU.eq(2)].copy()
    age_n=len(d)
    d=join_one_to_one(d,cys[['SEQN','SSCYSTAT']])
    d=join_one_to_one(d,btp[['SEQN','SSB2M','SSBTP','WTCYPEX6']])
    d=join_one_to_one(d,mort[['SEQN','elig','dead','cause','t_exam']])
    eligible=(d.elig.eq(1)&d.dead.isin([0,1])&d.t_exam.ge(0)&np.isfinite(d.t_exam)&d.HSSEX.isin([1,2])
              &d.WTCYPEX6.gt(0)&np.isfinite(d.WTCYPEX6))
    d=d.loc[eligible].copy().reset_index(drop=True)
    if d.empty or d[['SDPSTRA6','SDPPSU6']].isna().any().any() or not d.SDPPSU6.isin([1,2]).all():raise ValueError('Empty or invalid survey design')
    audit=dict(adult_source_n=len(adult),age40_75_n=age_n,eligible_protein_weight_n=len(d),adult_reader=adult_audit,
        mortality_reader=base.audit_mortality_reader(source/'mortality.dat'),
        cys_source_n=len(cys),btp_source_n=len(btp),cys_outside_adult_n=len(set(cys.SEQN)-set(adult.SEQN)),
        btp_outside_adult_n=len(set(btp.SEQN)-set(adult.SEQN)),weight_pairs_checked=int(mask.sum()),
        maximum_paired_weight_difference=float(differences.max()) if len(differences) else None,
        source_content_sha256=stable,source_manifest_sha256=base.sha(source/'manifest.json'),
        source_fingerprints={r['file']:r['sha256'] for r in manifest['files']})
    return d,gate,audit,manifest
