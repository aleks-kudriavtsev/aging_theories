"""New core definitions on the already-qualified public NHANES III sources."""
from __future__ import annotations
import json,re
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from ..cohort27 import importer as prior
from ..nfl22.analyze import sha
ROOT=Path(__file__).parent


def check_study():
    lock=json.loads((ROOT/'LOCK.json').read_text())
    for name,digest in lock['study_metadata_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Prespecified binary28 metadata changed')
    if sha(ROOT.parent/'cohort27/SOURCE_MANIFEST.json')!=lock['prior_source_manifest_sha256']:
        raise ValueError('Prior qualified source manifest changed')
    return lock


def read_core_checked(data_path,sas_path):
    fields=json.loads((ROOT/'CLINICAL_FIELDS.json').read_text())['positions_1_based_inclusive']
    sas=Path(sas_path).read_text(encoding='ascii')
    blocks=re.findall(r'(?ims)^[ \t]*INPUT\b(.*?);',sas)
    if len(blocks)!=1:raise ValueError('SAS INPUT unresolved')
    for name,(start,end) in fields.items():
        spans=re.findall(r'(?im)^[ \t]*'+name+r'[ \t]+(\d+)(?:[ \t]*-[ \t]*(\d+))?[ \t]*$',blocks[0])
        if [(int(a),int(b or a)) for a,b in spans]!=[(start,end)]:raise ValueError('Core INPUT definition mismatch: '+name)
    infile=re.findall(r'(?ims)^[ \t]*INFILE\b([^;]*);',sas)
    missover=any(re.search(r'\bMISSOVER\b',b,re.I) for b in infile)
    minimum=max(e for s,e in fields.values());lines=Path(data_path).read_text(encoding='ascii').splitlines();rows=[]
    for line in lines:
        if not ((minimum<=len(line)<=3346) if missover else len(line)==3346):raise ValueError('Core selected field truncation')
        row={}
        for name,(start,end) in fields.items():
            token=line[start-1:end].strip()
            if token in ('','.'):row[name]=np.nan
            elif token.isdigit():row[name]=float(token)
            else:raise ValueError('Non-numeric clinical source field')
        rows.append(row)
    manual=pd.DataFrame(rows)
    independent=pd.read_fwf(data_path,colspecs=[(s-1,e) for s,e in fields.values()],names=list(fields),header=None,na_values=['.'],dtype=float)
    pd.testing.assert_frame_equal(manual,independent,check_exact=True,check_dtype=False)
    if manual.SEQN.isna().any() or manual.SEQN.duplicated().any():raise ValueError('Invalid adult source IDs')
    return manual,dict(rows=len(manual),fields=len(fields),cell_comparisons=int(manual.size),exact_match=True,
        selected_fields_end=minimum,source_record_width_counts=dict(Counter(map(len,lines))),SAS_MISSOVER=missover)


def engineer(d):
    d=d.copy();never=d.HAR1.eq(2);current=d.HAR1.eq(1)&d.HAR3.eq(1);former=d.HAR1.eq(1)&d.HAR3.eq(2)
    known=never|current|former
    d['current']=np.where(known,current.astype(float),np.nan)
    d['former']=np.where(known,former.astype(float),np.nan)
    sbp=d.HAZMNK1R.where(np.isfinite(d.HAZMNK1R)&d.HAZMNK1R.gt(0)&d.HAZMNK1R.ne(888)&d.HAZNOK1R.isin(range(1,7)))
    age=(d.HSAGEIR-60)/10
    d['age_linear']=age;d['age_quadratic']=age**2;d['male']=d.HSSEX.map({1.:1.,2.:0.});d['sbp']=(sbp-130)/20
    return d

CORE=('age_linear','age_quadratic','male','current','former','sbp')
PROTEINS=('SSCYSTAT','SSB2M','SSBTP')


def load(source):
    check_study();source=Path(source)
    eligible,gate,audit,manifest=prior.load(source)
    prior_manifest=json.loads((ROOT.parent/'cohort27/SOURCE_MANIFEST.json').read_text())
    expected={r['file']:r['sha256'] for r in prior_manifest['files'] if r['file']!='index.htm'}
    actual={r['file']:r['sha256'] for r in manifest['files'] if r['file']!='index.htm'}
    if expected!=actual:raise ValueError('Qualified scientific source bytes changed before binary28')
    clinical,reader=read_core_checked(source/'adult.dat',source/'adult.sas')
    if len(clinical)!=20050:raise ValueError('New core reader changed official adult count')
    extra=[c for c in clinical if c not in eligible]+['SEQN']
    eligible=prior.join_one_to_one(eligible,clinical[extra]);eligible=engineer(eligible)
    eligible['weight']=eligible.WTCYPEX6
    eligible['stratum']='III:'+eligible.SDPSTRA6.astype(int).astype(str)
    eligible['SDMVPSU']=eligible.SDPPSU6
    eligible['cluster']=eligible.SDPSTRA6.astype(int).astype(str)+':'+eligible.SDPPSU6.astype(int).astype(str)
    protein_complete=np.isfinite(eligible[list(PROTEINS)]).all(axis=1).to_numpy()&eligible[list(PROTEINS)].gt(0).all(axis=1).to_numpy()
    if len(eligible)!=3899 or protein_complete.sum()!=3864:raise ValueError('Qualification27 domain mismatch')
    if (eligible.dead.eq(0)&eligible.t_exam.lt(96)).any():raise ValueError('Early censoring needs separately amended design')
    core_complete=np.isfinite(eligible[list(CORE)]).all(axis=1).to_numpy()
    common=protein_complete&core_complete
    audit.update(clinical_reader=reader,qualified_protein_n=int(protein_complete.sum()),core_complete_n=int(core_complete.sum()),
                 common_complete_n=int(common.sum()),core_missing_in_qualified_n=int((protein_complete&~core_complete).sum()),
                 missing_core={f:int((~np.isfinite(eligible[f])).sum()) for f in CORE},
                 model_origin='new_binary_logistic_no_frozen_model_coefficients')
    return eligible,common,audit,manifest
