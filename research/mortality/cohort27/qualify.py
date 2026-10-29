"""Real first-exam protein/long-followup qualification, no prediction metrics."""
from __future__ import annotations
import argparse,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
from ..nfl22.analyze import sha,dump
from .importer import load,ROOT
MARKERS=('SSCYSTAT','SSB2M','SSBTP')


def summarize_domain(d):
    w=d.WTCYPEX6.to_numpy(float)
    event=d.dead.eq(1)&d.t_exam.le(96)
    early=d.dead.eq(0)&d.t_exam.lt(96)
    return dict(n=len(d),observed_deaths8=int(event.sum()),early_censored8=int(early.sum()),
        horizon_years=8,known_status8=int((~early).sum()),
        minimum_survivor_followup_years=float(d.loc[d.dead.eq(0),'t_exam'].min()/12) if d.dead.eq(0).any() else None,
        effective_weighted_n=float(w.sum()**2/(w@w)) if len(w) else None,
        strata_n=int(d.SDPSTRA6.nunique()),psu_n=int(d[['SDPSTRA6','SDPPSU6']].drop_duplicates().shape[0]))


def run(source,out):
    out=Path(out)
    if out.exists():raise ValueError('New output directory required')
    d,gate,audit,manifest=load(source)
    available={f:np.isfinite(d[f]).to_numpy()&d[f].gt(0).to_numpy() for f in MARKERS}
    common=np.logical_and.reduce(list(available.values()))
    counts={f:{'positive_finite_n':int(available[f].sum()),'missing_n':int(d[f].isna().sum()),
                'nonpositive_finite_n':int((np.isfinite(d[f])&d[f].le(0)).sum()),
                'analytically_usable_above_LOQ_n':None} for f in MARKERS}
    result=dict(analysis='real_public_first_exam_protein_endpoint_qualification_not_model_validation',
        cohort='NHANESIII1988_1994',date='2026-10-29',plan_commit='a15532f',gate=gate,source_audit=audit,
        eligible_domain=summarize_domain(d),markers=counts,joint_three_reported_numeric_domain=summarize_domain(d.loc[common]),
        weighted_joint_numeric_coverage=float(np.average(common,weights=d.WTCYPEX6)),
        eight_year_AUC=None,eight_year_Brier=None,
        clinical_use_ready=False,participant_records_exported=False,
        model_evaluation_status='blocked_before_predictions_clinical_semantics_and_assay_transport_unresolved',
        plan_sha256=sha(ROOT/'ANALYSIS_PLAN.md'),lock_sha256=sha(ROOT/'LOCK.json'),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)):sha(p) for p in sorted(ROOT.parent.rglob('*.py'))},
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__))
    out.mkdir(parents=True);dump(out/'RESULTS.json',result);dump(out/'SOURCE_MANIFEST.json',manifest)
    print('COHORT27_RESULTS_BEGIN');print(json.dumps(result,indent=2,allow_nan=False));print('COHORT27_RESULTS_END')
    print('COHORT27_MANIFEST_BEGIN');print(json.dumps(manifest,indent=2,allow_nan=False));print('COHORT27_MANIFEST_END')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
