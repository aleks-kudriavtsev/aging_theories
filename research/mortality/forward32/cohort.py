"""Strict temporal/domain/native-weight qualification; aggregate diagnostics."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ..b2m31 import cohort as prior
from ..b2m31.features import CORE
from ..xpt_checked import read_xpt_checked
from ..binary28.model import psu_partition
from .retrieve import ROOT,PLAN_COMMIT,check_study,sha

PROTEINS=('SSCYST','SSB2M')


def counts(series,codes):
    return {str(k):int(series.eq(k).sum()) for k in codes}|{'missing':int(series.isna().sum())}


def native_period_frame(frame):
    """Use native W4early/W2late; never infer weights from late outcomes."""
    if not frame.cycle.isin([1999,2001,2003]).all():raise ValueError('Unsupported temporal cycle')
    d=frame.copy();early=d.cycle.isin([1999,2001])
    d['partition']=np.where(early,'train','late')
    d['weight6y_original']=d.weight
    d['weight']=np.where(early,d.WTSSCB4Y,d.WTSSCB2Y)
    return d


def _inventory(frame):
    return {str(int(st)):sorted(map(int,g.SDMVPSU.unique()))
            for st,g in frame.groupby('SDMVSTRA',sort=True)}


def _moments(values,weights):
    x=np.asarray(values,float);w=np.asarray(weights,float);ok=np.isfinite(x)&np.isfinite(w)&(w>0)
    if not ok.any():return dict(n=0,weighted_mean=None,weighted_SD=None)
    mean=float(np.average(x[ok],weights=w[ok]))
    return dict(n=int(ok.sum()),weighted_mean=mean,
                weighted_SD=float(np.sqrt(np.average((x[ok]-mean)**2,weights=w[ok]))))


def _sum(series):
    value=float(series.sum())
    return value if np.isfinite(value) else None


def _domain_summary(frame):
    y=frame.dead.eq(1)&frame.t_exam.le(96);w=frame.weight.to_numpy(float)
    native_ok=np.isfinite(w)&(w>0)
    finite_weights=bool(np.isfinite(w).all())
    wsum=float(w.sum()) if finite_weights else None
    wsquared=float(w@w) if finite_weights else None
    selected=np.where(frame.partition.eq('train'),frame.WTSSCB4Y,frame.WTSSCB2Y)
    return dict(n=len(frame),events8=int(y.sum()),nonevents8=int((~y).sum()),
        positive_weight_events8=int((y.to_numpy()&native_ok).sum()),
        positive_weight_nonevents8=int(((~y).to_numpy()&native_ok).sum()),
        event_bearing_psu_n=int(frame.loc[y&native_ok,'cluster'].nunique()),
        nonevent_bearing_psu_n=int(frame.loc[~y&native_ok,'cluster'].nunique()),
        strata_n=int(frame.stratum.nunique()),psu_n=int(frame.cluster.nunique()),
        cycle_counts={str(int(k)):int(v) for k,v in frame.cycle.value_counts().sort_index().items()},
        native_weight_sum=wsum,native_weight_squared_sum=wsquared,
        effective_n=float(w.sum()**2/(w@w)) if len(w) and np.isfinite(w).all() and (w@w)>0 else None,
        source_W4_sum=_sum(frame.WTSSCB4Y),source_W4_missing_n=int(frame.WTSSCB4Y.isna().sum()),
        source_W2_sum=_sum(frame.WTSSCB2Y),original_W6_sum=_sum(frame.weight6y_original),
        selected_native_source_max_abs_difference=float(np.max(np.abs(w-selected))) if len(w) and finite_weights else None,
        minimum_survivor_followup_months=float(frame.loc[frame.dead.eq(0),'t_exam'].min()) if frame.dead.eq(0).any() else None,
        pristine_counts=counts(frame.SSPRIS,[0,1]),
        protein_flag_counts={k:counts(frame[k],[0,1,2]) for k in ('SSCYSTL','SSB2ML')},
        outcome_free_feature_moments={k:_moments(frame[k],w) for k in CORE},
        outcome_free_log2_protein_moments={k:_moments(np.log2(frame[k].where(frame[k]>0)),w) for k in PROTEINS},
        per_cycle_weights={str(int(c)):dict(n=len(g),source_W4_sum=_sum(g.WTSSCB4Y),
            source_W2_sum=_sum(g.WTSSCB2Y),original_W6_sum=_sum(g.weight6y_original),
            native_weight_sum=_sum(g.weight) if np.isfinite(g.weight).all() else None) for c,g in frame.groupby('cycle',sort=True)})


def temporal_qualification(frame,common,full_design,lock):
    """Return an explicit blocked diagnostic rather than silently changing target."""
    common=np.asarray(common)
    if common.shape!=(len(frame),) or common.dtype!=bool:raise ValueError('Primary common mask invalid')
    if not isinstance(full_design,dict) or set(full_design)!=set(map(str,range(1,44))):
        raise ValueError('Full period masked-design reference is incomplete')
    d=native_period_frame(frame);blockers=[];periods={}
    for part,low,high in [('train',1,28),('late',29,43)]:
        base=d.loc[d.partition.eq(part)];target=d.loc[d.partition.eq(part).to_numpy()&common]
        expected={str(st):full_design[str(st)] for st in range(low,high+1)}
        actual=_inventory(base)
        if actual!=expected:blockers.append(part+':full_masked_design_inventory_lost')
        if _inventory(target)!=actual:blockers.append(part+':common_domain_lost_eligible_PSU_fold_equivalence')
        if any(len(psus)<2 for psus in actual.values()):blockers.append(part+':singleton_stratum')
        w=base.weight.to_numpy(float)
        if not np.isfinite(w).all() or (w<=0).any():blockers.append(part+':native_weight_unusable')
        if not len(target):blockers.append(part+':empty_primary_domain')
        if (base.dead.eq(0)&base.t_exam.lt(96)).any():blockers.append(part+':early_censoring')
        summary=_domain_summary(target)
        if summary['positive_weight_events8']<lock['minimum_positive_weight_event_rows_per_part']:
            blockers.append(part+':event_support_below_predeclared_minimum')
        if summary['positive_weight_nonevents8']<lock['minimum_positive_weight_nonevent_rows_per_part']:
            blockers.append(part+':nonevent_support_below_predeclared_minimum')
        periods[part]=dict(eligible=_domain_summary(base),common=summary,
            full_masked_design=expected,eligible_masked_design=actual,
            common_masked_design=_inventory(target),
            common_native_weighted_coverage=float(target.weight.sum()/base.weight.sum()) if np.isfinite(base.weight).all() and base.weight.sum()>0 else None,
            weight_variable='WTSSCB4Y' if part=='train' else 'WTSSCB2Y',
            estimand='received_complete_case_measurable_domain_of_1999_2002' if part=='train' else 'received_complete_case_measurable_domain_of_2003_2004')
    train=d.partition.eq('train');late=d.partition.eq('late')
    overlap=len(set(d.loc[train,'SEQN'])&set(d.loc[late,'SEQN']))
    psu_overlap=len(set(d.loc[train,'cluster'])&set(d.loc[late,'cluster']))
    if overlap or psu_overlap:blockers.append('temporal_participant_or_masked_PSU_overlap')
    early_clusters=d.loc[train,'cluster'].to_numpy()
    fold_mapping={}
    if len(set(early_clusters))<3:
        blockers.append('train:insufficient_eligible_clusters_for_inner_split')
    else:
        eligible_folds=psu_partition(early_clusters,lock['inner_salt'],folds=3)
        fold_mapping={key:int(eligible_folds[np.flatnonzero(early_clusters==key)[0]]) for key in sorted(set(early_clusters))}
    return d,dict(gate=dict(temporal_evaluation_allowed=not blockers,blockers=blockers,
        support_threshold_is_not_power_or_calibration_certificate=True),periods=periods,
        participant_overlap_n=overlap,masked_PSU_overlap_n=psu_overlap,
        eligible_derived_early_inner_fold_mapping=fold_mapping,
        common_eligible_PSU_fold_equivalence_required=True,
        no_six_year_pooled_score=True,no_late_recalibration=True,
        assay_batch_individually_identified=False,storage_duration_individually_identified=False)


def load(source):
    lock=check_study();source=Path(source)
    wrapper=json.loads((source/'forward32_manifest.json').read_text());base=json.loads((source/'manifest.json').read_text())
    if (wrapper.get('schema')!='forward32-retrieval/1' or wrapper.get('plan_commit')!=PLAN_COMMIT
            or wrapper.get('lock_sha256')!=sha(ROOT/'LOCK.json')
            or wrapper.get('base_manifest_sha256')!=sha(source/'manifest.json')
            or wrapper.get('source_fingerprints')!={r['file']:r['sha256'] for r in base['files']}):
        raise ValueError('Public source intake not bound to forward protocol')
    d,common,reported,old_audit,actual_base=prior.load(source)
    if actual_base!=base:raise ValueError('Base source manifest mismatch')
    d,temporal=temporal_qualification(d,common,old_audit['full_DEMO_masked_design'],lock)
    sscard=read_xpt_checked(source/'1999_SSCARD_A.xpt')
    bound=wrapper['source_fingerprints']['1999_SSCARD_A.xpt']
    if sscard.attrs['xpt_zero_audit']['sha256']!=bound:raise ValueError('Assay bytes changed before observed marginal audit')
    observed=dict(rows=len(sscard),columns=list(sscard.columns),
        flags={k:counts(sscard[k],[0,1,2]) for k in ('SSCYSTL','SSB2ML')},pristine=counts(sscard.SSPRIS,[0,1]),
        weights={k:dict(positive=int((np.isfinite(sscard[k])&sscard[k].gt(0)).sum()),
                       missing=int(sscard[k].isna().sum()),zero=int(sscard[k].eq(0).sum())) for k in ('WTSSCB4Y','WTSSCB2Y')},
        source_sha256=bound,individual_batch_or_assay_date_fields_provided=False)
    audit=dict(**temporal,prior31_source_audit=old_audit,observed_full_SSCARD_marginals=observed,
        forward_source_wrapper_sha256=sha(source/'forward32_manifest.json'),
        original31_common_n=int(common.sum()),prior31_reported_sensitivity_not_reanalysed=True)
    manifest=dict(schema='forward32-source-manifest/1',intake=wrapper,source=base)
    return d,common,audit,manifest
