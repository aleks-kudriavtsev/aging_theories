"""Retrospective endpoint feasibility across six previously studied NHANES cycles.

General examination sample, not selected for archived cardiac proteins. Export
aggregates only. No training, new external-validation claim, or restricted data.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from .xpt_checked import read_xpt_checked
from .validate_temporal08 import engineer_new,independent_mortality_reader,audit_mortality_reader
from .endpoint_contract11 import PUBLIC_CAUSES,TARGETS,classify_public,observability,fixed_horizon_states,precision_planning,mec_weight_1999_2010

CYCLES=(1999,2001,2003,2005,2007,2009)
FEATURES=('age','age_rcs','male','age_male','smoker_current','smoker_former','bmi','bmi_sq','sbp',
          'bp_treatment','diabetes_history','cvd_history','cancer_history','total_cholesterol','hdl',
          'hba1c','egfr_low','egfr_high','log_uacr','albumin','rdw','log_wbc','log_crp')


def fingerprint(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def read_verified(path: Path, manifests:dict,used:dict):
    if path.name not in manifests:raise ValueError('Scientific input missing from retrieved manifest: '+path.name)
    expected=manifests[path.name]
    if fingerprint(path)!=expected['sha256']:raise ValueError('Scientific input changed: '+path.name)
    used[path.name]={k:expected.get(k) for k in ('file','url','sha256','bytes')}
    return read_xpt_checked(path) if path.suffix=='.xpt' else independent_mortality_reader(path)


def cohort(sources: dict[int,Path]):
    manifest={};used={};flow=[];frames=[];reader_checks=[]
    for path in set(sources.values()):
        doc=json.loads((path/'manifest.json').read_text());rows=doc['files']
        for row in rows:
            name=row.get('file',row.get('name'))
            if row.get('sha256'):
                if name in manifest and row['sha256']!=manifest[name]['sha256']:raise ValueError('Conflicting manifest entries')
                manifest[name]={**row,'file':name}
    for year in CYCLES:
        root=sources[year];suff='' if year==1999 else '_'+chr(ord('B')+(year-2001)//2)
        demo=read_verified(root/f'{year}_DEMO{suff}.xpt',manifest,used)
        if demo.SEQN.duplicated().any():raise ValueError('Duplicate demographic IDs')
        ids=set(demo.SEQN)
        demo_cols=['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR']+(['WTMEC4YR'] if year<2003 else [])
        d=demo[demo.RIDAGEYR.between(40,79)][demo_cols].copy()
        age_n=len(d)
        names=[c+suff for c in ('BPX','BMX','SMQ','MCQ','BPQ','DIQ')]
        if year<=2003:
            names+=['LAB18' if year==1999 else 'L40'+suff]
            names+=[('LAB' if year==1999 else 'L')+str(k)+suff for k in (10,11,13,16,25)]
        else:names+=[c+suff for c in ('BIOPRO','TCHOL','HDL','GHB','ALB_CR','CBC','CRP')]
        for name in names:
            raw=read_verified(root/f'{year}_{name}.xpt',manifest,used)
            if not set(raw.SEQN)<=ids:raise ValueError('Unexpected IDs outside source cycle')
            d=d.merge(raw[['SEQN']+[c for c in raw if c not in d]],on='SEQN',how='left',validate='1:1')
        mortfile=root/f'{year}_mortality.dat';mort=read_verified(mortfile,manifest,used)
        reader_checks.append(audit_mortality_reader(mortfile))
        if set(mort.SEQN)!=ids:raise ValueError('Mortality and demographic identifiers differ')
        d=d.merge(mort,on='SEQN',how='left',validate='1:1')
        linked=d.elig.eq(1)&d.dead.isin([0,1])&d.t_exam.notna()&d.t_exam.ge(0)
        linked_n=int(linked.sum());d=d[linked&d.WTMEC2YR.gt(0)&np.isfinite(d.WTMEC2YR)&d.RIAGENDR.isin([1,2])].copy()
        # Harmonize only to the already reviewed common feature algebra.
        if year==1999:d['LBXSCR']=1.013*d.LBXSCR.where(d.LBXSCR>0)+.147
        if year==2001:d['LBXSCR']=d.LBDSCR
        if year<=2003:d['LBDHDD']=d['LBDHDL' if year<2003 else 'LBXHDD']
        d=engineer_new(d,year)
        d['complete_profile']=np.isfinite(d[list(FEATURES)]).all(axis=1)
        d['cycle']=year
        d['pooled_weight']=[mec_weight_1999_2010(year,float(w),float(v)) for w,v in
            zip(d.WTMEC2YR,d.WTMEC4YR if year<2003 else d.WTMEC2YR)]
        d['cause_state']=[classify_public(float(a),float(b),survey_start=year) for a,b in zip(d.dead,d.cause)]
        if d[['SDMVSTRA','SDMVPSU']].isna().any().any():raise ValueError('Missing survey design')
        flow.append({'cycle':year,'survey_n':len(demo),'age40_79_n':age_n,'linked_n':linked_n,
                     'eligible_examination_n':len(d),'complete_profile_n':int(d.complete_profile.sum()),
                     'minimum_alive_followup_years':float(d.loc[d.dead==0,'t_exam'].min()/12),
                     'zero_month_deaths':int(((d.dead==1)&(d.t_exam==0)).sum())})
        frames.append(d)
    all_data=pd.concat(frames,ignore_index=True)
    if all_data.SEQN.duplicated().any():raise ValueError('Participants duplicated across cycles')
    return all_data,flow,used,reader_checks


def table_audit(d:pd.DataFrame):
    events=[];intervals=[];completeness=[];checks=0
    for period,mask in [('1999_2010',np.ones(len(d),bool))]+[(str(y),d.cycle.eq(y).to_numpy()) for y in CYCLES]:
        sub=d.loc[mask].copy()
        for domain,base in [('eligible',sub),('complete_profile',sub[sub.complete_profile])]:
            for horizon in (1,5,8):
                if ((base.dead==0)&(base.t_exam<horizon*12)).any():raise ValueError('Unsupported binary horizon')
                dies=base.dead.eq(1)&base.t_exam.le(horizon*12)
                # Domain stratum summaries are transparent raw counts, not country probabilities.
                strata={'all':np.ones(len(base),bool),'female':base.RIAGENDR.eq(2).to_numpy(),
                  'male':base.RIAGENDR.eq(1).to_numpy(),'age40_59':base.RIDAGEYR.lt(60).to_numpy(),
                  'age60_79':base.RIDAGEYR.ge(60).to_numpy()}
                for label,stratum in strata.items():
                    b=base.loc[stratum];dead=dies.loc[b.index]
                    n=len(b);tot=int(dead.sum());w=b['pooled_weight' if period=='1999_2010' else 'WTMEC2YR'].to_numpy(float)
                    counted=0;weighted=0.
                    for code,info in PUBLIC_CAUSES.items():
                        is_event=dead&b.cause_state.eq(code);count=int(is_event.sum());counted+=count
                        fraction=float(np.sum(w*is_event.to_numpy())/np.sum(w)) if n else None
                        if fraction is not None:weighted+=fraction
                        events.append({'period':period,'domain':domain,'stratum':label,'horizon_years':horizon,
                          'cause_code':code,'cause_id':info['id'],'cause_label_ru':info['label_ru'],
                          'n':n,'all_deaths':tot,'cause_deaths':count,
                          'weighted_cause_fraction_in_analytic_sample':fraction,
                          'endpoint_type':'underlying_death','individual_risk':False})
                    if counted!=tot:raise ValueError('Deaths lost in endpoint partition')
                    if n and abs(weighted-float(np.average(dead,weights=w)))>1e-12:raise ValueError('Weighted partition mismatch')
                    checks+=2
            for lo,hi in [(0,1),(1,5),(5,8)]:
                for code in PUBLIC_CAUSES:
                    count=int((base.dead.eq(1)&base.cause_state.eq(code)&base.t_exam.le(hi*12)
                               &(base.t_exam.ge(0) if lo==0 else base.t_exam.gt(lo*12))).sum())
                    intervals.append({'period':period,'domain':domain,'interval_start':lo,'interval_end':hi,
                      'cause_code':code,'events':count,'zero_event_interval':count==0,
                      **{k:v for k,v in precision_planning(count).items() if k!='events'}})
        for feature in FEATURES:
            x=sub[feature].to_numpy(float)
            completeness.append({'period':period,'feature':feature,'eligible_n':len(sub),
              'finite_n':int(np.isfinite(x).sum()),'missing_n':int((~np.isfinite(x)).sum())})
    return events,intervals,completeness,checks


def run(old_source:Path,source08:Path,source09:Path,out:Path):
    if out.exists():raise ValueError('Use a new output directory; historical audit cannot be overwritten')
    roots={y:old_source if y<2005 else source08 if y<2009 else source09 for y in CYCLES}
    d,flow,used,reader=cohort(roots)
    ev,intervals,complete,checks=table_audit(d)
    # Previously frozen complete-case cohorts are an independent implementation check.
    if int(d.loc[d.cycle.isin([2005,2007]),'complete_profile'].sum())!=5287:
        raise ValueError('Previously reported 2005-2008 complete cohort not reproduced')
    if int(d.loc[d.cycle==2009,'complete_profile'].sum())!=3182:
        raise ValueError('Previously reported 2009-2010 complete cohort not reproduced')
    sub=d[d.complete_profile]
    result={'date':'2026-09-10','status':'retrospective_endpoint_feasibility_not_model_validation',
      'cycles':list(CYCLES),'age_range':[40,79],'followup_end':'2019-12-31',
      'new_independent_validation':False,'new_model_fit':False,'clinical_use_ready':False,
      'raw_participant_records_exported':False,'n_eligible':len(d),'n_complete_profile':len(sub),
      'deaths_5y_eligible':int((d.dead.eq(1)&d.t_exam.le(60)).sum()),
      'deaths_8y_eligible':int((d.dead.eq(1)&d.t_exam.le(96)).sum()),
      'deaths_5y_complete':int((sub.dead.eq(1)&sub.t_exam.le(60)).sum()),
      'deaths_8y_complete':int((sub.dead.eq(1)&sub.t_exam.le(96)).sum()),
      'data_files_verified':len(used),'flow':flow,'arithmetic_checks':checks,
      'independent_reader_cells_compared':sum(x['cell_comparisons'] for x in reader),
      'public_cause_groups':10,'unknown_cause_preserved':True,
      'general_examination_sample_not_protein_subsample':True,
      'weights_note':'Pooled: WTMEC4YR/3 for1999-2002, WTMEC2YR/6 for2003-2010. Per-cycle: WTMEC2YR. Descriptive analytic-sample fractions, not linkage-adjusted national risk',
      'future_unseen_cycles_not_retrieved':[2011,2013,2015,2017]}
    out.mkdir(parents=True)
    for name,rows in [('event_support.csv',ev),('interval_support.csv',intervals),('feature_completeness.csv',complete),
      ('endpoint_observability.csv',[{'survey_start':y,**observability(t,survey_start=y)} for y in (1999,2011,2015,2017) for t in TARGETS])]:
        pd.DataFrame(rows).to_csv(out/name,index=False)
    for name,value in [('results.json',result),('source_manifest.json',list(used.values())),('reader_audit.json',reader)]:
        (out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('old-source','source08','source09','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(a.old_source,a.source08,a.source09,a.out)
