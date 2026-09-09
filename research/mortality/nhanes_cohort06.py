"""Corrected-reader cohort, algebra copied from frozen nhanes_benchmark.make_cohort.
No legacy model/code is overwritten; only XPORT zero handling differs.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .nhanes_benchmark import read_mortality, yes_no
from .xpt_checked import read_xpt_checked as read_xpt

def make_cohort_checked(source: Path):
    frames=[]; flow=[]
    for yr,suff,bio in [(1999,'','LAB18'),(2001,'_B','L40_B'),(2003,'_C','L40_C')]:
        base=read_xpt(source/f'{yr}_DEMO{suff}.xpt')
        demo=base[['SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU']].copy()
        initial=len(demo); demo=demo[demo.RIDAGEYR.between(40,79)].copy()
        n_age=len(demo)
        labnames=([f'LAB{k}' for k in [10,11,13,16,25]] if yr==1999
                  else [f'L{k}{suff}' for k in [10,11,13,16,25]])
        for name in [f'{k}{suff}' for k in ['BPX','BMX','SMQ','MCQ','BPQ','DIQ']]+[bio]+labnames:
            d=read_xpt(source/f'{yr}_{name}.xpt')
            cols=['SEQN']+[c for c in d if c not in demo]
            demo=demo.merge(d[cols],on='SEQN',how='left',validate='1:1')
        mort=read_mortality(source/f'{yr}_mortality.dat')
        if set(base.SEQN)!=set(mort.SEQN): raise ValueError('Mortality/demographic key mismatch')
        demo=demo.merge(mort,on='SEQN',how='left',validate='1:1')
        demo=demo[demo.elig.eq(1)&demo.t_exam.notna()&demo.t_exam.ge(0)].copy()
        n_link=len(demo)
        # One specimen dataset spans all three survey cycles, despite suffix A.
        for name,cols in [('SSBNP_A',['WTSSCB2Y','WTSSCB4Y','SSBNP','SSBNPL','SSPRIS']),
                          ('SSTROP_A',['SSTNT','SSTNTLC','SSPRISTP']),
                          ('SSCARD_A',['SSCYST','SSCYSTL','SSB2M'])]:
            d=read_xpt(source/f'1999_{name}.xpt')
            demo=demo.merge(d[['SEQN']+cols],on='SEQN',how='left',validate='1:1')
        demo['weight']=demo['WTSSCB4Y' if yr<2003 else 'WTSSCB2Y']
        ok=demo.weight.gt(0)&demo[['SSBNP','SSTNT','SSCYST']].gt(0).all(axis=1)
        demo=demo[ok].copy(); demo['cycle']=yr
        flow.append({'cycle':yr,'survey_n':initial,'age40_79_n':n_age,
                     'eligible_with_exam_followup_n':n_link,'analysis_n':len(demo),
                     'excluded_cardio_assays_or_weight_n':n_link-len(demo)})
        scr=demo['LBDSCR' if yr==2001 else 'LBXSCR']
        if yr==1999: scr=1.013*scr+0.147  # CDC IDMS calibration, applied once
        demo['creatinine_standardized']=scr
        female=demo.RIAGENDR.eq(2)
        ratio=scr/np.where(female,.7,.9)
        egfr=142*np.minimum(ratio,1)**np.where(female,-.241,-.302)*np.maximum(ratio,1)**-1.2*.9938**demo.RIDAGEYR*np.where(female,1.012,1.)
        demo['egfr']=egfr
        demo['age']=(demo.RIDAGEYR-60)/10
        a=demo.RIDAGEYR
        # Restricted cubic spline, fixed knots (40,60,79); scaling only.
        demo['age_rcs']=(np.maximum(a-40,0)**3-(39/19)*np.maximum(a-60,0)**3+(20/19)*np.maximum(a-79,0)**3)/(39**2)
        demo['male']=demo.RIAGENDR.map({1.:1.,2.:0.})
        demo['age_male']=demo.age*demo.male
        known_never=demo.SMQ020.eq(2); known_smoker=demo.SMQ020.eq(1)&demo.SMQ040.isin([1,2,3])
        demo['smoker_current']=np.where(known_never,0,np.where(known_smoker,demo.SMQ040.isin([1,2]).astype(float),np.nan))
        demo['smoker_former']=np.where(known_never,0,np.where(known_smoker,demo.SMQ040.eq(3).astype(float),np.nan))
        demo['bmi']=(demo.BMXBMI-27)/5; demo['bmi_sq']=demo.bmi**2
        demo['sbp']=(demo[[f'BPXSY{i}' for i in [1,2,3]]].where(lambda x:x>0).mean(axis=1)-130)/20
        # Structural skips from never diagnosis/prescription are not unknown treatment.
        demo['bp_treatment']=yes_no(demo.BPQ050A)
        no_tx=demo.BPQ020.eq(2)|demo.BPQ040A.eq(2)
        demo.loc[no_tx,'bp_treatment']=0.
        demo['diabetes_history']=demo.DIQ010.map({1.:1.,2.:0.,3.:0.})
        cvd=demo[['MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']]
        demo['cvd_history']=np.where(cvd.eq(1).any(axis=1),1.,np.where(cvd.eq(2).all(axis=1),0.,np.nan))
        demo['cancer_history']=yes_no(demo.MCQ220)
        demo['total_cholesterol']=demo.LBXTC; demo['hdl']=demo['LBDHDL' if yr<2003 else 'LBXHDD']
        demo['hba1c']=demo.LBXGH; demo['albumin']=demo.LBXSAL*10 # g/L
        demo['egfr_low']=np.minimum(egfr,60)/15; demo['egfr_high']=np.maximum(egfr-60,0)/30
        # albumin ug/mL / creatinine mg/dL *100 = mg/g.
        demo['uacr']=100*demo.URXUMA/demo.URXUCR
        demo['log_uacr']=np.log2(demo.uacr.where(demo.uacr>0))
        demo['rdw']=demo.LBXRDW
        demo['log_wbc']=np.log2(demo.LBXWBCSI.where(demo.LBXWBCSI>0))
        demo['log_crp']=np.log2((demo.LBXCRP*10).where(demo.LBXCRP>0)) # mg/L
        demo['log_ntprobnp']=np.log2(demo.SSBNP)
        demo['log_troponin']=np.log2(demo.SSTNT)
        demo['log_cystatin']=np.log2(demo.SSCYST)
        demo['time_years']=demo.t_exam.clip(lower=.5)/12 # rounded zero-month deaths midpoint convention
        frames.append(demo)
    data=pd.concat(frames,ignore_index=True)
    if data.SEQN.duplicated().any(): raise ValueError('Duplicate persons across survey cycles')
    if ((data.dead==0)&(data.time_years<10)).any():
        raise ValueError('Incomplete 10-year ascertainment: use censoring-adjusted validation, not binary metrics')
    return data,flow
