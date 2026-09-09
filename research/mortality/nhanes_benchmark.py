"""Exploratory survey-weighted mortality benchmark; NOT a clinical calculator.

Train: NHANES 1999-2002; temporal holdout: 2003-2004; baseline ages 40-79.
Public-use 2019 LMF, exam origin. All-cause death (including infections/external
causes), not the strict noninfectious endpoint. Fixed 10-year administrative cap.
Input files are read locally; this module does not download or publish records.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize, brentq
from scipy.special import expit, logit
from sklearn.metrics import roc_auc_score

INTERVAL_ENDS = np.array([1., 5., 10.])
PENALTY = 1.0  # L2, sum weighted likelihood; subject weights normalized to mean 1
BASE = ['age','age_rcs','male','age_male']
CLINICAL = ['smoker_current','smoker_former','bmi','bmi_sq','sbp','bp_treatment',
            'diabetes_history','cvd_history','cancer_history']
ROUTINE = ['total_cholesterol','hdl','hba1c','egfr_low','egfr_high','log_uacr',
           'albumin','rdw','log_wbc','log_crp']
PANELS = {'M0_age_sex':BASE,
          'M1_routine':BASE+CLINICAL+ROUTINE,
          'M2_NTproBNP':BASE+CLINICAL+ROUTINE+['log_ntprobnp'],
          'M3_troponin':BASE+CLINICAL+ROUTINE+['log_ntprobnp','log_troponin'],
          'M4_cystatinC':BASE+CLINICAL+ROUTINE+['log_ntprobnp','log_troponin','log_cystatin']}
CONTINUOUS = set(BASE + ['bmi','bmi_sq','sbp'] + ROUTINE +
                 ['log_ntprobnp','log_troponin','log_cystatin']) - {'male'}
# Boolean/binary fields are not winsorized. No clinical decision cutoffs are used.


def read_xpt(path: Path) -> pd.DataFrame:
    if not path.read_bytes().startswith(b'HEADER RECORD'):
        raise ValueError(f'Not a SAS XPORT file: {path.name}')
    d=pd.read_sas(path,format='xport')
    if d.SEQN.isna().any() or d.SEQN.duplicated().any():
        raise ValueError(f'Nonunique/missing participant keys: {path.name}')
    return d


def read_mortality(path: Path) -> pd.DataFrame:
    """2019 public-use NHANES fixed-width layout; positions are 1-based in CDC.

    SEQN 1-6; ELIGSTAT 15; MORTSTAT 16; UCOD_LEADING 17-19;
    DIABETES 20; HYPERTEN 21; PERMTH_INT 43-45; PERMTH_EXM 46-48.
    Underlying causes and mortality conditions are NEVER prediction covariates.
    """
    lines=path.read_text().splitlines()
    if not lines or any(len(s) not in (46,47,48) for s in lines):
        raise ValueError('Unexpected NHANES public mortality record length')
    d=pd.read_fwf(path,colspecs=[(0,6),(14,15),(15,16),(16,19),(19,20),
                                (20,21),(42,45),(45,48)],
                  names=['SEQN','elig','dead','cause','death_diabetes',
                         'death_hypertension','t_int','t_exam'],na_values=['.'])
    if d.SEQN.duplicated().any() or not set(d.elig.unique()) <= {1,2,3}:
        raise ValueError('Invalid key or eligibility code')
    eligible=d.elig.eq(1)
    if not set(d.loc[eligible,'dead'].unique()) <= {0,1}:
        raise ValueError('Missing/invalid vital status in eligible records')
    if d.loc[~eligible,'dead'].notna().any():
        raise ValueError('Ineligible vital status must remain missing')
    return d


def yes_no(series):
    return series.map({1.:1.,2.:0.})


def weighted_quantile(x,w,q):
    x,w=np.asarray(x,float),np.asarray(w,float)
    keep=np.isfinite(x)&np.isfinite(w)&(w>0)
    if not keep.any(): raise ValueError('No usable training observations')
    order=np.argsort(x[keep],kind='stable'); xx=x[keep][order]; ww=w[keep][order]
    return float(np.interp(q,np.cumsum(ww)/ww.sum(),xx))


def make_cohort(source: Path):
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


class Preprocessor:
    def fit(self,data,features,weights):
        self.features=list(features); self.stats={}; self.columns=[]
        for f in features:
            x=data[f].to_numpy(float,copy=True); x[~np.isfinite(x)]=np.nan
            med=weighted_quantile(x,weights,.5)
            low,high=(weighted_quantile(x,weights,.005),weighted_quantile(x,weights,.995)) if f in CONTINUOUS else (-1e100,1e100)
            imputed=np.clip(np.where(np.isnan(x),med,x),low,high)
            mean=float(np.average(imputed,weights=weights)); sd=float(np.sqrt(np.average((imputed-mean)**2,weights=weights)))
            if sd<1e-10: sd=1.
            self.stats[f]={'median':med,'lower':low,'upper':high,'mean':mean,'sd':sd,'missing_indicator':bool(np.isnan(x).any())}
            self.columns.append(f)
            if self.stats[f]['missing_indicator']: self.columns.append(f+'__missing')
        return self

    def transform(self,data):
        columns=[]
        for f in self.features:
            st=self.stats[f]; x=data[f].to_numpy(float); missing=~np.isfinite(x)
            x=np.clip(np.where(missing,st['median'],x),st['lower'],st['upper'])
            columns.append((x-st['mean'])/st['sd'])
            if st['missing_indicator']: columns.append(missing.astype(float))
        z=np.column_stack(columns)
        if not np.isfinite(z).all(): raise ValueError('Nonfinite transformed predictors')
        return z


def interval_design(z,t,event):
    z=np.asarray(z,float); t=np.asarray(t,float); event=np.asarray(event,float)
    if (z.ndim!=2 or t.ndim!=1 or event.ndim!=1 or len(t)!=len(z) or len(event)!=len(t)
            or len(t)==0 or np.any(t<=0) or not np.isfinite(t).all()
            or not np.isfinite(z).all() or not np.isin(event,[0,1]).all()):
        raise ValueError('Invalid predictor matrix, follow-up or event coding')
    matrices=[]; exposure=[]; deaths=[]; ids=[]
    lower=0.
    for j,upper in enumerate(INTERVAL_ENDS):
        keep=t>lower; idx=np.flatnonzero(keep)
        dt=np.minimum(t[keep],upper)-lower
        a=np.zeros((len(idx),len(INTERVAL_ENDS))); a[:,j]=1
        matrices.append(np.column_stack([a,z[keep]])); exposure.append(dt)
        deaths.append(((event[keep]==1)&(t[keep]<=upper)).astype(float)); ids.append(idx)
        lower=upper
    return np.vstack(matrices),np.concatenate(exposure),np.concatenate(deaths),np.concatenate(ids)


def fit_piecewise(z,time,event,weight,penalty=PENALTY):
    x,exposure,deaths,ids=interval_design(z,time,event)
    w=np.asarray(weight,float)
    if w.ndim!=1 or len(w)!=len(z) or not np.isfinite(w).all() or np.any(w<=0):
        raise ValueError('Training weights must be finite and positive')
    if isinstance(penalty,bool) or not np.isfinite(penalty) or penalty<0:
        raise ValueError('Penalty must be finite and nonnegative')
    w=w/w.mean(); ew=w[ids]
    penalty_mask=np.r_[np.zeros(3),np.ones(z.shape[1])]
    initial=np.zeros(x.shape[1])
    for j in range(3):
        keep=x[:,j]==1
        initial[j]=np.log(max(float(np.sum(ew[keep]*deaths[keep])),.01)/np.sum(ew[keep]*exposure[keep]))
    def objective(beta):
        lp=x@beta
        if np.max(lp)>100: return 1e80,np.ones_like(beta)*1e40
        mu=exposure*np.exp(lp)
        loss=np.sum(ew*(mu-deaths*lp))+.5*penalty*np.sum(penalty_mask*beta**2)
        grad=x.T@(ew*(mu-deaths))+penalty*penalty_mask*beta
        return loss,grad
    result=minimize(objective,initial,jac=True,method='L-BFGS-B',options={'maxiter':2500,'ftol':1e-12,'gtol':1e-6,'maxls':40})
    if not result.success: raise RuntimeError(f'Optimizer failed: {result.message}')
    gradient=float(np.max(np.abs(objective(result.x)[1])))
    if gradient>.05: raise RuntimeError(f'Poor convergence: gradient {gradient}')
    return {'baseline_log_hazards':result.x[:3].tolist(),'coefficients':result.x[3:].tolist(),
            'penalty':penalty,'converged':True,'gradient_max_abs':gradient,'iterations':int(result.nit)}


def predict(z,model,horizon):
    if isinstance(horizon,bool) or not np.isfinite(horizon) or not 0<=horizon<=10:
        raise ValueError('Horizon must be in [0,10]; no 20-year extrapolation')
    z=np.asarray(z,float); coef=np.asarray(model['coefficients'],float)
    baseline=np.asarray(model['baseline_log_hazards'],float)
    if (z.ndim!=2 or coef.ndim!=1 or z.shape[1]!=len(coef) or baseline.shape!=(3,)
            or not np.isfinite(z).all() or not np.isfinite(coef).all()
            or not np.isfinite(baseline).all() or np.max(baseline)>700):
        raise ValueError('Invalid predictors or model coefficients')
    durations=np.maximum(0,np.minimum(INTERVAL_ENDS,horizon)-np.r_[0,INTERVAL_ENDS[:-1]])
    h0=float(durations@np.exp(model['baseline_log_hazards']))
    lp=np.asarray(z)@np.asarray(model['coefficients'])
    cumulative=h0*np.exp(np.minimum(lp,700))
    return -np.expm1(-cumulative)


def evaluate(y,p,w):
    y,p,w=np.asarray(y,float),np.asarray(p,float),np.asarray(w,float)
    if np.sum(w)<=0: raise ValueError('Empty weighted validation set')
    prev=float(np.average(y,weights=w)); pred=float(np.average(p,weights=w))
    auc=float(roc_auc_score(y,p,sample_weight=w)) if 0<prev<1 else None
    return {'n':len(y),'events':int(y.sum()),'weighted_observed':prev,
            'weighted_predicted':pred,'observed_expected':prev/pred if pred>0 else None,
            'auc':auc,'brier':float(np.average((y-p)**2,weights=w))}


def calibration(y,p,w):
    """Validation diagnostics only: fitted slope/intercept never alter predictions."""
    y=np.asarray(y,float); w=np.asarray(w,float); w=w/w.mean(); l=logit(np.clip(p,1e-9,1-1e-9))
    def fn(b):
        lp=b[0]+b[1]*l
        loss=np.sum(w*(np.logaddexp(0,lp)-y*lp))
        err=w*(expit(lp)-y)
        return loss,np.array([err.sum(),err@l])
    fit=minimize(fn,[0,1],jac=True,method='BFGS',options={'gtol':1e-5})
    cil=brentq(lambda a:float(np.sum(w*(expit(l+a)-y))),-40,40)
    return {'calibration_intercept':float(fit.x[0]),'calibration_slope':float(fit.x[1]),
            'calibration_in_the_large_logit':float(cil),'diagnostic_fit_success':bool(fit.success)}



def verify_sources(source: Path) -> dict:
    """Verify the complete retrieval manifest before reading any individual data."""
    manifest=json.loads((source/'manifest.json').read_text())
    checked=0
    for row in manifest['files']:
        name=row['file']
        if Path(name).name!=name: raise ValueError('Manifest path must be a plain filename')
        path=source/name
        if not path.is_file(): raise ValueError(f'Missing source: {name}')
        content=path.read_bytes()
        if len(content)!=row['bytes'] or hashlib.sha256(content).hexdigest()!=row['sha256']:
            raise ValueError(f'Source integrity failure: {name}')
        checked+=1
    return {'checked_files':checked,'errors_in_retrieval_manifest':len(manifest.get('errors',[])),
            'sha256_verified':True}


def research_predict(data: pd.DataFrame, model: dict, horizon: float, *,
                     country: str, acknowledge_research_only: bool=False) -> np.ndarray:
    """Replay a fitted model on engineered features; never a clinical decision API.

    Input transformations are documented in make_cohort and the model card.
    This deliberately requires explicit research acknowledgement and rejects
    non-US populations or ages outside 40-79. Missing lab values use training
    medians; age and published sex coding must not be missing/imputed.
    """
    if acknowledge_research_only is not True:
        raise ValueError('Explicit research-only acknowledgement required')
    if country!='US' or model.get('country')!='US' or model.get('clinical_use') is not False:
        raise ValueError('Model is a historical US research benchmark only')
    age=pd.to_numeric(data['RIDAGEYR'],errors='coerce')
    if not age.between(40,79).all() or not data['male'].isin([0.,1.]).all():
        raise ValueError('Age 40-79 and the published binary sex variable are required')
    if not np.allclose(data['age'],(age-60)/10) or not np.allclose(data['age_male'],data['age']*data['male']):
        raise ValueError('Inconsistent engineered age or sex interactions')
    pp=Preprocessor(); pp.features=model['features']; pp.stats=model['preprocessing'];pp.columns=model['columns']
    return predict(pp.transform(data),model,horizon)


def main(source: Path,out: Path,replicates=300):
    out.mkdir(parents=True,exist_ok=True)
    if isinstance(replicates,bool) or not isinstance(replicates,int) or replicates<100:
        raise ValueError('At least 100 conditional validation resamples required')
    source_audit=verify_sources(source)
    data,flow=make_cohort(source)
    train=data[data.cycle<2003].copy(); test=data[data.cycle==2003].copy()
    if set(train.SEQN)&set(test.SEQN): raise ValueError('Participant leakage')
    weights=train.weight.to_numpy(float); testw=test.weight.to_numpy(float)
    missing=[]
    for name,d in [('training',train),('temporal_validation',test)]:
        for f in PANELS['M4_cystatinC']:
            missing.append({'split':name,'feature':f,'n':len(d),'missing_n':int(d[f].isna().sum())})
    pd.DataFrame(missing).to_csv(out/'missingness.csv',index=False)
    coefficients={}; predictions={}; metrics=[]; subgroups=[]; calibration_bins=[]
    for name,features in PANELS.items():
        pp=Preprocessor().fit(train,features,weights)
        z=pp.transform(train); ztest=pp.transform(test)
        model=fit_piecewise(z,train.time_years.to_numpy(),train.dead.to_numpy(),weights)
        model.update({'columns':pp.columns,'preprocessing':pp.stats,'features':pp.features,
                      'endpoint':'all_cause_death','country':'US','baseline_age_range':[40,79],
                      'interval_ends_years':INTERVAL_ENDS.tolist(),'clinical_use':False})
        coefficients[name]=model; predictions[name]={}
        for horizon in [1,5,10]:
            pred=predict(ztest,model,horizon); predictions[name][horizon]=pred
            y=((test.dead==1)&(test.time_years<=horizon)).to_numpy(int)
            row={'model':name,'horizon_years':horizon,**evaluate(y,pred,testw)}
            row.update(calibration(y,pred,testw)); metrics.append(row)
            for subgroup,mask in [('female',test.male.eq(0)),('male',test.male.eq(1)),
                                  ('age40_59',test.RIDAGEYR.lt(60)),('age60_79',test.RIDAGEYR.ge(60))]:
                subgroups.append({'model':name,'horizon_years':horizon,'subgroup':subgroup,
                                  **evaluate(y[mask],pred[mask],testw[mask])})
            if horizon==10:
                # Weighted quantile bins within validation are descriptive only.
                edges=[weighted_quantile(pred,testw,q) for q in np.linspace(0,1,11)]
                bins=np.searchsorted(edges[1:-1],pred,side='right')
                for b in range(10):
                    mask=bins==b
                    if mask.any(): calibration_bins.append({'model':name,'decile':b+1,**evaluate(y[mask],pred[mask],testw[mask])})
        print(name,metrics[-1],flush=True)
    # Conditional uncertainty: Rao-Wu style resampling of m_h-1 PSUs per stratum,
    # with m_h/(m_h-1) multiplier. Validation only; models fixed, no training CI.
    rng=np.random.default_rng(20260909); rep_rows=[]; y10=((test.dead==1)&(test.time_years<=10)).to_numpy(int)
    strata=[]
    for sid,g in test.groupby('SDMVSTRA'):
        psus=sorted(g.SDMVPSU.unique())
        if len(psus)<2: raise ValueError('Cannot resample single-PSU stratum')
        strata.append((float(sid),psus))
    for b in range(replicates):
        mult=np.zeros(len(test))
        for sid,psus in strata:
            m=len(psus); chosen=rng.choice(psus,size=m-1,replace=True)
            for psu in chosen:
                mask=(test.SDMVSTRA.to_numpy()==sid)&(test.SDMVPSU.to_numpy()==psu)
                mult[mask]+=m/(m-1)
        wr=testw*mult
        values={name:evaluate(y10,predictions[name][10],wr) for name in PANELS}
        names=list(PANELS)
        for j,name in enumerate(names):
            row={'replicate':b,'model':name,**{key:values[name][key] for key in ['auc','brier','weighted_observed','weighted_predicted']}}
            if name not in ('M0_age_sex','M1_routine'):
                row['delta_auc_routine']=values[name]['auc']-values['M1_routine']['auc']
                row['delta_brier_routine']=values[name]['brier']-values['M1_routine']['brier']
            if j:
                prev=names[j-1];row['delta_auc_previous']=values[name]['auc']-values[prev]['auc'];row['delta_brier_previous']=values[name]['brier']-values[prev]['brier']
            rep_rows.append(row)
    rep=pd.DataFrame(rep_rows); intervals=[]
    for name,g in rep.groupby('model',sort=False):
        for key in ['auc','brier','weighted_observed','weighted_predicted','delta_auc_previous','delta_brier_previous','delta_auc_routine','delta_brier_routine']:
            v=g[key].dropna()
            if len(v): intervals.append({'model':name,'horizon_years':10,'metric':key,
                                        'p025':float(v.quantile(.025)),'p975':float(v.quantile(.975)),
                                        'replicates':len(v),'type':'conditional_validation_PSU_resampling_percentiles'})
    for name,rows in [('metrics',metrics),('subgroups',subgroups),('calibration_deciles',calibration_bins),('validation_intervals',intervals),('cohort_flow',flow)]:
        pd.DataFrame(rows).to_csv(out/f'{name}.csv',index=False)
    (out/'fitted_models.json').write_text(json.dumps(coefficients,ensure_ascii=False,indent=2,allow_nan=False))
    summaries=[]
    for label,d in [('training',train),('validation',test)]:
        row={'split':label,'n':len(d),'female_n':int(d.male.eq(0).sum()),
             'age_min':int(d.RIDAGEYR.min()),'age_max':int(d.RIDAGEYR.max()),
             'zero_month_deaths':int(((d.dead==1)&(d.t_exam==0)).sum()),
             'min_survivor_followup_years':float(d.loc[d.dead==0,'time_years'].min()),
             'weighted_age_mean':float(np.average(d.RIDAGEYR,weights=d.weight)),
             'weight_effective_n':float(d.weight.sum()**2/(d.weight**2).sum()),
             'ntprobnp_below_LOD_n':int(d.SSBNPL.eq(1).sum()),'ntprobnp_above_ULOQ_n':int(d.SSBNPL.eq(2).sum()),
             'troponin_initially_blinded_n':int(d.SSTNTLC.isin([1,2]).sum()),
             'nonpristine_ntprobnp_n':int(d.SSPRIS.eq(0).sum()),'nonpristine_troponin_n':int(d.SSPRISTP.eq(0).sum())}
        for h in [1,5,10]: row[f'deaths_{h}y']=int(((d.dead==1)&(d.time_years<=h)).sum())
        summaries.append(row)
    status={'version':'0.4','status':'exploratory_temporal_validation_not_clinical',
            'source_integrity':source_audit,
            'source_manifest_sha256':hashlib.sha256((source/'manifest.json').read_bytes()).hexdigest(),
            'cohort':summaries,'ridge_penalty':PENALTY,'test_model_refitting':False,
            'validation_strata':len(strata),'conditional_uncertainty_replicates':replicates,
            'limitations':['all-cause endpoint includes infection and external causes',
              'old US civilian noninstitutionalized sample; not calibrated for 2026, RU or DE',
              'public follow-up/causes partially perturbed; vital status not perturbed',
              'single temporal holdout; no independent geographic validation',
              'piecewise PH assumption and age-sex interaction limited to stated basis',
              'specimen and complete-assay selection may leave residual selection bias',
              'conditional intervals exclude training/model-selection uncertainty',
              'HDL assay change; no invented correction applied',
              'low troponin values are research-only; stored specimen assay transfer unvalidated'],
            'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}}
    (out/'status.json').write_text(json.dumps(status,indent=2,allow_nan=False))
    # Keep participant-level predictions in memory only, not in exported artifacts.
    return data,coefficients,predictions

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--replicates',type=int,default=300)
    args=parser.parse_args();main(args.source,args.out,args.replicates)
