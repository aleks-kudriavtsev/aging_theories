"""IPCW scoring for discrete-month administrative censoring, research only.

Weighted reverse Kaplan-Meier; deaths occur before censors tied at a month.
Censoring at the horizon is observed through that horizon. Censoring before it
has zero outcome weight, not a negative event label. No coefficient estimation.
Method: Graf et al., PMID10474158; Gerds & Schumacher, PMID17240660.
"""
import numpy as np
from .validate_temporal08 import weighted_auc


def outcome_weights(dead, cause, time, weights, horizon):
    d,c,t,w=map(lambda v:np.asarray(v,float),(dead,cause,time,weights))
    if (d.ndim!=1 or not len(d) or any(x.shape!=d.shape for x in (c,t,w))
        or not np.isin(d,[0,1]).all() or not np.isfinite(t).all() or np.any(t<0)
        or not np.isfinite(w).all() or np.any(w<0) or w.sum()<=0 or np.isinf(c).any()
        or type(horizon) not in (int,float) or not np.isfinite(horizon) or horizon<=0):
        raise ValueError('Invalid IPCW observations, weights or horizon')
    known=np.isfinite(c)
    if np.any((d==0)&known) or not np.isin(c[known],[1,2,10]).all():
        raise ValueError('Invalid three-category public cause partition')
    times,index=np.unique(t,return_inverse=True)
    total=np.bincount(index,weights=w,minlength=len(times))
    dw=np.bincount(index,weights=w*d,minlength=len(times))
    cw=total-dw
    risk=np.cumsum(total[::-1])[::-1]
    before=np.ones(len(times));g=1.
    for j,u in enumerate(times):
        before[j]=g
        if u>=horizon:break
        censor_risk=risk[j]-dw[j] # Events first at ties.
        if cw[j]>1e-12:
            if censor_risk<=0:raise ValueError('Invalid censoring risk set')
            g*=max(0.,1.-cw[j]/censor_risk)
    if g<=1e-12:raise ValueError('Censoring distribution lacks horizon support')
    event=(d==1)&(t<=horizon)
    status_known=event|(t>=horizon)
    factors=np.zeros(len(d))
    event_g=before[index[event]]
    if np.any(event_g<=0):raise ValueError('Unidentified event probability')
    factors[event]=1./event_g
    factors[~event&(t>=horizon)]=1./g
    cause_class=np.where(c==1,1,np.where(c==2,2,3))
    target=np.where(event,cause_class,0).astype(int) # Unknown statuses have zero factor.
    effective=w*factors
    if abs(effective.sum()/w.sum()-1)>1e-9:
        raise ArithmeticError('IPCW mass balance failed under discrete tie convention')
    return {'target':target,'effective_weights':effective,'ipcw_factors':factors,
            'known':status_known,'G_horizon_left':float(g),
            'early_censored_n':int((~status_known).sum()),'denominator':float(w.sum())}


def brier_states(target, probabilities, effective_weights, denominator):
    p=np.asarray(probabilities,float);y=np.asarray(target,int);w=np.asarray(effective_weights,float)
    if (p.ndim!=2 or len(p)!=len(y) or w.shape!=y.shape or not np.isfinite(p).all()
        or np.any((p<0)|(p>1)) or not np.allclose(p.sum(1),1,atol=1e-10,rtol=0)
        or np.any(y<0) or np.any(y>=p.shape[1]) or not np.isfinite(w).all() or np.any(w<0)
        or not np.isfinite(denominator) or denominator<=0):
        raise ValueError('Invalid Brier inputs')
    loss=np.sum((p-np.eye(p.shape[1])[y])**2,axis=1)
    return float(w@loss/denominator)


def binary_scores(target,risk,weights,outcomes,k=0):
    y=target==k if k else target>0
    p=np.asarray(risk,float);w=np.asarray(weights,float);ew=outcomes['effective_weights'];den=w.sum()
    if p.shape!=y.shape or not np.isfinite(p).all() or np.any((p<0)|(p>1)):
        raise ValueError('Invalid risk vector')
    valid=(ew>0)
    auc=weighted_auc(y,p,ew) if np.unique(y[valid]).size==2 else None
    observed=float(ew@y/den);predicted=float(w@p/den)
    return {'n':len(y),'events':int(y.sum()),'early_censored_n':outcomes['early_censored_n'],
      'auc':auc,'brier':float(ew@((y-p)**2)/den),'observed':observed,'predicted':predicted,
      'observed_expected':observed/predicted if predicted>0 else None,'bias_pp':100*(predicted-observed),
      'G_horizon_left':outcomes['G_horizon_left'],'metric_weighting':'survey_times_IPCW'}


def aalen_johansen_reference(dead,cause,time,weights,horizon):
    """Independent weighted product-limit CIF; all deaths retained as events."""
    d,c,t,w=map(lambda a:np.asarray(a,float),(dead,cause,time,weights))
    survival=1.;cif=np.zeros(3)
    for u in sorted(set(t[t<=horizon])):
        n=w[t>=u].sum()
        if n<=0:continue
        step=(t==u)&(d==1);counts=np.array([w[step&(c==1)].sum(),w[step&(c==2)].sum(),w[step&~np.isin(c,[1,2])].sum()])
        cif+=survival*counts/n
        survival*=1-counts.sum()/n
    return np.r_[survival,cif]
