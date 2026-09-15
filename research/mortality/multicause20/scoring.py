"""IPCW scoring for arbitrary mutually exclusive causes, survey-domain bootstrap."""
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.special import expit,logit

def outcomes(time,dead,cause,w,horizon,k):
    t=np.asarray(time,float);d=np.asarray(dead,int);c=np.asarray(cause,int);w=np.asarray(w,float)
    if len(t)==0 or np.any(w<0) or w.sum()<=0 or np.any(t<=0) or not np.isfinite(t).all() or not np.isin(d,[0,1]).all():raise ValueError('Invalid observations')
    if np.any((d==1)&((c<1)|(c>k))):raise ValueError('All deaths need a cause including explicit residual')
    times,ind=np.unique(t,return_inverse=True)
    total=np.bincount(ind,weights=w);dw=np.bincount(ind,weights=w*d);cw=total-dw;risk=np.cumsum(total[::-1])[::-1]
    before=np.ones(len(times));g=1.
    for j,u in enumerate(times):
        before[j]=g
        if u>=horizon:break
        denom=risk[j]-dw[j]
        if cw[j]>1e-10:
            if denom<=0:raise ValueError('No censoring support')
            g*=max(0.,1-cw[j]/denom)
    if g<1e-6:raise ValueError('Censoring positivity failure')
    ev=(d==1)&(t<=horizon);known=ev|(t>=horizon);factor=np.zeros(len(t))
    factor[ev]=1/before[ind[ev]];factor[~ev&(t>=horizon)]=1/g
    y=np.where(ev,c,0);ew=w*factor
    if abs(ew.sum()/w.sum()-1)>1e-8:raise ArithmeticError('IPCW mass balance')
    return y,ew,{'early_censored_n':int((~known).sum()),'G_horizon':float(g)}

def aalen_johansen(time,dead,cause,w,h,k):
    t=np.asarray(time);d=np.asarray(dead);c=np.asarray(cause);w=np.asarray(w);s=1.;cif=np.zeros(k)
    for u in np.unique(t[(t<=h)&(d==1)]):
        risk=w[t>=u].sum();ev=(t==u)&(d==1)
        by=np.bincount(c[ev],weights=w[ev],minlength=k+1)[1:]
        if risk>0:cif+=s*by/risk;s*=1-by.sum()/risk
    return np.r_[s,cif]

def auc(y,p,w):
    order=np.argsort(p,kind='stable');ps=np.asarray(p)[order];ys=np.asarray(y)[order];ws=np.asarray(w)[order]
    start=np.r_[0,np.flatnonzero(np.diff(ps))+1];pos=np.add.reduceat(ws*ys,start);neg=np.add.reduceat(ws*(1-ys),start)
    if pos.sum()<=0 or neg.sum()<=0:return None
    return float((pos*(np.cumsum(neg)-neg+.5*neg)).sum()/(pos.sum()*neg.sum()))

def metrics(y,p,w,ew):
    den=w.sum();loss=(p*p).sum(1)-2*p[np.arange(len(y)),y]+1
    rows=[{'outcome':'multistate','brier':float(ew@loss/den)}]
    for j in range(p.shape[1]):
        yes=y>0 if j==0 else y==j;risk=1-p[:,0] if j==0 else p[:,j]
        obs=float(ew@yes/den);pred=float(w@risk/den)
        rows.append({'outcome':'all_cause' if j==0 else j,'brier':float(ew@((yes-risk)**2)/den),
            'auc':auc(yes,risk,ew),'observed':obs,'predicted':pred,'observed_expected':obs/pred if pred>0 else None,
            'mean_bias_pp':100*(pred-obs),'events':int(yes.sum())})
    return rows

def calibration(y,risk,ew):
    if np.unique(y[ew>0]).size<2:return {'status':'insufficient_events'}
    w=ew/ew.mean();x=logit(np.clip(risk,1e-10,1-1e-10));X=np.c_[np.ones(len(x)),x]
    def fun(b):
        lp=X@b;return float(w@(np.logaddexp(0,lp)-y*lp)),X.T@(w*(expit(lp)-y))
    r=minimize(fun,[0.,1.],jac=True,method='BFGS',options={'maxiter':1000,'gtol':1e-5})
    if np.max(abs(fun(r.x)[1]))>1e-3:return {'status':'not_converged'}
    cil=brentq(lambda a:w@(expit(x+a)-y),-40,40)
    return {'status':'computed','intercept':float(r.x[0]),'slope':float(r.x[1]),'calibration_in_large':float(cil)}

def resample(frame,rng):
    st=frame.stratum.to_numpy();psu=frame.psu.to_numpy();mult=np.zeros(len(frame))
    for s in np.unique(st):
        idx=np.flatnonzero(st==s);units=np.unique(psu[idx]);m=len(units)
        if m<2:raise ValueError('Singleton design stratum')
        sampled=rng.choice(units,m-1,replace=True)
        for u in units:mult[idx[psu[idx]==u]]=np.sum(sampled==u)*m/(m-1)
    return frame.weight.to_numpy(float)*mult
