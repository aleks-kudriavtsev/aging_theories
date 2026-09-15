"""Penalized cause-specific piecewise hazards, with same-training partial pooling.

Sex-stratified baseline hazards; interval ends 5/10/20 years. All probabilities
are jointly integrated. No literature hazard ratios become fitted coefficients.
"""
import numpy as np
from scipy.optimize import minimize
ENDS=np.array([5.,10.,20.]);STARTS=np.r_[0.,ENDS[:-1]]

def preprocessing(frame,features):
    x=frame[features].to_numpy(float)
    if not np.isfinite(x).all():raise ValueError('Complete observed profiles required')
    lo,hi=np.quantile(x,[.005,.995],axis=0);x=np.clip(x,lo,hi)
    mean=x.mean(axis=0);sd=x.std(axis=0);sd[sd==0]=1.
    return dict(features=list(features),lower=lo.tolist(),upper=hi.tolist(),mean=mean.tolist(),sd=sd.tolist())

def design(frame,pp):
    x=frame[pp['features']].to_numpy(float)
    if not np.isfinite(x).all():raise ValueError('Missing or nonfinite model input')
    return (np.clip(x,pp['lower'],pp['upper'])-pp['mean'])/pp['sd']

def fit_one(z,time,sex,event,weight,penalty,beta_target=None,alpha_target=None):
    n,p=z.shape;sex=np.asarray(sex,int)-1;time=np.asarray(time,float);e=np.asarray(event,bool)
    w=np.asarray(weight,float)
    if not np.isin(sex,[0,1]).all() or np.any(time<=0) or np.any(w<=0) or not np.isfinite(z).all():raise ValueError('Invalid training inputs')
    duration=np.maximum(0,np.minimum(time[:,None],ENDS)-STARTS)
    exposure=np.zeros((n,6));exposure[np.arange(n)[:,None],sex[:,None]*3+np.arange(3)]=duration
    y=np.zeros((n,6));h=np.minimum(np.searchsorted(ENDS,time,side='left'),2)
    idx=np.flatnonzero(e&(time<=20));y[idx,sex[idx]*3+h[idx]]=1.
    if beta_target is None:beta_target=np.zeros(p)
    if alpha_target is None:
        alpha_target=np.log(((w[:,None]*y).sum(0)+.5)/((w[:,None]*exposure).sum(0)+.5))
    target=np.asarray(alpha_target);bt=np.asarray(beta_target)
    def objective(theta):
        alpha=theta[:6];beta=theta[6:];lp=z@beta;rate=np.exp(alpha);mult=np.exp(lp)
        expected=mult[:,None]*exposure*rate
        residual=w[:,None]*(expected-y)
        val=np.sum(w*(expected.sum(1)-y@alpha-y.sum(1)*lp))+.125*np.sum((alpha-target)**2)+penalty/2*np.sum((beta-bt)**2)
        grad=np.r_[residual.sum(0)+.25*(alpha-target),z.T@residual.sum(1)+penalty*(beta-bt)]
        return float(val),grad
    initial=np.r_[np.clip(target,-25,5),np.clip(bt,-5,5)]
    bounds=[(-25,5)]*6+[(-5,5)]*p
    result=minimize(objective,initial,jac=True,method='L-BFGS-B',bounds=bounds,
        options={'maxiter':1500,'ftol':1e-15,'gtol':1e-7,'maxls':50})
    gradient=objective(result.x)[1]
    projected=gradient.copy()
    for j,(low,high) in enumerate(bounds):
        if (result.x[j]<=low+1e-8 and gradient[j]>0) or (result.x[j]>=high-1e-8 and gradient[j]<0):projected[j]=0
    pg=float(np.max(abs(projected)))
    if pg>1e-3:raise ValueError('Optimizer convergence failure: '+str(pg))
    return {'baseline_log_hazards':result.x[:6].reshape(2,3).tolist(),'coefficients':result.x[6:].tolist(),
        'raw_events':int(y.sum()),'penalty':penalty,'projected_gradient_max':pg,'iterations':int(result.nit),
        'boundary_parameters':int(sum(abs(result.x[j]-low)<1e-6 or abs(result.x[j]-high)<1e-6 for j,(low,high) in enumerate(bounds)))}

def fit(frame,features,causes,penalty):
    pp=preprocessing(frame,features);z=design(frame,pp);w=frame.weight.to_numpy(float);w=w/w.mean()
    time=frame.time.to_numpy(float);sex=frame.sex.to_numpy(int);dead=frame.dead.to_numpy(bool)
    overall=fit_one(z,time,sex,dead,w,penalty)
    models={};total=float(w@(dead&(time<=20)))
    for c in causes:
        event=frame.model_cause.eq(c).to_numpy()&dead
        mask=sex==1 if c=='prostate_cancer' else np.ones(len(frame),bool)
        share=(float(w@(event&(time<=20)))+.5)/(total+.5*len(causes))
        anchor=np.asarray(overall['baseline_log_hazards']).ravel()+np.log(share)
        models[c]=fit_one(z[mask],time[mask],sex[mask],event[mask],w[mask],penalty,overall['coefficients'],anchor)
        models[c]['female_structural_zero']=c=='prostate_cancer'
    return {'preprocessing':pp,'causes':models,'cause_order':list(causes),'interval_ends':ENDS.tolist(),
        'allcause_pooling_target':overall,'training_n':len(frame),'penalty':penalty,
        'clinical_use_ready':False,'training_country':'historical_US','age_range':[25,74]}

def predict(frame,model,horizon):
    if not isinstance(horizon,(float,int)) or not 0<horizon<=20:raise ValueError('Unsupported horizon')
    z=design(frame,model['preprocessing']);sidx=frame.sex.to_numpy(int)-1
    if not np.isin(sidx,[0,1]).all():raise ValueError('Unknown sex')
    rates=[]
    for c in model['cause_order']:
        m=model['causes'][c];r=np.exp((z@np.array(m['coefficients']))[:,None]+np.array(m['baseline_log_hazards'])[sidx])
        if m['female_structural_zero']:r[sidx==1]=0
        rates.append(r)
    rates=np.stack(rates,axis=1);s=np.ones(len(frame));cif=np.zeros((len(frame),len(rates[0])))
    for j,(lo,hi) in enumerate(zip(STARTS,ENDS)):
        dt=max(0.,min(horizon,hi)-lo);total=rates[:,:,j].sum(1);mass=-np.expm1(-total*dt)
        cif+=s[:,None]*mass[:,None]*rates[:,:,j]/total[:,None];s*=np.exp(-total*dt)
    p=np.column_stack([s,cif])
    if not np.isfinite(p).all() or np.any(p<0) or not np.allclose(p.sum(1),1,atol=1e-12):raise ArithmeticError('Probability integrity failure')
    return p
