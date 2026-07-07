import itertools
import numpy as np
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from cdl.stats import r2_det
from cdl.blackbox import tune_svr, svr_pred

L_STAR=10; N_OUTER=5

def build_terms(d,fn,allow_nonlinear=False,g=1e-6):
    T=[dict(name=fn[i],cost=1,fn=(lambda X,i=i:X[:,i])) for i in range(d)]
    if allow_nonlinear:
        T+=[dict(name=f"1/{fn[i]}",cost=3,fn=(lambda X,i=i:1.0/(X[:,i]+g))) for i in range(d)]
        T+=[dict(name=f"{fn[i]}^2",cost=3,fn=(lambda X,i=i:X[:,i]**2)) for i in range(d)]
    T+=[dict(name=f"{fn[i]}*{fn[j]}",cost=3,fn=(lambda X,i=i,j=j:X[:,i]*X[:,j]))
        for i,j in itertools.combinations(range(d),2)]
    if allow_nonlinear:
        T+=[dict(name=f"{fn[i]}/{fn[j]}",cost=3,fn=(lambda X,i=i,j=j:X[:,i]/(X[:,j]+g)))
            for i,j in itertools.permutations(range(d),2)]
    return T
def _design(X,ts):
    M=np.column_stack([t["fn"](X) for t in ts]) if ts else np.empty((X.shape[0],0))
    return np.where(np.isfinite(M),M,0.0)
def greedy_sr_fit(Xtr,ytr,terms,budget=L_STAR):
    sel=[]; rem=list(terms)
    nodes=lambda s:sum(t["cost"] for t in s)+len(s)
    best=r2_det(ytr,np.full_like(ytr,ytr.mean(),dtype=float))
    while True:
        cb=None; cbr=best
        for t in rem:
            tr=sel+[t]
            if nodes(tr)>budget: continue
            M=_design(Xtr,tr); mu=M.mean(0); sd=M.std(0); sd=np.where(sd<1e-12,1,sd)
            try: lr=LinearRegression().fit((M-mu)/sd,ytr); r=r2_det(ytr,lr.predict((M-mu)/sd))
            except Exception: continue
            if r>cbr+1e-9: cbr=r; cb=t
        if cb is None: break
        sel.append(cb); rem.remove(cb); best=cbr
        if nodes(sel)>=budget: break
    if not sel: return [],None,None,float(ytr.mean())
    M=_design(Xtr,sel); mu=M.mean(0); sd=M.std(0); sd=np.where(sd<1e-12,1,sd)
    return sel,mu,sd,LinearRegression().fit((M-mu)/sd,ytr)
def greedy_sr_pred(Xte,sel,mu,sd,lr):
    if not sel: return np.full(Xte.shape[0],lr,dtype=float)
    return lr.predict((_design(Xte,sel)-mu)/sd)
def screen(Xtr,ytr,k):
    sc=np.array([abs(stats.spearmanr(Xtr[:,j],ytr).correlation) if np.std(Xtr[:,j])>1e-12 else 0
                 for j in range(Xtr.shape[1])]); sc=np.nan_to_num(sc)
    return np.argsort(-sc)[:k]
def run_anchor(X,y,feat,allow_nonlinear=False,topk=40,fast=False):
    n,d=X.shape; nf=min(N_OUTER,n) if n>=N_OUTER else max(2,n//2)
    kf=KFold(nf,shuffle=True,random_state=0); osr=np.full(n,np.nan); obb=np.full(n,np.nan)
    for a,b in kf.split(X):
        Xtr,Xte,ytr=X[a],X[b],y[a]
        if topk and d>topk:
            kp=screen(Xtr,ytr,topk); Xtr2,Xte2=Xtr[:,kp],Xte[:,kp]; fn=[feat[i] for i in kp]
        else: Xtr2,Xte2,fn=Xtr,Xte,feat
        sel,mu,sd,lr=greedy_sr_fit(Xtr2,ytr,build_terms(Xtr2.shape[1],fn,allow_nonlinear))
        osr[b]=greedy_sr_pred(Xte2,sel,mu,sd,lr)
        obb[b]=svr_pred(Xtr2,ytr,Xte2,tune_svr(Xtr2,ytr,fast=fast))
    return float(r2_det(y,osr)), float(r2_det(y,obb))
