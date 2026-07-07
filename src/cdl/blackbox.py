import numpy as np
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from cdl.stats import r2_det

N_INNER=3

def tune_svr(Xtr,ytr,seed=0,fast=False):
    Cs=[1.,100.] if fast else [1.,10.,100.,1000.]; Gs=["scale",.1] if fast else ["scale",.01,.1,1.]
    Es=[.1] if fast else [.01,.1,.2]
    inner=KFold(min(N_INNER,max(2,len(ytr)//3)),shuffle=True,random_state=seed)
    best=None; bs=-np.inf
    for C in Cs:
        for gm in Gs:
            for e in Es:
                sc=[]
                for a,b in inner.split(Xtr):
                    xs=StandardScaler().fit(Xtr[a]); ym=ytr[a].mean(); ys=ytr[a].std() or 1
                    m=SVR(kernel="rbf",C=C,gamma=gm,epsilon=e).fit(xs.transform(Xtr[a]),(ytr[a]-ym)/ys)
                    sc.append(r2_det(ytr[b],m.predict(xs.transform(Xtr[b]))*ys+ym))
                s=np.nanmean(sc)
                if s>bs: bs=s; best=(C,gm,e)
    return dict(C=best[0],gamma=best[1],epsilon=best[2])
def svr_pred(Xtr,ytr,Xte,p):
    xs=StandardScaler().fit(Xtr); ym=ytr.mean(); ys=ytr.std() or 1
    m=SVR(kernel="rbf",C=p["C"],gamma=p["gamma"],epsilon=p["epsilon"]).fit(xs.transform(Xtr),(ytr-ym)/ys)
    return m.predict(xs.transform(Xte))*ys+ym
