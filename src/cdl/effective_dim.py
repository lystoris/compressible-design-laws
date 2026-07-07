import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor

def effective_d(X,y,groups,seed=0):
    rf=RandomForestRegressor(n_estimators=300,random_state=seed,n_jobs=-1).fit(X,y)
    imp=rf.feature_importances_.astype(float)
    g=pd.Series(imp).groupby(pd.Series(groups)).sum()
    s=g.values; s=s/s.sum() if s.sum()>0 else s
    return float(1.0/np.sum(s**2)) if s.sum()>0 else float(g.size)
