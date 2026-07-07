import numpy as np
import pandas as pd
from scipy import stats

def r2_det(yt, yp):
    yt = np.asarray(yt, float); yp = np.asarray(yp, float)
    ssr = np.sum((yt - yp) ** 2); sst = np.sum((yt - yt.mean()) ** 2)
    return np.nan if sst <= 0 else 1.0 - ssr / sst

def partial_spearman(a, b, c):
    ra, rb, rc = [stats.rankdata(v) for v in (a, b, c)]
    Z = np.column_stack([np.ones_like(rc), rc])
    res = lambda u: u - Z @ np.linalg.lstsq(Z, u, rcond=None)[0]
    return float(np.corrcoef(res(ra), res(rb))[0, 1])

def eta2(values, groups):
    s = pd.Series(np.asarray(values, float)); g = pd.Series(list(groups))
    grand = s.mean()
    ss_between = s.groupby(g).apply(lambda x: len(x) * (x.mean() - grand) ** 2).sum()
    ss_total = float(np.sum((s - grand) ** 2))
    return float(ss_between / ss_total) if ss_total > 0 else 0.0
