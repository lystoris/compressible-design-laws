#!/usr/bin/env python3
"""SI-A calibration: does the independent RF participation-ratio effective-d estimator
(used throughout Fig 3/4) recover the TRUE effective dimensionality?

Standalone check: generate synthetic landscapes with a KNOWN number of active design
variables (d_eff) embedded among d_nom total variables (the rest inert decoys), across
three functional families, then apply the identical estimator used in the population panel
(cdl.effective_dim.effective_d / ../paper/v2 round-08 panel.py::effective_d) and record
estimated vs set. Fully synthetic -- no dependency on results/ or data/.

Run with the system Python that has scikit-learn: /usr/bin/python3 scripts/build_calibration_data.py
Writes: figures/data/figS_calibration.csv
"""
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "figures", "data", "figS_calibration.csv")

# --- the estimator, identical to cdl.effective_dim.effective_d / round-08 panel.py --------
def effective_d(X, y, groups, seed=0):
    rf = RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=-1).fit(X, y)
    imp = rf.feature_importances_.astype(float)
    g = pd.Series(imp).groupby(pd.Series(groups)).sum()
    s = g.values
    s = s / s.sum() if s.sum() > 0 else s
    return float(1.0 / np.sum(s ** 2)) if s.sum() > 0 else float(g.size)

# --- landscape generator: d_eff active vars among d_nom, rest are inert decoys ----------
D_NOM = 12
D_EFF_SET = [2, 3, 5, 9]        # true effective d (must be <= D_NOM)
FAMILIES = ["additive", "michaelis_menten", "random_gp"]
N = 500                          # designs per landscape
REPS = 25                        # replicate landscapes per (family, d_eff)
NOISE = 0.10                     # Gaussian noise as fraction of signal SD

def make_y(family, Xa, rng):
    d_eff = Xa.shape[1]
    if family == "additive":
        return Xa @ np.ones(d_eff)                       # equal-weight linear
    if family == "michaelis_menten":
        K = 0.5
        return (Xa / (K + Xa)).sum(axis=1)               # saturating (Xa in (0,1))
    if family == "random_gp":
        M = 40                                           # random Fourier features
        omega = rng.normal(size=(d_eff, M)) * 1.5
        b = rng.uniform(0, 2 * np.pi, size=M)
        a = rng.normal(size=M)
        return np.cos(Xa @ omega + b) @ a                # smooth random function
    raise ValueError(family)


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    for family in FAMILIES:
        for d_eff in D_EFF_SET:
            for rep in range(REPS):
                seed = hash((family, d_eff, rep)) % (2 ** 31)
                rng = np.random.default_rng(seed)
                X = rng.random((N, D_NOM))                # U(0,1); cols >= d_eff are decoys
                y = make_y(family, X[:, :d_eff], rng)
                y = y + rng.normal(0, NOISE * (y.std() + 1e-9), N)
                groups = np.arange(D_NOM)                 # each variable is its own group
                est = effective_d(X, y, groups, seed=seed % 10000)
                rows.append(dict(family=family, d_nom=D_NOM, d_eff_set=d_eff,
                                 rep=rep, eff_d_est=round(est, 3)))
            sub = [r["eff_d_est"] for r in rows if r["family"] == family and r["d_eff_set"] == d_eff]
            print(f"{family:16s} set={d_eff:2d}  est median={np.median(sub):.2f} "
                  f"[{np.min(sub):.2f}, {np.max(sub):.2f}]  (n_decoys={D_NOM - d_eff})")

    pd.DataFrame(rows).to_csv(OUT, index=False)
    print("wrote", OUT, "rows:", len(rows))


if __name__ == "__main__":
    main()
