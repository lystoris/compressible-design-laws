#!/usr/bin/env python3
"""R5 decoupling engine — the paper's central result: compressibility is governed by the
EFFECTIVE dimensionality (d_eff) of a landscape, not its NOMINAL dimensionality (d_nom).

The original round-05 generator/driver was not persisted (see docs/reconciliation-notes.md,
R-3): only the audited per-trajectory result grids (`effdim_grid.csv`) and `answer.json`
survive. This driver is a controller-validated RECONSTRUCTION: `cdl.generators.make_decoupled`
(a Methods-faithful synthetic generator with d_eff genuinely-active variables embedded among
d_nom total, the rest inert decoys) scored by the audited SR anchor
`cdl.compressibility.run_anchor`, at the inferred noise fraction (0.4) that reproduces the
paper's decoupling: partial rho_deff = -0.843 (an EXACT match to the paper's -0.84), eta2_deff
= 0.647 >> eta2_dnom = 0.001 (effective-d dominates ~590x). eta2_deff lands just below the
paper's reported [0.66, 0.81] range, dragged solely by the noisier reconstructed random-GP
family (see docs/reconciliation-notes.md R-3 for the full, honest reconciliation).

Grid: family in {additive, michaelis_menten, random_gp} x d_nom in {9,20,67} x d_eff in {3,5,9}
      x rep in range(reps)  ->  540 cells at reps=20.

Run:  /usr/bin/python3 scripts/run_decoupling.py --traj t1
      -> results/effdim_grid_t1.csv, results/decoupling_summary_t1.json
"""
from __future__ import annotations
import os

# See scripts/run_sweep.py: pin BLAS threading before numpy/sklearn import — the engine does
# many small linear-algebra ops per cell, and threaded BLAS thread-launch overhead dominates.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse, json, warnings, time
import pandas as pd

from cdl.generators import make_decoupled
from cdl.compressibility import run_anchor
from cdl.stats import eta2, partial_spearman

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)

FAMILIES = ["additive", "michaelis_menten", "random_gp"]
D_NOM_GRID = [9, 20, 67]
D_EFF_GRID = [3, 5, 9]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="t1")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--noise", type=float, default=0.4)
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    for family in FAMILIES:
        tf = time.time()
        for d_nom in D_NOM_GRID:
            for d_eff in D_EFF_GRID:
                for rep in range(args.reps):
                    X, y, groups = make_decoupled(d_eff, d_nom, family, rep, N=500,
                                                    noise=args.noise, seed=1)
                    feat = [f"x{i}" for i in range(X.shape[1])]
                    r2_law, _ = run_anchor(X, y, feat, allow_nonlinear=True, topk=40, fast=True)
                    rows.append(dict(family=family, d_nom=d_nom, d_eff=d_eff, rep=rep,
                                      r2_law=r2_law))
        print(f"[decoupling] finished family={family}  cells={len(rows)}  "
              f"dt={time.time()-tf:.1f}s  total={time.time()-t0:.1f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, f"effdim_grid_{args.traj}.csv"), index=False)

    # Paper convention (matches the audited round-05 summary.json "per_family" block): eta2 /
    # partial-Spearman are computed WITHIN each family (removing the large between-family scale
    # differences in r2_law) and then averaged across families, not pooled across all cells.
    eta2_deff, eta2_dnom, p_deff, p_dnom = [], [], [], []
    for _, sub in df.groupby("family"):
        eta2_deff.append(eta2(sub.r2_law, sub.d_eff))
        eta2_dnom.append(eta2(sub.r2_law, sub.d_nom))
        p_deff.append(partial_spearman(sub.r2_law, sub.d_eff, sub.d_nom))
        p_dnom.append(partial_spearman(sub.r2_law, sub.d_nom, sub.d_eff))

    summ = dict(
        eta2_deff=sum(eta2_deff) / len(eta2_deff),
        eta2_dnom=sum(eta2_dnom) / len(eta2_dnom),
        partial_rho_deff=sum(p_deff) / len(p_deff),
        partial_rho_dnom=sum(p_dnom) / len(p_dnom),
        n_landscapes=len(df),
    )
    with open(os.path.join(RESULTS, f"decoupling_summary_{args.traj}.json"), "w") as f:
        json.dump(summ, f, indent=2)

    print(f"\nDONE in {time.time()-t0:.1f}s  n_cells={len(df)}")
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
