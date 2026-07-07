#!/usr/bin/env python3
"""Synthetic 405-cell sweep driver (Fig 3A/B) — compressibility Lambda vs dimensionality d,
epistasis k, noise sigma, sample size N, and encoding.

Ported VERBATIM (cell ordering + seeds) from round-02 trajectory t3's audited driver, `main()`'s
sweep factorial:
research-loop/sr-compressible-design-laws/rounds/round-02/runs/round-02-t3/work/analysis.py:482-564

Scores each cell with the DEPTH-budgeted symbolic-regression engine under nested CV
(cdl.depth_sr.lambda_depth_nestedcv) against a KernelRidge(RBF) black-box ceiling
(cdl.depth_sr.krr_cv_r2), then evaluates top-1/top-5 selection regret + Spearman on an
independent 400-design selection pool (cdl.depth_sr.selection_fitness).

Run:  /usr/bin/python3 scripts/run_sweep.py --traj t3   -> results/sweep_grid.csv (405 rows)
"""
from __future__ import annotations
import os

# The engine does thousands of small linear-algebra ops (per-fold least-squares, small KRR
# fits) inside Python-level loops; threaded BLAS spawns worker threads for each tiny op and
# thread-launch overhead dominates, making the sweep ~10x slower than single-threaded BLAS.
# Must be set before numpy/sklearn are imported.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse, warnings, time
import numpy as np
import pandas as pd
from scipy import stats

from cdl.generators import make_generator
from cdl.depth_sr import (
    krr_cv_r2, krr_fit_predict, lambda_depth_nestedcv, selection_fitness,
    screen_keep, standardize_cols, build_feature_bank, greedy_law, predict_law,
)

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)

# Seed pinned to 0 (matches the sandbox launcher that pinned round-02-t3's RNG); trajectory
# tag only labels the output, per the source's convention across t1/t2/t3 drivers.
RNG = np.random.RandomState(0)

D_GRID = [3, 5, 9, 20, 67]
K_GRID = [1, 2, 3]
SIG_GRID = [0.0, 0.1, 0.3]
N_GRID = [30, 100, 300]
ENC_GRID = ["integer", "continuous", "binary"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="t3")
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    cell_id = 0
    for d in D_GRID:
        d_screen = min(6, d)
        td = time.time()
        for k in K_GRID:
            for enc in ENC_GRID:
                sample_X, true_f = make_generator(d, k, enc, 1000 + len(rows))
                Xpool = sample_X(400, seed=777)               # independent selection design pool
                ypool_true = true_f(Xpool)
                sig_signal = ypool_true.std() if ypool_true.std() > 1e-9 else 1.0
                for sigma in SIG_GRID:
                    for N in N_GRID:
                        cell_id += 1
                        Xtr = sample_X(N, seed=2000 + cell_id)
                        ytr = true_f(Xtr) + RNG.normal(0, sigma * sig_signal, size=N)

                        # black-box ceiling (KernelRidge RBF, alpha tuned, N capped)
                        r2_bb = krr_cv_r2(Xtr, ytr, n_splits=4)

                        # law compressibility (t3: DEPTH + nested CV)
                        r2_law, depth_sel, terms_sel = lambda_depth_nestedcv(
                            Xtr, ytr, depth_grid=(1, 2, 3), max_terms=4,
                            max_vars=d_screen, outer_splits=4, inner_splits=3)

                        denom = max(r2_bb, 0.05)
                        Lambda = float(np.clip(max(r2_law, 0.0) / denom, 0.0, 1.0))

                        # final law on ALL training data (at selected depth) for selection
                        keep_idx = screen_keep(Xtr, ytr, d_screen)
                        Xtr_k = Xtr[:, keep_idx]
                        Xtr_sk, muk, sdk = standardize_cols(Xtr_k)
                        feats_final = build_feature_bank(Xtr_sk, int(round(depth_sel)) or 1, d_screen)
                        law_final = greedy_law(Xtr_sk, ytr, feats_final, 4)

                        def law_predict(Xq, keep_idx=keep_idx, muk=muk, sdk=sdk,
                                        law=law_final, ymean=ytr.mean()):
                            Xq = np.asarray(Xq, float)[:, keep_idx]
                            if law is None:
                                return np.full(Xq.shape[0], ymean)
                            return predict_law(law, (Xq - muk) / sdk)

                        sel_law = selection_fitness(law_predict, Xpool, ypool_true)

                        # black-box selection (fit once on the full training set, capped)
                        bb_pred_pool, _ = krr_fit_predict(Xtr, ytr, Xpool, alpha=None, seed=0)
                        sel_bb = dict(
                            top1_regret=float(max(0.0, (ypool_true.max() - ypool_true[int(np.argmax(bb_pred_pool))]) /
                                                  (ypool_true.max() if abs(ypool_true.max()) > 1e-9 else 1.0))),
                            spearman=float(stats.spearmanr(bb_pred_pool, ypool_true).correlation))
                        order_bb = np.argsort(bb_pred_pool)[::-1][:5]
                        sel_bb["top5_regret"] = float(max(0.0, (ypool_true.max() - ypool_true[order_bb].max()) /
                                                          (ypool_true.max() if abs(ypool_true.max()) > 1e-9 else 1.0)))

                        rows.append(dict(
                            cell=cell_id, d=d, k=k, sigma=sigma, N=N, encoding=enc,
                            cluster=f"d{d}_k{k}_s{sigma}",          # (d,k,sigma) cluster id
                            r2_bb=r2_bb, r2_law=r2_law, Lambda=Lambda,
                            depth_sel=depth_sel, terms_sel=terms_sel,
                            law_top1_regret=sel_law["top1_regret"],
                            law_top5_regret=sel_law["top5_regret"],
                            law_spearman=sel_law["spearman"],
                            bb_top1_regret=sel_bb["top1_regret"],
                            bb_top5_regret=sel_bb["top5_regret"],
                            bb_spearman=sel_bb["spearman"]))
        print(f"[sweep] finished d={d}  cells={len(rows)}  dt={time.time()-td:.1f}s  total={time.time()-t0:.1f}s",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, "sweep_grid.csv"), index=False)

    sel_fit = 1.0 - df["law_top1_regret"].values
    rho_primary = float(stats.spearmanr(df["Lambda"].values, sel_fit).correlation)
    print(f"\nDONE in {time.time()-t0:.1f}s  n_cells={len(df)}  rho(Lambda, 1-top1_regret)={rho_primary:.3f}")


if __name__ == "__main__":
    main()
