#!/usr/bin/env python3
"""
van Lent kinetic-simulator driver -- reconstructs Fig 6A (compressibility vs
training size N) and Fig 6B (selection payoff: does the interpretable law
out-select a tuned black-box ceiling?) from round-04/round-07.

Ports round-04's audited LOG-FLUX OMP engine (cdl.sim_engine, a near-verbatim
port of research-loop round-04-t1/work/analysis.py) and drives it two ways:

  Fig 6A (sim_Ncurve.csv): a SINGLE deterministic trajectory across
    N in {50,100,200,500,1000} using one fixed-seed rng that advances across
    the N loop -- exactly round-04 main()'s N-curve construction. One
    eval_N() call per N; DETERMINISTIC (reproduces the golden closely).

  Fig 6B (selection_regime.csv): an N x seed grid (seed in range(--seeds)).
    Each (N, seed) draws a fresh N-design training set (seed unique to that
    cell), fits the law + black-box ceiling, and ranks the FULL 279,936-design
    space to record the single-pick top-1 regret and top-100 overlap for both
    models. round-07's exact driver was never persisted, so this is a
    reconstruction -- exact per-cell values are not expected to match, but the
    aggregate payoff direction (law_ov >> bb_ov, especially at N>=100) is the
    paper's Fig 6B claim.

The law is fit to log(flux) (the multiplicative kinetic scale); predictions
back-transform via exp(clip(...)); r2_det is always scored on the flux scale.

Run:  /usr/bin/python3 scripts/run_simulator.py --seeds 20
"""
from __future__ import annotations
import os
import time
import warnings
import argparse

import numpy as np
import pandas as pd
from numpy.random import default_rng

from cdl.simulator import load_simulator
from cdl.sim_engine import eval_N, SEED, NS

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")

NCURVE_COLS = [
    "N", "r2_law", "r2_blackbox", "lambda_ratio",
    "top1_regret_law", "top1_regret_law_pick", "best_in_plateau_regret_law",
    "n_tied_best_law", "top100_overlap_law",
    "top1_regret_blackbox", "top1_regret_blackbox_pick", "n_tied_best_blackbox",
    "top100_overlap_blackbox", "Lbudget", "bb_rank_universe",
]
SEL_COLS = ["N", "seed", "law_top1", "bb_top1", "law_ov", "bb_ov"]


def load_data(data_dir):
    path = os.path.join(data_dir, "vanlent-simulator.csv")
    X, flux, names = load_simulator(path)
    logflux = np.log(flux)
    n_total = X.shape[0]
    true_best = float(flux.max())
    true_top100_set = set(np.argpartition(-flux, 100)[:100].tolist())
    return X, flux, logflux, true_best, true_top100_set, n_total


def run_ncurve(X, flux, logflux, true_best, true_top100_set):
    """Fig 6A: single deterministic trajectory across N (fixed SEED). One rng
    is created once and advances across the NS loop -- this matches round-04
    main()'s N-curve construction exactly (NOT a fresh rng per N)."""
    rng = default_rng(SEED)
    rows = []
    for N in NS:
        res = eval_N(N, X, flux, logflux, true_best, true_top100_set, rng)
        rows.append({k: res[k] for k in NCURVE_COLS})
        print(f"[Ncurve N={N:4d}] r2_law={res['r2_law']:.3f} r2_bb={res['r2_blackbox']:.3f} "
              f"lam={res['lambda_ratio']:.3f} top100_ov_law={res['top100_overlap_law']} "
              f"top100_ov_bb={res['top100_overlap_blackbox']} Lbudget={res['Lbudget']}",
              flush=True)
    return pd.DataFrame(rows, columns=NCURVE_COLS)


def run_selection_regime(X, flux, logflux, true_best, true_top100_set, n_seeds):
    """Fig 6B: N x seed grid. Each (N, seed) draws N training designs with a
    seed unique to that cell (independent of loop order), fits law + bb via
    eval_N, and records the selection payoff: single-pick top-1 regret and
    top-100 overlap over the full ground-truth space, for both models."""
    rows = []
    for N in NS:
        t0 = time.time()
        for seed in range(n_seeds):
            rng = default_rng((int(seed), int(N)))
            res = eval_N(N, X, flux, logflux, true_best, true_top100_set, rng)
            rows.append(dict(
                N=int(N), seed=int(seed),
                law_top1=float(res["top1_regret_law_pick"]),
                bb_top1=float(res["top1_regret_blackbox_pick"]),
                law_ov=int(res["top100_overlap_law"]),
                bb_ov=int(res["top100_overlap_blackbox"]),
            ))
        print(f"[selection N={N:4d}] done ({n_seeds} seeds) t={time.time() - t0:.0f}s",
              flush=True)
    return pd.DataFrame(rows, columns=SEL_COLS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data", "cleaned"),
                     help="directory containing vanlent-simulator.csv")
    ap.add_argument("--seeds", type=int, default=20,
                     help="seeds per N for the Fig 6B selection-payoff grid")
    args = ap.parse_args()

    os.makedirs(RESULTS, exist_ok=True)
    X, flux, logflux, true_best, true_top100_set, n_total = load_data(args.data_dir)
    print(f"[data] {n_total} designs (6^7={6 ** 7}); true_best flux={true_best:.4f}",
          flush=True)

    ncurve = run_ncurve(X, flux, logflux, true_best, true_top100_set)
    ncurve.to_csv(os.path.join(RESULTS, "sim_Ncurve.csv"), index=False)
    print("[done] wrote results/sim_Ncurve.csv", flush=True)

    sel = run_selection_regime(X, flux, logflux, true_best, true_top100_set, args.seeds)
    sel.to_csv(os.path.join(RESULTS, "selection_regime.csv"), index=False)
    print("[done] wrote results/selection_regime.csv", flush=True)

    print("\n== selection payoff medians (law_ov vs bb_ov) by N ==")
    print(sel.groupby("N")[["law_ov", "bb_ov"]].median())


if __name__ == "__main__":
    main()
