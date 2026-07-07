#!/usr/bin/env python3
"""
van Lent kinetic-simulator selection driver -- reconstructs round-07's
selection payoff (Fig 6A/6B): does a bounded-complexity sparse-basis SR law
(fit on log(flux), the true multiplicative kinetic scale) SELECT better
designs than a tuned RandomForest ceiling, over the full 279,936-design
ground truth?

Ports round-04's audited LOG-FLUX OMP engine (cdl.sim_engine, verbatim port
of research-loop round-04-t1/work/analysis.py) and drives it over a grid of
training sizes N x seeds, exactly as round-07's selection-regime sweep did.
Fitting the law on RAW flux (instead of log(flux)) does NOT reproduce this
result -- the true law flux ~ EA*EC/EG is linear in log-space, so the OMP
sparse-basis fit must operate on the log target.

Run:  /usr/bin/python3 scripts/run_simulator.py --seeds 20
"""
from __future__ import annotations
import os, json, warnings, argparse, time
import numpy as np
import pandas as pd

from cdl.simulator import load_simulator, DEFAULT_PATH
from cdl.sim_engine import (
    r2_det, fit_omp_pareto, omp_predict_log, fit_rf_tuned, selection_metrics,
    MAX_TERMS,
)

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)

NS = [50, 100, 200, 500, 1000]
N_TEST = 4000          # disjoint held-out test rows for r2_law / r2_bb


def eval_draw(N, seed, X, flux, logflux, true_best, true_top100_set):
    """One (N, seed) training-subset draw: fit law (log-flux OMP) + bb (tuned
    RF on raw flux), score held-out R2_det, and rank the FULL 279,936-design
    space for selection (top-100 overlap, top-1 regret)."""
    n_total = X.shape[0]
    rng = np.random.RandomState(seed)
    tr_idx = rng.choice(n_total, N, replace=False)
    mask = np.ones(n_total, dtype=bool)
    mask[tr_idx] = False
    remaining = np.flatnonzero(mask)
    te_idx = rng.choice(remaining, size=min(N_TEST, remaining.size), replace=False)

    Xtr, Xte = X[tr_idx], X[te_idx]
    ytr_log = logflux[tr_idx]
    ytr_flux = flux[tr_idx]
    yte_flux = flux[te_idx]

    # ---- law: sparse-basis OMP fit on log(flux) ----
    models, maxL = fit_omp_pareto(Xtr, ytr_log, max_terms=MAX_TERMS)
    law_model = models[maxL]
    yp_log_te = omp_predict_log(law_model, Xte)
    yp_flux_te = np.exp(np.clip(yp_log_te, -50, 50))
    r2_law = r2_det(yte_flux, yp_flux_te)

    # ---- bb: tuned RandomForest fit directly on raw flux ----
    rf, _ = fit_rf_tuned(Xtr, ytr_flux)
    yte_bb_flux = rf.predict(Xte)
    r2_bb = r2_det(yte_flux, yte_bb_flux)

    # ---- selection over the FULL ground-truth space (exact) ----
    law_pred_full = omp_predict_log(law_model, X)     # log scale; monotone -> ranking valid
    bb_pred_full = rf.predict(X)                        # raw flux scale

    law_sel = selection_metrics(law_pred_full, flux, true_best, true_top100_set,
                                tie_rng=np.random.default_rng(seed))
    bb_sel = selection_metrics(bb_pred_full, flux, true_best, true_top100_set,
                               tie_rng=np.random.default_rng(seed))

    return dict(
        N=int(N), seed=int(seed),
        law_top1=float(law_sel["pick_top1_regret"]),
        bb_top1=float(bb_sel["pick_top1_regret"]),
        law_ov=int(law_sel["top100_overlap"]),
        bb_ov=int(bb_sel["top100_overlap"]),
        r2_law=float(r2_law), r2_bb=float(r2_bb),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--data", default=DEFAULT_PATH)
    args = ap.parse_args()

    X, flux, names = load_simulator(args.data)
    n_total = X.shape[0]
    true_best = float(flux.max())
    true_top100_set = set(np.argpartition(-flux, 100)[:100].tolist())
    logflux = np.log(flux)
    print(f"[data] {n_total} designs; true_best flux={true_best:.4f}", flush=True)

    sel_rows, ncurve_rows = [], []
    for N in NS:
        t0 = time.time()
        for seed in range(args.seeds):
            r = eval_draw(N, seed, X, flux, logflux, true_best, true_top100_set)
            sel_rows.append(dict(N=r["N"], seed=r["seed"], law_top1=r["law_top1"],
                                 bb_top1=r["bb_top1"], law_ov=r["law_ov"], bb_ov=r["bb_ov"]))
            ncurve_rows.append(dict(N=r["N"], seed=r["seed"], r2_law=r["r2_law"], r2_bb=r["r2_bb"]))
        print(f"[N={N:4d}] done ({args.seeds} seeds) t={time.time()-t0:.0f}s", flush=True)

    sel_df = pd.DataFrame(sel_rows)
    sel_df.to_csv(os.path.join(RESULTS, "selection_regime.csv"), index=False)

    nc_df = pd.DataFrame(ncurve_rows)
    ncurve_summary = nc_df.groupby("N", as_index=False)[["r2_law", "r2_bb"]].mean()
    ncurve_summary.to_csv(os.path.join(RESULTS, "sim_Ncurve.csv"), index=False)

    per_N = {}
    for N, g in sel_df.groupby("N"):
        per_N[str(int(N))] = dict(
            law_top1=round(float(g["law_top1"].mean()), 3),
            bb_top1=round(float(g["bb_top1"].mean()), 3),
            law_top100=round(float(g["law_ov"].mean()), 1),
            bb_top100=round(float(g["bb_ov"].mean()), 1),
        )
    summary = dict(per_N=per_N, n_seeds=args.seeds,
                   law_wins_region_all_N=all(per_N[k]["law_top100"] > per_N[k]["bb_top100"] for k in per_N))
    with open(os.path.join(RESULTS, "selection_regime_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print("\n== per-N selection summary ==")
    print(json.dumps(summary, indent=2))
    print("[done] wrote results/selection_regime.csv, results/selection_regime_summary.json, "
          "results/sim_Ncurve.csv")


if __name__ == "__main__":
    main()
