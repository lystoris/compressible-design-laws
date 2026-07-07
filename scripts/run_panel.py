#!/usr/bin/env python3
"""
Curated-panel driver — reads a cleaned-data manifest, encodes each dataset to a numeric
design matrix (via cdl.encoding.load), scores compressibility with the audited SR anchor
(cdl.compressibility.run_anchor) plus an independent RF participation-ratio effective-d
(cdl.effective_dim.effective_d), and fits the population law (compressibility vs eff-d |
nom-d) via cdl.stats.partial_spearman.

Ported from panel.py:159-205 (research-loop/sr-compressible-design-laws/rounds/round-08).

Run:  /usr/bin/python3 scripts/run_panel.py --smoke                      # keystone + a few
      /usr/bin/python3 scripts/run_panel.py --traj t3 --data-dir data/cleaned  # full panel
"""
from __future__ import annotations
import os, json, warnings, argparse, time
import numpy as np, pandas as pd
from scipy import stats

from cdl.encoding import load
from cdl.compressibility import run_anchor
from cdl.effective_dim import effective_d
from cdl.stats import partial_spearman

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="t3")
    ap.add_argument("--smoke", action="store_true")
    # 2000 = the subsample cap used for the published panel_t3 golden / reported partials (paper Methods says 1,500 — stale; see reconciliation note).
    ap.add_argument("--max-n", type=int, default=2000)
    ap.add_argument("--only", default="")
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data", "cleaned"))
    args = ap.parse_args()
    seed = {"t1": 1, "t2": 2, "t3": 3}.get(args.traj, 0)
    data_dir = args.data_dir

    man = pd.read_csv(os.path.join(data_dir, "manifest.csv"))
    man = man[man.status.isin(["clean", "needs-review"])]
    if args.only:
        ids = args.only.split(",")
    elif args.smoke:
        ids = ["poelwijk", "gb1", "avgfp", "pca-paperA", "beta-carotene"]
    else:
        ids = man.id.tolist()

    rows = []
    for did in ids:
        r = man[man.id == did]
        if r.empty:
            print(f"[skip] {did} not in manifest")
            continue
        enc = r.iloc[0]["encoding"]
        grp = r.iloc[0]["group"]
        t0 = time.time()
        try:
            X, y, feat, groups, nomd = load(os.path.join(data_dir, f"{did}.csv"), enc,
                                             max_n=args.max_n, seed=seed)
            effd = effective_d(X, y, groups, seed=seed)
            # nonlinear (reciprocal/ratio) terms only make sense for CONTINUOUS design vars;
            # on binary/one-hot features 1/x and x/x' are degenerate. Gate on feature TYPE, uniformly.
            frac_binary = float(np.mean([(np.unique(X[:, j]).size <= 2) for j in range(X.shape[1])]))
            allow_nl = frac_binary < 0.5
            r2l, r2b = run_anchor(X, y, feat, allow_nonlinear=allow_nl, topk=40, fast=args.smoke)
            rows.append(dict(id=did, group=grp, nom_d=nomd, eff_d=round(effd, 2),
                              r2_law=round(r2l, 3), r2_bb=round(r2b, 3), N=len(y)))
            print(f"  {did:18s} grp={grp:9s} nom_d={nomd:4d} eff_d={effd:6.2f} "
                  f"r2_law={r2l:+.3f} r2_bb={r2b:+.3f} ({time.time()-t0:.0f}s)", flush=True)
        except Exception as e:
            print(f"  [err] {did}: {e}", flush=True)

    if not rows:
        return
    df = pd.DataFrame(rows)
    tag = "smoke" if args.smoke else args.traj
    df.to_csv(os.path.join(RESULTS, f"panel_{tag}.csv"), index=False)

    if len(df) >= 4:
        c = df.r2_law.values
        e = df.eff_d.values
        n = df.nom_d.values.astype(float)
        summ = dict(traj=tag, n=len(df),
                    partial_effd=partial_spearman(c, e, n),
                    partial_nomd=partial_spearman(c, n, e),
                    rho_effd=float(stats.spearmanr(c, e).correlation),
                    rho_nomd=float(stats.spearmanr(c, n).correlation))
        summ["effd_beats_nomd"] = abs(summ["partial_effd"]) > abs(summ["partial_nomd"]) and summ["partial_effd"] < 0
        json.dump(summ, open(os.path.join(RESULTS, f"summary_{tag}.json"), "w"), indent=2)
        print("\n== population ==")
        print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
