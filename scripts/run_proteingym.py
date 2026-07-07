#!/usr/bin/env python3
"""
ProteinGym population run — the anti-cherry-pick scale result.
Iterates the 69 COMBINATORIAL (multi-mutant) ProteinGym assays, encodes each as per-position
amino-acid one-hot (group=position), and reuses the audited engine (cdl.compressibility.run_anchor)
to score compressibility (r2_law / r2_bb) + INDEPENDENT effective-d (cdl.effective_dim.effective_d,
RF participation ratio). Then fits the population law: partial Spearman(compress, eff-d | nom-d)
vs (compress, nom-d | eff-d).

`load_assay` is its OWN loader (does not use cdl.encoding.load): it reads only
["mutant","DMS_score"], builds per-position amino-acid one-hot (group=position) via
cdl.encoding.parse_aa, drops zero-variance columns, and returns None when <2 groups or <40 rows.

Ported from proteingym_run.py (research-loop/sr-compressible-design-laws/rounds/round-08).

Run:  /usr/bin/python3 scripts/run_proteingym.py --traj t3 --data-dir data/proteingym
"""
from __future__ import annotations
import os, json, warnings, argparse, time
import numpy as np, pandas as pd
from scipy import stats

from cdl.encoding import parse_aa
from cdl.compressibility import run_anchor
from cdl.effective_dim import effective_d
from cdl.stats import partial_spearman

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)


def load_assay(did, data_dir, max_n, seed):
    pg = os.path.join(data_dir, "DMS_ProteinGym_substitutions")
    df = pd.read_csv(os.path.join(pg, did + ".csv"), usecols=["mutant", "DMS_score"]).dropna(subset=["DMS_score"])
    if len(df) > max_n:
        df = df.sample(max_n, random_state=seed)
    y = df["DMS_score"].values.astype(float)
    aas = df["mutant"].apply(parse_aa)
    pos = sorted({p for d in aas for p in d})
    cols = []
    names = []
    groups = []
    for p in pos:
        letters = sorted({d.get(p, "_") for d in aas})
        if len(letters) < 2:
            continue
        for L in letters[1:]:
            cols.append([1.0 if d.get(p, "_") == L else 0.0 for d in aas])
            names.append(f"p{p}_{L}")
            groups.append(f"p{p}")
    if not cols:
        return None
    X = np.array(cols, float).T
    keep = [i for i in range(X.shape[1]) if np.std(X[:, i]) > 1e-12]
    X = X[:, keep]
    names = [names[i] for i in keep]
    groups = [groups[i] for i in keep]
    if X.shape[1] < 2 or len(y) < 40:
        return None
    return X, y, names, groups, len(set(groups))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="t3")
    ap.add_argument("--max-n", type=int, default=1500)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data", "proteingym"))
    args = ap.parse_args()
    seed = {"t1": 1, "t2": 2, "t3": 3}.get(args.traj, 0)
    data_dir = args.data_dir
    pg = os.path.join(data_dir, "DMS_ProteinGym_substitutions")

    ref = pd.read_csv(os.path.join(data_dir, "proteingym_ref.csv"))
    combo = ref[ref["includes_multiple_mutants"].astype(str).str.lower().isin(["true", "1", "yes"])]
    ids = combo["DMS_id"].tolist()
    if args.limit:
        ids = ids[: args.limit]

    rows = []
    for k, did in enumerate(ids):
        if not os.path.exists(os.path.join(pg, did + ".csv")):
            continue
        t0 = time.time()
        try:
            L = load_assay(did, data_dir, args.max_n, seed)
            if L is None:
                print(f"[skip] {did}")
                continue
            X, y, feat, groups, nomd = L
            effd = effective_d(X, y, groups, seed=seed)
            # proteingym is one-hot/binary encoded => nonlinear (reciprocal/ratio) terms are degenerate
            r2l, r2b = run_anchor(X, y, feat, allow_nonlinear=False, topk=40, fast=False)
            rows.append(dict(id=did, nom_d=nomd, eff_d=round(effd, 2),
                              r2_law=round(r2l, 3), r2_bb=round(r2b, 3), N=len(y)))
            print(f"[{k+1}/{len(ids)}] {did:34s} nom_d={nomd:4d} eff_d={effd:6.2f} "
                  f"r2_law={r2l:+.3f} r2_bb={r2b:+.3f} ({time.time()-t0:.0f}s)", flush=True)
        except Exception as e:
            print(f"[err] {did}: {e}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, f"proteingym_{args.traj}.csv"), index=False)
    if len(df) < 4:
        return

    c = df.r2_law.values
    e = df.eff_d.values
    n = df.nom_d.values.astype(float)
    summ = dict(traj=args.traj, n_assays=len(df),
                partial_effd=partial_spearman(c, e, n),
                partial_nomd=partial_spearman(c, n, e),
                rho_effd=float(stats.spearmanr(c, e).correlation),
                rho_nomd=float(stats.spearmanr(c, n).correlation),
                n_compress_0p7=int((c >= 0.7).sum()))
    summ["effd_beats_nomd"] = abs(summ["partial_effd"]) > abs(summ["partial_nomd"]) and summ["partial_effd"] < 0
    json.dump(summ, open(os.path.join(RESULTS, f"proteingym_summary_{args.traj}.json"), "w"), indent=2)
    print("\n== ProteinGym population (" + args.traj + ") ==")
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
