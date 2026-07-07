#!/usr/bin/env python3
"""Generate the paper tables from THIS REPO's reproduced results/.

Ported from paper/v2/tables/make_tables.py (which read a research-loop round-08 dir).

NON-DESTRUCTIVE BY DESIGN (same rationale as scripts/build_figure_data.py): Table 1 and
Table S1/S2 report the median [min, max] over the THREE analysis trajectories (t1/t2/t3),
which a default reproduction run does not produce (this repo reproduces trajectory t3).
The committed tables/*.csv|md are the PAPER-CANONICAL tables. This script regenerates them
only when panel_t{1,2,3}.csv AND proteingym_t{1,2,3}.csv are ALL present in results/;
otherwise it PRESERVES the committed canonical tables and prints how to regenerate.

To fully regenerate: run `run_panel --traj t1|t2|t3` and `run_proteingym --traj t1|t2|t3`,
then re-run this script.

Note: the Table S2 columns are one audited engine re-run under three random SEEDS (t1/t2/t3 =
seeds 1/2/3), i.e. stability to resampling -- NOT three different SR engines. (A stale legend
in an earlier manuscript draft mislabelled them as "sparse-OMP / enumerative / depth" engines;
that legend belongs to the synthetic sweep/decoupling, not the curated panel -- flagged for the
manuscript correction, MANIFEST F-A.)
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
MAN = os.path.join(ROOT, "data", "cleaned", "manifest.csv")
TRAJ = ["t1", "t2", "t3"]


def med_range(vals):
    v = np.array(vals, float)
    return f"{np.median(v):.2f} [{v.min():.2f}, {v.max():.2f}]"


SYS = {
    "beta-carotene": ("Y. lipolytica / β-carotene", "3 integer gene copies"),
    "tryptophan": ("S. cerevisiae / tryptophan", "5 genes × 6 promoters"),
    "limonene": ("E. coli / limonene", "9 pathway protein levels"),
    "isoprenol": ("P. putida / isoprenol", "binary multiplex-CRISPRi (≤67)"),
    "jervis-RBS": ("E. coli / limonene", "combinatorial RBS (one-hot)"),
    "vanlent-pCA": ("S. cerevisiae / p-coumaric acid", "19 genes × promoter strength"),
    "vanlent-simulator": ("in silico kinetic / flux", "7 enzyme levels (full 6^7)"),
    "poelwijk": ("eqFP611 / brightness", "13 binary positions (complete 2^13)"),
    "gb1": ("protein G B1 / fitness", "4 sites × 20 aa"),
    "avgfp": ("avGFP / brightness", "~233 positions"),
    "aav": ("AAV2 capsid / viability", "28-position segment"),
    "smanski-nif": ("E. coli (nif) / nitrogenase", "refactored cluster part-slots"),
    "pca-paperA": ("S. cerevisiae / p-coumaric acid", "6 categorical pathway factors"),
}


def _have_all_trajectories():
    need = [f"panel_{t}.csv" for t in TRAJ] + [f"proteingym_{t}.csv" for t in TRAJ]
    return all(os.path.exists(os.path.join(RESULTS, n)) for n in need)


def to_md(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def main():
    if not _have_all_trajectories():
        print("[preserve] tables/*.csv|md are canonical (median over t1/t2/t3).")
        print("           results/ has only a subset of trajectories, so the committed")
        print("           paper-canonical tables are kept unchanged.")
        print("           To regenerate: run run_panel + run_proteingym for --traj t1, t2, t3,")
        print("           then re-run this script.")
        return

    man = pd.read_csv(MAN).set_index("id")
    cur = {t: pd.read_csv(os.path.join(RESULTS, f"panel_{t}.csv")).set_index("id") for t in TRAJ}
    ids = cur["t3"].index.tolist()
    rows = []
    for i in ids:
        eff = [cur[t].loc[i, "eff_d"] for t in TRAJ]
        r2 = [cur[t].loc[i, "r2_law"] for t in TRAJ]
        bb = [cur[t].loc[i, "r2_bb"] for t in TRAJ]
        sysp, enc = SYS.get(i, (i, ""))
        comp = "yes" if np.median(r2) >= 0.7 else ("partial" if np.median(r2) >= 0.4 else "no")
        rows.append(dict(Dataset=i, System=sysp, Encoding=enc,
                         **{"Nominal d": int(cur["t3"].loc[i, "nom_d"]),
                            "Effective d": med_range(eff), "N": int(cur["t3"].loc[i, "N"]),
                            "Compressibility R2det": med_range(r2), "Black-box ceiling": med_range(bb),
                            "Compresses": comp,
                            "Domain": man.loc[i, "group"] if i in man.index else ""}))
    t1 = pd.DataFrame(rows).sort_values("Compressibility R2det", ascending=False)
    t1.to_csv(os.path.join(HERE, "Table1_panel.csv"), index=False)

    s2 = []
    for i in ids:
        d = dict(Dataset=i)
        for t in TRAJ:
            d[f"R2det_{t}"] = round(cur[t].loc[i, "r2_law"], 3)
            d[f"effd_{t}"] = round(cur[t].loc[i, "eff_d"], 2)
        s2.append(d)
    pd.DataFrame(s2).to_csv(os.path.join(HERE, "TableS2_trajectories.csv"), index=False)

    pg = {t: pd.read_csv(os.path.join(RESULTS, f"proteingym_{t}.csv")).set_index("id") for t in TRAJ}
    pids = pg["t3"].index.tolist()
    prows = []
    for i in pids:
        eff = [pg[t].loc[i, "eff_d"] for t in TRAJ]
        r2 = [pg[t].loc[i, "r2_law"] for t in TRAJ]
        prows.append(dict(Assay=i, **{"Nominal d": int(pg["t3"].loc[i, "nom_d"]),
                                       "Effective d (median)": round(np.median(eff), 2),
                                       "Compressibility R2det (median)": round(np.median(r2), 3),
                                       "Study": "Tsuboyama_2023" if "Tsuboyama" in i
                                       else i.split("_")[1] if "_" in i else ""}))
    pd.DataFrame(prows).sort_values("Compressibility R2det (median)", ascending=False).to_csv(
        os.path.join(HERE, "TableS1_proteingym.csv"), index=False)

    open(os.path.join(HERE, "Table1_panel.md"), "w").write(
        "# Table 1. The compressibility panel\n\n" + to_md(t1) + "\n")

    print(f"WROTE Table1_panel (.csv/.md), TableS1_proteingym.csv, TableS2_trajectories.csv "
          f"({len(t1)} datasets, {len(prows)} assays)")


if __name__ == "__main__":
    main()
