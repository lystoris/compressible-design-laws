#!/usr/bin/env python3
"""Assemble tidy CSVs for the paper figures from THIS REPO's reproduced results/.

Ported from paper/v2/figures/data/prep_data.py. The original read from a deep
research-loop round-08 work directory (`RL`/`R8`); this version reads only from this
repo's own `results/` (everything `scripts/run_*.py` reproduces/reconstructs -- see
docs/reconciliation-notes.md for how each results/ file was produced and any deltas
from the original round runs).

Run: /usr/bin/python3 scripts/build_figure_data.py
Writes into figures/data/.

------------------------------------------------------------------------------------
PROVENANCE -- every figures/data/*.csv is either:

(A) RECOMPUTED HERE, each run, from results/*.csv|json (this repo's own reproduced or
    controller-validated reconstructed outputs):
      fig3_population.csv    <- results/panel_t3.csv + results/proteingym_t3.csv
      fig3_partials.csv      <- results/summary_t3.json + results/proteingym_summary_t3.json
      fig3_pseudorep.csv     <- results/proteingym_t3.csv (Tsuboyama pseudoreplication groups)
      fig2_effdim_pooled.csv <- results/effdim_grid_t1.csv
      fig2_eta2.csv          <- results/effdim_grid_t1.csv (per-family eta2 via cdl.stats.eta2)
      fig5_selection.csv     <- results/selection_regime_summary.json
      fig5_enzyme_sens.csv   <- cdl.simulator.enzyme_sensitivity on data/cleaned/vanlent-simulator.csv

    NOTE on trajectory coverage: this repo reproduces trajectory t3 for the curated/ProteinGym
    panels and t1 for the decoupling grid (see docs/reconciliation-notes.md) -- it does NOT
    have t1/t2/t3 versions of every intermediate the way the original round-08 assembly did.
    fig3_partials.csv therefore carries one (stream, axis) point per available trajectory
    (t3 only) instead of the original's three; it is NOT padded with fabricated t1/t2 rows.
    Likewise fig2_effdim_pooled.csv / fig2_eta2.csv carry t1 only (this repo's reconstructed
    decoupling grid; see docs/reconciliation-notes.md R-3), not a t1/t2/t3 pool.

(B) a RECORDED AUDITED ARTIFACT: committed once, verbatim, in figures/data/, with its
    provenance documented below -- read (never re-typed as Python constants) by this
    script's callers (the R figure scripts) directly. This script does not write these
    files; it only checks they are present, so a missing one fails loudly instead of
    silently producing a blank/wrong panel.

      fig4_poelwijk_probe.csv -- round-08 Poelwijk link-transform probe: raw R2_law 0.590,
                                  log 0.606, rank-normal 0.531, interpretable 0.318 (r2_bb
                                  alongside). THE FIX: the original prep_data.py hand-typed
                                  these four numbers inline as `dict(transform=..., r2_law=...)`
                                  literals despite its own module docstring claiming "no
                                  hand-typed stats" -- the honesty gap this task closes. They
                                  are now a single committed CSV with this comment as their
                                  provenance record (round-08 Poelwijk link-probe), regenerated
                                  by a future `scripts/run_panel.py --only poelwijk --probe`
                                  mode if one is added; until then this is the recorded value.
      anchor_table.csv         -- round-03-t3 audited anchor table (Fig 2a: absolute R2det vs
                                  ratio-to-black-box). Round-03's generating driver was not
                                  persisted for this repo (only the per-trajectory result CSVs
                                  survive -- same pattern as R-3 in docs/reconciliation-notes.md)
                                  -- copied in verbatim from the audited round-03-t3
                                  work/results/anchor_table.csv.
      betacarotene_oof.csv     -- out-of-fold beta-carotene law fit (Fig 5b: titre ~ C + C/B +
                                  A*C). Copied verbatim from paper/v2/figures/data/; the fitting
                                  notebook is not part of this repo's scripts/.
      figS_calibration.csv     -- RF participation-ratio estimator calibration sweep (Fig S2).
                                  Copied verbatim from paper/v2/figures/data/.
                                  scripts/build_calibration_data.py (ported from
                                  prep_calibration.py) can regenerate an equivalent sweep from
                                  scratch, but its cells are reps of random synthetic
                                  landscapes and won't be bit-identical to the committed one
                                  actually plotted in the paper, so the recorded CSV is kept
                                  (see that script's docstring).

These four files are committed directly in figures/data/ and are NOT written by this
script; main() only checks they are present.
------------------------------------------------------------------------------------
"""
from __future__ import annotations
import os
import json

import pandas as pd
from scipy import stats

from cdl.stats import eta2
from cdl.simulator import load_simulator, enzyme_sensitivity

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
FIG_DATA = os.path.join(ROOT, "figures", "data")
CLEANED = os.path.join(ROOT, "data", "cleaned")

PROTEIN_CUR = {"poelwijk", "gb1", "avgfp", "aav"}

RECORDED_ARTIFACTS = [
    "fig4_poelwijk_probe.csv",
    "anchor_table.csv",
    "betacarotene_oof.csv",
    "figS_calibration.csv",
]


def domain(idv, src):
    if src == "proteingym":
        return "protein"
    if idv == "vanlent-simulator":
        return "simulator"
    return "protein" if idv in PROTEIN_CUR else "metabolic"


def build_population():
    """Fig 4a/b/c/d (population): curated panel_t3 + ProteinGym proteingym_t3, joined + labelled."""
    cur = pd.read_csv(os.path.join(RESULTS, "panel_t3.csv"))
    cur["source"] = "curated"
    pg = pd.read_csv(os.path.join(RESULTS, "proteingym_t3.csv"))
    pg["source"] = "proteingym"
    pg["group"] = "protein-dms"
    pop = pd.concat([cur, pg], ignore_index=True)
    pop["domain"] = [domain(i, s) for i, s in zip(pop["id"], pop["source"])]
    pop["compresses"] = (pop["r2_law"] >= 0.7).map({True: "full", False: "partial-or-none"})
    out = os.path.join(FIG_DATA, "fig3_population.csv")
    pop.to_csv(out, index=False)
    return pop, out


def build_partials():
    """Fig 4c (partial-correlation forest): curated + ProteinGym, trajectory t3 only
    (this repo reproduces one trajectory per stream; see module docstring)."""
    rows = []
    for stream, fname in [
        ("curated", "summary_t3.json"),
        ("ProteinGym", "proteingym_summary_t3.json"),
    ]:
        s = json.load(open(os.path.join(RESULTS, fname)))
        rows.append(dict(stream=stream, traj=s["traj"], axis="effective d", partial=s["partial_effd"]))
        rows.append(dict(stream=stream, traj=s["traj"], axis="nominal d", partial=s["partial_nomd"]))
    out = os.path.join(FIG_DATA, "fig3_partials.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    return out


def build_pseudorep():
    """Fig 4d: ProteinGym pseudoreplication-corrected rho(compress, effective d), recomputed
    from results/proteingym_t3.csv (Tsuboyama-2023 studies vs. 15 other independent studies)."""
    d = pd.read_csv(os.path.join(RESULTS, "proteingym_t3.csv"))
    d["is_tsu"] = d["id"].str.contains("Tsuboyama")
    d["study"] = d["id"].str.extract(r"_([A-Za-z0-9]+_\d{4})")

    def rho(sub):
        return stats.spearmanr(sub.r2_law, sub.eff_d).correlation

    g = d.groupby(d["study"].where(~d.is_tsu, "Tsuboyama_2023")).agg(
        r2_law=("r2_law", "median"), eff_d=("eff_d", "median")
    )
    prow = [
        dict(level="all 69 (pooled)", n=len(d), rho=rho(d)),
        dict(level="Tsuboyama within-protocol", n=int(d.is_tsu.sum()), rho=rho(d[d.is_tsu])),
        dict(level="15 other studies", n=int((~d.is_tsu).sum()), rho=rho(d[~d.is_tsu])),
        dict(level="study-collapsed", n=len(g), rho=stats.spearmanr(g.r2_law, g.eff_d).correlation),
    ]
    out = os.path.join(FIG_DATA, "fig3_pseudorep.csv")
    pd.DataFrame(prow).to_csv(out, index=False)
    return out


def build_effdim_pooled():
    """Fig 3c (decoupling boxplots): results/effdim_grid_t1.csv, this repo's reconstructed
    trajectory (see docs/reconciliation-notes.md R-3). Only t1 is available."""
    e1 = pd.read_csv(os.path.join(RESULTS, "effdim_grid_t1.csv"))
    e1 = e1.copy()
    e1["traj"] = "t1"
    out = os.path.join(FIG_DATA, "fig2_effdim_pooled.csv")
    e1.to_csv(out, index=False)
    return e1, out


def build_eta2(effdim_pooled):
    """Fig 3d (eta2 effective vs nominal per family): recomputed from the same t1 grid via
    cdl.stats.eta2, per-family (paper convention -- see tests/test_decoupling_golden.py)."""
    rows = []
    for fam, sub in effdim_pooled.groupby("family"):
        rows.append(dict(
            traj="t1", family=fam,
            eta2_eff=eta2(sub.r2_law, sub.d_eff),
            eta2_nom=eta2(sub.r2_law, sub.d_nom),
        ))
    out = os.path.join(FIG_DATA, "fig2_eta2.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    return out


def build_selection():
    """Fig 6b (selection payoff): law vs black-box top-100/top-1 recovery vs training size N."""
    sj = json.load(open(os.path.join(RESULTS, "selection_regime_summary.json")))
    srows = [
        dict(N=int(N), law_top100=v["law_top100"], bb_top100=v["bb_top100"],
             law_top1=v["law_top1"], bb_top1=v["bb_top1"])
        for N, v in sj["per_N"].items()
    ]
    out = os.path.join(FIG_DATA, "fig5_selection.csv")
    pd.DataFrame(srows).sort_values("N").to_csv(out, index=False)
    return out


def build_enzyme_sensitivity():
    """Fig 6c (enzyme sensitivity): |Spearman| of each simulator enzyme knob with flux, via
    cdl.simulator.enzyme_sensitivity over the full 279,936-row enumerated table (deterministic;
    the original prep_data.py subsampled 5000 rows with a fixed seed -- the full table is small
    enough to use directly here, so no sampling is needed)."""
    sim_path = os.path.join(CLEANED, "vanlent-simulator.csv")
    X, flux, _names = load_simulator(sim_path)
    sens = enzyme_sensitivity(X, flux)
    sens_df = pd.DataFrame(
        [dict(enzyme=name, abs_spearman=val) for name, val in sens.items()]
    ).sort_values("abs_spearman", ascending=False)
    out = os.path.join(FIG_DATA, "fig5_enzyme_sens.csv")
    sens_df.to_csv(out, index=False)
    return out


def verify_recorded_artifacts():
    """The four RECORDED AUDITED ARTIFACTS (see module docstring) must already be committed
    in figures/data/ -- this script never writes or fabricates them."""
    missing = [f for f in RECORDED_ARTIFACTS if not os.path.exists(os.path.join(FIG_DATA, f))]
    if missing:
        raise FileNotFoundError(
            "Missing recorded audited artifact(s) in figures/data/: " + ", ".join(missing) +
            ". These are committed source artifacts (see the provenance block at the top of "
            "this script) -- restore them from their documented provenance rather than "
            "fabricating replacement values."
        )
    probe = pd.read_csv(os.path.join(FIG_DATA, "fig4_poelwijk_probe.csv"))
    need = {"transform", "r2_law", "r2_bb"}
    assert need <= set(probe.columns) and len(probe) == 4, \
        "fig4_poelwijk_probe.csv is missing/malformed (expected columns %s, 4 rows)" % need
    return RECORDED_ARTIFACTS


def main():
    os.makedirs(FIG_DATA, exist_ok=True)
    pop, pop_path = build_population()
    partials_path = build_partials()
    pseudorep_path = build_pseudorep()
    effdim_pooled, effdim_path = build_effdim_pooled()
    eta2_path = build_eta2(effdim_pooled)
    sel_path = build_selection()
    sens_path = build_enzyme_sensitivity()
    recorded = verify_recorded_artifacts()

    print("WROTE (recomputed from results/ each run):")
    for p in (pop_path, partials_path, pseudorep_path, effdim_path, eta2_path, sel_path, sens_path):
        print(" ", os.path.relpath(p, ROOT))
    print("\nPRESENT (recorded audited artifacts, not written by this script):")
    for f in recorded:
        print(" ", os.path.join("figures", "data", f))

    print("\nfig3_population (n):", len(pop), "| domains:", pop["domain"].value_counts().to_dict())


if __name__ == "__main__":
    main()
