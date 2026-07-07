#!/usr/bin/env python3
"""Assemble tidy CSVs for the figures from THIS repo's `results/` (+ `data/cleaned/`).

Run: /usr/bin/python3 scripts/build_figure_data.py
Writes into figures/data/. Ported from ../paper/v2/figures/data/prep_data.py, repointed from
the paper's multi-trajectory research-loop roots (RL/R8) to this repo's single committed
trajectory of results (results/*_t3.csv, results/*_t3.json, results/effdim_grid_t1.csv, ...).

WHAT THIS SCRIPT REGENERATES from results/ + data/cleaned/ committed to this repo:
  - figures/data/fig3_population.csv    <- results/panel_t3.csv + results/proteingym_t3.csv
  - figures/data/fig5_selection.csv     <- results/selection_regime_summary.json
  - figures/data/fig5_enzyme_sens.csv   <- data/cleaned/vanlent-simulator.csv

WHAT IS SHIPPED PRE-COMPUTED under figures/data/ and NOT regenerated here, because the
paper's figure needs MULTIPLE trajectories or a round this repo did not reconstruct (this
script refuses to fabricate them from a single trajectory or from thin air):
  - fig3_partials.csv    needs curated summary_t1/t2/t3.json + proteingym_summary_t1/t2/t3.json
                          (run_panel.py --traj t1/t2/t3, run_proteingym.py --traj t1/t2/t3);
                          this repo's results/ only carries the t3 pair.
  - fig3_pseudorep.csv   needs proteingym_t1.csv specifically (run_proteingym.py --traj t1);
                          this repo's results/ only carries proteingym_t3.csv.
  - fig2_effdim_pooled.csv, fig2_eta2.csv
                          need effdim_grid.csv for t1/t2/t3 (run_decoupling.py --traj t1/t2/t3)
                          with a per-family eta2 breakdown; this repo's results/ only has
                          effdim_grid_t1.csv (+ a whole-grid decoupling_summary_t1.json, which
                          is not broken down per-family), so the paper's pooled/eta2 figure
                          data is committed as-is.
  - betacarotene_oof.csv out-of-fold predictions from the audited beta-carotene law fit; the
                          paper's own prep_data.py does not produce this file either (it ships
                          it directly under figures/data/) -- committed as-is here too.
  - figures/data/anchor_table.csv (read directly by figures/R/fig2.R, not by this script)
                          round-03 was NOT reconstructed for this repo (see
                          docs/reconciliation-notes.md); this is the committed, audited
                          round-03-t3 anchor_table.csv, copied in verbatim.

figures/data/figS_calibration.csv (Fig S2) is fully synthetic (no results/ dependency) and is
regenerated separately by scripts/build_calibration_data.py.

DO NOT expand this script to recompute fig3_partials.csv / fig3_pseudorep.csv /
fig2_effdim_pooled.csv / fig2_eta2.csv from the single t3/t1 trajectory in results/ -- that
would silently drop the paper's multi-trajectory evidence in those panels (e.g. Fig 4c's
partial-correlation forest is supposed to show 3 trajectories per stream). If you have
genuinely regenerated results/{summary,proteingym_summary}_{t1,t2}.json and
results/effdim_grid_{t2,t3}.csv, extend this script deliberately and update this comment --
don't let an automated edit quietly re-add that logic.
"""
from __future__ import annotations
import os
import json

import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
FIG_DATA = os.path.join(ROOT, "figures", "data")
CLEANED = os.path.join(ROOT, "data", "cleaned")

METABOLIC = {"beta-carotene", "tryptophan", "limonene", "isoprenol", "jervis-RBS",
             "vanlent-pCA", "pca-paperA", "smanski-nif"}
PROTEIN_CUR = {"poelwijk", "gb1", "avgfp", "aav"}


def domain(idv, src):
    if src == "proteingym":
        return "protein"
    if idv == "vanlent-simulator":
        return "simulator"
    return "protein" if idv in PROTEIN_CUR else "metabolic"


def build_population():
    """Fig 4a/b (population): curated panel_t3 + ProteinGym proteingym_t3, joined + labelled."""
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
    """Fig 6c (enzyme sensitivity): |Spearman| of each simulator enzyme knob with flux."""
    sim = pd.read_csv(os.path.join(CLEANED, "vanlent-simulator.csv"))
    sim = sim.sample(min(5000, len(sim)), random_state=0)
    enz = [c for c in sim.columns if c != "phenotype"]
    sens = [
        dict(enzyme=c, abs_spearman=abs(stats.spearmanr(sim[c], sim["phenotype"]).correlation))
        for c in enz
    ]
    out = os.path.join(FIG_DATA, "fig5_enzyme_sens.csv")
    pd.DataFrame(sens).sort_values("abs_spearman", ascending=False).to_csv(out, index=False)
    return out


def verify_poelwijk_probe():
    """Fig 5c (Poelwijk link-transform probe): de-hardcoded -- previously 4 typed
    dict(transform=..., r2_law=<const>) literals in prep_data.py; now just a read-through
    check of the recorded audited Poelwijk link-probe outputs (Fig 5C); see
    docs/reconciliation-notes.md. Not re-derived here: the probe needs the round-03 Poelwijk
    driver, which was not reconstructed for this repo (same caveat as anchor_table.csv)."""
    path = os.path.join(FIG_DATA, "fig4_poelwijk_probe.csv")
    probe = pd.read_csv(path)
    need = {"transform", "r2_law", "r2_bb"}
    assert need <= set(probe.columns) and len(probe) == 4, \
        f"fig4_poelwijk_probe.csv at {path} is missing/malformed (expected columns {need}, 4 rows)"
    return path


def main():
    os.makedirs(FIG_DATA, exist_ok=True)
    pop, pop_path = build_population()
    sel_path = build_selection()
    sens_path = build_enzyme_sensitivity()
    probe_path = verify_poelwijk_probe()

    print("WROTE:")
    for p in (pop_path, sel_path, sens_path):
        print(" ", os.path.relpath(p, ROOT))
    print("VERIFIED (pass-through, not regenerated):", os.path.relpath(probe_path, ROOT))

    print("\nfig3_population (n):", len(pop), "| domains:", pop["domain"].value_counts().to_dict())
    print(
        "\nPre-computed under figures/data/ (need multi-trajectory runs or an "
        "unreconstructed round -- see module docstring):\n"
        "  fig3_partials.csv, fig3_pseudorep.csv, fig2_effdim_pooled.csv, fig2_eta2.csv,\n"
        "  betacarotene_oof.csv, anchor_table.csv"
    )


if __name__ == "__main__":
    main()
