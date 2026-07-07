#!/usr/bin/env python3
"""Assemble tidy CSVs for the paper figures from THIS REPO's reproduced results/.

Ported from paper/v2/figures/data/prep_data.py. The original read from a deep
research-loop round-08 work directory; this version reads only from this repo's own
`results/`.

Run: /usr/bin/python3 scripts/build_figure_data.py
Writes into figures/data/.

------------------------------------------------------------------------------------
NON-DESTRUCTIVE BY DESIGN
------------------------------------------------------------------------------------
The committed `figures/data/*.csv` are the PAPER-CANONICAL figure inputs, so the
committed/rendered figures match the published paper. Several panels need data from
ALL THREE analysis trajectories (t1/t2/t3), which a default reproduction run does not
produce (this repo reproduces trajectory t3 for the curated/ProteinGym panels and t1
for the decoupling grid -- see docs/reconciliation-notes.md). To avoid silently
replacing a complete canonical CSV with a partial one, each builder below only writes
when its FULL required inputs are present in results/; otherwise it PRESERVES the
committed canonical file and prints a [preserve] line explaining what to run to
regenerate it.

Coverage after a default run (t3 panel + t1 decoupling + simulator):
  REGENERATED  fig3_population.csv   (panel_t3 + proteingym_t3 -- matches paper)
               fig5_enzyme_sens.csv  (full simulator table -- matches paper)
  PRESERVED    fig3_partials.csv     (needs summary_t{1,2,3} + proteingym_summary_t{1,2,3})
               fig3_pseudorep.csv    (paper used proteingym_t1)
               fig2_effdim_pooled.csv / fig2_eta2.csv (need effdim_grid_t{1,2,3})
               fig5_selection.csv    (paper's audited round-07; run_simulator reproduces
                                      the CLAIM -- verified by tests/test_simulator_golden.py
                                      -- but per-cell values differ; see R-3-style note)
  RECORDED     fig4_poelwijk_probe.csv, anchor_table.csv, betacarotene_oof.csv,
               figS_calibration.csv  (committed audited artifacts; never written here)

To FULLY regenerate the multi-trajectory panels from source, first run all trajectories
(`run_panel --traj t1|t2|t3`, `run_proteingym --traj t1|t2|t3`, `run_decoupling --traj
t1|t2|t3`) then re-run this script.

The Poelwijk probe (fig4_poelwijk_probe.csv) closes a real honesty gap: prep_data.py
hand-typed its four numbers inline as dict literals despite claiming "no hand-typed
stats"; they are now a committed CSV read (never re-typed), with provenance recorded.
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

_wrote, _preserved = [], []


def _have(*names):
    """True iff every named results/ file exists."""
    return all(os.path.exists(os.path.join(RESULTS, n)) for n in names)


def _preserve(fig_name, needed):
    _preserved.append((fig_name, needed))
    print(f"  [preserve] figures/data/{fig_name} (canonical) -- to regenerate, run: {needed}")


def domain(idv, src):
    if src == "proteingym":
        return "protein"
    if idv == "vanlent-simulator":
        return "simulator"
    return "protein" if idv in PROTEIN_CUR else "metabolic"


def build_population():
    """Fig 4a/b (population scatter): curated panel_t3 + ProteinGym proteingym_t3 (t3 == the
    paper's population trajectory, so this matches the canonical)."""
    if not _have("panel_t3.csv", "proteingym_t3.csv"):
        _preserve("fig3_population.csv", "run_panel --traj t3 && run_proteingym --traj t3")
        return None
    cur = pd.read_csv(os.path.join(RESULTS, "panel_t3.csv")); cur["source"] = "curated"
    pg = pd.read_csv(os.path.join(RESULTS, "proteingym_t3.csv"))
    pg["source"] = "proteingym"; pg["group"] = "protein-dms"
    pop = pd.concat([cur, pg], ignore_index=True)
    pop["domain"] = [domain(i, s) for i, s in zip(pop["id"], pop["source"])]
    pop["compresses"] = (pop["r2_law"] >= 0.7).map({True: "full", False: "partial-or-none"})
    pop.to_csv(os.path.join(FIG_DATA, "fig3_population.csv"), index=False)
    _wrote.append("fig3_population.csv")
    return pop


def build_partials():
    """Fig 4c (partial-correlation forest): needs all three trajectories per stream."""
    need = [f"summary_{t}.json" for t in ("t1", "t2", "t3")] + \
           [f"proteingym_summary_{t}.json" for t in ("t1", "t2", "t3")]
    if not _have(*need):
        _preserve("fig3_partials.csv", "run_panel/run_proteingym for --traj t1,t2,t3")
        return
    rows = []
    for stream, prefix in [("curated", "summary"), ("ProteinGym", "proteingym_summary")]:
        for t in ("t1", "t2", "t3"):
            s = json.load(open(os.path.join(RESULTS, f"{prefix}_{t}.json")))
            rows.append(dict(stream=stream, traj=t, axis="effective d", partial=s["partial_effd"]))
            rows.append(dict(stream=stream, traj=t, axis="nominal d", partial=s["partial_nomd"]))
    pd.DataFrame(rows).to_csv(os.path.join(FIG_DATA, "fig3_partials.csv"), index=False)
    _wrote.append("fig3_partials.csv")


def build_pseudorep():
    """Fig 4d: ProteinGym pseudoreplication-corrected rho. The paper's canonical uses trajectory
    t1; regenerate only if proteingym_t1 is present (else preserve canonical)."""
    if not _have("proteingym_t1.csv"):
        _preserve("fig3_pseudorep.csv", "run_proteingym --traj t1")
        return
    d = pd.read_csv(os.path.join(RESULTS, "proteingym_t1.csv"))
    d["is_tsu"] = d["id"].str.contains("Tsuboyama")
    d["study"] = d["id"].str.extract(r"_([A-Za-z0-9]+_\d{4})")

    def rho(sub):
        return stats.spearmanr(sub.r2_law, sub.eff_d).correlation

    g = d.groupby(d["study"].where(~d.is_tsu, "Tsuboyama_2023")).agg(
        r2_law=("r2_law", "median"), eff_d=("eff_d", "median"))
    prow = [
        dict(level="all 69 (pooled)", n=len(d), rho=rho(d)),
        dict(level="Tsuboyama within-protocol", n=int(d.is_tsu.sum()), rho=rho(d[d.is_tsu])),
        dict(level="15 other studies", n=int((~d.is_tsu).sum()), rho=rho(d[~d.is_tsu])),
        dict(level="study-collapsed", n=len(g), rho=stats.spearmanr(g.r2_law, g.eff_d).correlation),
    ]
    pd.DataFrame(prow).to_csv(os.path.join(FIG_DATA, "fig3_pseudorep.csv"), index=False)
    _wrote.append("fig3_pseudorep.csv")


def build_effdim_and_eta2():
    """Fig 3c/3d (decoupling boxplots + eta2): need effdim_grid_t{1,2,3}."""
    need = [f"effdim_grid_{t}.csv" for t in ("t1", "t2", "t3")]
    if not _have(*need):
        _preserve("fig2_effdim_pooled.csv", "run_decoupling --traj t1,t2,t3")
        _preserve("fig2_eta2.csv", "run_decoupling --traj t1,t2,t3")
        return
    pooled = []
    for t in ("t1", "t2", "t3"):
        e = pd.read_csv(os.path.join(RESULTS, f"effdim_grid_{t}.csv")); e["traj"] = t
        pooled.append(e)
    pooled = pd.concat(pooled, ignore_index=True)
    pooled.to_csv(os.path.join(FIG_DATA, "fig2_effdim_pooled.csv"), index=False)
    rows = []
    for (t, fam), sub in pooled.groupby(["traj", "family"]):
        rows.append(dict(traj=t, family=fam,
                         eta2_eff=eta2(sub.r2_law, sub.d_eff),
                         eta2_nom=eta2(sub.r2_law, sub.d_nom)))
    pd.DataFrame(rows).to_csv(os.path.join(FIG_DATA, "fig2_eta2.csv"), index=False)
    _wrote.extend(["fig2_effdim_pooled.csv", "fig2_eta2.csv"])


def build_selection():
    """Fig 6b (selection payoff). The committed canonical is the paper's audited round-07 output.
    run_simulator.py reproduces the CLAIM (law out-selects bb at every N; verified by
    tests/test_simulator_golden.py) but the round-07 driver was not persisted, so per-cell values
    differ -- the canonical is preserved for figure fidelity. Set env CDL_REGEN_SELECTION=1 to
    overwrite it with this repo's reconstruction."""
    if os.environ.get("CDL_REGEN_SELECTION") != "1" or not _have("selection_regime_summary.json"):
        _preserve("fig5_selection.csv",
                  "run_simulator (then set CDL_REGEN_SELECTION=1 to plot the reconstruction)")
        return
    sj = json.load(open(os.path.join(RESULTS, "selection_regime_summary.json")))
    srows = [dict(N=int(N), law_top100=v["law_top100"], bb_top100=v["bb_top100"],
                  law_top1=v["law_top1"], bb_top1=v["bb_top1"]) for N, v in sj["per_N"].items()]
    pd.DataFrame(srows).sort_values("N").to_csv(os.path.join(FIG_DATA, "fig5_selection.csv"), index=False)
    _wrote.append("fig5_selection.csv")


def build_enzyme_sensitivity():
    """Fig 6c (enzyme sensitivity): |Spearman| of each simulator enzyme with flux over the full
    enumerated table (deterministic; matches the paper canonical to 2 dp)."""
    sim_path = os.path.join(CLEANED, "vanlent-simulator.csv")
    if not os.path.exists(sim_path):
        _preserve("fig5_enzyme_sens.csv", "fetch data/cleaned/vanlent-simulator.csv (make data)")
        return
    X, flux, _names = load_simulator(sim_path)
    sens = enzyme_sensitivity(X, flux)
    pd.DataFrame([dict(enzyme=n, abs_spearman=v) for n, v in sens.items()]) \
        .sort_values("abs_spearman", ascending=False) \
        .to_csv(os.path.join(FIG_DATA, "fig5_enzyme_sens.csv"), index=False)
    _wrote.append("fig5_enzyme_sens.csv")


def verify_recorded_artifacts():
    """The four RECORDED AUDITED ARTIFACTS must already be committed in figures/data/."""
    missing = [f for f in RECORDED_ARTIFACTS if not os.path.exists(os.path.join(FIG_DATA, f))]
    if missing:
        raise FileNotFoundError(
            "Missing recorded audited artifact(s) in figures/data/: " + ", ".join(missing) +
            ". These are committed source artifacts (see the header) -- restore them from their "
            "documented provenance rather than fabricating replacement values.")
    probe = pd.read_csv(os.path.join(FIG_DATA, "fig4_poelwijk_probe.csv"))
    assert {"transform", "r2_law", "r2_bb"} <= set(probe.columns) and len(probe) == 4, \
        "fig4_poelwijk_probe.csv is missing/malformed (expected transform,r2_law,r2_bb; 4 rows)"


def main():
    os.makedirs(FIG_DATA, exist_ok=True)
    print("Building figure data from results/ (non-destructive: partial inputs preserve canonical)\n")
    build_population()
    build_partials()
    build_pseudorep()
    build_effdim_and_eta2()
    build_selection()
    build_enzyme_sensitivity()
    verify_recorded_artifacts()

    print("\nREGENERATED from results/:", ", ".join(_wrote) or "(none)")
    print("PRESERVED (canonical; need full-trajectory runs):",
          ", ".join(f for f, _ in _preserved) or "(none)")
    print("RECORDED artifacts present:", ", ".join(RECORDED_ARTIFACTS))


if __name__ == "__main__":
    main()
