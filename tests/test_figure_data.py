"""Smoke test for the figure-data build: runs scripts/build_figure_data.py against the
committed results/ and checks every expected figures/data/*.csv was (re)written and is
non-empty. Fast (only reads results/ + data/cleaned/vanlent-simulator.csv), not marked slow.
"""
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
FIG_DATA = os.path.join(ROOT, "figures", "data")

# recomputed fresh, each run, from results/ (see scripts/build_figure_data.py provenance docstring)
RECOMPUTED = [
    "fig3_population.csv",
    "fig3_partials.csv",
    "fig3_pseudorep.csv",
    "fig2_effdim_pooled.csv",
    "fig2_eta2.csv",
    "fig5_selection.csv",
    "fig5_enzyme_sens.csv",
]

# recorded audited artifacts: committed verbatim, only checked (never written) by the build script
RECORDED = [
    "fig4_poelwijk_probe.csv",
    "anchor_table.csv",
    "betacarotene_oof.csv",
    "figS_calibration.csv",
]


def test_build_figure_data_produces_all_csvs():
    subprocess.run(
        [sys.executable, os.path.join(ROOT, "scripts", "build_figure_data.py")],
        check=True, cwd=ROOT,
    )
    for fname in RECOMPUTED + RECORDED:
        path = os.path.join(FIG_DATA, fname)
        assert os.path.exists(path), f"missing {path}"
        df = pd.read_csv(path)
        assert len(df) > 0, f"{path} is empty"
        assert len(df.columns) > 0, f"{path} has no columns"


def test_fig4_poelwijk_probe_has_no_hand_typed_stand_in():
    # The de-hardcoding fix: fig4_poelwijk_probe.csv must be a real, labelled 4-row table
    # (not regenerated inline from Python literals -- see build_figure_data.py docstring).
    df = pd.read_csv(os.path.join(FIG_DATA, "fig4_poelwijk_probe.csv"))
    assert len(df) == 4
    assert set(df["transform"]) == {"raw", "log", "rank-normal", "interpretable (binary, linear+interaction)"}


def test_fig3_population_matches_curated_plus_proteingym():
    pop = pd.read_csv(os.path.join(FIG_DATA, "fig3_population.csv"))
    cur = pd.read_csv(os.path.join(ROOT, "results", "panel_t3.csv"))
    pg = pd.read_csv(os.path.join(ROOT, "results", "proteingym_t3.csv"))
    assert len(pop) == len(cur) + len(pg)
    assert set(pop["domain"]) <= {"metabolic", "protein", "simulator"}
