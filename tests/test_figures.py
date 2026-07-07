import os
import subprocess
import sys

import pandas as pd
import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
FIG_DATA = os.path.join(ROOT, "figures", "data")

# figure-data files build_figure_data.py regenerates from this repo's committed results/
REGENERATED = ["fig3_population.csv", "fig5_selection.csv", "fig5_enzyme_sens.csv"]

# figure-data files shipped pre-computed (need multi-trajectory runs or an unreconstructed
# round -- see the module docstring in scripts/build_figure_data.py); build_figure_data.py
# must leave these alone, so we only check they exist and are non-empty.
PRECOMPUTED = [
    "fig3_partials.csv", "fig3_pseudorep.csv", "fig2_effdim_pooled.csv", "fig2_eta2.csv",
    "betacarotene_oof.csv", "anchor_table.csv", "fig4_poelwijk_probe.csv",
]

REQUIRED = os.path.join(ROOT, "results", "panel_t3.csv")


@pytest.mark.skipif(not os.path.exists(REQUIRED),
                    reason="results/ not present (run scripts/run_panel.py etc. first)")
def test_build_figure_data_runs_and_writes_expected_csvs():
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "build_figure_data.py")],
                    check=True, cwd=ROOT)
    for name in REGENERATED + PRECOMPUTED:
        path = os.path.join(FIG_DATA, name)
        assert os.path.exists(path), f"missing figure-data file: {name}"
        df = pd.read_csv(path)
        assert len(df) > 0, f"empty figure-data file: {name}"


def test_precomputed_figure_data_present():
    # These ship committed regardless of whether results/ is present (see docstring in
    # scripts/build_figure_data.py for why they are not regenerated from a single trajectory).
    for name in PRECOMPUTED:
        path = os.path.join(FIG_DATA, name)
        assert os.path.exists(path), f"missing figure-data file: {name}"
        df = pd.read_csv(path)
        assert len(df) > 0, f"empty figure-data file: {name}"


def test_r_scripts_present_and_repointed():
    r_dir = os.path.join(ROOT, "figures", "R")
    expected = ["fig2.R", "fig3.R", "fig4.R", "fig5.R", "fig6.R",
                "figS_calibration.R", "figS_sweep_phase.R", "make_all.R", "theme_pub.R"]
    for name in expected:
        path = os.path.join(r_dir, name)
        assert os.path.exists(path), f"missing R script: {name}"
    # no leftover references to the source repo's research-loop / paper/v2 paths
    for name in expected:
        text = open(os.path.join(r_dir, name)).read()
        assert "research-loop" not in text, f"{name} still references research-loop"
        assert "paper/v2" not in text, f"{name} still references paper/v2"
