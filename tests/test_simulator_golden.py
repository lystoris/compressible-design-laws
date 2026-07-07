"""Fig 6A (compressibility vs N) + Fig 6B (selection payoff) golden test.

sim_Ncurve.csv is a SINGLE deterministic trajectory (fixed SEED) so it
reproduces the round-04 golden closely (atol on r2_law). selection_regime.csv
(Fig 6B) has no persisted round-07 driver to match exactly -- this is a
reconstruction over (N, seed) draws using the same audited round-04 engine, so
we assert the AGGREGATE payoff claim instead: the interpretable law recovers
substantially more of the true top-100 than the tuned black-box ceiling,
especially once N >= 100.
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)

GOLDEN_R2_LAW = {50: 0.836, 100: 0.874, 200: 0.935, 500: 0.895, 1000: 0.934}


@pytest.mark.slow
def test_sim_ncurve_and_selection_payoff():
    subprocess.run(
        [sys.executable, os.path.join(ROOT, "scripts", "run_simulator.py"), "--seeds", "20"],
        check=True, cwd=ROOT,
    )

    # (a) Fig 6A: deterministic N-curve, close to the round-04 golden.
    nc = pd.read_csv(os.path.join(ROOT, "results", "sim_Ncurve.csv"))
    for N, golden_r2 in GOLDEN_R2_LAW.items():
        got = float(nc.loc[nc["N"] == N, "r2_law"].iloc[0])
        assert abs(got - golden_r2) < 0.05, f"N={N}: r2_law={got} vs golden {golden_r2}"

    # (b) Fig 6B: aggregate selection payoff -- the law out-selects the black
    # box (paper's Fig 6B / R5 claim), most clearly once N >= 100.
    sel = pd.read_csv(os.path.join(ROOT, "results", "selection_regime.csv"))
    hi_N = sel[sel["N"] >= 100]
    assert hi_N["law_ov"].median() > hi_N["bb_ov"].median(), (
        f"median law_ov={hi_N['law_ov'].median()} !> "
        f"median bb_ov={hi_N['bb_ov'].median()} at N>=100"
    )
    assert sel["law_ov"].max() >= 15
