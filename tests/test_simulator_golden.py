"""R7 simulator selection payoff (Fig 6B) — reconstruction test.

The round-07 selection driver was not persisted as code; this is a reconstruction using the
ported round-04 log-flux OMP engine (`cdl.sim_engine`) + `cdl.selection`. It reproduces the
CLAIM decisively: the interpretable law recovers more of the top-producing region than the tuned
black box at EVERY training size N. Exact per-cell values differ from the golden (different
sampling RNG), so we assert the direction + rough magnitude, not exact golden numbers.

Controller-confirmed reconstruction (20 seeds):
  law_top100  ≈ {50:22, 100:33, 200:34, 500:34, 1000:32}
  bb_top100   ≈ {50: 5, 100: 6, 200:15, 500:21, 1000:30}
  golden (round-07) law {16.2,21.6,29.0,34.2,34.9} vs bb {6.5,6.9,11.2,18.1,26.1}
"""
import os
import subprocess
import sys
import json

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)


@pytest.mark.slow
def test_law_out_selects_blackbox_at_every_N():
    subprocess.run(
        [sys.executable, os.path.join(ROOT, "scripts", "run_simulator.py"), "--seeds", "20"],
        check=True, cwd=ROOT,
    )
    summ = json.load(open(os.path.join(ROOT, "results", "selection_regime_summary.json")))
    per_N = summ.get("per_N", summ)
    for N in ["50", "100", "200", "500", "1000"]:
        law = per_N[N]["law_top100"]
        bb = per_N[N]["bb_top100"]
        # the central R5 claim: the interpretable law recovers more of the top-100 than the black box
        assert law > bb, f"N={N}: law_top100={law} !> bb_top100={bb}"
    # and the law reaches a substantial fraction of the top region at the largest N
    assert per_N["1000"]["law_top100"] > 25
