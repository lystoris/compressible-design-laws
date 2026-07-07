import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")

# verify.py reads the full reproduced canon from results/; that directory isn't
# committed until Task 15, so skip cleanly if it's absent (e.g. a fresh checkout).
_REQUIRED = [
    "panel_t3.csv", "summary_t3.json", "proteingym_t3.csv",
    "proteingym_summary_t3.json", "sweep_grid.csv", "effdim_grid_t1.csv",
    "decoupling_summary_t1.json", "sim_Ncurve.csv", "selection_regime.csv",
]


@pytest.mark.skipif(
    not all(os.path.exists(os.path.join(RESULTS, f)) for f in _REQUIRED),
    reason="results/ canon not present (generate via scripts/run_*.py or fetch release artifacts)",
)
def test_verify_exits_zero_on_reproduced_canon():
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify.py")],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"verify.py failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
    )
    assert "PASS" in proc.stdout
    assert "FAIL" not in proc.stdout
