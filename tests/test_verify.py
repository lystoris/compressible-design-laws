"""verify.py must exit 0 against the committed canonical results/ (the reproduction gate)."""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "results", "panel_t3.csv")),
    reason="results/ not present -- run scripts/run_*.py or `make reproduce` first",
)
def test_verify_passes():
    r = subprocess.run([sys.executable, os.path.join(ROOT, "verify.py")],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, f"verify.py failed:\n{r.stdout}\n{r.stderr}"
