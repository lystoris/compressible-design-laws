import os, subprocess, sys
import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

HERE = os.path.dirname(__file__); ROOT = os.path.dirname(HERE)


@pytest.mark.slow
def test_sweep_median_falloff():
    """Reproduce the round-02-t3 audited 405-cell sweep and check the Fig 3A/B headline:
    median r2_law falloff with d, and Spearman(Lambda-driving r2_law, selection fitness)."""
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "run_sweep.py"), "--traj", "t3"],
                    check=True, cwd=ROOT)
    g = pd.read_csv(os.path.join(ROOT, "results", "sweep_grid.csv"))
    assert len(g) == 405

    med = g.groupby("d")["r2_law"].median()
    for d, exp in [(3, 0.98), (9, 0.75), (20, 0.44), (67, 0.045)]:
        assert abs(med.loc[d] - exp) < 0.06, (d, med.loc[d])

    rho = spearmanr(g["r2_law"], 1 - g["law_top1_regret"]).correlation
    assert abs(rho - 0.71) < 0.08

    golden = pd.read_csv(os.path.join(HERE, "fixtures", "golden", "sweep_grid.csv"))
    assert len(golden) == 405
    gmed = golden.groupby("d")["r2_law"].median()
    assert np.allclose(med.values, gmed.reindex(med.index).values, atol=0.08)
