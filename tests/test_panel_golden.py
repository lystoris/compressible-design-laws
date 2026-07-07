import os, subprocess, sys, pandas as pd, numpy as np, pytest
HERE = os.path.dirname(__file__); ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "cleaned")

@pytest.mark.slow  # full 13-dataset panel run (~12 min); excluded from the fast suite
@pytest.mark.skipif(not os.path.exists(os.path.join(DATA,"manifest.csv")),
                    reason="cleaned data not present (fetch from Zenodo)")
def test_panel_t3_matches_golden():
    subprocess.run([sys.executable, os.path.join(ROOT,"scripts","run_panel.py"),
                    "--traj","t3","--data-dir",DATA], check=True, cwd=ROOT)
    got = pd.read_csv(os.path.join(ROOT,"results","panel_t3.csv")).set_index("id").sort_index()
    exp = pd.read_csv(os.path.join(HERE,"fixtures/golden/panel_t3.csv")).set_index("id").sort_index()
    common = got.index.intersection(exp.index)
    assert len(common) >= 13
    # seed-locked engine => tight reproduction of compressibility + effective-d
    assert np.allclose(got.loc[common,"r2_law"], exp.loc[common,"r2_law"], atol=0.02)
    assert np.allclose(got.loc[common,"eff_d"], exp.loc[common,"eff_d"], atol=0.15)
