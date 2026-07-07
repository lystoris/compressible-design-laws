import os, subprocess, sys, pandas as pd, numpy as np, pytest
HERE = os.path.dirname(__file__); ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "proteingym")

@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "proteingym_ref.csv")),
                    reason="proteingym data not present (fetch DMS_ProteinGym_substitutions)")
def test_proteingym_t3_matches_golden():
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "run_proteingym.py"),
                    "--traj", "t3", "--data-dir", DATA], check=True, cwd=ROOT)
    got = pd.read_csv(os.path.join(ROOT, "results", "proteingym_t3.csv")).set_index("id").sort_index()
    exp = pd.read_csv(os.path.join(HERE, "fixtures/golden/proteingym_t3.csv")).set_index("id").sort_index()
    common = got.index.intersection(exp.index)
    assert len(common) >= 69
    # seed-locked engine => tight reproduction of compressibility + effective-d
    assert np.allclose(got.loc[common, "r2_law"], exp.loc[common, "r2_law"], atol=0.03)
    assert np.allclose(got.loc[common, "eff_d"], exp.loc[common, "eff_d"], atol=0.2)
