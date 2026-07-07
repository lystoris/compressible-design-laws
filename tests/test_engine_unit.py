import numpy as np
from cdl.compressibility import run_anchor
from cdl.effective_dim import effective_d

def test_run_anchor_recovers_linear_signal():
    rng = np.random.default_rng(0); X = rng.normal(size=(200,3))
    y = 2*X[:,0] - X[:,1]                     # simple linear law
    r2_law, r2_bb = run_anchor(X, y, ["a","b","c"], allow_nonlinear=True, topk=40)
    assert r2_law > 0.9

def test_effective_d_counts_active_groups():
    rng = np.random.default_rng(0); X = rng.normal(size=(300,6))
    y = X[:,0] + X[:,1]                        # 2 active, 4 decoys
    groups = [f"g{i}" for i in range(6)]
    assert 1.5 < effective_d(X, y, groups) < 3.0
