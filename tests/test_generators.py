import numpy as np
from numpy.linalg import lstsq

from cdl.generators import make_sweep_cell, make_decoupled


def test_sweep_shapes():
    X, y, feat = make_sweep_cell(d=5, k=2, sigma=0.1, N=100, encoding="continuous", seed=1)
    assert X.shape == (100, 5) and len(y) == 100 and len(feat) == 5


def test_decoupled_only_active_drive_signal():
    X, y, groups = make_decoupled(d_eff=3, d_nom=9, family="additive", rep=0, N=500, noise=0.0, seed=0)
    assert X.shape == (500, 9)
    # signal must be recoverable from the first 3 columns alone
    for cols, thresh in [(slice(0, 3), 0.9), (slice(3, 9), 0.2)]:
        A = np.column_stack([X[:, cols], np.ones(500)])
        pred = A @ lstsq(A, y, rcond=None)[0]
        r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
        assert (r2 > thresh) if cols == slice(0, 3) else (r2 < thresh)


def test_decoupled_is_deterministic():
    a = make_decoupled(3, 9, "additive", rep=0, N=100, noise=0.1, seed=0)
    b = make_decoupled(3, 9, "additive", rep=0, N=100, noise=0.1, seed=0)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and a[2] == b[2]
