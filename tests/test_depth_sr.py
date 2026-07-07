import numpy as np

from cdl.depth_sr import r2_det


def test_r2_det_degenerate_target_returns_zero_not_nan():
    # Degenerate (zero-variance) target -> 0.0, NOT nan. A nan here would make every
    # `nan > best` comparison False and silently freeze the depth/alpha grid search.
    val = r2_det([2, 2, 2], [1, 2, 3])
    assert val == 0.0
    assert not np.isnan(val)


def test_r2_det_perfect_fit_returns_one():
    assert r2_det([1, 2, 3], [1, 2, 3]) == 1.0
