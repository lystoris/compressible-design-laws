import numpy as np
from cdl.stats import r2_det, partial_spearman, eta2

def test_r2_det_perfect():
    assert r2_det([1,2,3],[1,2,3]) == 1.0

def test_r2_det_zero_variance_is_nan():
    assert np.isnan(r2_det([2,2,2],[1,2,3]))

def test_partial_spearman_removes_confound():
    rng = np.random.default_rng(0)
    c = rng.normal(size=200); a = c + 0.01*rng.normal(size=200); b = c + 0.01*rng.normal(size=200)
    # a and b correlate only through c -> partial ~ 0
    assert abs(partial_spearman(a, b, c)) < 0.2

def test_eta2_all_between():
    vals = [1,1,1,5,5,5]; grp = ["x","x","x","y","y","y"]
    assert eta2(vals, grp) > 0.99
