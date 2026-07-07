import os
import numpy as np
import pytest
from cdl.selection import top1_regret, topk_overlap
from cdl.simulator import load_simulator, enzyme_sensitivity

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
DATA_PATH = os.path.join(ROOT, "data", "cleaned", "vanlent-simulator.csv")


def test_top1_regret_zero_when_argmax_matches():
    yt = np.array([1, 2, 3, 10])
    yp = np.array([0, 0, 0, 5])
    assert top1_regret(yt, yp) == 0.0


def test_top1_regret_positive_when_wrong():
    yt = np.array([1, 2, 3, 10])
    yp = np.array([5, 0, 0, 0])  # picks index 0 (yt=1), not the true best (yt=10)
    regret = top1_regret(yt, yp)
    assert regret > 0.0
    assert regret == pytest.approx((10 - 1) / 10)


def test_topk_overlap_full():
    yt = np.arange(200)
    yp = np.arange(200)
    assert topk_overlap(yt, yp, k=100) == 100


def test_topk_overlap_partial():
    # True top-100 (by yt) is indices [100..199].
    yt = np.arange(200)
    # Predicted top-100 is exactly indices [50..149]: those get distinct
    # nonnegative scores, everything else is pinned below zero so there is
    # no ambiguity about which 100 indices are "predicted top-100".
    yp = np.full(200, -1.0)
    yp[50:150] = np.arange(100)
    # Overlap of {100..199} and {50..149} is {100..149} -> 50 indices.
    assert topk_overlap(yt, yp, k=100) == 50


@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(DATA_PATH), reason="vanlent-simulator.csv not present")
def test_load_simulator_shape():
    X, flux, names = load_simulator(DATA_PATH)
    assert X.shape == (279936, 7)
    assert len(flux) == 279936
    assert names == ["EA", "EB", "EC", "ED", "EE", "EF", "EG"]

    sens = enzyme_sensitivity(X, flux)
    ranked = sorted(sens, key=sens.get, reverse=True)
    assert ranked[0] == "EG"
    assert sens["EG"] > 0.7
    assert 0.3 < sens["EA"] < 0.5
    for name in ("EB", "ED", "EE", "EF"):
        assert sens[name] < 0.1, (name, sens[name])
