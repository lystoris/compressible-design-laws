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


def test_topk_overlap_counts_intersection():
    yt = np.arange(200)
    yp = np.arange(200)
    assert topk_overlap(yt, yp, k=100) == 100


@pytest.mark.skipif(not os.path.exists(DATA_PATH), reason="vanlent-simulator.csv not present")
def test_load_simulator_shape():
    X, flux, enzyme_names = load_simulator(DATA_PATH)
    assert X.shape == (279936, 7)
    assert len(flux) == 279936
    assert enzyme_names == ["EA", "EB", "EC", "ED", "EE", "EF", "EG"]


@pytest.mark.skipif(not os.path.exists(DATA_PATH), reason="vanlent-simulator.csv not present")
def test_enzyme_sensitivity_reproduces_fig6c():
    X, flux, enzyme_names = load_simulator(DATA_PATH)
    sens = enzyme_sensitivity(X, flux, enzyme_names)
    expected = {"EG": 0.80, "EA": 0.40, "EC": 0.24}
    for name, exp in expected.items():
        assert abs(sens[name] - exp) < 0.03, (name, sens[name])
    for name in ("EB", "ED", "EE", "EF"):
        assert sens[name] < 0.10, (name, sens[name])
