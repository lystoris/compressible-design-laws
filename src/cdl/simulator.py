import os
import numpy as np
import pandas as pd
from scipy import stats

ENZYME_NAMES = ["EA", "EB", "EC", "ED", "EE", "EF", "EG"]
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_PATH = os.path.join(_REPO_ROOT, "data", "cleaned", "vanlent-simulator.csv")


def load_simulator(path=DEFAULT_PATH):
    """Load the van Lent kinetic-simulator ground-truth table: the full
    6^7 = 279,936 enzyme-level design enumeration with table-driven flux
    (NOT a computed kinetic model). Returns (X, flux, enzyme_names)."""
    df = pd.read_csv(path)
    X = df[ENZYME_NAMES].to_numpy(float)
    flux = df["phenotype"].to_numpy(float)
    return X, flux, list(ENZYME_NAMES)


def enzyme_sensitivity(X, flux, enzyme_names):
    """Absolute Spearman correlation of each enzyme column with flux, over
    all rows. Reproduces Fig 6C."""
    X = np.asarray(X, float); flux = np.asarray(flux, float)
    return {
        name: float(abs(stats.spearmanr(X[:, i], flux).correlation))
        for i, name in enumerate(enzyme_names)
    }
