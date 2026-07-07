import numpy as np


def safe_div(a, b):
    b = np.asarray(b, float)
    b = np.where(np.abs(b) < 1e-6, 1e-6, b)
    return np.asarray(a, float) / b


def _make_generator(d, k, encoding, gen_seed):
    """Return (sample_X, true_f). Ported from round-02 analysis.py make_generator:
    equal-magnitude linear term per var + Michaelis-Menten saturating terms on d//3
    vars + pairwise (k>=2) and three-way (k>=3) interactions."""
    rs = np.random.RandomState(gen_seed)
    w = rs.uniform(0.5, 1.5, size=d) * rs.choice([-1, 1], size=d)
    n_sat = max(1, d // 3)
    sat_idx = rs.choice(d, size=n_sat, replace=False)
    Vmax = rs.uniform(1.0, 3.0, size=n_sat)
    Km = rs.uniform(0.5, 2.0, size=n_sat)
    pairs = []
    if k >= 2:
        for _ in range(max(1, d // 2)):
            i, j = rs.choice(d, size=2, replace=False)
            pairs.append((i, j, rs.uniform(0.4, 1.0) * rs.choice([-1, 1])))
    triples = []
    if k >= 3:
        for _ in range(max(1, d // 3)):
            i, j, l = rs.choice(d, size=3, replace=False)
            triples.append((i, j, l, rs.uniform(0.3, 0.8) * rs.choice([-1, 1])))

    def encode(Xraw):
        if encoding == "integer":
            return np.round(Xraw * 5).astype(float)     # copy-number 0..5
        if encoding == "binary":
            return (Xraw > 0.5).astype(float)
        return Xraw.astype(float)                         # continuous 0..1

    def sample_X(N, seed):
        r = np.random.RandomState(seed)
        return encode(r.uniform(0, 1, size=(N, d)))

    def true_f(X):
        X = np.asarray(X, float)
        Xn = X / 5.0 if encoding == "integer" else X
        out = Xn @ w
        for m, j in enumerate(sat_idx):
            out += Vmax[m] * safe_div(Xn[:, j], Km[m] + Xn[:, j])
        for (i, j, c) in pairs:
            out += c * Xn[:, i] * Xn[:, j]
        for (i, j, l, c) in triples:
            out += c * Xn[:, i] * Xn[:, j] * Xn[:, l]
        return out

    return sample_X, true_f


def make_sweep_cell(d, k, sigma, N, encoding, seed):
    """One 405-sweep cell: linear + saturating + interaction landscape, sampled N times
    and noised at sigma*std(y_clean). Returns (X, y, feat)."""
    sample_X, true_f = _make_generator(d, k, encoding, gen_seed=seed)
    X = sample_X(N, seed)
    y_clean = true_f(X)
    rng = np.random.default_rng(seed)
    y = y_clean + sigma * np.std(y_clean) * rng.standard_normal(N)
    feat = [f"x{i}" for i in range(d)]
    return X, y, feat


def make_decoupled(d_eff, d_nom, family, rep, N=500, noise=0.1, seed=0):
    """d_eff genuinely active variables among d_nom total (rest are inert decoys).
    family in {"additive", "michaelis_menten", "random_gp"}. Returns (X, y, groups)."""
    combo = (int(seed), int(d_eff), int(d_nom), int(rep), hash(family) & 0xFFFF)
    rng = np.random.default_rng(combo)
    X = rng.uniform(0, 1, size=(N, d_nom))
    Xa = X[:, :d_eff]

    if family == "additive":
        w = rng.uniform(0.5, 1.5, size=d_eff) * rng.choice([-1, 1], size=d_eff)
        signal = Xa @ w
    elif family == "michaelis_menten":
        Vmax = rng.uniform(1.0, 3.0, size=d_eff)
        Km = rng.uniform(0.5, 2.0, size=d_eff)
        signal = np.sum(Vmax * safe_div(Xa, Km + Xa), axis=1)
    elif family == "random_gp":
        n_feat = 40
        Wf = rng.normal(size=(d_eff, n_feat))
        b = rng.uniform(0, 2 * np.pi, n_feat)
        beta = rng.normal(size=n_feat)
        signal = np.sqrt(2.0 / n_feat) * np.cos(Xa @ Wf + b) @ beta
    else:
        raise ValueError(f"unknown family: {family!r}")

    y = signal + noise * np.std(signal) * rng.standard_normal(N)
    groups = [f"g{i}" for i in range(d_nom)]
    return X, y, groups
