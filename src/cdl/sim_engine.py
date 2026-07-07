"""LOG-FLUX sparse-basis SR engine for the van Lent kinetic-simulator landscape.

Ported (near-verbatim, warnings-filter and figure/CLI code stripped) from
round-04's audited engine:
  research-loop/sr-compressible-design-laws/rounds/round-04/runs/round-04-t1/work/analysis.py

The kinetics are strongly multiplicative (flux ~ EA*EC/EG-type terms), so the
law is fit on TARGET = log(flux) with a library built on BOTH raw enzyme
levels and log-enzyme levels; OrthogonalMatchingPursuit then selects <=10
terms. Fitting directly on raw flux does NOT reproduce the round-04/07
result -- this module is a LIBRARY (no global warnings filter; callers that
want quiet sklearn/scipy output should filter warnings themselves).
"""
import numpy as np
from numpy.random import default_rng
from sklearn.linear_model import OrthogonalMatchingPursuit
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score

SEED = 0
EPS = 1e-3
MAX_TERMS = 10
COLS = ["EA", "EB", "EC", "ED", "EE", "EF", "EG"]


# --------------------------------------------------------------------------
# R2_det = coefficient of determination (NOT correlation^2).
# --------------------------------------------------------------------------
def r2_det(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot <= 1e-12:
        return 0.0
    return 1.0 - ss_res / ss_tot


# --------------------------------------------------------------------------
# SPARSE-BASIS feature library.
# Built on a feature block that stacks RAW enzyme levels and LOG enzyme levels
# (kinetics are multiplicative/saturating). Returns design matrix + names.
# 'saturating' family 1/(1+z), x/(1+z) uses z = train-fold min-max scaled value
# in [0,1] (leakage-free; offsets passed in).  Library size is bounded because
# d=7 is small, so we can afford the full deg-3 cross on the RAW levels plus the
# linear/pairwise on log levels.
# --------------------------------------------------------------------------
def build_library(Xraw, sat_lo, sat_hi):
    """Xraw: (n,7) raw enzyme levels (>0). Returns (Phi, names) WITHOUT an
    explicit intercept column (OMP fits its own intercept)."""
    n, p = Xraw.shape
    Xlog = np.log(Xraw)                      # multiplicative -> additive
    cols, names = [], []

    # ---- RAW block ----
    # linear
    for i in range(p):
        cols.append(Xraw[:, i]); names.append(f"{COLS[i]}")
    # squares
    for i in range(p):
        cols.append(Xraw[:, i] ** 2); names.append(f"{COLS[i]}^2")
    # pairwise products
    for i in range(p):
        for j in range(i + 1, p):
            cols.append(Xraw[:, i] * Xraw[:, j]); names.append(f"{COLS[i]}*{COLS[j]}")
    # degree-3 monomials: cubes + x_i^2*x_j  (d=7 -> affordable)
    for i in range(p):
        cols.append(Xraw[:, i] ** 3); names.append(f"{COLS[i]}^3")
    for i in range(p):
        for j in range(p):
            if i != j:
                cols.append((Xraw[:, i] ** 2) * Xraw[:, j])
                names.append(f"{COLS[i]}^2*{COLS[j]}")
    # reciprocal / saturating on raw min-max scaled value
    rng_sc = np.where((sat_hi - sat_lo) > 1e-12, (sat_hi - sat_lo), 1.0)
    Z = np.clip((Xraw - sat_lo) / rng_sc, 0.0, 1.0)
    for i in range(p):
        cols.append(1.0 / (1.0 + Z[:, i])); names.append(f"1/(1+{COLS[i]})")
    for i in range(p):
        for j in range(p):
            if i != j:
                cols.append(Xraw[:, i] / (1.0 + Z[:, j]))
                names.append(f"{COLS[i]}/(1+{COLS[j]})")
    # also a plain reciprocal 1/x (Michaelis-Menten draining term, esp. EG)
    for i in range(p):
        cols.append(1.0 / Xraw[:, i]); names.append(f"1/{COLS[i]}")

    # ---- LOG block (linear + pairwise in log -> captures power-law kinetics) ----
    for i in range(p):
        cols.append(Xlog[:, i]); names.append(f"log{COLS[i]}")
    for i in range(p):
        for j in range(i + 1, p):
            cols.append(Xlog[:, i] * Xlog[:, j]); names.append(f"log{COLS[i]}*log{COLS[j]}")

    Phi = np.column_stack(cols)
    Phi = np.nan_to_num(Phi, nan=0.0, posinf=0.0, neginf=0.0)
    return Phi, names


def fit_omp_pareto(Xtr_raw, ytr, max_terms=MAX_TERMS):
    """Fit the sparse-basis SR via OMP on TRAIN (leakage-free: library scaling +
    saturating offsets use TRAIN stats only). Returns a dict of per-L models for
    L=1..max_terms, each a closure-ready parameter bundle. ytr is the (log-scale)
    target. Complexity = number of selected non-intercept terms."""
    sat_lo = Xtr_raw.min(axis=0)
    sat_hi = Xtr_raw.max(axis=0)
    Phi_tr, names = build_library(Xtr_raw, sat_lo, sat_hi)
    # standardise every library column using TRAIN stats
    mu = Phi_tr.mean(axis=0)
    sd = Phi_tr.std(axis=0)
    sd = np.where(sd > 1e-12, sd, 1.0)
    Phi_tr_s = (Phi_tr - mu) / sd
    ymu = float(ytr.mean())
    yc = ytr - ymu

    n_feat = Phi_tr_s.shape[1]
    maxL = max(1, min(max_terms, n_feat, len(ytr) - 1))
    models = {}
    for L in range(1, maxL + 1):
        omp = OrthogonalMatchingPursuit(n_nonzero_coefs=L, fit_intercept=False)
        omp.fit(Phi_tr_s, yc)
        models[L] = dict(coef=omp.coef_.copy(), ymu=ymu, mu=mu, sd=sd,
                         sat_lo=sat_lo, sat_hi=sat_hi, names=names)
    return models, maxL


def omp_predict_log(model, Xnew_raw):
    """Predict the (log-scale) target for new raw designs (leakage-free)."""
    Phi, _ = build_library(Xnew_raw, model["sat_lo"], model["sat_hi"])
    Phi_s = (Phi - model["mu"]) / model["sd"]
    return Phi_s @ model["coef"] + model["ymu"]


def law_terms(model):
    """Human-readable active terms (name, standardized coef) sorted by |coef|."""
    coef = model["coef"]
    nz = np.where(np.abs(coef) > 1e-10)[0]
    terms = [(model["names"][k], float(coef[k])) for k in nz]
    terms.sort(key=lambda t: -abs(t[1]))
    return terms


# --------------------------------------------------------------------------
# RandomForest ceiling, tuned by a small inner CV (genuinely-strong denominator,
# per the Round-2 lesson that a loose ceiling distorts lambda).
# --------------------------------------------------------------------------
def fit_rf_tuned(Xtr, ytr, inner_seed=SEED):
    grid = [dict(max_depth=md, max_features=mf)
            for md in (6, 12, None)
            for mf in (1.0, "sqrt")]
    inner = KFold(n_splits=3, shuffle=True, random_state=inner_seed)
    best = (-np.inf, grid[0])
    for params in grid:
        sc = []
        for itr, iva in inner.split(Xtr):
            rf = RandomForestRegressor(n_estimators=120, random_state=inner_seed,
                                       n_jobs=-1, **params)
            rf.fit(Xtr[itr], ytr[itr])
            sc.append(r2_score(ytr[iva], rf.predict(Xtr[iva])))
        m = float(np.mean(sc))
        if m > best[0]:
            best = (m, params)
    _, params = best
    rf = RandomForestRegressor(n_estimators=400, random_state=inner_seed,
                               n_jobs=-1, **params)
    rf.fit(Xtr, ytr)
    return rf, dict(n_estimators=400, **params)


# --------------------------------------------------------------------------
# Selection metrics over the FULL ground-truth space (exact).
# pred_full: predicted score (any monotone scale) for every one of the 279,936
# designs. flux_full: true flux_G.
#
# TIE HANDLING (load-bearing): the sparse law contains saturating basis terms
# 1/(1+Z) and 1/x with Z = train min-max-scaled value CLIPPED to [0,1]. On the
# extreme corner of the design grid these clip to constants, so the law assigns
# the SAME maximal predicted score to a PLATEAU of several designs it genuinely
# cannot tell apart. A naive np.argmax then tie-breaks to the lowest INDEX (an
# arbitrary plateau member), which froze top1_regret at a fixed artefactual value
# across all N. The honest metric for a model that says "these M designs are all
# my best" is the EXPECTED regret over that tied-best plateau (you would pick one
# at random), reported alongside a deterministic random-tiebreak single pick.
# Top-100 overlap uses a deterministic tiny jitter so the 100-cut is well-defined.
# --------------------------------------------------------------------------
def selection_metrics(pred_full, flux_full, true_best, true_top100_set, tie_rng):
    mx = pred_full.max()
    tied = np.where(pred_full >= mx - 1e-9)[0]          # the tied-best plateau
    n_tied = int(tied.size)
    # expected (mean) top-1 regret over the plateau the law calls "best"
    mean_top1_regret = float(np.mean((true_best - flux_full[tied]) / true_best))
    # a single random-tiebreak pick (reproducible via tie_rng)
    pick = int(tied[tie_rng.integers(0, n_tied)])
    pick_top1_regret = float((true_best - flux_full[pick]) / true_best)
    # best-case over the plateau (does the true optimum SIT in the plateau?)
    best_in_plateau_regret = float(np.min((true_best - flux_full[tied]) / true_best))
    # top-100 by predicted score, ties broken by deterministic jitter
    jitter = tie_rng.uniform(-1e-9, 1e-9, size=pred_full.shape[0])
    score = pred_full + jitter
    pred_top100 = set(np.argpartition(-score, 100)[:100].tolist())
    overlap = len(pred_top100 & true_top100_set)
    return dict(mean_top1_regret=mean_top1_regret, pick_top1_regret=pick_top1_regret,
                best_in_plateau_regret=best_in_plateau_regret, n_tied_best=n_tied,
                top100_overlap=int(overlap), pick=pick)
