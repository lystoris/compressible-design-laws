"""DEPTH-budgeted symbolic-regression engine (compressibility scorer) + KRR black-box ceiling.

Ported VERBATIM from round-02 trajectory t3's audited engine:
research-loop/sr-compressible-design-laws/rounds/round-02/runs/round-02-t3/work/analysis.py
(functions safe_div/safe_inv/standardize_cols/screen_keep, build_feature_bank, greedy_law/
predict_law, the KernelRidge(RBF) black-box ceiling (_median_gamma/_cap_subsample/
krr_fit_predict/krr_cv_r2), lambda_depth_nestedcv (the compressibility scorer, DEPTH budget +
nested k-fold), and selection_fitness). r2_det is the round-02 engine's OWN local version
(threshold 1e-12, returns 0.0 on degenerate variance) — deliberately NOT cdl.stats.r2_det
(which is faithful to panel.py and returns np.nan on degenerate variance). The distinction
matters: a NaN would make every `NaN > best` comparison False and silently freeze the
depth/alpha grid search on its first candidate, so the engine keeps its own 0.0 semantics.

A "law" is a sparse linear combination of basis-feature trees built by a recursive grammar
whose budget is expression DEPTH. Leaves (depth 1) are raw vars + saturating transforms of a
single var (Michaelis-Menten x/(1+|x|), 1/x, x^2). Depth-2 nodes combine two leaves with
{+,-,*,/}; depth-3 nodes combine a depth-2 node with a leaf. Complexity charged = DEPTH of the
deepest selected tree (node count kept only as a diagnostic).

This module is a LIBRARY: no global warnings filters here.
"""
import itertools
import numpy as np
from scipy import stats

from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import KFold

# ----------------------------------------------------------------------------------------
# 0. Helpers
# ----------------------------------------------------------------------------------------


def r2_det(y_true, y_pred):
    """Coefficient of determination on the held-out fold. NOT correlation^2.

    Ported VERBATIM from round-02 analysis.py (~line 56): returns 0.0 (NOT np.nan) when the
    target variance is degenerate (<= 1e-12). This 0.0 semantics is load-bearing — the engine's
    grid searches compare `r2 > best`, and a NaN there would freeze selection on the first
    candidate. Distinct on purpose from cdl.stats.r2_det (panel.py-faithful, returns np.nan)."""
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot <= 1e-12:
        return 0.0
    return 1.0 - ss_res / ss_tot


def safe_div(a, b):
    b = np.asarray(b, float)
    b = np.where(np.abs(b) < 1e-6, 1e-6, b)
    return np.asarray(a, float) / b


def safe_inv(a):
    return safe_div(np.ones_like(np.asarray(a, float)), a)


def standardize_cols(X):
    X = np.asarray(X, float)
    mu = X.mean(0)
    sd = X.std(0)
    sd = np.where(sd < 1e-9, 1.0, sd)
    return (X - mu) / sd, mu, sd


def screen_keep(X, y, max_vars):
    """Top-max_vars variables by |Spearman| with y (variable pre-screen for high d).
    Fit on TRAIN data only by the caller."""
    d = X.shape[1]
    if d <= max_vars:
        return np.arange(d)
    scores = np.zeros(d)
    for j in range(d):
        try:
            s = stats.spearmanr(X[:, j], y).correlation
            scores[j] = 0.0 if not np.isfinite(s) else abs(s)
        except Exception:
            scores[j] = 0.0
    return np.sort(np.argsort(scores)[::-1][:max_vars])

# ----------------------------------------------------------------------------------------
# 1. DEPTH-budgeted symbolic regression (compact recursive tree search)
# ----------------------------------------------------------------------------------------

UNARY = {
    "inv": safe_inv,
    "mm":  lambda x: safe_div(x, 1.0 + np.abs(x)),   # Michaelis-Menten-like saturation
    "sq":  lambda x: np.asarray(x, float) ** 2,
}
BINARY = {
    "+": lambda a, b: a + b,
    "-": lambda a, b: a - b,
    "*": lambda a, b: a * b,
    "/": lambda a, b: safe_div(a, b),
}


def build_feature_bank(Xtr, depth_budget, max_vars):
    """Return list of (name, depth, col_train, builder_fn). builder_fn maps a standardized X
    matrix -> a column, so the IDENTICAL transform is applied to held-out folds (no leakage).
    Bounded enumeration keeps high-d / depth-3 tractable. Xtr is already standardized."""
    n, d = Xtr.shape
    feats = []
    norm_cols = []   # normalized columns, for vectorized near-duplicate rejection

    def add(name, depth, col, builder):
        col = np.asarray(col, float)
        if not np.all(np.isfinite(col)):
            return
        s = col.std()
        if s < 1e-9:
            return
        cs = (col - col.mean()) / s
        if norm_cols:
            M = np.asarray(norm_cols)
            if np.max(np.abs(M @ cs)) / n > 0.999:   # ~identical to an existing feature
                return
        norm_cols.append(cs)
        feats.append((name, depth, col, builder))

    # depth 1: raw vars + unary transforms
    for j in range(d):
        add(f"x{j}", 1, Xtr[:, j], (lambda X, j=j: X[:, j]))
        for uname, uf in UNARY.items():
            add(f"{uname}(x{j})", 1, uf(Xtr[:, j]), (lambda X, j=j, uf=uf: uf(X[:, j])))

    leaves = list(feats)

    # depth 2: leaf (op) leaf, plus squared leaves
    if depth_budget >= 2:
        maxpairs = 140
        cnt = 0
        for (na, _da, ca, ba), (nb, _db, cb, bb) in itertools.combinations(leaves, 2):
            if cnt >= maxpairs:
                break
            for opname, opf in BINARY.items():
                add(f"({na}{opname}{nb})", 2, opf(ca, cb),
                    (lambda X, ba=ba, bb=bb, opf=opf: opf(ba(X), bb(X))))
            cnt += 1
        for (na, _da, ca, ba) in leaves:
            add(f"({na}*{na})", 2, ca * ca, (lambda X, ba=ba: ba(X) * ba(X)))

    depth2 = [f for f in feats if f[1] == 2]

    # depth 3: depth-2 node (op) leaf  (bounded; tight caps -- high-d is incompressible anyway)
    if depth_budget >= 3:
        maxpairs3 = 70
        cnt = 0
        leaf_subset = leaves[:max(3, min(len(leaves), 2 * max_vars))]
        for (na, _da, ca, ba) in depth2[:35]:
            stop = False
            for (nb, _db, cb, bb) in leaf_subset:
                for opname, opf in (("+", BINARY["+"]), ("-", BINARY["-"]), ("*", BINARY["*"])):
                    add(f"({na}{opname}{nb})", 3, opf(ca, cb),
                        (lambda X, ba=ba, bb=bb, opf=opf: opf(ba(X), bb(X))))
                cnt += 1
                if cnt >= maxpairs3:
                    stop = True
                    break
            if stop:
                break
    return feats


def greedy_law(Xtr, ytr, feats, max_terms):
    """Forward-selection of up to max_terms features + intercept (least squares refit each step)."""
    n = Xtr.shape[0]
    ytr = np.asarray(ytr, float)
    cols = [np.asarray(f[2], float) for f in feats]
    norm_cols = []
    for c in cols:
        s = c.std()
        norm_cols.append((c - c.mean()) / (s if s > 1e-9 else 1.0))
    chosen = []
    residual = ytr - ytr.mean()
    available = list(range(len(feats)))
    for _ in range(max_terms):
        rn = np.linalg.norm(residual)
        if rn < 1e-9:
            break
        best_j, best_corr = None, 0.0
        for j in available:
            cn = np.linalg.norm(norm_cols[j])
            if cn < 1e-12:
                continue
            cc = abs(float(np.dot(norm_cols[j], residual)) / (cn * rn + 1e-12))
            if cc > best_corr:
                best_corr, best_j = cc, j
        if best_j is None or best_corr < 1e-4:
            break
        chosen.append(best_j)
        available.remove(best_j)
        A = np.column_stack([cols[j] for j in chosen] + [np.ones(n)])
        coef, *_ = np.linalg.lstsq(A, ytr, rcond=None)
        residual = ytr - A @ coef
    if not chosen:
        return None
    A = np.column_stack([cols[j] for j in chosen] + [np.ones(n)])
    coef, *_ = np.linalg.lstsq(A, ytr, rcond=None)
    return dict(builders=[feats[j][3] for j in chosen], coef=coef,
                depth=max(feats[j][1] for j in chosen), n_terms=len(chosen))


def predict_law(law, X):
    cols = [b(X) for b in law["builders"]]
    A = np.column_stack(cols + [np.ones(np.asarray(X).shape[0])])
    return A @ law["coef"]

# ----------------------------------------------------------------------------------------
# 2. Black-box ceiling -- KernelRidge (RBF), closed-form, alpha tuned by small inner CV.
#    (Replaces GaussianProcessRegressor, which stalled the previous attempt.)
# ----------------------------------------------------------------------------------------
BB_CAP = 150          # cap black-box training subsample for speed (recorded in caveats)
ALPHA_GRID = (1e-3, 1e-2, 1e-1, 1.0)


def _median_gamma(Xs):
    """RBF gamma from the median pairwise distance heuristic (on standardized X)."""
    m = Xs.shape[0]
    if m > 120:
        idx = np.random.RandomState(0).choice(m, 120, replace=False)
        Xsub = Xs[idx]
    else:
        Xsub = Xs
    D = euclidean_distances(Xsub, Xsub)
    med = np.median(D[D > 0]) if np.any(D > 0) else 1.0
    med = med if med > 1e-6 else 1.0
    return 1.0 / (2.0 * med ** 2)


def _cap_subsample(X, y, cap, seed):
    n = len(y)
    if n <= cap:
        return X, y
    idx = np.random.RandomState(seed).choice(n, cap, replace=False)
    return X[idx], y[idx]


def krr_fit_predict(Xtr, ytr, Xte, alpha=None, cap=BB_CAP, seed=0):
    """Standardize on train, RBF KernelRidge. If alpha is None, pick it by a quick 3-fold
    inner CV on the (capped) training data. Returns test predictions on Xte."""
    Xtr, ytr = _cap_subsample(np.asarray(Xtr, float), np.asarray(ytr, float), cap, seed)
    Xs, mu, sd = standardize_cols(Xtr)
    Xte_s = (np.asarray(Xte, float) - mu) / sd
    ymu = ytr.mean(); ystd = ytr.std() if ytr.std() > 1e-9 else 1.0
    yz = (ytr - ymu) / ystd
    gamma = _median_gamma(Xs)

    if alpha is None:
        n = len(yz)
        kf = KFold(n_splits=min(3, max(2, n // 8)), shuffle=True, random_state=0)
        best_alpha, best_r2 = ALPHA_GRID[0], -1e18
        for a in ALPHA_GRID:
            pr = np.zeros(n)
            for itr, ite in kf.split(Xs):
                kr = KernelRidge(alpha=a, kernel="rbf", gamma=gamma)
                kr.fit(Xs[itr], yz[itr])
                pr[ite] = kr.predict(Xs[ite])
            r2 = r2_det(yz, pr)
            if r2 > best_r2:
                best_r2, best_alpha = r2, a
        alpha = best_alpha

    kr = KernelRidge(alpha=alpha, kernel="rbf", gamma=gamma)
    kr.fit(Xs, yz)
    return kr.predict(Xte_s) * ystd + ymu, alpha


def krr_cv_r2(X, y, n_splits=4):
    """Outer k-fold test R2_det for the KernelRidge ceiling; alpha tuned per outer-train fold."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n = len(y)
    k = min(n_splits, max(2, n // 8))
    kf = KFold(n_splits=k, shuffle=True, random_state=0)
    preds = np.zeros(n)
    for tr, te in kf.split(X):
        preds[te], _ = krr_fit_predict(X[tr], y[tr], X[te], alpha=None, seed=0)
    return r2_det(y, preds)

# ----------------------------------------------------------------------------------------
# 3. Compressibility Lambda via NESTED k-fold  (t3: DEPTH budget, outer 4 / inner 3)
# ----------------------------------------------------------------------------------------


def lambda_depth_nestedcv(X, y, depth_grid=(1, 2, 3), max_terms=4,
                          max_vars=6, outer_splits=4, inner_splits=3):
    """Outer k-fold = honest test R2_det of the best law. Inner k-fold inside each outer-train
    fold selects the DEPTH budget maximizing inner-val R2_det (tie-break toward shallower);
    refit at that depth on the full outer-train fold, predict the outer-test fold. Variable
    pre-screen (train-only) keeps high-d tractable.
    Returns (R2_law_test, median_selected_depth, median_n_terms)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n = len(y)
    k_out = min(outer_splits, max(2, n // 8))
    kf_out = KFold(n_splits=k_out, shuffle=True, random_state=0)
    preds = np.zeros(n)
    chosen_depths, chosen_terms = [], []

    for tr, te in kf_out.split(X):
        Xtr_full, ytr = X[tr], y[tr]
        Xte = X[te]
        keep = screen_keep(Xtr_full, ytr, max_vars)
        Xtr = Xtr_full[:, keep]; Xte_k = Xte[:, keep]
        Xtr_s, mu, sd = standardize_cols(Xtr)
        Xte_s = (Xte_k - mu) / sd

        k_in = min(inner_splits, max(2, len(tr) // 8))
        kf_in = KFold(n_splits=k_in, shuffle=True, random_state=1)
        best_depth, best_inner = depth_grid[0], -1e18
        for depth in depth_grid:
            inner_pred = np.zeros(len(tr)); ok = True
            for itr, ite in kf_in.split(Xtr_s):
                feats = build_feature_bank(Xtr_s[itr], depth, max_vars)
                if not feats:
                    ok = False; break
                law = greedy_law(Xtr_s[itr], ytr[itr], feats, max_terms)
                inner_pred[ite] = (ytr[itr].mean() if law is None
                                   else predict_law(law, Xtr_s[ite]))
            if not ok:
                continue
            r2_in_adj = r2_det(ytr, inner_pred) - 0.01 * depth   # prefer compression on ties
            if r2_in_adj > best_inner:
                best_inner, best_depth = r2_in_adj, depth

        feats = build_feature_bank(Xtr_s, best_depth, max_vars)
        law = greedy_law(Xtr_s, ytr, feats, max_terms)
        if law is None:
            preds[te] = ytr.mean(); chosen_depths.append(best_depth); chosen_terms.append(0)
        else:
            preds[te] = predict_law(law, Xte_s)
            chosen_depths.append(law["depth"]); chosen_terms.append(law["n_terms"])

    return r2_det(y, preds), float(np.median(chosen_depths)), float(np.median(chosen_terms))

# ----------------------------------------------------------------------------------------
# 4. Selection fitness over the full design space (uses the GROUND-TRUTH map)
# ----------------------------------------------------------------------------------------


def selection_fitness(model_predict_fn, X_design, y_true_design):
    y_true = np.asarray(y_true_design, float)
    y_pred = np.asarray(model_predict_fn(X_design), float)
    true_best = y_true.max()
    denom = true_best if abs(true_best) > 1e-9 else 1.0
    pick = int(np.argmax(y_pred))
    top1_regret = float(max(0.0, (true_best - y_true[pick]) / denom))
    try:
        sp = stats.spearmanr(y_pred, y_true).correlation
        sp = 0.0 if not np.isfinite(sp) else float(sp)
    except Exception:
        sp = 0.0
    order = np.argsort(y_pred)[::-1][:5]
    top5_regret = float(max(0.0, (true_best - y_true[order].max()) / denom))
    return dict(top1_regret=top1_regret, top5_regret=top5_regret, spearman=sp)
