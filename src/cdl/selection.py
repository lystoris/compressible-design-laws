import numpy as np


def top1_regret(y_true, y_pred):
    """Regret of picking the argmax of y_pred, normalized by the true best
    value (round-04 definition: NOT range-normalized).
    0 = picked design matches the true best. Ties in y_pred break to the
    first argmax (deterministic, via np.argmax)."""
    yt = np.asarray(y_true, float); yp = np.asarray(y_pred, float)
    pick = int(np.argmax(yp))
    true_best = float(yt.max())
    return float((true_best - yt[pick]) / true_best)


def topk_overlap(y_true, y_pred, k=100):
    """Count of indices in the predicted top-k (by y_pred) that are also in
    the true top-k (by y_true)."""
    yt = np.asarray(y_true, float); yp = np.asarray(y_pred, float)
    k = min(k, yt.shape[0])
    true_topk = set(np.argpartition(-yt, k - 1)[:k].tolist()) if k > 0 else set()
    pred_topk = set(np.argpartition(-yp, k - 1)[:k].tolist()) if k > 0 else set()
    return len(true_topk & pred_topk)
