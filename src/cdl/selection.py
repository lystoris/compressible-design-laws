import numpy as np


def top1_regret(y_true, y_pred):
    """Regret of picking the argmax of y_pred, normalized by the true range.
    0 = picked design matches the true best; 1 = picked the worst design."""
    yt = np.asarray(y_true, float); yp = np.asarray(y_pred, float)
    denom = yt.max() - yt.min()
    if denom <= 0:
        return 0.0
    return float((yt.max() - yt[np.argmax(yp)]) / denom)


def topk_overlap(y_true, y_pred, k=100):
    """Count of indices in the predicted top-k (argsort of y_pred) that are
    also in the true top-k (argsort of y_true)."""
    yt = np.asarray(y_true, float); yp = np.asarray(y_pred, float)
    k = min(k, yt.shape[0])
    true_topk = set(np.argsort(-yt)[:k].tolist())
    pred_topk = set(np.argsort(-yp)[:k].tolist())
    return len(true_topk & pred_topk)
