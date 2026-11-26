"""Property models with calibrated uncertainty.

Random-forest ensemble over Morgan + physchem features. Uncertainty is the spread
across trees, scaled by a calibration factor fitted on out-of-bag residuals — an
uncalibrated sigma makes active learning chase noise.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

from .chem import feature_vector
from .library import Compound


class PropertyModel:
    def __init__(self, n_estimators: int = 300, seed: int = 0, min_samples_leaf: int = 1):
        self.rf = RandomForestRegressor(
            n_estimators=n_estimators, min_samples_leaf=min_samples_leaf,
            max_features="sqrt", n_jobs=-1, random_state=seed, oob_score=False,
        )
        self.calibration = 1.0
        self.cv_r2: float | None = None
        self.cv_mae: float | None = None
        self.n_train = 0

    # ---------------------------------------------------------------- fit
    def fit(self, compounds: list[Compound], target: str = "measured") -> "PropertyModel":
        X = np.array([feature_vector(c.smiles) for c in compounds])
        y = np.array([getattr(c, target) for c in compounds], dtype=float)
        self.n_train = len(y)
        if self.n_train >= 8:
            self._cross_validate(X, y)
        self.rf.fit(X, y)
        return self

    def _cross_validate(self, X, y, folds: int = 5):
        """Honest internal estimate + the sigma calibration factor."""
        folds = min(folds, max(2, len(y) // 4))
        kf = KFold(n_splits=folds, shuffle=True, random_state=0)
        preds = np.zeros_like(y)
        sigmas = np.ones_like(y)
        for train_idx, test_idx in kf.split(X):
            m = RandomForestRegressor(
                n_estimators=120, max_features="sqrt", n_jobs=-1, random_state=0
            ).fit(X[train_idx], y[train_idx])
            per_tree = np.stack([t.predict(X[test_idx]) for t in m.estimators_])
            preds[test_idx] = per_tree.mean(axis=0)
            sigmas[test_idx] = per_tree.std(axis=0)
        resid = np.abs(preds - y)
        ss_tot = float(((y - y.mean()) ** 2).sum())
        self.cv_r2 = round(1 - float(((y - preds) ** 2).sum()) / ss_tot, 3) if ss_tot else None
        self.cv_mae = round(float(resid.mean()), 3)
        # scale predicted sigma so its mean matches observed mean absolute error
        mean_sigma = float(sigmas.mean())
        self.calibration = round(float(resid.mean() / mean_sigma), 3) if mean_sigma > 1e-6 else 1.0
        self.calibration = float(np.clip(self.calibration, 0.2, 5.0))

    # ------------------------------------------------------------ predict
    def predict(self, compounds: list[Compound]) -> tuple[np.ndarray, np.ndarray]:
        X = np.array([feature_vector(c.smiles) for c in compounds])
        per_tree = np.stack([t.predict(X) for t in self.rf.estimators_])
        mu = per_tree.mean(axis=0)
        sigma = per_tree.std(axis=0) * self.calibration
        return mu, np.maximum(sigma, 1e-3)

    def apply(self, compounds: list[Compound]) -> None:
        if not compounds:
            return
        mu, sigma = self.predict(compounds)
        for c, m, s in zip(compounds, mu, sigma):
            c.predicted = round(float(m), 3)
            c.uncertainty = round(float(s), 3)

    def report(self) -> dict:
        return {
            "n_train": self.n_train, "cv_r2": self.cv_r2, "cv_mae": self.cv_mae,
            "sigma_calibration": self.calibration, "n_estimators": self.rf.n_estimators,
        }


def regression_metrics(pred: list[float], obs: list[float]) -> dict:
    if len(pred) < 2:
        return {"n": len(pred), "r2": None, "mae": None, "rmse": None, "spearman": None}
    p, o = np.array(pred, dtype=float), np.array(obs, dtype=float)
    ss_tot = float(((o - o.mean()) ** 2).sum())
    r2 = round(1 - float(((o - p) ** 2).sum()) / ss_tot, 3) if ss_tot > 1e-9 else None
    # A constant prediction has no ranking, so the correlation is undefined
    # rather than zero. With mid-ranks its rank vector really is constant, so
    # this guard now fires instead of being bypassed by arbitrary tie order.
    rank_p, rank_o = _rank(p), _rank(o)
    sd = rank_p.std() * rank_o.std()
    spearman = round(float(((rank_p - rank_p.mean()) * (rank_o - rank_o.mean())).mean() / sd), 3) if sd > 1e-9 else None
    return {
        "n": len(p), "r2": r2,
        "mae": round(float(np.abs(p - o).mean()), 3),
        "rmse": round(float(np.sqrt(((p - o) ** 2).mean())), 3),
        "spearman": spearman,
    }


def _rank(a: np.ndarray) -> np.ndarray:
    """Mid-ranks, which is what Spearman is defined on.

    This used to assign distinct ranks in argsort order, so tied values were
    ordered arbitrarily. For ordinary predictions that is a rounding-level
    error, but it makes a *constant* predictor look like a ranking: its
    arbitrary ranks have real variance, the "is it constant" guard below never
    fires, and the correlation comes back as a small random number instead of
    undefined. Any benchmark carrying a constant-predictor floor would have
    scored that floor as noise around zero rather than as no ranking at all.
    """
    order = a.argsort(kind="mergesort")
    r = np.empty(len(a), dtype=float)
    r[order] = np.arange(len(a), dtype=float)
    # average the ranks within each run of equal values
    srt = a[order]
    i = 0
    while i < len(srt):
        j = i + 1
        while j < len(srt) and srt[j] == srt[i]:
            j += 1
        if j - i > 1:
            r[order[i:j]] = (i + j - 1) / 2.0
        i = j
    return r
