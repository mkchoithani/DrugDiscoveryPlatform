"""Property models with calibrated uncertainty.

Random-forest ensemble over Morgan + physchem features. Uncertainty is the spread
across trees, scaled by a calibration factor fitted on out-of-bag residuals — an
uncalibrated sigma makes active learning chase noise.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

from .chem import feature_vector
from .library import Compound

N_PHYS = 12          # feature_vector = fingerprint ++ 12 physchem descriptors


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


class DescriptorRidge:
    """Ridge on the twelve physchem descriptors alone. No fingerprint.

    Section 18 measured this against the forest on aqueous solubility and the
    forest did not win: the ranking was a tie (d-rho -0.001, CI [-0.033, +0.032])
    and the absolute error was worse (MAE 1.053 against 0.731, R2 0.612 against
    0.821). The mechanism is not a tuning failure. Solubility spans twelve log
    units, a tree cannot predict outside the range of its training leaves, and so
    held-out scaffolds in the tails get pulled toward the mean. A linear model
    extrapolates, which is exactly what that endpoint needs.

    Sigma is the spread across ridges fitted to bootstrap resamples, scaled by the
    same out-of-fold calibration the forest uses. A linear model has no ensemble
    to disagree with itself, so without the resampling there would be no
    uncertainty at all -- and every acquisition function except greedy spends part
    of its batch on sigma.
    """

    def __init__(self, n_boot: int = 60, seed: int = 0,
                 alphas: tuple[float, ...] = tuple(np.logspace(-3, 3, 13))):
        self.n_boot, self.seed, self.alphas = n_boot, seed, list(alphas)
        self.ridge: RidgeCV | None = None
        self.members: list[RidgeCV] = []
        self.calibration = 1.0
        self.cv_r2: float | None = None
        self.cv_mae: float | None = None
        self.n_train = 0

    @staticmethod
    def _phys(compounds: list[Compound]) -> np.ndarray:
        return np.array([feature_vector(c.smiles)[-N_PHYS:] for c in compounds], dtype=float)

    def fit(self, compounds: list[Compound], target: str = "measured") -> "DescriptorRidge":
        X = self._phys(compounds)
        y = np.array([getattr(c, target) for c in compounds], dtype=float)
        self.n_train = len(y)
        if self.n_train >= 8:
            self._cross_validate(X, y)
        self.ridge = RidgeCV(alphas=self.alphas).fit(X, y)
        rng = np.random.default_rng(self.seed)
        self.members = []
        for _ in range(self.n_boot):
            idx = rng.integers(0, len(y), len(y))
            if len(np.unique(y[idx])) < 2:        # a degenerate resample fits nothing
                continue
            self.members.append(RidgeCV(alphas=self.alphas).fit(X[idx], y[idx]))
        return self

    def _cross_validate(self, X, y, folds: int = 5):
        folds = min(folds, max(2, len(y) // 4))
        kf = KFold(n_splits=folds, shuffle=True, random_state=0)
        preds = np.zeros_like(y)
        for train_idx, test_idx in kf.split(X):
            m = RidgeCV(alphas=self.alphas).fit(X[train_idx], y[train_idx])
            preds[test_idx] = m.predict(X[test_idx])
        resid = np.abs(preds - y)
        ss_tot = float(((y - y.mean()) ** 2).sum())
        self.cv_r2 = round(1 - float(((y - preds) ** 2).sum()) / ss_tot, 3) if ss_tot else None
        self.cv_mae = round(float(resid.mean()), 3)
        self._mean_resid = float(resid.mean())

    def predict(self, compounds: list[Compound]) -> tuple[np.ndarray, np.ndarray]:
        X = self._phys(compounds)
        mu = self.ridge.predict(X)
        if self.members:
            spread = np.stack([m.predict(X) for m in self.members]).std(axis=0)
        else:
            spread = np.zeros(len(X))
        mean_spread = float(spread.mean())
        if getattr(self, "_mean_resid", None) and mean_spread > 1e-6:
            self.calibration = float(np.clip(self._mean_resid / mean_spread, 0.2, 50.0))
        return mu, np.maximum(spread * self.calibration, 1e-3)

    def apply(self, compounds: list[Compound]) -> None:
        if not compounds:
            return
        mu, sigma = self.predict(compounds)
        for c, m, sg in zip(compounds, mu, sigma):
            c.predicted = round(float(m), 3)
            c.uncertainty = round(float(sg), 3)

    def report(self) -> dict:
        return {
            "kind": "descriptor_ridge", "n_train": self.n_train, "cv_r2": self.cv_r2,
            "cv_mae": self.cv_mae, "sigma_calibration": round(self.calibration, 3),
            "n_bootstrap": len(self.members),
            "alpha": round(float(self.ridge.alpha_), 4) if self.ridge is not None else None,
        }


def _grouped_oof_mae(model_factory, compounds: list[Compound], target: str,
                     folds: int = 5) -> float | None:
    """Out-of-fold MAE with whole scaffolds held out of each fold.

    Random folds are the wrong instrument for this choice and section 18 proved
    it the hard way. On aqueous solubility, random cross-validation says the
    forest is better (MAE 0.624 against the ridge's 0.729) and a selector reading
    that picks the forest every time. Hold scaffolds out instead and the order
    reverses (1.053 against 0.731), which is also what the held-out scaffold
    split independently reports. Random folds leak: near-identical analogues land
    on both sides, the forest interpolates between them, and the score describes
    a regime no production model is ever in.

    This is the same argument section 1 makes about reporting the random split as
    if it were generalisation. It applies with equal force to choosing a model.
    """
    groups = [c.scaffold_id for c in compounds]
    by_group: dict[str, list[int]] = {}
    for i, g in enumerate(groups):
        by_group.setdefault(g, []).append(i)
    if len(by_group) < 2:
        return None
    # greedy balanced assignment: largest scaffolds first into the emptiest fold
    order = sorted(by_group.values(), key=len, reverse=True)
    buckets: list[list[int]] = [[] for _ in range(min(folds, len(by_group)))]
    for idx in order:
        buckets.sort(key=len)
        buckets[0].extend(idx)

    y = np.array([getattr(c, target) for c in compounds], dtype=float)
    errs: list[float] = []
    for held in buckets:
        if not held or len(held) == len(compounds):
            continue
        mask = np.ones(len(compounds), dtype=bool)
        mask[held] = False
        train = [c for c, keep in zip(compounds, mask) if keep]
        test = [compounds[i] for i in held]
        if len(train) < 8 or len({getattr(c, target) for c in train}) < 2:
            continue
        mu, _ = model_factory().fit(train, target=target).predict(test)
        errs.extend(np.abs(mu - y[held]).tolist())
    return round(float(np.mean(errs)), 3) if errs else None


def choose_model(compounds: list[Compound], target: str = "measured",
                 folds: int = 5, seed: int = 0) -> tuple[object, dict]:
    """Fit both and keep the one that measures better with scaffolds held out.

    Section 18 produced opposite answers on two ADME endpoints -- the forest wins
    logD outright and loses solubility to twelve descriptors -- so "which model"
    is not a property of this platform, it is a property of the endpoint. Writing
    `if endpoint == "solubility"` would be a guess dressed as a decision, and the
    next endpoint would arrive with no rule at all.

    Two things this gets right only because getting them wrong was measured first:

    1. The comparison is on MAE, not rho. Solubility was a *tie* on ranking and a
       loss on error, so a selector reading rho keeps the forest and is wrong in
       exactly the case this exists for.
    2. The folds hold out whole scaffolds. The first version used each model's own
       random-fold CV and picked the forest 5 times out of 5 on solubility -- the
       honest answer to the question it was asked, and the wrong question.

    Held-out test data is never touched: the folds are cut inside the training
    split only.
    """
    f_mae = _grouped_oof_mae(lambda: PropertyModel(seed=seed), compounds, target, folds)
    r_mae = _grouped_oof_mae(lambda: DescriptorRidge(seed=seed), compounds, target, folds)
    if f_mae is None or r_mae is None:
        forest = PropertyModel(seed=seed).fit(compounds, target=target)
        return forest, {"chosen": "forest", "scoring": "scaffold-grouped out-of-fold MAE",
                        "reason": "too few scaffolds to compare — defaulted to the forest",
                        "forest_mae": f_mae, "ridge_mae": r_mae}
    if r_mae < f_mae:
        model, chosen = DescriptorRidge(seed=seed).fit(compounds, target=target), "ridge"
    else:
        model, chosen = PropertyModel(seed=seed).fit(compounds, target=target), "forest"
    return model, {
        "chosen": chosen,
        "scoring": "scaffold-grouped out-of-fold MAE",
        "forest_mae": f_mae,
        "ridge_mae": r_mae,
        "margin": round(abs(f_mae - r_mae), 3),
        "reason": f"{chosen} has the lower MAE with whole scaffolds held out",
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
