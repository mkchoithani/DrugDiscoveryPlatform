"""Linear model on an RMT-denoised chemical correlation matrix.

Lee et al., PNAS 2019, 116, 3373 (DOI 10.1073/pnas.1810847116). The regime that
paper addresses is exactly ours: thousands of fingerprint features, tens to
hundreds of measured compounds. In that regime the sample correlation matrix is
mostly noise — with p features and n molecules, random data alone produces
correlations as large as the real ones. Random matrix theory says precisely how
large: the eigenvalues of a pure-noise correlation matrix fall inside the
Marchenko-Pastur bulk, whose upper edge is sigma^2 (1 + sqrt(p/n))^2. Eigenvalues
above that edge carry chemical signal; everything inside the bulk is undersampling.

So: keep the eigenvalues above the edge, replace the rest by their common average
(eigenvalue clipping, Bouchaud and Potters), and solve the linear problem against
that cleaned matrix instead of the raw one. The result is generalised least
squares with the noise directions flattened rather than inverted — which is what
makes plain least squares explode here.

Same fit / predict / apply / report surface as PropertyModel and ChempropModel,
so model_benchmark.py can put all three on identical splits.

Interpretable by construction: `top_features()` returns the signal eigenvectors
and the fitted weights, which are fingerprint bits and physchem descriptors, not
hidden units.
"""
from __future__ import annotations

import numpy as np
from sklearn.model_selection import KFold

from .chem import feature_vector
from .library import Compound
from .models import regression_metrics


def marchenko_pastur_edge(n: int, p: int, variance: float = 1.0) -> float:
    """Upper edge of the noise bulk for an n x p standardised data matrix."""
    if n <= 0 or p <= 0:
        raise ValueError("n and p must be positive")
    return float(variance * (1.0 + np.sqrt(p / n)) ** 2)


def clip_spectrum(eigenvalues: np.ndarray, n: int, p: int) -> tuple[np.ndarray, float, float]:
    """Split a correlation spectrum into signal and noise.

    `eigenvalues` holds the non-zero eigenvalues of the sample correlation matrix,
    descending; when p > n there are p - len(eigenvalues) further zeros, and they
    count towards the trace.

    The noise variance is 1, not fitted. For a *correlation* matrix the diagonal is
    all ones, so the trace is exactly p and the mean eigenvalue is exactly 1 — which
    is what the null variance means. Fitting it self-consistently (set the edge from
    a variance, re-estimate the variance from whatever falls below the edge, repeat)
    is the textbook recipe for a covariance matrix and it DIVERGES here: with p > n
    the p - n zero eigenvalues pull each estimate downwards, which lowers the edge,
    which admits more eigenvalues as signal, which lowers it further. It ran to
    edge = 0 with every eigenvalue called signal and a noise floor of 1e-8, giving
    the noise directions a weight of 1e8. That is where the first version's
    predictions went.

    Returns (signal eigenvalues, the flat noise level, the bulk edge).
    """
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    if n <= 1:
        raise ValueError("need at least two rows to estimate a correlation matrix")
    edge = marchenko_pastur_edge(n, p)
    keep = eigenvalues[eigenvalues > edge]
    remaining = p - len(keep)
    # The flat level absorbs everything below the edge, zeros included, so that the
    # cleaned matrix has the same trace as the raw one.
    noise_level = ((float(p) - float(keep.sum())) / remaining) if remaining > 0 else 1.0
    return keep, max(noise_level, 1e-6), edge


class RmtModel:
    """Generalised least squares against a spectrum-cleaned correlation matrix."""

    def __init__(self, mode: str = "project", shrinkage: float = 0.0, seed: int = 0,
                 max_components: int | None = None, clean: bool = True):
        # Two defensible readings of the paper, both measured in model_benchmark.py:
        #
        #   "clip"    invert the cleaned correlation matrix literally. This is
        #             generalised least squares: it WHITENS, so a signal direction
        #             with eigenvalue 90 is damped by 1/90 while a noise direction
        #             at the flat level keeps weight ~1. Right for a discriminant,
        #             wrong for ranking — it throws away the structure RMT just found.
        #   "project" solve inside the signal subspace only and discard the noise
        #             complement: principal component regression on the components
        #             RMT selected. This is what "clean the correlations" should mean
        #             when the aim is to predict.
        if mode not in ("clip", "project"):
            raise ValueError("mode must be 'clip' or 'project'")
        self.mode = mode
        # clean=False is the EXPERIMENTAL CONTROL, not a second method. It takes the
        # identical path -- same standardisation, same solve, same predictive variance
        # -- and only skips the Marchenko-Pastur clipping, keeping every non-zero
        # eigenvalue as if it were signal. That isolates the paper's contribution:
        # any difference between clean=True and clean=False is the cleaning and
        # nothing else.
        #
        # With p > n the raw correlation matrix is singular, so the control cannot
        # run at shrinkage 0 -- its smallest kept eigenvalues are ~1e-16 and the
        # weights explode. It therefore needs a ridge, which is exactly the ad-hoc
        # regularisation RMT is proposed as a principled replacement for. The
        # comparison is only fair if that ridge is chosen properly, so the benchmark
        # picks it by cross-validation on the training split and never by hand.
        self.clean = clean
        # shrinkage is a ridge added to the cleaned matrix; 0 means pure RMT,
        # and the benchmark reports whether it is needed at all.
        self.shrinkage = shrinkage
        self.seed = seed
        self.max_components = max_components
        self.mean_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None
        self.keep_: np.ndarray | None = None       # column mask of usable features
        self.weights_: np.ndarray | None = None
        self.intercept_ = 0.0
        self.n_signal = 0
        self.noise_level = 1.0
        self.edge = 0.0
        self.residual_var = 1.0
        self.calibration = 1.0
        self.cv_r2: float | None = None
        self.cv_mae: float | None = None
        self.n_train = 0
        self.n_features = 0
        self._eigvecs: np.ndarray | None = None    # signal directions, for inspection
        self._eigvals: np.ndarray | None = None
        # (signal directions, their eigenvalues, the flat noise level) — the three
        # pieces the cleaned inverse needs. Set by every path through _fit_matrix,
        # including the degenerate one, so _leverage can never find it missing.
        self._inv_cache: tuple = (np.zeros((1, 0)), np.zeros(0), 1.0)

    # -------------------------------------------------------------- internals
    def _standardise(self, X: np.ndarray) -> np.ndarray:
        return (X[:, self.keep_] - self.mean_) / self.scale_

    def _fit_matrix(self, X: np.ndarray, y: np.ndarray) -> None:
        """The whole method, on raw feature and value arrays."""
        n_raw, _ = X.shape
        # A feature that never varies in this training set carries no information
        # and would divide by zero on standardisation.
        spread = X.std(axis=0)
        self.keep_ = spread > 1e-9
        if not self.clean and self.shrinkage <= 0:
            raise ValueError(
                "the uncleaned control is singular when p > n: it needs a positive "
                "shrinkage, chosen by cross-validation on the training split")
        if not self.keep_.any():
            self.mean_ = np.zeros(1); self.scale_ = np.ones(1)
            self.keep_ = np.zeros(X.shape[1], dtype=bool); self.keep_[0] = True
            self.weights_ = np.zeros(1); self.intercept_ = float(y.mean())
            self.residual_var = float(y.var()) or 1.0
            self.n_features, self.n_signal = 1, 0
            self._inv_cache = (np.zeros((1, 0)), np.zeros(0), 1.0)
            return

        Xk = X[:, self.keep_]
        self.mean_ = Xk.mean(axis=0)
        self.scale_ = Xk.std(axis=0)
        Z = (Xk - self.mean_) / self.scale_
        n, p = Z.shape
        self.n_features = p
        self.intercept_ = float(y.mean())
        yc = y - self.intercept_

        # Eigenvalues of C = Z^T Z / n, from the SVD of Z: cheaper than forming
        # the p x p matrix when p >> n, and it gives the eigenvectors directly.
        _, s, Vt = np.linalg.svd(Z, full_matrices=False)
        eigvals = (s ** 2) / n
        if self.clean:
            signal, noise_level, edge = clip_spectrum(eigvals, n, p)
        else:
            # No clipping: every direction the data spans is treated as signal.
            # The raw correlation matrix has eigenvalue exactly 0 on the p - n
            # directions the data does not span, so the ridged inverse weights those
            # by 1/shrinkage -- NOT by 1/(1 + shrinkage), which is what a noise level
            # of 1 would give. Setting the level to 0 makes the floor exactly the
            # ridge, so this really is (C_raw + shrinkage I)^-1 and not a cleaned
            # matrix wearing the control's name.
            signal, edge = eigvals, 0.0
            noise_level = 0.0
            if self.shrinkage <= 0:
                raise ValueError(
                    "the uncleaned control is singular when p > n: it needs a "
                    "positive shrinkage, chosen by cross-validation on the "
                    "training split")
        if self.max_components is not None:
            signal = signal[: self.max_components]
        k = len(signal)
        self.n_signal, self.noise_level, self.edge = k, noise_level, edge

        V = Vt[:k].T                                     # p x k signal directions
        self._eigvecs, self._eigvals = V, signal

        # C_clean = V diag(signal) V^T + noise_level (I - V V^T), so its inverse is
        # V diag(1/signal) V^T + (1/noise_level)(I - V V^T) -- no p x p inverse needed.
        lam = signal + self.shrinkage
        floor = noise_level + self.shrinkage
        b = Z.T @ yc / n                                  # feature-value covariance
        if self.mode == "project":
            # Least squares restricted to the signal subspace. The scores Z V have
            # variance lam, so the coefficient on each component is (V^T b) / lam,
            # and the noise complement is dropped rather than inverted.
            self.weights_ = V @ ((V.T @ b) / lam) if k else np.zeros(p)
        elif k:
            proj = V.T @ b
            self.weights_ = V @ (proj / lam) + (b - V @ proj) / floor
        else:
            self.weights_ = b / floor                     # no signal: ridge-like

        resid = yc - Z @ self.weights_
        dof = max(n - k - 1, 1)
        self.residual_var = max(float((resid ** 2).sum()) / dof, 1e-6)
        # "project" has no noise complement, so its predictive variance must not
        # borrow one: an infinite floor sends that term to zero.
        self._inv_cache = (V, lam, np.inf if self.mode == "project" else floor)

    def _leverage(self, Z: np.ndarray) -> np.ndarray:
        """x^T C_clean^-1 x, the Bayesian predictive variance term."""
        V, lam, floor = self._inv_cache
        if V.size:
            proj = Z @ V
            quad = (proj ** 2 / lam).sum(axis=1)
            if np.isfinite(floor):
                quad = quad + ((Z ** 2).sum(axis=1) - (proj ** 2).sum(axis=1)) / floor
        else:
            quad = ((Z ** 2).sum(axis=1) / floor if np.isfinite(floor)
                    else np.zeros(len(Z)))
        return np.maximum(quad, 0.0)

    # -------------------------------------------------------------------- fit
    def fit(self, compounds: list[Compound], target: str = "measured") -> "RmtModel":
        X = np.array([feature_vector(c.smiles) for c in compounds], dtype=float)
        y = np.array([getattr(c, target) for c in compounds], dtype=float)
        self.n_train = len(y)
        if self.n_train >= 8:
            self._cross_validate(X, y)
        self._fit_matrix(X, y)
        return self

    def fit_matrix_only(self, X: np.ndarray, y: np.ndarray) -> "RmtModel":
        """Fit from feature and value arrays that are already built.

        `fit` re-featurises every compound from SMILES. A screening benchmark
        scores the same tens of thousands of molecules with several arms, so it
        featurises once and hands the matrix to each; this is that entry point.
        No cross-validation here -- the caller holds out its own split.
        """
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        self.n_train = len(y)
        self._fit_matrix(X, y)
        return self

    def _cross_validate(self, X: np.ndarray, y: np.ndarray, folds: int = 5) -> None:
        """Out-of-fold estimate, and the sigma calibration factor."""
        folds = min(folds, max(2, len(y) // 4))
        kf = KFold(n_splits=folds, shuffle=True, random_state=self.seed)
        preds, sigmas = np.zeros_like(y), np.ones_like(y)
        for train_idx, test_idx in kf.split(X):
            fold = RmtModel(mode=self.mode, shrinkage=self.shrinkage, seed=self.seed,
                            max_components=self.max_components, clean=self.clean)
            fold._fit_matrix(X[train_idx], y[train_idx])
            mu, sigma = fold._predict_matrix(X[test_idx])
            preds[test_idx], sigmas[test_idx] = mu, sigma
        resid = np.abs(preds - y)
        ss_tot = float(((y - y.mean()) ** 2).sum())
        self.cv_r2 = round(1 - float(((y - preds) ** 2).sum()) / ss_tot, 3) if ss_tot else None
        self.cv_mae = round(float(resid.mean()), 3)
        mean_sigma = float(sigmas.mean())
        self.calibration = (float(np.clip(resid.mean() / mean_sigma, 0.2, 5.0))
                            if mean_sigma > 1e-6 else 1.0)

    # ---------------------------------------------------------------- predict
    def _predict_matrix(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        Z = self._standardise(np.asarray(X, dtype=float))
        mu = Z @ self.weights_ + self.intercept_
        n = max(self.n_train, 1)
        var = self.residual_var * (1.0 + self._leverage(Z) / n)
        return mu, np.sqrt(var) * self.calibration

    def predict(self, compounds: list[Compound]) -> tuple[np.ndarray, np.ndarray]:
        X = np.array([feature_vector(c.smiles) for c in compounds], dtype=float)
        mu, sigma = self._predict_matrix(X)
        if len(mu) != len(compounds):        # never return a silently short answer
            raise RuntimeError(
                f"rmt returned {len(mu)} predictions for {len(compounds)} compounds")
        return mu.astype(float), np.maximum(sigma.astype(float), 1e-3)

    def apply(self, compounds: list[Compound]) -> None:
        if not compounds:
            return
        mu, sigma = self.predict(compounds)
        for c, m, s in zip(compounds, mu, sigma):
            c.predicted = round(float(m), 3)
            c.uncertainty = round(float(s), 3)

    # ------------------------------------------------------------ inspection
    def top_features(self, k: int = 10) -> list[tuple[int, float]]:
        """The feature indices this model leans on hardest, with their weights.

        Indices are into feature_vector's output: 0..2047 are Morgan bits, the
        last twelve are the scaled physchem block.
        """
        if self.weights_ is None:
            return []
        cols = np.flatnonzero(self.keep_)
        order = np.argsort(np.abs(self.weights_))[::-1][:k]
        return [(int(cols[i]), round(float(self.weights_[i]), 4)) for i in order]

    def report(self) -> dict:
        return {
            "model": f"rmt-linear-{self.mode}" + ("" if self.clean else "-uncleaned"),
            "n_train": self.n_train,
            "n_features": self.n_features, "n_signal_eigenvalues": self.n_signal,
            "mp_edge": round(self.edge, 3), "noise_level": round(self.noise_level, 4),
            "cv_r2": self.cv_r2, "cv_mae": self.cv_mae,
            "sigma_calibration": round(self.calibration, 3),
            "shrinkage": self.shrinkage,
        }


def ols_baseline_metrics(X: np.ndarray, y: np.ndarray, Xt: np.ndarray,
                         yt: np.ndarray) -> dict:
    """Plain least squares on the same features, for the contrast the paper draws."""
    Xd = np.hstack([X, np.ones((len(X), 1))])
    w, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    pred = np.hstack([Xt, np.ones((len(Xt), 1))]) @ w
    return regression_metrics([float(v) for v in pred], [float(v) for v in yt])
