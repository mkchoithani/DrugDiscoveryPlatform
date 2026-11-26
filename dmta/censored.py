"""Use the compounds that were too weak to measure.

A potency reported as "> 10 uM" is not missing data. It is a measurement: this
compound is definitely bad. Dropping it throws away the most confident statement the
assay made, and it biases the training set towards actives — the model never sees
what failure looks like.

This is not a corner case in real data. In the Open Source Malaria file, **51 of 123
compounds** are reported as "> X". Until now `dmta/ingest.py` dropped every one.

## The direction trap

Censoring direction flips under the p-scale conversion, and getting it backwards
inverts the dataset:

| as reported | means | on the p-scale (pIC50) |
|---|---|---|
| `IC50 > 10 uM` | weaker than 10 uM | pIC50 **<** 5.0 — left-censored |
| `IC50 < 1 nM` | more potent than 1 nM | pIC50 **>** 9.0 — right-censored |
| `pIC50 > 8` | already log, more potent | pIC50 **>** 8 — right-censored, no flip |

So the flip happens exactly when the reported units are a concentration, because
pIC50 = -log10(molar) is monotonically *decreasing* in concentration. When the column
is already on a log scale, the inequality stands as written.

## The model

Tobit-style expectation-maximisation around whatever regressor is already in use, so
the incumbent random forest keeps its job:

1. Fit on the fully observed compounds alone.
2. Predict the censored ones. A left-censored compound known to be below `b` gets the
   truncated expectation of its predictive distribution restricted to (-inf, b] —
   which is the model's own guess when that guess already respects the bound, and
   pulled below `b` when it does not.
3. Refit on observed plus imputed values, and repeat.

Each censored compound therefore contributes what is actually known about it — "no
higher than this" — rather than either nothing at all or a fictitious point value at
the bound.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .library import Compound
from .models import PropertyModel


@dataclass
class Censoring:
    """One censored observation, already on the p-scale."""
    bound: float
    direction: str        # "left" = true value is BELOW bound; "right" = above

    def __post_init__(self):
        if self.direction not in ("left", "right"):
            raise ValueError("direction must be 'left' or 'right'")


def censoring_on_p_scale(operator: str, bound_on_p_scale: float,
                         reported_in_log_units: bool) -> Censoring:
    """Translate a reported inequality into censoring on the p-scale.

    `operator` is as the source file wrote it (">", ">=", "<", "<="), and
    `bound_on_p_scale` is that row's numeric bound already converted. The flip is
    applied here, in one place, so no caller has to remember which way it goes.
    """
    if operator not in (">", ">=", "<", "<="):
        raise ValueError(f"not a censoring operator: {operator!r}")
    reported_greater = operator.startswith(">")
    # On a concentration scale, "greater" means weaker, and weaker is a LOWER
    # p-value. On a log scale the inequality already points the right way.
    more_potent = reported_greater if reported_in_log_units else not reported_greater
    return Censoring(bound=bound_on_p_scale,
                     direction="right" if more_potent else "left")


def _truncated_expectation(mu: np.ndarray, sigma: np.ndarray, bound: np.ndarray,
                           left: np.ndarray) -> np.ndarray:
    """E[X | X <= b] for left-censored rows, E[X | X >= b] for right-censored ones.

    Normal truncation: the mean shifts by sigma * phi(z) / Phi(z). When the model
    already predicts inside the allowed region the correction is small; when it
    predicts outside, this drags the value back across the bound, which is the whole
    point — a compound known to be weak must not be imputed as potent.
    """
    sigma = np.maximum(sigma, 1e-6)
    z = (bound - mu) / sigma
    # phi/Phi, guarding the far tail where Phi underflows to zero
    pdf = np.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)
    cdf = 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))
    lower = mu - sigma * pdf / np.maximum(cdf, 1e-8)          # E[X | X <= b]
    upper = mu + sigma * pdf / np.maximum(1.0 - cdf, 1e-8)    # E[X | X >= b]
    out = np.where(left, lower, upper)
    # Never let the imputation cross the bound: that is the one thing we know.
    out = np.where(left, np.minimum(out, bound), np.maximum(out, bound))
    return out


class CensoredModel:
    """Tobit-style EM around any model with fit / predict / apply."""

    def __init__(self, base=None, iterations: int = 5, factory=None):
        self.factory = factory or (lambda: PropertyModel(n_estimators=300))
        self.model = base or self.factory()
        self.iterations = iterations
        self.n_observed = 0
        self.n_censored = 0
        self.imputations_: np.ndarray | None = None

    def fit(self, observed: list[Compound], censored: list[Compound] | None = None,
            censoring: list[Censoring] | None = None,
            target: str = "measured") -> "CensoredModel":
        censored = censored or []
        censoring = censoring or []
        if len(censored) != len(censoring):
            raise ValueError(f"{len(censored)} censored compounds but "
                             f"{len(censoring)} censoring records")
        self.n_observed, self.n_censored = len(observed), len(censored)
        if not observed:
            raise ValueError("need at least some fully observed measurements")

        self.model = self.factory().fit(observed, target=target)
        if not censored:
            return self

        bound = np.array([c.bound for c in censoring], dtype=float)
        left = np.array([c.direction == "left" for c in censoring])
        # Start each censored compound at its bound, then let the model pull it into
        # the allowed region. Starting AT the bound is the naive practice this is
        # meant to improve on, so the first iteration is deliberately that baseline.
        imputed = bound.copy()

        for _ in range(self.iterations):
            previous = imputed.copy()
            saved = [getattr(c, target, None) for c in censored]
            for c, v in zip(censored, imputed):
                setattr(c, target, float(v))
            try:
                self.model = self.factory().fit(observed + censored, target=target)
                mu, sigma = self.model.predict(censored)
            finally:
                for c, v in zip(censored, saved):
                    setattr(c, target, v)
            imputed = _truncated_expectation(np.asarray(mu), np.asarray(sigma),
                                             bound, left)
            if np.max(np.abs(imputed - previous)) < 1e-3:
                break

        # Final fit with the settled imputations in place.
        saved = [getattr(c, target, None) for c in censored]
        for c, v in zip(censored, imputed):
            setattr(c, target, float(v))
        try:
            self.model = self.factory().fit(observed + censored, target=target)
        finally:
            for c, v in zip(censored, saved):
                setattr(c, target, v)
        self.imputations_ = imputed
        return self

    def predict(self, compounds: list[Compound]):
        return self.model.predict(compounds)

    def apply(self, compounds: list[Compound]) -> None:
        self.model.apply(compounds)

    def report(self) -> dict:
        base = self.model.report() if hasattr(self.model, "report") else {}
        return {**base, "model": "censored-em", "n_observed": self.n_observed,
                "n_censored": self.n_censored, "em_iterations": self.iterations,
                "mean_imputed": (round(float(np.mean(self.imputations_)), 3)
                                 if self.imputations_ is not None else None)}
