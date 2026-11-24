"""Batch selection: acquisition functions plus the constraints a real lab imposes.

A batch is not the top-N list. It must be diverse (no point making ten near-identical
analogues), it should reuse intermediates (cheaper and faster to make), and it has to
fit a reagent budget.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from .chem import bulk_tanimoto, fingerprint
from .library import Compound

ACQUISITIONS = ("ucb", "greedy", "expected_improvement", "thompson", "random", "uncertainty")


def acquisition_values(pool: list[Compound], kind: str, beta: float, best_so_far: float,
                       rng: np.random.Generator) -> np.ndarray:
    """All acquisitions return a value on the design-score scale (0–1).

    Optimism is applied to the *uncertain* quantity — predicted potency — and then
    passed back through the same multi-parameter profile. A compound the model is
    unsure about is worth making only if an optimistic potency would still clear the
    rest of the profile; this is what stops exploration from buying unmakeable bricks.
    """
    from .scoring import score as mpo_score

    score = np.array([c.score or 0.0 for c in pool])
    mu = np.array([c.predicted if c.predicted is not None else 0.0 for c in pool])
    sigma = np.array([c.uncertainty if c.uncertainty is not None else 1.0 for c in pool])

    if kind == "greedy":
        return score
    if kind == "random":
        return rng.random(len(pool))
    if kind == "uncertainty":
        return sigma / max(1e-6, sigma.max())

    if kind == "thompson":
        optimistic = rng.normal(mu, sigma)
    elif kind == "expected_improvement":
        z = (mu - best_so_far) / np.maximum(sigma, 1e-6)
        ei = (mu - best_so_far) * _norm_cdf(z) + sigma * _norm_pdf(z)
        optimistic = mu + ei           # EI expressed as a potency uplift
    else:                               # ucb
        optimistic = mu + beta * sigma

    return _rescore_with_potency(pool, optimistic, mpo_score)


def _rescore_with_potency(pool: list[Compound], potency: np.ndarray, mpo_score) -> np.ndarray:
    """Re-run the MPO profile with a substituted potency, then restore state."""
    out = np.empty(len(pool))
    for i, c in enumerate(pool):
        real_pred, real_score, real_obj = c.predicted, c.score, c.objectives
        c.predicted = float(potency[i])
        out[i] = mpo_score(c)
        c.predicted, c.score, c.objectives = real_pred, real_score, real_obj
    return out


def _norm_cdf(z):
    return 0.5 * (1 + np.vectorize(math.erf)(z / math.sqrt(2)))


def _norm_pdf(z):
    return np.exp(-0.5 * z ** 2) / math.sqrt(2 * math.pi)


def select_batch(
    library: list[Compound],
    size: int = 16,
    acquisition: str = "ucb",
    beta: float = 0.25,
    diversity: float = 0.15,
    intermediate_bonus: float = 0.03,
    budget: float | None = None,
    max_per_scaffold: int | None = None,
    seed: int = 0,
) -> list[Compound]:
    """Greedy sequential selection with a diversity penalty and lab constraints."""
    rng = np.random.default_rng(seed)
    pool = [c for c in library if c.measured is None]
    if not pool:
        return []
    measured = [c.measured for c in library if c.measured is not None]
    best_so_far = max(measured) if measured else 0.0

    base = acquisition_values(pool, acquisition, beta, best_so_far, rng)
    fps = np.stack([fingerprint(c.smiles) for c in pool])

    picked: list[Compound] = []
    picked_idx: list[int] = []
    max_sim = np.zeros(len(pool))
    chosen_intermediates: set[str] = set()
    per_scaffold: dict[str, int] = defaultdict(int)
    spend = 0.0

    for _ in range(size):
        bonus = np.array([
            intermediate_bonus if c.intermediate_key in chosen_intermediates else 0.0
            for c in pool
        ])
        value = base + bonus - diversity * max_sim
        value[picked_idx] = -np.inf
        for i, c in enumerate(pool):
            if max_per_scaffold is not None and per_scaffold[c.scaffold_id] >= max_per_scaffold:
                value[i] = -np.inf
            elif budget is not None and spend + c.cost > budget:
                value[i] = -np.inf
        if not np.isfinite(value).any():
            break

        i = int(np.argmax(value))
        c = pool[i]
        picked.append(c)
        picked_idx.append(i)
        chosen_intermediates.add(c.intermediate_key)
        per_scaffold[c.scaffold_id] += 1
        spend += c.cost
        max_sim = np.maximum(max_sim, bulk_tanimoto(fps[i], fps))

    return picked


def diversity_report(batch: list[Compound]) -> dict:
    if len(batch) < 2:
        return {"mean_pairwise_tanimoto": None, "n_scaffolds": len(batch), "n_intermediates": len(batch)}
    fps = np.stack([fingerprint(c.smiles) for c in batch])
    sims = []
    for i in range(len(batch)):
        s = bulk_tanimoto(fps[i], fps)
        sims.extend(s[j] for j in range(len(batch)) if j != i)
    return {
        "mean_pairwise_tanimoto": round(float(np.mean(sims)), 3),
        "max_pairwise_tanimoto": round(float(np.max(sims)), 3),
        "n_scaffolds": len({c.scaffold_id for c in batch}),
        "n_intermediates": len({c.intermediate_key for c in batch}),
    }
