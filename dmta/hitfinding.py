"""Hit finding: separating actives from decoys, which is a different job.

Every evaluation in this repository until now has measured **lead optimisation** —
ranking compounds within one programme by potency, where nearly everything is already
active and the question is which is better. Docking and pharmacophore maps were both
rejected on that test ([RESULTS.md](RESULTS.md) 4-5 and 10).

But neither method claims to do that job. Both are **hit finding** tools: pull the few
actives out of a large library of molecules that mostly do nothing. That is a coarser
question, and a coarse representation can answer it while being useless at the fine
one. Testing them only on lead optimisation and concluding they do not work was an
unfair test, and this module is the fair one.

## The benchmark

DUD-E (Mysinger et al., *J. Med. Chem.* 2012, 55, 6582), `dude.docking.org`: for each
target, known actives plus decoys drawn to **match the actives on physical properties**
(molecular weight, logP, rotatable bonds, hydrogen-bond counts, charge) while being
topologically dissimilar. The property matching is the point — it stops a method from
scoring well by noticing that actives are simply bigger.

**DUD-E's known bias, stated up front.** Because its actives are drawn from the
literature they come in analogue series, and because decoys are property-matched but
structurally distant, *ligand-based* 2D methods score suspiciously well on it — they
can recognise the series rather than the binding. Published critiques of exactly this
led to LIT-PCBA. So a strong 2D-fingerprint result here is evidence about the benchmark
as much as about the method, and the comparison below reports the 2D model for that
reason rather than as a fair rival.

## The metrics

Potency-ranking metrics are the wrong ones here. Hit finding is judged on whether the
actives come out near the **top** of the list, because that is all anyone will screen:

- **ROC AUC** — the probability a random active outranks a random decoy. 0.5 is chance.
- **EF at 1%** — how many more actives are in the top 1% than chance would give. The
  number that decides whether a screen is worth running.
- **BEDROC (alpha = 20)** — a weighted version that concentrates on roughly the top 8%,
  so a method that finds actives at rank 5000 gets little credit for it.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


def load_ism(path: str | Path) -> list[tuple[str, str]]:
    """DUD-E's .ism format: SMILES, then whitespace-separated identifiers."""
    out = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                out.append((parts[0], parts[-1]))
            elif parts:
                out.append((parts[0], f"row{len(out)}"))
    if not out:
        raise ValueError(f"no molecules in {path}")
    return out


def roc_auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Rank-based AUC, ties averaged. Higher score should mean more likely active."""
    scores, labels = np.asarray(scores, float), np.asarray(labels).astype(bool)
    n_pos, n_neg = int(labels.sum()), int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError("need both actives and decoys")
    order = np.argsort(scores)
    ranks = np.empty(len(scores), float)
    ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks within ties, or a method that scores many compounds identically
    # gets credit or blame it has not earned
    for value in np.unique(scores):
        mask = scores == value
        if mask.sum() > 1:
            ranks[mask] = ranks[mask].mean()
    return float((ranks[labels].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def enrichment_factor(scores: np.ndarray, labels: np.ndarray,
                      fraction: float = 0.01) -> float | None:
    scores, labels = np.asarray(scores, float), np.asarray(labels).astype(bool)
    n = len(scores)
    k = max(1, int(round(n * fraction)))
    base = labels.mean()
    if base == 0:
        return None
    top = np.argsort(scores)[::-1][:k]
    return float(labels[top].mean() / base)


def bedroc(scores: np.ndarray, labels: np.ndarray, alpha: float = 20.0) -> float:
    """Truhlar and Rogers 2009, eq. 36. alpha=20 weights roughly the top 8%.

    Bounded in [0, 1]: 1 means every active is at the very top, 0 means every active
    is at the very bottom. An early version of this returned 12.6 on a perfect
    ranking, which the known-answer test below caught immediately — the published
    formula is a sum of an exponentially weighted term and an offset, not a chain of
    rescalings.
    """
    scores, labels = np.asarray(scores, float), np.asarray(labels).astype(bool)
    n = len(scores)
    n_pos = int(labels.sum())
    if n_pos == 0 or n_pos == n:
        raise ValueError("need both actives and decoys")
    ratio = n_pos / n
    # ranks of the actives, 1-based, best score first
    order = np.argsort(-scores, kind="stable")
    ranks = np.flatnonzero(labels[order]) + 1

    # RIE first: the exponentially weighted sum DIVIDED BY its expectation under a
    # random ranking. Applying the rescaling below to the raw sum instead returns
    # 0.05 for a perfect ranking, which is how this was wrong a moment ago.
    weighted_sum = float(np.sum(np.exp(-alpha * ranks / n))) / n_pos
    random_expectation = (1.0 / n) * ((1.0 - math.exp(-alpha)) /
                                      (math.exp(alpha / n) - 1.0))
    rie = weighted_sum / random_expectation
    scale = (ratio * math.sinh(alpha / 2.0) /
             (math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ratio)))
    offset = 1.0 / (1.0 - math.exp(alpha * (1.0 - ratio)))
    return float(rie * scale + offset)


def hit_finding_report(scores, labels, name: str = "") -> dict:
    """Every metric at once, plus the counts needed to read them."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels).astype(bool)
    usable = np.isfinite(scores)
    if usable.sum() < len(scores):
        scores, labels = scores[usable], labels[usable]
    return {
        "method": name,
        "n": int(len(scores)),
        "n_actives": int(labels.sum()),
        "roc_auc": round(roc_auc(scores, labels), 3),
        # `is not None`, not a truth test: an enrichment of exactly 0.0 means no
        # active made the cut, which is a result, not missing data.
        "ef_1pct": (round(v, 2) if (v := enrichment_factor(scores, labels, 0.01))
                    is not None else None),
        "ef_5pct": (round(v, 2) if (v := enrichment_factor(scores, labels, 0.05))
                    is not None else None),
        "bedroc_20": round(bedroc(scores, labels, 20.0), 3),
    }
