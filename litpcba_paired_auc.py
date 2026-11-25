"""Is one method's AUC really above another's? A paired bootstrap on LIT-PCBA.

[RESULTS.md](RESULTS.md) 15 needed this. Docking reached AUC 0.622 on LIT-PCBA MAPK1
against molecular size's 0.577, which reads as a win. Every method here ranks the
*same* compounds, so the two AUC estimates are correlated and an unpaired error bar
(Hanley-McNeil, say) overstates the uncertainty of their *difference*. Resampling the
compounds and recomputing both AUCs on each resample does not.

Actives and inactives are resampled separately so every resample holds the real class
ratio; with 77 actives against 2,000 inactives, pooled resampling would let the active
count wander and inflate the spread.

Docking scores come from the pose cache written by the main benchmark, so this runs in
seconds rather than re-docking for 87 minutes. It therefore only works after
`litpcba_benchmark.py ... --dock` has populated that cache for the same receptor.

    python3 litpcba_paired_auc.py <litpcba_target_dir> [receptor] [max_inactives]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from dmta.docking import apply_scores, dock, engine_available
from dmta.hitfinding import roc_auc
from dmta.models import PropertyModel
from dmta.structure import prepare
from litpcba_benchmark import read_split

RESAMPLES = 4000


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    directory = Path(sys.argv[1])
    receptor = sys.argv[2] if len(sys.argv) > 2 else "4QTE"
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 2000

    if not engine_available():
        sys.exit("smina not found — this script reads the cache it writes")

    test = (read_split(directory, "active_V.smi", True, None)
            + read_split(directory, "inactive_V.smi", False, cap))
    train = (read_split(directory, "active_T.smi", True, None)
             + read_split(directory, "inactive_T.smi", False, cap))
    labels = np.array([c.truth for c in test], float)
    print(f"test {len(test)} compounds, {int(labels.sum())} actives "
          f"({labels.mean()*100:.2f}%)", flush=True)

    n = apply_scores(test, dock(test, prepare(receptor), verbose=False))
    print(f"docking scores available for {n}/{len(test)} "
          f"(from the cache; a long wait here means it is not populated)", flush=True)

    for c in train:
        c.measured = c.truth
    try:
        mu, _ = PropertyModel(n_estimators=300).fit(train).predict(test)
    finally:
        for c in train:
            c.measured = None

    columns = {
        # lower kcal/mol is a better fit, so negate for higher = better
        "docking": np.array([-c.docking_score if c.docking_score is not None else np.nan
                             for c in test], float),
        "size": np.array([c.desc["heavy_atoms"] for c in test], float),
        "forest": np.asarray(mu, float),
    }
    # One compound set for every method, or the comparison is not paired.
    usable = np.ones(len(test), bool)
    for v in columns.values():
        usable &= np.isfinite(v)
    print(f"compounds scored by every method: {int(usable.sum())}\n")
    y = labels[usable]
    columns = {k: v[usable] for k, v in columns.items()}
    for name, v in columns.items():
        print(f"  {name:<9} AUC {roc_auc(v, y):.4f}")

    rng = np.random.default_rng(0)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    pairs = [("docking", "size"), ("docking", "forest"), ("forest", "size")]
    diffs = {p: np.empty(RESAMPLES) for p in pairs}
    for b in range(RESAMPLES):
        idx = np.concatenate([rng.choice(pos, pos.size, replace=True),
                              rng.choice(neg, neg.size, replace=True)])
        yb = y[idx]
        auc = {k: roc_auc(v[idx], yb) for k, v in columns.items()}
        for a, c in pairs:
            diffs[(a, c)][b] = auc[a] - auc[c]

    print(f"\nstratified paired bootstrap, {RESAMPLES} resamples "
          f"({pos.size} actives and {neg.size} inactives resampled separately)")
    head = f"{'comparison':<20}{'mean dAUC':>11}{'95% CI':>22}{'P(d<=0)':>10}"
    print(head); print("-" * len(head))
    for (a, c), v in diffs.items():
        lo, hi = np.percentile(v, [2.5, 97.5])
        print(f"{a + ' - ' + c:<20}{v.mean():>+11.4f}"
              f"   [{lo:>+.4f}, {hi:>+.4f}]{(v <= 0).mean():>10.4f}")
        if lo <= 0 <= hi:
            print(f"{'':<20}  the interval contains zero: this gap is not evidence "
                  f"that {a} ranks better than {c}")


if __name__ == "__main__":
    main()
