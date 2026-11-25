"""Chemprop in the regime its own paper validated: a learning curve, not one number.

[RESULTS.md](RESULTS.md) 7 ran Chemprop against the random forest on individual
programmes of 175-884 compounds and reported that it lost. That is the small-data
regime, and it is the regime where fingerprints plus a forest are *known* to be
competitive — Yang et al. and the Chemprop papers since report the D-MPNN's advantage
appearing on datasets an order of magnitude larger. So section 7 tested the model
outside the conditions its authors claimed, and then used the word "rejected".

This is the test we owed the paper. One target's entire cleaned set, one fixed scaffold
split, and the *training set* subsampled upwards. If Chemprop's advantage is a data-size
effect, the paired gap shrinks and crosses as n grows, and the crossover point is the
number that should decide which model a programme uses. If the gap is flat, size was
never the explanation and section 7's conclusion holds for any n we will ever have.

The test set is held fixed across every training size, so the curves are comparable:
only the amount of training data changes.

    python3 chemprop_scale.py [target] [standard] [sizes]

    python3 chemprop_scale.py CHEMBL203 IC50 250,500,1000,2000,4000,0
                                                            # 0 means "all of it"
"""
from __future__ import annotations

import sys
import time

import numpy as np

from dmta.chembl import load_compounds
from dmta.models import PropertyModel, regression_metrics
from dmta.validate import scaffold_split

SIZES = [250, 500, 1000, 2000, 4000, 0]       # 0 = every training compound


def fit_and_score(model, train, test) -> tuple[dict, float]:
    for c in train:
        c.measured = c.truth
    t0 = time.time()
    try:
        model.fit(train)
        mu, _ = model.predict(test)
    finally:
        for c in train:
            c.measured = None
    m = regression_metrics([float(x) for x in mu], [c.truth for c in test])
    return m, time.time() - t0


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "CHEMBL203"
    standard = sys.argv[2] if len(sys.argv) > 2 else "IC50"
    sizes = ([int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else SIZES)

    from dmta.chemprop_model import ChempropModel, chemprop_available
    if not chemprop_available():
        sys.exit("chemprop is not installed — this script has nothing to compare")

    compounds = load_compounds(target, standard)
    train_all, test = scaffold_split(compounds, test_frac=0.2)
    print(f"{target} {standard}: {len(compounds)} compounds "
          f"-> {len(train_all)} train / {len(test)} test (scaffold split, fixed)\n")

    rng = np.random.default_rng(0)
    header = (f"{'n train':>9}{'forest rho':>13}{'chemprop rho':>15}"
              f"{'paired diff':>13}{'forest MAE':>12}{'cp MAE':>9}"
              f"{'forest s':>10}{'cp s':>8}")
    print(header); print("-" * len(header))

    rows = []
    for size in sizes:
        n = len(train_all) if size in (0, None) or size > len(train_all) else size
        idx = rng.choice(len(train_all), size=n, replace=False)
        train = [train_all[i] for i in sorted(idx)]

        f_m, f_s = fit_and_score(PropertyModel(n_estimators=300), train, test)
        c_m, c_s = fit_and_score(ChempropModel(epochs=40), train, test)
        if f_m["spearman"] is None or c_m["spearman"] is None:
            print(f"{n:>9}   a model returned no ranking — skipped")
            continue
        diff = c_m["spearman"] - f_m["spearman"]
        rows.append((n, f_m["spearman"], c_m["spearman"], diff))
        print(f"{n:>9}{f_m['spearman']:>13.3f}{c_m['spearman']:>15.3f}"
              f"{diff:>+13.3f}{f_m['mae']:>12.3f}{c_m['mae']:>9.3f}"
              f"{f_s:>10.1f}{c_s:>8.1f}", flush=True)
        if n == len(train_all):
            break                      # asking for more than exists twice is wasted

    if len(rows) < 3:
        sys.exit("\nnot enough sizes completed to read a trend")

    print("\n" + "=" * len(header))
    ns = np.array([r[0] for r in rows], float)
    diffs = np.array([r[3] for r in rows], float)
    slope = float(np.polyfit(np.log10(ns), diffs, 1)[0])
    first, last = diffs[0], diffs[-1]
    print(f"  paired difference (chemprop - forest) at n={int(ns[0])}: {first:+.3f}")
    print(f"  at n={int(ns[-1])}: {last:+.3f}")
    print(f"  trend: {slope:+.3f} of rho per 10x more training data")

    if last > 0 and first <= 0:
        crossover = next((int(n) for n, d in zip(ns, diffs) if d > 0), None)
        print(f"\n  -> IT CROSSES OVER, at about n = {crossover}. Section 7 measured the"
              f"\n     small-data regime and the result does not generalise upwards."
              f"\n     A programme with more than ~{crossover} measured compounds should"
              f"\n     use chemprop; below that the forest is free and just as good.")
    elif slope > 0.02:
        print(f"\n  -> The gap is closing as data grows ({slope:+.3f} per decade) but has"
              f"\n     not crossed by n={int(ns[-1])}. Size is part of the explanation."
              f"\n     Re-test if a programme ever reaches several times this much data.")
    else:
        print(f"\n  -> The gap is FLAT in n ({slope:+.3f} per decade). Data size was not"
              f"\n     the explanation, and section 7's conclusion holds at every size"
              f"\n     this target can reach. That is a stronger result than section 7's.")

    print("\n  Stated fairly to chemprop, exactly as section 7 did: default"
          "\n  hyperparameters, 40 epochs, no validation split and no early stopping."
          "\n  A tuned D-MPNN would do better than this curve shows.")


if __name__ == "__main__":
    main()
