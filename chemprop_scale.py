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

import json
import sys
import time
from pathlib import Path

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
    sizes = ([int(x) for x in sys.argv[3].split(",")]
             if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else SIZES)
    save_to = None
    if "--save" in sys.argv:
        i = sys.argv.index("--save")
        if i + 1 >= len(sys.argv):
            sys.exit("--save needs a filename")
        save_to = Path(sys.argv[i + 1])

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

    if save_to is not None:
        save_to.write_text(json.dumps(
            {"target": target, "standard": standard,
             "n_test": len(test), "n_train_available": len(train_all),
             "curve": [{"n": r[0], "forest": r[1], "chemprop": r[2], "diff": r[3]}
                       for r in rows]}, indent=1))
        print(f"\n  curve -> {save_to}")

    print("\n" + "=" * len(header))
    for line in read_trend([r[0] for r in rows], [r[3] for r in rows]):
        print(line)

    print("\n  Stated fairly to chemprop, exactly as section 7 did: default"
          "\n  hyperparameters, 40 epochs, no validation split and no early stopping."
          "\n  A tuned D-MPNN would do better than this curve shows.")


def read_trend(ns_in, diffs_in, plateau_from: int = 3) -> list[str]:
    """Read the curve, allowing for the fact that it may have two regimes.

    The first version of this fitted ONE log-linear slope through every point and
    reported it. On the real curve that gave +0.125 of rho per decade, which
    extrapolates to a crossover at n = 7,111 -- while the measured difference at
    n = 7,173 was -0.033, nowhere near zero. The slope was dominated by the n = 250
    point, where the gap is -0.233.

    That is the same error as comparing a saturated enrichment factor (section 14):
    a summary statistic used where it is not valid. Fit a straight line through a
    curve with a knee in it and the line describes the knee, not where you are now.

    So the last `plateau_from` points get their own slope. If the two disagree in
    size or sign, the recent one is what describes the current regime, and the
    overall slope is reported only as history.
    """
    import numpy as np
    ns, diffs = np.asarray(ns_in, float), np.asarray(diffs_in, float)
    overall = float(np.polyfit(np.log10(ns), diffs, 1)[0])
    out = [f"  paired difference (chemprop - forest) at n={int(ns[0])}: {diffs[0]:+.3f}",
           f"  at n={int(ns[-1])}: {diffs[-1]:+.3f}",
           f"  slope over the whole curve: {overall:+.3f} of rho per 10x data"]

    if len(ns) < plateau_from + 1:
        out.append("  (too few sizes to look for a second regime)")
        recent, level, spread = overall, float(diffs[-1]), 0.0
    else:
        tail_n, tail_d = ns[-plateau_from:], diffs[-plateau_from:]
        recent = float(np.polyfit(np.log10(tail_n), tail_d, 1)[0])
        level, spread = float(tail_d.mean()), float(tail_d.std(ddof=1))
        out.append(f"  slope over the last {plateau_from} sizes "
                   f"(n={int(tail_n[0])}-{int(tail_n[-1])}): {recent:+.3f} per 10x")
        out.append(f"  level over those sizes: {level:+.3f} (spread {spread:.3f})")

    flat = abs(recent) < 0.03
    two_regime = abs(overall - recent) > 0.05

    if diffs[-1] > 0:
        out.append("\n  -> IT HAS CROSSED. chemprop is ahead at the largest size tested.")
    elif flat and two_regime:
        out.append(
            f"\n  -> TWO REGIMES, AND IT HAS PLATEAUED. The gap closes steeply while data"
            f"\n     is scarce, then settles at {level:+.3f} and stops moving"
            f"\n     ({recent:+.3f} per decade over the last {plateau_from} sizes)."
            f"\n     The {overall:+.3f} whole-curve slope is driven by the smallest size and"
            f"\n     must NOT be extrapolated: it would predict a crossover that the last"
            f"\n     {plateau_from} points rule out."
            f"\n     Data size explains the large gap at small n and not the residual one."
            f"\n     There is no n in reach where switching pays.")
    elif flat:
        out.append(f"\n  -> FLAT at {level:+.3f}. Data size was never the explanation.")
    elif recent > 0:
        out.append(
            f"\n  -> STILL CLOSING at {recent:+.3f} per decade, not yet crossed."
            f"\n     Re-test if a programme reaches several times this much data.")
    else:
        out.append(f"\n  -> WIDENING at {recent:+.3f} per decade. More data is not the answer.")
    return out


if __name__ == "__main__":
    main()
