"""Does a model really beat the forest? A paired test over programmes.

[RESULTS.md](RESULTS.md) 15 bootstrapped a paired AUC difference and found that
docking's apparent edge over molecular size contained zero. The model comparisons in
sections 7 and 8 never got the same treatment: they were decided on a median gap and a
2 s.e. screen, and then summarised with the word "rejected". That is a stronger word
than a median supports, so this script applies the section 15 standard to them.

Three tests, because each answers a different question and they can disagree:

  - **Bootstrap over programmes.** Resample the programmes with replacement and
    recompute the mean paired difference. The programme is the unit of resampling, not
    the compound: a model that wins on big programmes and loses on small ones has not
    won, and resampling compounds would hide that. Gives a confidence interval.
  - **Wilcoxon signed-rank.** Uses the *size* of each programme's difference, not just
    its direction. This is the one that matters when a few programmes move a lot.
  - **Sign test.** Counts wins only. Deliberately crude, and reported because a win
    record is what gets quoted ("11 of 18") and it is worth seeing how little that says.

The verdict is deliberately three-valued. "No evidence of a difference" is not
"rejected", and the distinction is the whole point of writing this.

    python3 model_paired_test.py <results.json from model_benchmark.py --save>
"""
from __future__ import annotations

import json
import sys
from math import comb
from pathlib import Path

import numpy as np

RESAMPLES = 20000
BASELINE = "forest"
METRIC = "spearman"


def sign_test(wins: int, n: int) -> float:
    """Two-sided P of a record at least this lopsided, under a fair coin."""
    k = max(wins, n - wins)
    tail = sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def signed_rank_p(diff: np.ndarray, rng: np.random.Generator,
                  resamples: int = RESAMPLES) -> float:
    """Wilcoxon signed-rank, P by permutation of the signs.

    Exact enough and avoids a scipy dependency. Ties in |diff| get average ranks,
    and zero differences are dropped, as the test requires.
    """
    d = diff[diff != 0]
    if d.size == 0:
        return 1.0
    order = np.argsort(np.abs(d))
    ranks = np.empty(d.size, float)
    ranks[order] = np.arange(1, d.size + 1)
    # average ranks within ties of |d|
    absd = np.abs(d)
    for v in np.unique(absd):
        m = absd == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    observed = float(np.sum(ranks * np.sign(d)))
    flips = rng.choice([-1.0, 1.0], size=(resamples, d.size))
    null = flips @ ranks
    return float((np.abs(null) >= abs(observed)).mean())


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    payload = json.loads(Path(sys.argv[1]).read_text())
    names = [n for n in payload["models"] if n != BASELINE]
    progs = payload["programmes"]
    if not names:
        sys.exit(f"nothing to compare: the file holds only {BASELINE}")

    rng = np.random.default_rng(0)
    print(f"{len(progs)} programmes, paired against the {BASELINE}, "
          f"metric = {METRIC}\n")

    for name in names:
        pairs = [(p["scores"][name][METRIC], p["scores"][BASELINE][METRIC], p["n"],
                  p["document"]) for p in progs
                 if p["scores"][name].get(METRIC) is not None
                 and p["scores"][BASELINE].get(METRIC) is not None]
        if len(pairs) < 4:
            print(f"{name}: only {len(pairs)} usable programmes — not testing\n")
            continue
        diff = np.array([a - b for a, b, _, _ in pairs], float)
        wins = int((diff > 0).sum())

        idx = rng.integers(0, diff.size, size=(RESAMPLES, diff.size))
        boot = diff[idx].mean(axis=1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        p_sr = signed_rank_p(diff, rng)
        p_sign = sign_test(wins, diff.size)

        print(f"{name} vs {BASELINE}   ({diff.size} programmes)")
        print(f"  mean paired difference   {diff.mean():+.4f}")
        print(f"  median                   {float(np.median(diff)):+.4f}")
        print(f"  95% CI (bootstrap)       [{lo:+.4f}, {hi:+.4f}]")
        print(f"  Wilcoxon signed-rank     P = {p_sr:.4f}")
        print(f"  sign test ({wins}/{diff.size} wins)".ljust(27) + f"P = {p_sign:.4f}")

        contains_zero = lo <= 0 <= hi
        if contains_zero:
            # The interval's own width says how big a real effect could still hide here.
            bound = max(abs(lo), abs(hi))
            print(f"  -> NO EVIDENCE OF A DIFFERENCE. The interval contains zero, so this"
                  f"\n     is a tie, not a rejection. It does not show the gap is zero:"
                  f"\n     an effect up to {bound:.3f} of {METRIC} is still inside the interval."
                  f"\n     Keeping the {BASELINE} is a cost decision, not a measured win.")
        elif diff.mean() > 0:
            print(f"  -> BETTER than the {BASELINE}, and the interval excludes zero.")
        else:
            print(f"  -> WORSE than the {BASELINE}, and the interval excludes zero. "
                  f"A real rejection.")

        worst = sorted(pairs, key=lambda t: t[0] - t[1])[:2]
        best = sorted(pairs, key=lambda t: t[1] - t[0])[:2]
        print("     largest losses: " + ", ".join(
            f"{d} ({a - b:+.3f}, n={n})" for a, b, n, d in worst))
        print("     largest wins:   " + ", ".join(
            f"{d} ({a - b:+.3f}, n={n})" for a, b, n, d in best))
        print()


if __name__ == "__main__":
    main()
