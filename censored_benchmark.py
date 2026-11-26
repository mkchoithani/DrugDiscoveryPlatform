"""Does using the too-weak-to-measure compounds actually help?

A known-answer test. Real programmes where every potency IS known, artificially
censored: every training compound below a bound is replaced by "weaker than this",
which is what the assay would have reported. The held-out set keeps its true values,
so the question has a definite answer.

Four ways to handle it, on identical splits:

  oracle    train on the true values. Not available in real life; the ceiling.
  drop      discard the censored compounds. What dmta/ingest.py did until now.
  bound     use the bound as if it were a real measurement. Common, and wrong:
            it asserts a precise value the assay never reported.
  censored  Tobit-style EM — each censored compound contributes "no higher
            than this" (dmta/censored.py).

The censoring fraction defaults to 0.4, because that is roughly real: 51 of 123
Open Source Malaria compounds are reported as "> X".

    python3 censored_benchmark.py [n_programmes] [censor_fraction]
"""
from __future__ import annotations

import statistics
import sys
import time

import numpy as np

from dmta.censored import CensoredModel, Censoring
from dmta.chembl import programmes
from dmta.models import PropertyModel, regression_metrics
from dmta.validate import scaffold_split

TARGETS = [("CHEMBL203", "IC50"), ("CHEMBL279", "IC50"), ("CHEMBL1824", "IC50"),
           ("CHEMBL244", "IC50"), ("CHEMBL205", "Ki"), ("CHEMBL233", "Ki")]


def evaluate(predictions, test) -> dict:
    truth = [c.truth for c in test]
    m = regression_metrics([float(v) for v in predictions], truth)
    threshold = float(np.percentile(truth, 90))
    active = np.array(truth) >= threshold
    k = max(1, len(test) // 10)
    top = np.argsort(predictions)[::-1][:k]
    base = float(active.mean())
    m["enrichment"] = round(float(active[top].mean()) / base, 2) if base else None
    return m


def arms(train, test, bound):
    """Train each way on the same split and score against held-out truth."""
    observed = [c for c in train if c.truth >= bound]
    censored = [c for c in train if c.truth < bound]
    out = {}

    def with_measured(compounds, values, fn):
        saved = [c.measured for c in compounds]
        for c, v in zip(compounds, values):
            c.measured = v
        try:
            return fn()
        finally:
            for c, v in zip(compounds, saved):
                c.measured = v

    # oracle: every true value
    out["oracle"] = with_measured(
        train, [c.truth for c in train],
        lambda: evaluate(PropertyModel(n_estimators=300).fit(train).predict(test)[0], test))

    # drop: observed only
    if len(observed) >= 20:
        out["drop"] = with_measured(
            observed, [c.truth for c in observed],
            lambda: evaluate(PropertyModel(n_estimators=300).fit(observed).predict(test)[0], test))

    # bound: censored compounds asserted to sit exactly at the bound
    out["bound"] = with_measured(
        train, [c.truth if c.truth >= bound else bound for c in train],
        lambda: evaluate(PropertyModel(n_estimators=300).fit(train).predict(test)[0], test))

    # censored: EM, contributing "no higher than the bound"
    if censored and len(observed) >= 20:
        records = [Censoring(bound=bound, direction="left") for _ in censored]
        def run():
            model = CensoredModel(iterations=4).fit(observed, censored, records)
            return evaluate(model.predict(test)[0], test)
        out["censored"] = with_measured(observed, [c.truth for c in observed], run)
    return out, len(observed), len(censored)


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    fraction = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4

    pool = []
    for target, std in TARGETS:
        try:
            pool += [p for p in programmes(target, std, min_size=80)
                     if p["spread"] >= 1.5]
        except Exception as exc:
            print(f"  skipping {target}: {exc}")
    pool.sort(key=lambda p: -p["n"])
    pool = pool[:limit]
    print(f"{len(pool)} programmes, censoring the weakest {fraction:.0%} of each "
          f"training set\n")

    names = ["oracle", "drop", "bound", "censored"]
    header = f"{'paper':<16}{'n':>5}{'cens':>6}" + "".join(f"{n:>10}" for n in names)
    print(header); print("-" * len(header))

    rows, t0 = [], time.time()
    for prog in pool:
        train, test = scaffold_split(prog["compounds"], test_frac=0.25)
        if len(train) < 50 or len(test) < 10:
            continue
        bound = float(np.percentile([c.truth for c in train], fraction * 100))
        result, n_obs, n_cens = arms(train, test, bound)
        if "censored" not in result or "drop" not in result:
            continue
        rows.append(result)
        print(f"{prog['document']:<16}{prog['n']:>5}{n_cens:>6}"
              + "".join(f"{(result[n]['spearman'] if result.get(n) and result[n]['spearman'] is not None else float('nan')):>10.3f}"
                        for n in names), flush=True)

    if not rows:
        sys.exit("no programme produced a usable split")

    print(f"\n{len(rows)} programmes in {time.time() - t0:.0f}s\n")
    print(f"  {'arm':<12}{'median rho':>12}{'median MAE':>12}{'median enrich':>15}")
    for n in names:
        rho = [r[n]["spearman"] for r in rows if r.get(n) and r[n]["spearman"] is not None]
        mae = [r[n]["mae"] for r in rows if r.get(n) and r[n]["mae"] is not None]
        enr = [r[n]["enrichment"] for r in rows if r.get(n) and r[n]["enrichment"] is not None]
        print(f"  {n:<12}{statistics.median(rho):>12.3f}{statistics.median(mae):>12.3f}"
              f"{statistics.median(enr):>15.2f}")

    # Paired against dropping them, which is the behaviour being replaced.
    print("\n  paired against 'drop', programme by programme:")
    for n in ("bound", "censored", "oracle"):
        paired = [r[n]["spearman"] - r["drop"]["spearman"] for r in rows
                  if r.get(n) and r[n]["spearman"] is not None
                  and r["drop"]["spearman"] is not None]
        if len(paired) < 2:
            continue
        med = statistics.median(paired)
        spread = statistics.pstdev(paired)
        stderr = spread / (len(paired) ** 0.5)
        wins = sum(1 for d in paired if d > 0)
        verdict = ("no reliable difference" if abs(med) < 2 * stderr
                   else ("BETTER" if med > 0 else "WORSE"))
        print(f"    {n:<10} median {med:+.3f}  2 s.e. {2 * stderr:.3f}  "
              f"wins {wins}/{len(paired)}  -> {verdict}")

    print("\n" + "=" * 68)
    cen = [r["censored"]["spearman"] - r["drop"]["spearman"] for r in rows
           if r["censored"]["spearman"] is not None and r["drop"]["spearman"] is not None]
    bnd = [r["bound"]["spearman"] - r["drop"]["spearman"] for r in rows
           if r["bound"]["spearman"] is not None and r["drop"]["spearman"] is not None]
    mc, mb = statistics.median(cen), statistics.median(bnd)
    se = statistics.pstdev(cen) / (len(cen) ** 0.5)
    if mc > 2 * se and mc > mb:
        print(f"  VERDICT: keeping the censored compounds helps, by {mc:+.3f} rho.")
        print("  Wire it into ingest.py and stop dropping them.")
    elif mb > 2 * se and mb >= mc:
        print(f"  VERDICT: the naive bound substitution wins ({mb:+.3f} vs {mc:+.3f}).")
        print("  Use the censored compounds, but the EM machinery is not earning its")
        print("  keep — prefer the simpler treatment.")
    else:
        print(f"  VERDICT: no reliable gain ({mc:+.3f} rho, 2 s.e. {2 * se:.3f}).")
        print("  Dropping censored compounds costs less than it appeared to. Keep")
        print("  them in ingest for provenance, but do not claim a modelling win.")


if __name__ == "__main__":
    main()
