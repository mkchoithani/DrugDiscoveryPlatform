"""Compare two models on identical data: does the neural network earn its keep?

Paired comparison across many real programmes. Both models see exactly the same
scaffold split of exactly the same compounds, so the difference between them is the
model and nothing else. One dataset could not settle this; dozens can.

    python3 model_benchmark.py [n_programmes] [target ...]
"""
from __future__ import annotations

import statistics
import sys
import time

from dmta.chembl import programmes
from dmta.models import PropertyModel, regression_metrics
from dmta.validate import _spearman, scaffold_split

TARGETS = [("CHEMBL203", "IC50"), ("CHEMBL279", "IC50"), ("CHEMBL1824", "IC50"),
           ("CHEMBL244", "IC50"), ("CHEMBL205", "Ki"), ("CHEMBL233", "Ki")]
MIN_SPREAD = 1.5


def evaluate_model(model, train, test) -> dict:
    """Train on `train`, score `test`. Identical protocol for every model."""
    for c in train:
        c.measured = c.truth
    try:
        model.fit(train)
        mu, sigma = model.predict(test)
    finally:
        for c in train:
            c.measured = None
    truth = [c.truth for c in test]
    m = regression_metrics([float(x) for x in mu], truth)

    import numpy as np
    k = max(1, len(test) // 10)
    top = np.argsort(mu)[::-1][:k]
    threshold = float(np.percentile([c.truth for c in train], 90))
    active = np.array(truth) >= threshold
    base = float(active.mean())
    m["enrichment"] = round(float(active[top].mean()) / base, 2) if base else None
    m["sigma_ranks_error"] = _spearman(sigma, np.abs(np.array(mu) - np.array(truth)))
    return m


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    try:
        from dmta.chemprop_model import ChempropModel, chemprop_available
    except ImportError:
        print("chemprop_model not importable — is chemprop installed?")
        raise

    if not chemprop_available():
        sys.exit("chemprop is not installed: pip install chemprop")

    pool = []
    for target, std in TARGETS:
        try:
            pool += [(target, p) for p in programmes(target, std, min_size=60)
                     if p["spread"] >= MIN_SPREAD]
        except Exception as exc:
            print(f"  skipping {target}: {exc}")
    pool.sort(key=lambda tp: -tp[1]["n"])
    pool = pool[:limit]
    print(f"comparing on {len(pool)} programmes "
          f"from {len({t for t, _ in pool})} targets\n")

    header = (f"{'paper':<16}{'n':>5}  {'forest rho':>11}{'chemprop rho':>14}"
              f"{'forest MAE':>12}{'chemprop MAE':>14}{'winner':>10}")
    print(header); print("-" * len(header))

    rows, t0 = [], time.time()
    for target, prog in pool:
        train, test = scaffold_split(prog["compounds"], test_frac=0.25)
        if len(train) < 30 or len(test) < 10:
            continue
        rf = evaluate_model(PropertyModel(n_estimators=300), train, test)
        cp = evaluate_model(ChempropModel(epochs=40), train, test)
        better = ("chemprop" if (cp["spearman"] or -9) > (rf["spearman"] or -9) else "forest")
        rows.append({"document": prog["document"], "n": prog["n"], "rf": rf, "cp": cp,
                     "winner": better})
        print(f"{prog['document']:<16}{prog['n']:>5}  {_f(rf['spearman']):>11}"
              f"{_f(cp['spearman']):>14}{_f(rf['mae']):>12}{_f(cp['mae']):>14}"
              f"{better:>10}", flush=True)

    if not rows:
        sys.exit("no programme produced a usable split")

    rf_rho = [r["rf"]["spearman"] for r in rows if r["rf"]["spearman"] is not None]
    cp_rho = [r["cp"]["spearman"] for r in rows if r["cp"]["spearman"] is not None]
    rf_mae = [r["rf"]["mae"] for r in rows if r["rf"]["mae"] is not None]
    cp_mae = [r["cp"]["mae"] for r in rows if r["cp"]["mae"] is not None]
    wins = sum(1 for r in rows if r["winner"] == "chemprop")

    print(f"\n{len(rows)} programmes in {time.time() - t0:.0f}s\n")
    print(f"  median rho   forest {statistics.median(rf_rho):.3f}   "
          f"chemprop {statistics.median(cp_rho):.3f}")
    print(f"  median MAE   forest {statistics.median(rf_mae):.3f}   "
          f"chemprop {statistics.median(cp_mae):.3f}")
    print(f"  chemprop ranked better in {wins}/{len(rows)} programmes "
          f"({wins / len(rows):.0%})")

    # Paired differences say more than two medians: is the gap consistent?
    paired = [(r["cp"]["spearman"] - r["rf"]["spearman"]) for r in rows
              if r["cp"]["spearman"] is not None and r["rf"]["spearman"] is not None]
    med = statistics.median(paired)
    print(f"\n  median paired difference in rho: {med:+.3f}")
    if len(paired) > 5:
        spread = statistics.pstdev(paired)
        print(f"  spread of that difference: {spread:.3f}")
        if abs(med) < spread / (len(paired) ** 0.5) * 2:
            print("\n  VERDICT: no reliable difference. The gap is smaller than the noise")
            print("  between programmes. Keep the random forest — it trains in seconds,")
            print("  needs no GPU, and the neural network buys nothing measurable here.")
        elif med > 0:
            print(f"\n  VERDICT: chemprop is better by {med:+.3f} rho, consistently.")
            print("  Worth the extra training cost; switch the default model.")
        else:
            print(f"\n  VERDICT: chemprop is WORSE by {med:.3f} rho on this data.")
            print("  Keep the random forest.")


def _f(v, dp=3) -> str:
    return "—" if v is None else f"{v:.{dp}f}"


if __name__ == "__main__":
    main()
