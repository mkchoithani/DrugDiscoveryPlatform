"""Run the evaluation across many real programmes, not one.

A single good result is an anecdote. This splits a target's data into single-paper
campaigns — one team, one series — and evaluates each, so the answer is a
distribution: in how many programmes does the loop actually help, and in which
does it not?

    python3 programme_sweep.py [target] [min_size] [max_programmes]
"""
from __future__ import annotations

import statistics
import sys
import time

from dmta.chembl import programmes
from dmta.evaluate import run_loop
from dmta.ingest import IngestReport
from dmta.validate import evaluate as split_evaluate
from dmta.validate import scaffold_split

MIN_SPREAD = 1.5          # below this there is nothing for a model to learn


def assess(prog: dict) -> dict:
    compounds = prog["compounds"]
    row = {"document": prog["document"], "n": prog["n"], "scaffolds": prog["scaffolds"],
           "spread": prog["spread"], "year": prog["year"]}
    try:
        train, test = scaffold_split(compounds, test_frac=0.25)
        if len(train) < 20 or len(test) < 8:
            row["note"] = "too few compounds to hold chemotypes out"
            return row
        m = split_evaluate(train, test, n_estimators=150)
        row["rho"] = m["spearman"]
        row["enrichment"] = m["enrichment"]
    except Exception as exc:
        row["note"] = f"split failed: {exc}"
        return row

    loop = run_loop(compounds, top_n=max(5, prog["n"] // 20), objective="potency")
    row["lift"] = loop["decile_lift"]
    row["speedup"] = loop["speedup_to_best"]
    row["tested"] = loop["tested_total"]
    row["found"] = f"{loop['decile_found']}/{loop['decile_size']}"
    return row


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "CHEMBL203"
    min_size = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10

    found = programmes(target, min_size=min_size)
    usable = [p for p in found if p["spread"] >= MIN_SPREAD]
    skipped = [p for p in found if p["spread"] < MIN_SPREAD]
    print(f"{target}: {len(found)} programmes of {min_size}+ compounds; "
          f"{len(usable)} have enough potency spread to learn from "
          f"({len(skipped)} skipped, spread < {MIN_SPREAD} log)\n")

    header = (f"{'paper':<16}{'n':>5}{'chemo':>7}{'spread':>8}"
              f"{'rho':>7}{'enrich':>8}{'lift':>7}{'speedup':>9}{'found':>10}")
    print(header); print("-" * len(header))

    rows, t0 = [], time.time()
    for prog in usable[:limit]:
        row = assess(prog)
        rows.append(row)
        if "note" in row:
            print(f"{row['document']:<16}{row['n']:>5}{row['scaffolds']:>7}"
                  f"{row['spread']:>8}   {row['note']}")
            continue
        print(f"{row['document']:<16}{row['n']:>5}{row['scaffolds']:>7}{row['spread']:>8}"
              f"{_f(row['rho']):>7}{_f(row['enrichment'], 2):>8}{_f(row['lift'], 2):>7}"
              f"{_f(row['speedup'], 2):>9}{row['found']:>10}", flush=True)

    done = [r for r in rows if "note" not in r]
    if not done:
        print("\nno programme could be evaluated")
        return

    rhos = [r["rho"] for r in done if r["rho"] is not None]
    lifts = [r["lift"] for r in done if r["lift"] is not None]
    speedups = [r["speedup"] for r in done if r["speedup"] is not None]
    helped = sum(1 for r in done if (r["lift"] or 0) >= 1.3 or (r["speedup"] or 0) >= 1.5)

    print(f"\n{len(done)} programmes evaluated in {time.time() - t0:.0f}s")
    print(f"  model ranks unseen chemotypes: median rho {statistics.median(rhos):.3f}"
          f"  (range {min(rhos):.2f} to {max(rhos):.2f})")
    if lifts:
        print(f"  loop enrichment over chance:   median {statistics.median(lifts):.2f}x"
              f"  (range {min(lifts):.2f} to {max(lifts):.2f})")
    if speedups:
        print(f"  speed to the best compound:    median {statistics.median(speedups):.2f}x")
    print(f"\n  the loop helped in {helped} of {len(done)} programmes "
          f"({helped / len(done):.0%})")
    print("  'helped' means 1.3x or better enrichment, or reaching the best compound "
          "at least 1.5x sooner than random ordering.")


def _f(v, dp=3) -> str:
    return "—" if v is None else f"{v:.{dp}f}"


if __name__ == "__main__":
    main()
