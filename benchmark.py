"""Does the loop actually beat screening? Compare acquisition strategies on the same landscape.

  python3 benchmark.py [cycles] [batch_size] [seeds...]
"""
from __future__ import annotations

import statistics
import sys
import time

from dmta.assay import SimulatedAssay
from dmta.campaign import Campaign, CycleSettings

STRATEGIES = ["ucb", "greedy", "expected_improvement", "uncertainty", "random"]


def run(strategy: str, cycles: int, batch: int, seed: int) -> dict:
    c = Campaign(
        name=f"bench-{strategy}-{seed}",
        settings=CycleSettings(batch_size=batch, acquisition=strategy, beta=1.0,
                               diversity=0.0 if strategy == "random" else 0.15,
                               oracle_seed=seed, assay_noise=0.2),
    )
    oracle: SimulatedAssay = c.assay  # type: ignore[assignment]
    target_set = oracle.top_fraction(0.05)
    for _ in range(cycles):
        c.run_cycle()
    tested = c.measured
    found = sum(1 for x in tested if x.id in target_set)
    prospective = [h["prospective"] for h in c.history if h["prospective"]]
    return {
        "strategy": strategy,
        "tested": len(tested),
        "best": round(max(x.measured for x in tested), 2),
        "ceiling": round(oracle.ceiling(), 2),
        "top5pct_found": found,
        "top5pct_total": len(target_set),
        "mean_r2": round(statistics.mean(p["r2"] for p in prospective if p["r2"] is not None), 3) if prospective else None,
        "mean_spearman": round(statistics.mean(p["spearman"] for p in prospective if p["spearman"] is not None), 3) if prospective else None,
        "steps": sum(h["plan"]["shared_steps"] for h in c.history),
        "naive_steps": sum(h["plan"]["naive_steps"] for h in c.history),
        "cost": round(sum(h["plan"]["reagent_cost"] for h in c.history)),
    }


def main():
    cycles = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    batch = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    seeds = [int(s) for s in sys.argv[3:]] or [42, 7, 101]

    print(f"{cycles} cycles x {batch} compounds, {len(seeds)} landscapes (seeds {seeds})\n")
    header = f"{'strategy':<22}{'best':>7}{'ceiling':>9}{'top5% found':>13}{'R2':>7}{'rho':>7}{'steps':>8}{'$':>8}"
    print(header)
    print("-" * len(header))

    rows = {}
    t0 = time.time()
    for s in STRATEGIES:
        runs = [run(s, cycles, batch, seed) for seed in seeds]
        rows[s] = runs
        best = statistics.mean(r["best"] for r in runs)
        ceil = statistics.mean(r["ceiling"] for r in runs)
        found = statistics.mean(r["top5pct_found"] for r in runs)
        total = runs[0]["top5pct_total"]
        r2s = [r["mean_r2"] for r in runs if r["mean_r2"] is not None]
        rhos = [r["mean_spearman"] for r in runs if r["mean_spearman"] is not None]
        print(f"{s:<22}{best:>7.2f}{ceil:>9.2f}"
              f"{f'{found:.1f} / {total}':>13}"
              f"{(statistics.mean(r2s) if r2s else float('nan')):>7.2f}"
              f"{(statistics.mean(rhos) if rhos else float('nan')):>7.2f}"
              f"{statistics.mean(r['steps'] for r in runs):>8.0f}"
              f"{statistics.mean(r['cost'] for r in runs):>8.0f}")

    ucb = statistics.mean(r["top5pct_found"] for r in rows["ucb"])
    rnd = statistics.mean(r["top5pct_found"] for r in rows["random"])
    naive = statistics.mean(r["naive_steps"] for r in rows["ucb"])
    steps = statistics.mean(r["steps"] for r in rows["ucb"])
    tested = rows["ucb"][0]["tested"]
    print(f"\nUCB enrichment over random screening: {ucb / max(rnd, 0.5):.2f}x "
          f"({ucb:.1f} vs {rnd:.1f} of the top 5%, from {tested} compounds made)")
    print(f"Shared intermediates: {steps:.0f} reactions instead of {naive:.0f} "
          f"({(1 - steps / naive) * 100:.0f}% fewer)")
    print(f"\nwall clock {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
