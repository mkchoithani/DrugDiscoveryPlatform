"""Retrospective active learning on real ChEMBL data.

The loop may only see a compound's measured potency once it decides to make it.
That makes this a fair test of the selection strategy against real SAR, including
all its activity cliffs and assay noise.

    python3 chembl_benchmark.py [target] [cycles] [batch] [standard_type]
"""
from __future__ import annotations

import statistics
import sys
import time

from dmta.campaign import Campaign, CycleSettings
from dmta.chembl import dataset_summary, load_compounds
from dmta.scoring import CHEMBL_PROFILE

STRATEGIES = ["ucb", "greedy", "expected_improvement", "random"]
SEEDS = [1, 2]


def run(library, strategy: str, cycles: int, batch: int, seed: int) -> dict:
    lib = [c for c in library]
    for c in lib:                       # fresh campaign state each run
        c.measured = c.predicted = c.uncertainty = c.score = None
        c.cycle_made = None
    c = Campaign(name=f"chembl-{strategy}-{seed}", library=lib,
                 profile=list(CHEMBL_PROFILE),
                 settings=CycleSettings(batch_size=batch, acquisition=strategy,
                                        diversity=0.0 if strategy == "random" else 0.15,
                                        beta=1.0, oracle_seed=seed, hit_threshold=8.0))
    target = c.assay.top_fraction(0.05)
    for _ in range(cycles):
        c.run_cycle()
    tested = c.measured
    rho = [h["prospective"]["spearman"] for h in c.history
           if h["prospective"] and h["prospective"]["spearman"] is not None]
    return {
        "tested": len(tested),
        "best": max(x.measured for x in tested),
        "actives": sum(1 for x in tested if x.measured >= 8.0),
        "top5pct": sum(1 for x in tested if x.id in target),
        "top5pct_total": len(target),
        "late_rho": round(statistics.mean(rho[-3:]), 3) if len(rho) >= 3 else None,
        "mean_mw": round(statistics.mean(x.desc["mw"] for x in tested), 1),
        "mean_qed": round(statistics.mean(x.desc["qed"] for x in tested), 3),
    }


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "CHEMBL203"
    cycles = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    batch = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    std_type = sys.argv[4] if len(sys.argv) > 4 else "IC50"

    library = load_compounds(target, std_type, verbose=True)
    info = dataset_summary(library)
    print(f"\n{target} ({std_type}): {info['compounds']} compounds, {info['scaffolds']} scaffolds, "
          f"{info['actives_ge_8']} actives (pChEMBL ≥ 8), "
          f"published {info['year_range'][0]}–{info['year_range'][1]}")
    print(f"{cycles} cycles x {batch} compounds = {cycles * batch} made "
          f"({cycles * batch / info['compounds']:.1%} of the library)\n")

    header = (f"{'strategy':<22}{'best':>7}{'actives≥8':>11}{'top 5% found':>14}"
              f"{'late ρ':>9}{'mean MW':>9}{'mean QED':>10}")
    print(header); print("-" * len(header))
    results = {}
    t0 = time.time()
    for s in STRATEGIES:
        runs = [run(library, s, cycles, batch, seed) for seed in SEEDS]
        results[s] = runs
        best = statistics.mean(r["best"] for r in runs)
        act = statistics.mean(r["actives"] for r in runs)
        top = statistics.mean(r["top5pct"] for r in runs)
        rho = [r["late_rho"] for r in runs if r["late_rho"] is not None]
        print(f"{s:<22}{best:>7.2f}{f'{act:.1f} / {runs[0]["tested"]}':>11}"
              f"{f'{top:.1f} / {runs[0]["top5pct_total"]}':>14}"
              f"{(statistics.mean(rho) if rho else float('nan')):>9.2f}"
              f"{statistics.mean(r['mean_mw'] for r in runs):>9.0f}"
              f"{statistics.mean(r['mean_qed'] for r in runs):>10.2f}")

    base_rate = info["actives_ge_8"] / info["compounds"]
    made = results["ucb"][0]["tested"]
    for s in STRATEGIES:
        act = statistics.mean(r["actives"] for r in results[s]) / made
        print(f"\n{s:<22} hit rate {act:.1%}  ({act / base_rate:.2f}x the library's {base_rate:.1%})")
    print(f"\nwall clock {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
