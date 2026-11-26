"""Audit our own headline result against PostEra's DEL Simulator landscape.

RESULTS.md 2 reports that UCB selection finds 2.75x as many of the top 5% as random
screening. That was measured on a landscape we wrote ourselves, and a selection
strategy can look good on one landscape for reasons that have nothing to do with
the method. This runs the identical comparison on two landscapes:

  ours    dmta/assay.py's SimulatedAssay over our enumerated 924-compound library
  theirs  the DEL Simulator's affinities over a size-matched sample of its library,
          anchored on 2,660 real soluble-epoxide-hydrolase binders

Everything else is held fixed: same cycles, same batch size, same acquisition
functions, same scoring profile, same assay noise, same number of landscape
instances. If the ordering of acquisition functions survives the change of
landscape, the result is about the method. If it does not, our headline number is
a property of our own simulator and RESULTS.md has to say so.

The profile is POTENCY_LED deliberately. The metric is potency-only ("did you find
the top 5%"), and scoring against a profile where potency is a third of the weight
penalises the loop for obeying its instructions — a mistake this repository already
made once. The DEL compounds also have no enumerated route, so route objectives
could not apply to them anyway.

    python3 del_benchmark.py [cycles] [batch] [del_data_dir] [seeds...]
"""
from __future__ import annotations

import statistics
import sys
import time
from pathlib import Path

from dmta.assay import SimulatedAssay
from dmta.campaign import Campaign, CycleSettings, fresh_library
from dmta.del_landscape import NoisyReplayAssay, landscape_summary, load_del_landscape
from dmta.scoring import POTENCY_LED_PROFILE

STRATEGIES = ["ucb", "greedy", "expected_improvement", "uncertainty", "random"]
DEFAULT_DEL_DIR = "/tmp/delsim/del-simulator/data/experiments/example"


def our_landscape(seed: int):
    library = fresh_library()
    return library, SimulatedAssay(library, seed=seed, noise=0.2)


def their_landscape(data_dir: Path, seed: int, n: int, sampling: str = "uniform"):
    library = load_del_landscape(data_dir / "library.csv",
                                 data_dir / "affinity_seh_512.csv", n=n, seed=seed,
                                 sampling=sampling)
    return library, NoisyReplayAssay(library, noise=0.2, seed=seed)


def run(strategy: str, library, assay, cycles: int, batch: int, seed: int) -> dict:
    c = Campaign(
        name=f"del-bench-{strategy}-{seed}",
        settings=CycleSettings(batch_size=batch, acquisition=strategy, beta=1.0,
                               diversity=0.0 if strategy == "random" else 0.15,
                               oracle_seed=seed, assay_noise=0.2),
        profile=list(POTENCY_LED_PROFILE),
        library=library, assay=assay,
    )
    target_set = assay.top_fraction(0.05)
    for _ in range(cycles):
        c.run_cycle()
    tested = c.measured
    found = sum(1 for x in tested if x.id in target_set)
    prospective = [h["prospective"] for h in c.history if h["prospective"]]
    rhos = [p["spearman"] for p in prospective if p and p.get("spearman") is not None]
    return {
        "found": found, "target_total": len(target_set), "tested": len(tested),
        "best": max(x.measured for x in tested),
        "ceiling": assay.ceiling(),
        "rho": statistics.mean(rhos) if rhos else None,
    }


def report(label: str, library_for, cycles: int, batch: int, seeds: list[int]) -> dict:
    print(f"\n--- {label} ---")
    first_library, first_assay = library_for(seeds[0])
    summary = landscape_summary_of(first_library, first_assay)
    print("  landscape: " + ", ".join(f"{k} {v}" for k, v in summary.items()))
    header = (f"  {'strategy':<22}{'top5% found':>13}{'best':>8}{'ceiling':>9}{'rho':>7}")
    print(header); print("  " + "-" * (len(header) - 2))

    out = {}
    for s in STRATEGIES:
        runs = []
        for seed in seeds:
            library, assay = library_for(seed)
            runs.append(run(s, library, assay, cycles, batch, seed))
        found = statistics.mean(r["found"] for r in runs)
        total = statistics.mean(r["target_total"] for r in runs)
        rhos = [r["rho"] for r in runs if r["rho"] is not None]
        out[s] = found
        print(f"  {s:<22}{f'{found:.1f} / {total:.0f}':>13}"
              f"{statistics.mean(r['best'] for r in runs):>8.2f}"
              f"{statistics.mean(r['ceiling'] for r in runs):>9.2f}"
              f"{(statistics.mean(rhos) if rhos else float('nan')):>7.2f}", flush=True)
    ratio = out["ucb"] / max(out["random"], 0.5)
    print(f"  UCB over random: {ratio:.2f}x  ({out['ucb']:.1f} vs {out['random']:.1f})")
    return {"found": out, "ratio": ratio}


def landscape_summary_of(library, assay) -> dict:
    s = landscape_summary(library, assay)
    return {"n": s["n"], "median": s["median"], "sd": s["sd"],
            "top5%≥": s["top5pct_threshold"], "scaffolds": s["scaffolds"]}


def main():
    cycles = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    batch = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    data_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(DEFAULT_DEL_DIR)
    seeds = [int(s) for s in sys.argv[4:]] or [42, 7, 101]

    missing = [f for f in ("library.csv", "affinity_seh_512.csv")
               if not (data_dir / f).exists()]
    if missing:
        sys.exit(f"missing {', '.join(missing)} in {data_dir}\n"
                 "Generate them in the del-simulator repository first — see "
                 "dmta/del_landscape.py for the two commands.")

    n = len(fresh_library())            # size-match their sample to our library
    print(f"{cycles} cycles x {batch} compounds, {len(seeds)} landscape instances "
          f"per arm (seeds {seeds})")
    print(f"identical protocol on both; their library sampled to n={n} to match ours")

    t0 = time.time()
    ours = report("our landscape (dmta SimulatedAssay, 924 cmpds / 5 scaffolds)",
                  our_landscape, cycles, batch, seeds)
    theirs = report("their landscape, UNIFORM sample (every compound its own chemotype)",
                    lambda seed: their_landscape(data_dir, seed, n, "uniform"),
                    cycles, batch, seeds)
    dense = report("their landscape, DENSE sample (products sharing two of three cycles)",
                   lambda seed: their_landscape(data_dir, seed, n, "dense"),
                   cycles, batch, seeds)

    print("\n" + "=" * 62)
    print(f"  UCB over random, our landscape:            {ours['ratio']:.2f}x")
    print(f"  UCB over random, their landscape uniform:  {theirs['ratio']:.2f}x")
    print(f"  UCB over random, their landscape dense:    {dense['ratio']:.2f}x")

    # The ordering matters more than the ratio: a method that ranks strategies the
    # same way on an independently authored landscape is a method, not an artefact.
    def ranking(found: dict) -> list[str]:
        return sorted(found, key=lambda k: -found[k])
    ours_rank = ranking(ours["found"])
    print(f"\n  our ranking:            {' > '.join(ours_rank)}")
    print(f"  their uniform ranking:  {' > '.join(ranking(theirs['found']))}")
    print(f"  their dense ranking:    {' > '.join(ranking(dense['found']))}")

    print()
    if theirs["ratio"] > 1.0 and dense["ratio"] > 1.0:
        print("  VERDICT: the result survives an independently authored landscape,")
        print("  whatever the shape of the sample. It is about the method.")
    elif dense["ratio"] > 1.0:
        print("  VERDICT: the result survives, but only where the library has")
        print("  near-neighbours to generalise from. On their landscape the loop beats")
        print("  random screening in a DENSE sample and not in a uniform one — so what")
        print("  the 2.75x really depends on is library structure, not our simulator.")
        print("  RESULTS.md must state that precondition: active learning needs a")
        print("  library whose members resemble each other, and a uniform draw over a")
        print("  large combinatorial space is not one.")
    elif theirs["ratio"] > 1.0:
        print("  VERDICT: survives in the uniform sample but not the dense one — the")
        print("  opposite of the expected failure, and worth investigating before any")
        print("  of this is relied on.")
    else:
        print("  VERDICT: DOES NOT SURVIVE. The loop does not beat random screening on")
        print("  an independently authored landscape in either sample shape. The")
        print("  headline number in RESULTS.md is a property of our own simulator and")
        print("  must be reported as such.")
    print(f"\n  wall clock {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
