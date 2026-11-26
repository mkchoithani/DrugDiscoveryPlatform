"""Evaluate a customer's dataset: can a model learn it, and would the loop have helped?

Two questions, deliberately separate:

  1. Can a model rank these compounds?  (random / scaffold / time splits)
  2. Replaying history, how many experiments would the loop have needed to reach
     the best compounds the programme actually found?

Question 2 is the one that sells; question 1 is the one that says whether to believe it.
Both are answered against results the customer already has, so every claim is checkable.
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field

import numpy as np

from .campaign import Campaign, CycleSettings
from .cost import Ledger
from .ingest import IngestReport
from .library import Compound
from .scoring import CHEMBL_PROFILE, POTENCY_LED_PROFILE
from .validate import evaluate as split_evaluate
from .validate import random_split, scaffold_split, time_split


@dataclass
class Evaluation:
    ingest: IngestReport
    splits: dict = field(default_factory=dict)
    loop: dict = field(default_factory=dict)
    verdict: dict = field(default_factory=dict)
    settings: dict = field(default_factory=dict)
    # What the evaluation cost to produce. Default "local", so the headline figure is the
    # marginal cash cost on hardware already owned -- which is the honest number for a run
    # on your own machine, and is zero. Pass a rented rate to see what it would cost on
    # someone else's. See dmta/cost.py.
    cost: Ledger = field(default_factory=Ledger)
    cost_rate: str = "local"

    def as_dict(self) -> dict:
        return {"ingest": self.ingest.as_dict(), "splits": self.splits,
                "loop": self.loop, "verdict": self.verdict, "settings": self.settings,
                "cost": self.cost.as_dict(self.cost_rate)}


def run_splits(compounds: list[Compound], n_estimators: int = 300, verbose: bool = False) -> dict:
    out: dict[str, dict] = {}
    tr, te = random_split(compounds)
    out["random"] = split_evaluate(tr, te, n_estimators)
    try:
        tr, te = scaffold_split(compounds)
        out["scaffold"] = split_evaluate(tr, te, n_estimators) if te and tr else {
            "error": "not enough distinct scaffolds to hold any out"}
    except Exception as exc:
        out["scaffold"] = {"error": str(exc)}
    try:
        tr, te, cutoff = time_split(compounds)
        out["time"] = split_evaluate(tr, te, n_estimators)
        out["time"]["cutoff_year"] = cutoff
    except Exception as exc:
        out["time"] = {"error": str(exc)}
    if verbose:
        for name, m in out.items():
            note = m.get("error") or f"rho {m.get('spearman')}, enrichment {m.get('enrichment')}x"
            print(f"  {name:<10} {note}", flush=True)
    return out


def _random_baseline_curve(truth: np.ndarray, batch: int, cycles: int, draws: int = 200,
                           seed: int = 7) -> list[float]:
    """Expected best-found after each batch if compounds were picked at random."""
    rng = np.random.default_rng(seed)
    n = len(truth)
    curve = []
    for c in range(1, cycles + 1):
        k = min(n, batch * c)
        best = [truth[rng.choice(n, k, replace=False)].max() for _ in range(draws)]
        curve.append(round(float(np.mean(best)), 3))
    return curve


def run_loop(compounds: list[Compound], batch_size: int | None = None,
             cycles: int | None = None, top_n: int = 10, seed: int = 1,
             objective: str = "potency", verbose: bool = False) -> dict:
    """Replay the programme: the loop may only see a result once it picks that compound.

    `objective` must match how success is measured. The replay is scored on potency,
    so it runs potency-led by default; the balanced profile also weighs MW, QED and
    synthetic accessibility, which legitimately steers away from potent-but-heavy
    compounds and would be marked down by a potency-only yardstick.
    """
    n = len(compounds)
    batch_size = batch_size or max(8, min(32, round(n * 0.04)))
    cycles = cycles or max(4, min(12, (n // 2) // batch_size))

    truth = np.array([c.truth for c in compounds])
    best_possible = float(truth.max())
    target_ids = {c.id for c in sorted(compounds, key=lambda c: -c.truth)[:top_n]}

    fresh = []
    for c in compounds:                       # a detached copy; the loop must start blind
        import copy
        d = copy.copy(c)
        d.measured = d.predicted = d.uncertainty = d.score = None
        d.cycle_made = None
        fresh.append(d)

    profile = POTENCY_LED_PROFILE if objective == "potency" else CHEMBL_PROFILE
    campaign = Campaign(name="evaluation", library=fresh, profile=list(profile),
                        settings=CycleSettings(batch_size=batch_size, acquisition="greedy",
                                               diversity=0.15, beta=0.25, oracle_seed=seed,
                                               hit_threshold=float(np.percentile(truth, 90))))
    history = []
    found_curve, best_curve = [], []
    for _ in range(cycles):
        s = campaign.run_cycle()
        tested = campaign.measured
        found = sum(1 for c in tested if c.id in target_ids)
        history.append(s)
        found_curve.append(found)
        best_curve.append(round(max(c.measured for c in tested), 3))
        if verbose:
            print(f"    cycle {s['cycle']}: {s['tested_total']:>4} tested, "
                  f"best {s['best_measured']}, {found}/{top_n} of the top set found",
                  flush=True)

    tested_total = len(campaign.measured)
    random_curve = _random_baseline_curve(truth, batch_size, cycles)

    # how many experiments to reach the programme's own best compound
    reached_best = next((i + 1 for i, b in enumerate(best_curve) if b >= best_possible), None)
    experiments_for_best = reached_best * batch_size if reached_best else None

    rho = [h["prospective"]["spearman"] for h in history
           if h["prospective"] and h["prospective"]["spearman"] is not None]

    # Top-decile enrichment: a steadier measure than "did it find these exact 10",
    # which is dominated by where a handful of outliers happen to sit.
    decile_cut = float(np.percentile(truth, 90))
    decile_ids = {c.id for c in compounds if c.truth >= decile_cut}
    decile_found = sum(1 for c in campaign.measured if c.id in decile_ids)
    decile_expected = len(decile_ids) * tested_total / n

    # Random ordering reaches the single best compound at position (n+1)/2 on average.
    random_experiments_for_best = round((n + 1) / 2)

    return {
        "objective": objective,
        "library_size": n,
        "batch_size": batch_size,
        "cycles": cycles,
        "tested_total": tested_total,
        "fraction_tested": round(tested_total / n, 3),
        "top_n": top_n,
        "found_of_top": found_curve[-1] if found_curve else 0,
        "found_curve": found_curve,
        "best_curve": best_curve,
        "random_curve": random_curve,
        "best_possible": round(best_possible, 3),
        "best_found": best_curve[-1] if best_curve else None,
        "random_expected_best": random_curve[-1] if random_curve else None,
        "experiments_for_best": experiments_for_best,
        "experiments_available": n,
        "late_rho": round(statistics.mean(rho[-3:]), 3) if len(rho) >= 3 else None,
        "random_found_expected": round(len(target_ids) * tested_total / n, 2),
        "decile_size": len(decile_ids),
        "decile_found": decile_found,
        "decile_expected": round(decile_expected, 2),
        "decile_lift": round(decile_found / decile_expected, 2) if decile_expected else None,
        "random_experiments_for_best": random_experiments_for_best,
        "speedup_to_best": (round(random_experiments_for_best / experiments_for_best, 2)
                            if experiments_for_best else None),
    }


def verdict(splits: dict, loop: dict) -> dict:
    """A plain recommendation that must agree with every number it is based on.

    The model can rank well and the replay still disappoint, or the reverse. When the
    two disagree the verdict says so rather than quoting whichever is flattering.
    """
    scaffold = splits.get("scaffold", {})
    rho = scaffold.get("spearman")
    lift = loop.get("decile_lift")
    speedup = loop.get("speedup_to_best")
    model_ok = rho is not None and rho >= 0.5
    loop_ok = (lift or 0) >= 1.3 or (speedup or 0) >= 1.5

    def loop_sentence() -> str:
        bits = [f"replaying the programme, the loop reached the best compound after "
                f"{loop['experiments_for_best']} experiments against about "
                f"{loop['random_experiments_for_best']} for random ordering"
                if loop.get("experiments_for_best") else
                f"the loop did not reach the best compound within {loop['tested_total']} experiments"]
        bits.append(f"and found {loop['decile_found']} of the top {loop['decile_size']} "
                    f"compounds against {loop['decile_expected']} expected by chance")
        return ", ".join(bits)

    if rho is None:
        headline, detail = "inconclusive", "The scaffold split could not be run on this dataset."
    elif model_ok and loop_ok:
        headline = "the loop should help"
        detail = (f"A model ranks chemotypes it has never seen at rho {rho}, and {loop_sentence()}.")
    elif model_ok and not loop_ok:
        headline = "the model works, the replay does not yet show a saving"
        detail = (f"Ranking of unseen chemotypes is strong (rho {rho}), but {loop_sentence()}. "
                  "That usually means the best compounds sit in chemotypes the loop never "
                  "explored. Worth one prospective cycle, not a commitment.")
    elif not model_ok and loop_ok:
        headline = "promising, but the model is the weak link"
        detail = (f"The replay looks good — {loop_sentence()} — yet the model ranks unseen "
                  f"chemotypes only at rho {rho}. The saving may not survive on new chemistry.")
    elif rho >= 0.3:
        headline = "promising, worth a prospective trial"
        detail = (f"Ranking of unseen chemotypes is moderate (rho {rho}), and {loop_sentence()}. "
                  "Worth one real cycle before committing.")
    else:
        headline = "not yet — the data does not support it"
        detail = (f"A model cannot rank chemotypes it has not seen on this dataset "
                  f"(rho {rho}), and {loop_sentence()}. More data, or a narrower series, "
                  "is needed first.")

    time_split = splits.get("time", {})
    caveat = None
    if "error" not in time_split and time_split.get("spearman") is not None:
        if time_split["spearman"] < 0.3:
            caveat = (f"The time split is weak (rho {time_split['spearman']}): on this "
                      "programme the model does not extrapolate to later chemistry. Expect it "
                      "to optimise within a series rather than invent new ones.")
    return {"headline": headline, "detail": detail, "lift_vs_random": lift, "caveat": caveat}


def evaluate_dataset(compounds: list[Compound], ingest: IngestReport,
                     batch_size: int | None = None, cycles: int | None = None,
                     top_n: int = 10, n_estimators: int = 300,
                     objective: str = "potency", verbose: bool = False,
                     cost_rate: str = "local") -> Evaluation:
    led = Ledger()
    if verbose:
        print(f"  modelling {len(compounds)} compounds…", flush=True)
    with led.step("three splits, model fitted and scored on each"):
        splits = run_splits(compounds, n_estimators=n_estimators)
    if verbose:
        print("  replaying the programme…", flush=True)
    with led.step("retrospective loop replay"):
        loop = run_loop(compounds, batch_size=batch_size, cycles=cycles, top_n=top_n,
                        objective=objective, verbose=verbose)
    # The denominators. "verdicts" is 1 on purpose: one run answers one question, and
    # that is the number a buyer is actually paying for. compounds_scored is reported
    # beside it because it is the figure that flatters, and hiding it would be its own
    # kind of dishonesty -- it is simply not the headline.
    led.count(verdicts=1, compounds_scored=len(compounds),
              models_fitted=3 + loop.get("cycles", 0))
    return Evaluation(
        ingest=ingest, splits=splits, loop=loop, verdict=verdict(splits, loop),
        settings={"n_estimators": n_estimators, "top_n": top_n,
                  "batch_size": loop["batch_size"], "cycles": loop["cycles"]},
        cost=led, cost_rate=cost_rate,
    )
