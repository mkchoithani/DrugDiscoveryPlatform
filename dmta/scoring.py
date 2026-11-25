"""Multi-parameter optimisation: desirability curves over predicted and computed properties.

The profile is data, not code — a project can be re-weighted without touching the engine.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Callable, Literal

from .library import Compound


def ramp_up(x, lo, hi):      # 0 below lo, 1 above hi
    if x >= hi: return 1.0
    if x <= lo: return 0.0
    return (x - lo) / (hi - lo)


def ramp_down(x, lo, hi):    # 1 below lo, 0 above hi
    if x <= lo: return 1.0
    if x >= hi: return 0.0
    return (hi - x) / (hi - lo)


def window(x, lo, hi, soft):
    if lo <= x <= hi: return 1.0
    d = lo - x if x < lo else x - hi
    return max(0.0, 1.0 - d / soft)


@dataclass
class Objective:
    key: str
    label: str
    weight: float
    source: Literal["predicted", "descriptor", "route", "docking"]
    field: str
    kind: Literal["up", "down", "window"]
    lo: float
    hi: float
    soft: float = 1.0
    critical: bool = False   # if 0, the whole compound is killed

    def value(self, c: Compound) -> float:
        if self.source == "docking":
            x = c.docking_score            # None until the compound has been docked
        elif self.source == "predicted":
            x = c.predicted if c.predicted is not None else (c.measured if c.measured is not None else 6.0)
        elif self.source == "descriptor":
            x = c.desc.get(self.field)
        else:
            x = getattr(c, self.field)
        if x is None:
            return 0.5
        if self.kind == "up":
            return ramp_up(x, self.lo, self.hi)
        if self.kind == "down":
            return ramp_down(x, self.lo, self.hi)
        return window(x, self.lo, self.hi, self.soft)


# A lead-optimisation profile for an oral small molecule. Edit the numbers, not the code.
DEFAULT_PROFILE: list[Objective] = [
    Objective("potency", "Potency pIC50", 3.0, "predicted", "pIC50", "up", 6.0, 8.0),
    Objective("mw", "MW ≤ 500", 1.0, "descriptor", "mw", "down", 450, 600),
    Objective("clogp", "cLogP 1–4", 1.5, "descriptor", "clogp", "window", 1.0, 4.0, soft=2.0),
    Objective("tpsa", "TPSA 40–120", 1.0, "descriptor", "tpsa", "window", 40, 120, soft=40),
    Objective("hbd", "HBD ≤ 3", 0.5, "descriptor", "hbd", "down", 3, 6),
    Objective("rotb", "Rotatable bonds ≤ 8", 0.5, "descriptor", "rotb", "down", 8, 14),
    Objective("qed", "QED", 1.0, "descriptor", "qed", "up", 0.4, 0.8),
    Objective("sa", "Synthetic accessibility", 0.5, "descriptor", "sa_score", "down", 3.5, 6.5),
    Objective("route", "Route confidence", 2.0, "route", "route_confidence", "up", 0.4, 0.85),
    Objective("cost", "Reagent cost", 0.5, "route", "cost", "down", 60, 200),
    Objective("lead", "Lead time", 0.5, "route", "lead_days", "down", 7, 21),
]

# For real compounds pulled from the literature: no enumerated route, so
# synthesisability leans on the Ertl SA score instead of route confidence.
CHEMBL_PROFILE: list[Objective] = [
    Objective("potency", "Potency pIC50", 3.0, "predicted", "pIC50", "up", 6.0, 9.0),
    Objective("mw", "MW ≤ 500", 1.0, "descriptor", "mw", "down", 500, 650),
    Objective("clogp", "cLogP 1–5", 1.0, "descriptor", "clogp", "window", 1.0, 5.0, soft=2.0),
    Objective("tpsa", "TPSA 40–140", 0.75, "descriptor", "tpsa", "window", 40, 140, soft=40),
    Objective("hbd", "HBD ≤ 3", 0.5, "descriptor", "hbd", "down", 3, 6),
    Objective("rotb", "Rotatable bonds ≤ 10", 0.5, "descriptor", "rotb", "down", 10, 16),
    Objective("qed", "QED", 1.0, "descriptor", "qed", "up", 0.3, 0.7),
    Objective("sa", "Synthetic accessibility", 1.0, "descriptor", "sa_score", "down", 3.5, 6.0),
]

# Docking as a ninth objective. NOT enabled by default: turn it on for a target
# only after `dock_validate.py` shows the score ranks known compounds on that
# structure. An unvalidated docking score is a confident-looking random number.
DOCKING_OBJECTIVE = Objective(
    "docking", "Docking score (kcal/mol)", 2.0, "docking", "docking_score",
    "down", -10.0, -6.0)   # -10 or better scores 1.0, -6 or worse scores 0.0

CHEMBL_DOCKING_PROFILE: list[Objective] = CHEMBL_PROFILE + [DOCKING_OBJECTIVE]

# For a retrospective replay scored on potency: the objective must match the
# question being asked. Judging a loop on potency while it optimises drug-likeness
# marks it down for obeying its instructions. Light property guards remain so it
# cannot simply chase the heaviest molecule in the set.
POTENCY_LED_PROFILE: list[Objective] = [
    Objective("potency", "Potency (p-scale)", 8.0, "predicted", "pIC50", "up", 5.0, 9.0),
    Objective("mw", "MW ≤ 650", 0.5, "descriptor", "mw", "down", 650, 850),
    Objective("sa", "Synthetic accessibility", 0.5, "descriptor", "sa_score", "down", 5.0, 7.5),
]

PROFILES = {
    "oral_lead_opt": DEFAULT_PROFILE,
    "potency_led": POTENCY_LED_PROFILE,
    "chembl": CHEMBL_PROFILE,
    "chembl_with_docking": CHEMBL_DOCKING_PROFILE,
}


def score(c: Compound, profile: list[Objective] | None = None) -> float:
    profile = profile or DEFAULT_PROFILE
    parts, num, den = {}, 0.0, 0.0
    for o in profile:
        v = max(0.0, min(1.0, o.value(c)))
        parts[o.key] = round(v, 3)
        if o.critical and v == 0.0:
            c.objectives, c.score = parts, 0.0
            return 0.0
        num += o.weight * v
        den += o.weight
    c.objectives = parts
    c.score = round(num / den, 4)
    return c.score


def score_all(lib: list[Compound], profile: list[Objective] | None = None) -> None:
    for c in lib:
        score(c, profile)


def profile_as_dicts(profile: list[Objective] | None = None) -> list[dict]:
    return [asdict(o) for o in (profile or DEFAULT_PROFILE)]
