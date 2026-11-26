"""What an answer cost, reported next to the answer.

Every number this platform produces has a price, and until now that price appeared only
in prose: "a few hours of CPU", "96× less compute". Those are the strongest claims we
have for a customer who cannot afford a physics stack, and they were footnotes.

THE RULE THIS MODULE EXISTS TO ENFORCE
--------------------------------------
Time is measured. Money is assumed. They must never be printed as if they were the same
kind of number.

A wall-clock second is a fact: a clock was read twice. A dollar is a fact multiplied by a
rental rate that depends on who is buying, where, and when -- a rate this process did not
observe. So:

  * `measured` always carries seconds, and is always populated.
  * `assumed` carries the rate and where it came from, or is absent.
  * `estimate` exists only when a rate was supplied. **No rate means no money**, not a
    default of zero -- a silent zero is the most expensive kind of wrong number, because
    it reads as "free" rather than as "unknown".

THE DENOMINATOR TRAP
--------------------
"Cost per compound scored" is the flattering denominator: scoring is vectorised and the
library is large, so it is always a fraction of a cent, and it would be true and useless.
The number that decides anything is the cost of the whole run that produced a verdict,
because that is what a buyer actually spends to learn something. Both are reported, and
`headline_denominator` says which one to quote.

THE COMPARISON THAT MATTERS
---------------------------
Compute cost is meaningless in the abstract and obvious in chemistry units. A run that
costs less than one thousandth of a single synthesised compound is not "cheap compute";
it is a decision that was free compared with the decision it informs. The anchors below
are published per-compound costs, carried as a RANGE with sources, because the honest
state of that literature is a range.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Rate:
    """A price per hour. An assumption, and it travels with its provenance."""
    usd_per_hour: float
    label: str
    source: str
    gpu: bool = False       # does this rate rent a GPU the run may not have used?


# Rented compute. The GPU figures are from the Runpod catalogue read recorded in
# ROADMAP.md on 2026-11-24 (Community cloud); they move, and a stale rate quoted as
# current is exactly the failure this module is built to prevent, hence the date in
# every source string.
RATES = {
    "local": Rate(0.0, "this machine",
                  "hardware already owned: the marginal cash cost of a run is zero"),
    "cpu_cloud": Rate(0.05, "a rented CPU core-hour",
                      "order-of-magnitude placeholder; replace with your own invoice"),
    "rtx4090": Rate(0.34, "RTX 4090, Community cloud",
                    "Runpod catalogue read, ROADMAP.md, 2026-11-24", gpu=True),
    "a100_80": Rate(1.19, "A100 80GB PCIe, Community cloud",
                    "Runpod catalogue read, ROADMAP.md, 2026-11-24", gpu=True),
    "h100_80": Rate(2.69, "H100 80GB HBM3, Community cloud",
                    "Runpod catalogue read, ROADMAP.md, 2026-11-24", gpu=True),
}

# What one compound costs to make and test. A RANGE, because that is what the published
# evidence supports -- two self-reported project figures, not an industry benchmark.
# Quoting a single number here would manufacture a precision the sources do not have.
COMPOUND_COST_USD = (565.0, 1000.0)
COMPOUND_COST_SOURCES = (
    "€565 per compound including assay, single project, ChemRxiv 2025-01 "
    "(DMTA cycle challenges, a unified AI-driven system)",
    "~US$1,000 per compound for custom synthesis, arXiv 2026 (D2ML, "
    "direct-to-biology screening)",
)


@dataclass
class Step:
    name: str
    wall_s: float
    cpu_s: float


@dataclass
class Ledger:
    """Accumulates what a run spent. Use it as a context manager per step.

        led = Ledger()
        with led.step("fit the model"):
            ...
    """
    steps: list[Step] = field(default_factory=list)
    units: dict[str, int] = field(default_factory=dict)

    def step(self, name: str) -> "_StepTimer":
        return _StepTimer(self, name)

    def count(self, **units: int) -> None:
        """Record denominators: compounds scored, models fitted, verdicts produced."""
        self.units.update(units)

    @property
    def wall_s(self) -> float:
        return sum(s.wall_s for s in self.steps)

    @property
    def cpu_s(self) -> float:
        return sum(s.cpu_s for s in self.steps)

    # -------------------------------------------------------------- reporting

    def as_dict(self, rate: Rate | str | None = None) -> dict:
        """Measured time always; money only if a rate was supplied.

        Cost uses WALL time, not CPU time, because a rented machine bills by the clock
        whether or not every core was busy. CPU time is reported beside it so that a run
        which is parallel, or which is mostly waiting on the network, is visible as such
        rather than hidden inside one number.
        """
        out: dict = {
            "measured": {
                "wall_s": round(self.wall_s, 3),
                "cpu_s": round(self.cpu_s, 3),
                "steps": [{"name": s.name, "wall_s": round(s.wall_s, 3),
                           "cpu_s": round(s.cpu_s, 3)} for s in self.steps],
            },
            "units": dict(self.units),
            "estimate": None,
            "assumed": None,
        }
        if rate is None:
            out["note"] = ("No rate supplied, so no cost is reported. Time is measured; "
                           "money is an assumption and is not invented here.")
            return out

        r = RATES[rate] if isinstance(rate, str) else rate
        usd = self.wall_s / 3600.0 * r.usd_per_hour
        out["assumed"] = {"usd_per_hour": r.usd_per_hour, "label": r.label,
                          "source": r.source}
        out["estimate"] = {
            "usd": usd,
            "basis": "wall-clock seconds × the assumed hourly rate",
            "per_unit": {k: (usd / v if v else None) for k, v in self.units.items()},
        }
        out["compounds_not_made"] = compound_equivalent(usd)
        if r.gpu:
            # The platform's own result is that the random forest ties Chemprop at a
            # fraction of the compute, and the forest never touches a GPU. Pricing a
            # CPU run at a GPU rate therefore overstates the bill -- in the direction
            # that makes our own cost argument look WORSE, which is the safe direction,
            # but a reader should not have to work that out.
            out["warning"] = (
                f"Priced at a GPU rate ({r.label}), but this run used no GPU. The figure "
                f"is what renting that machine would have cost for the elapsed time, not "
                f"what the work required. For CPU-only work the 'cpu_cloud' or 'local' "
                f"rate is the honest one.")
        return out

    def summary(self, rate: Rate | str | None = None,
                headline_denominator: str = "verdicts") -> str:
        """One paragraph for a human, with the assumption stated in it, not under it."""
        d = self.as_dict(rate)
        m = d["measured"]
        n = len(self.steps)
        line = (f"{m['wall_s']:.1f} s wall / {m['cpu_s']:.1f} s CPU"
                f" across {n} step{'' if n == 1 else 's'}")
        if d["estimate"] is None:
            return line + ". No hourly rate was supplied, so no cost is shown."

        r, e = d["assumed"], d["estimate"]
        line += (f". At {r['usd_per_hour']:.2f} USD/hour ({r['label']}) that is an"
                 f" estimated {_money(e['usd'])}")
        if d.get("warning"):
            line += f" [{d['warning'].split('.')[0]}]"
        per = e["per_unit"].get(headline_denominator)
        if per is not None:
            line += f", {_money(per)} per {headline_denominator.rstrip('s')}"
        lo, hi = d["compounds_not_made"]
        lo_c, hi_c = min(COMPOUND_COST_USD), max(COMPOUND_COST_USD)
        if hi <= 0:
            line += (f". One compound made and tested costs ${lo_c:,.0f}–${hi_c:,.0f};"
                     " this run added nothing to that, because the machine was already"
                     " bought.")
        else:
            line += (f". One compound made and tested costs ${lo_c:,.0f}–${hi_c:,.0f},"
                     f" so this run cost {_fraction(lo, hi)} of one compound.")
        return line


def compound_equivalent(usd: float) -> tuple[float, float]:
    """This run, priced in compounds not made. Returns (low, high) fractions.

    The low end divides by the HIGH per-compound cost, so the pair is ordered the way a
    reader expects: the smallest claim first. Reporting only the flattering end of this
    range would be the same error the per-compound denominator makes.
    """
    hi_cost, lo_cost = max(COMPOUND_COST_USD), min(COMPOUND_COST_USD)
    return (usd / hi_cost, usd / lo_cost)


def _money(usd: float) -> str:
    if usd == 0:
        return "$0.00 (hardware already owned)"
    if usd < 0.01:
        return f"${usd:.5f}"
    return f"${usd:,.2f}"


def _fraction(lo: float, hi: float) -> str:
    """Express a fraction of one compound the way a chemist would say it."""
    if hi <= 0:
        return "nothing measurable"
    if hi < 0.001:
        return f"under one thousandth ({lo:.2g}–{hi:.2g})"
    return f"{lo:.3g}–{hi:.3g}"


class _StepTimer:
    def __init__(self, ledger: Ledger, name: str):
        self.ledger, self.name = ledger, name

    def __enter__(self):
        self._w, self._c = time.perf_counter(), time.process_time()
        return self

    def __exit__(self, *exc):
        self.ledger.steps.append(Step(
            self.name,
            time.perf_counter() - self._w,
            time.process_time() - self._c,
        ))
        return False
