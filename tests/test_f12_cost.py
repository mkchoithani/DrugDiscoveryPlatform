"""Tests for cost reporting — F12.

A cost number is easy to produce and easy to make dishonest, and every failure mode here
is one that would still render a plausible-looking page. These pin the constraints rather
than the arithmetic:

  * no rate must mean NO money, never a silent zero;
  * the assumed rate must travel with its source, so a stale price cannot be quoted as
    current;
  * the flattering denominator must not become the headline;
  * a GPU rate on a run that used no GPU must say so.

In their own file, like test_f9_profile.py, because sibling sessions are editing
test_platform.py.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dmta.cost import (COMPOUND_COST_SOURCES, COMPOUND_COST_USD,   # noqa: E402
                       RATES, Ledger, Rate, compound_equivalent)


def _ledger() -> Ledger:
    led = Ledger()
    with led.step("a step"):
        pass
    led.count(verdicts=1, compounds_scored=9000)
    return led


def test_no_rate_means_no_money_rather_than_a_silent_zero():
    """The most expensive wrong number is a zero that means "unknown".

    A report showing $0.00 because nobody configured a price reads as "free" and will be
    repeated by whoever saw it. Absence has to be visible.
    """
    d = _ledger().as_dict(None)
    assert d["estimate"] is None
    assert d["assumed"] is None
    assert "no cost is reported" in d["note"].lower()
    assert "0" not in d["note"][:40], "the note must not lead with a number"


def test_time_is_always_measured_even_without_a_rate():
    """Seconds are a fact and do not depend on anyone's price list."""
    d = _ledger().as_dict(None)
    assert d["measured"]["wall_s"] >= 0.0
    assert d["measured"]["cpu_s"] >= 0.0
    assert len(d["measured"]["steps"]) == 1


def test_every_rate_carries_a_source():
    """A price with no provenance cannot be audited or spotted as stale.

    The GPU rates came from a catalogue read on a particular day; without the date in the
    source string, a figure that has since moved still looks authoritative.
    """
    for name, r in RATES.items():
        assert r.source and len(r.source) > 20, f"{name} has no usable source"
        assert isinstance(r.usd_per_hour, float)


def test_the_estimate_is_labelled_as_resting_on_an_assumption():
    d = _ledger().as_dict("a100_80")
    assert d["assumed"]["source"]
    assert "assumed" in d["estimate"]["basis"] or "×" in d["estimate"]["basis"]
    # the measured block must not acquire a money field
    assert not any("usd" in k for k in d["measured"])


def test_a_gpu_rate_on_a_cpu_run_says_so():
    """Our own argument is that this needs no GPU.

    Pricing a CPU run at a GPU rate overstates the bill. That direction is the safe one --
    it makes our cost claim look worse, not better -- but a reader should not have to
    deduce it.
    """
    assert "no GPU" in _ledger().as_dict("h100_80")["warning"]
    assert "warning" not in _ledger().as_dict("cpu_cloud")
    assert "warning" not in _ledger().as_dict("local")


def test_cost_uses_wall_time_not_cpu_time():
    """A rented machine bills by the clock, busy or not.

    Using CPU seconds would understate a run that waits on the network and overstate a
    parallel one. Both are reported; only wall-clock is charged.
    """
    led = Ledger()
    from dmta.cost import Step
    led.steps.append(Step("busy", wall_s=3600.0, cpu_s=14400.0))   # 4 cores, one hour
    led.count(verdicts=1)
    d = led.as_dict(Rate(1.0, "one dollar an hour", "a test fixture, not a real price"))
    assert d["estimate"]["usd"] == pytest.approx(1.0), (
        "one wall-clock hour at $1/hour must cost $1, regardless of core count")


def test_the_flattering_denominator_is_reported_but_is_not_the_headline():
    """Cost per compound scored is true and decides nothing.

    Scoring is vectorised, so that figure is always a fraction of a cent and would make
    any platform look good. The headline denominator is the verdict: one run answers one
    question, and that is what is being bought.
    """
    led = _ledger()
    d = led.as_dict("a100_80")
    assert "compounds_scored" in d["estimate"]["per_unit"], "it must still be reported"
    assert d["estimate"]["per_unit"]["compounds_scored"] < \
        d["estimate"]["per_unit"]["verdicts"], "scoring must be the cheaper denominator"
    assert "per verdict" in led.summary("a100_80", headline_denominator="verdicts")


def test_the_per_compound_anchor_is_a_range_with_sources():
    """Two single-project self-reports are not an industry benchmark.

    Collapsing them to one number would manufacture a precision the literature does not
    have, and that number would then be quoted back at us.
    """
    assert len(COMPOUND_COST_USD) == 2
    assert COMPOUND_COST_USD[0] != COMPOUND_COST_USD[1]
    assert len(COMPOUND_COST_SOURCES) == 2
    for src in COMPOUND_COST_SOURCES:
        assert any(ch.isdigit() for ch in src), "a source without a year is not a source"


def test_the_compound_equivalent_reports_the_modest_end_first():
    """Quoting only the flattering end of a range is the same error twice."""
    lo, hi = compound_equivalent(100.0)
    assert lo < hi
    assert lo == pytest.approx(100.0 / max(COMPOUND_COST_USD))
    assert hi == pytest.approx(100.0 / min(COMPOUND_COST_USD))


def test_the_local_rate_is_zero_and_says_why():
    """Zero here means something different from the zero this module refuses to print.

    On hardware already owned the marginal cash cost really is nil, and the label has to
    carry that reason so it cannot be confused with "no rate configured".
    """
    r = RATES["local"]
    assert r.usd_per_hour == 0.0
    assert "already owned" in r.source
    assert "hardware already owned" in _ledger().summary("local")
