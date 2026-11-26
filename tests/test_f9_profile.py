"""Tests for the antibacterial profile — section 19.

In its own file rather than appended to test_platform.py, because two other sessions were
editing that file when these were written.

What these pin is not "the profile scores antibacterials well". Of course it does; it was
fitted to them. What they pin are the three properties that make the fit honest, each of
which a plausible future edit would silently break:

  * the profile must still RANK — a profile that scores a whole class 1.000 is useless,
    and that is exactly how the first hand-written attempt failed;
  * QED must stay out, because it is fitted to oral drugs and would re-import the bias
    the profile exists to remove;
  * the rotatable-bond term must stay STRICTER than the default, because the one piece of
    the literature in here is a Gram-negative entry guideline, and a future "let's relax
    this" would quietly delete it.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dmta.library import Compound                                      # noqa: E402
from dmta.scoring import (ANTIBACTERIAL_PROFILE, DEFAULT_PROFILE,      # noqa: E402
                          Objective, PROFILES, score)


def _compound(**desc) -> Compound:
    """A compound carrying only descriptors — enough for the property objectives."""
    full = {"mw": 400.0, "clogp": 0.5, "tpsa": 140.0, "hbd": 3, "hba": 6, "rotb": 4,
            "rings": 2, "aromatic_rings": 1, "fsp3": 0.4, "heavy_atoms": 28,
            "formal_charge": 0, "qed": 0.3, "sa_score": 4.0, "stereocenters": 2}
    full.update(desc)
    c = Compound(id="t", smiles="C", inchikey="K", scaffold_id="s",
                 intermediate_key="s", intermediate_smiles="", route=[], desc=full,
                 alerts=[], overall_yield=1.0, route_confidence=1.0, cost=0.0,
                 lead_days=0)
    c.predicted = 7.0
    return c


def test_the_antibacterial_profile_is_registered():
    assert PROFILES["antibacterial"] is ANTIBACTERIAL_PROFILE


def test_qed_is_absent_from_the_antibacterial_profile():
    """QED is trained on oral drugs.

    Scoring a class that is deliberately not oral-drug-like by how much it resembles
    oral drugs assumes the conclusion. It was dropped rather than down-weighted, and a
    re-added QED term would quietly restore the bias section 19 measured.
    """
    assert "qed" not in {o.key for o in ANTIBACTERIAL_PROFILE}
    assert "qed" in {o.key for o in DEFAULT_PROFILE}, "the default is the comparison"


def test_the_rotatable_bond_term_is_stricter_than_the_default():
    """The eNTRy guideline is the one term here taken from the literature.

    Everything else in this profile is wider than the default, so a future edit that
    relaxes the whole thing "for consistency" would delete the only part that encodes
    how a compound gets into a Gram-negative cell.
    """
    ab = next(o for o in ANTIBACTERIAL_PROFILE if o.key == "rotb")
    df = next(o for o in DEFAULT_PROFILE if o.key == "rotb")
    assert ab.lo < df.lo, f"antibacterial rotB tolerance {ab.lo} is not stricter than {df.lo}"
    assert ab.lo == 5, "the eNTRy guideline is 5 rotatable bonds"


def test_a_polar_antibacterial_scores_better_under_the_antibacterial_profile():
    """Ciprofloxacin-like properties: polar, low cLogP, high TPSA.

    This is the defect in one compound. The default profile's cLogP window starts at 1.0
    and its TPSA window ends at 120, so a compound in the middle of the marketed
    antibacterial distribution is marked down by both.
    """
    c = _compound(mw=331.3, clogp=-0.72, tpsa=74.6, hbd=2, rotb=3, qed=0.68, sa_score=2.9)
    default = score(c, DEFAULT_PROFILE)
    abx = score(c, ANTIBACTERIAL_PROFILE)
    assert abx > default, f"antibacterial profile {abx} did not beat default {default}"


def test_the_antibacterial_profile_still_ranks():
    """The failure mode that killed the first attempt, pinned.

    A hand-widened profile scored the median marketed antibacterial exactly 1.000. A
    profile that gives a whole class full marks cannot order it, and ordering is the only
    thing a profile does. Three compounds spread across the class must come out strictly
    ordered, not tied at the top.
    """
    middle = _compound(mw=400.0, clogp=0.5, tpsa=140.0, hbd=3, rotb=4, sa_score=4.0)
    edge = _compound(mw=560.0, clogp=2.6, tpsa=215.0, hbd=5, rotb=7, sa_score=5.0)
    far = _compound(mw=900.0, clogp=-6.0, tpsa=400.0, hbd=12, rotb=15, sa_score=7.5)
    s_mid = score(middle, ANTIBACTERIAL_PROFILE)
    s_edge = score(edge, ANTIBACTERIAL_PROFILE)
    s_far = score(far, ANTIBACTERIAL_PROFILE)
    assert s_mid > s_edge > s_far, f"not ordered: {s_mid}, {s_edge}, {s_far}"
    assert s_mid < 1.0, "a profile that scores its own class 1.000 has no gradient left"


def test_the_antibacterial_profile_does_not_simply_approve_of_everything():
    """Looser is not the same as class-specific.

    The first attempt lifted the control arm's median too, which meant it was not a
    repair but a lower bar. A compound far outside the antibacterial class must still be
    rejected by the antibacterial profile.
    """
    greasy = _compound(mw=520.0, clogp=6.8, tpsa=35.0, hbd=0, rotb=12, sa_score=3.0)
    assert score(greasy, ANTIBACTERIAL_PROFILE) < 0.6


def test_potency_still_carries_weight():
    """The profile must not become a pure property filter.

    A class-specific property profile that forgot to optimise activity would select
    well-behaved molecules that do nothing.
    """
    potency = next(o for o in ANTIBACTERIAL_PROFILE if o.key == "potency")
    assert potency.weight >= max(o.weight for o in ANTIBACTERIAL_PROFILE if o.key != "potency")


def test_ciprofloxacin_is_not_marked_down_by_the_antibacterial_profile():
    """The regression that caught the section 19 error, pinned with the real drug.

    Ciprofloxacin's real descriptors. An earlier fitted profile used the class's 25th-75th
    percentile windows and scored it 0.836, because its TPSA of 74.6 falls below that
    lower edge -- fitting to the middle 50% of a class excludes the edges of that same
    class. The fitted distribution cannot tell you ciprofloxacin matters; only chemistry
    can, which is why this test names a drug instead of a percentile.
    """
    cipro = _compound(mw=331.35, clogp=-0.72, tpsa=74.57, hbd=2, hba=5,
                      rotb=3, qed=0.68, sa_score=2.9)
    # The claim is about the property guards, so measure those and not the full score --
    # the full score also carries a potency term, and a helper's placeholder potency
    # would otherwise decide whether this test passes.
    guards = [o for o in ANTIBACTERIAL_PROFILE if o.key != "potency"]
    s = sum(o.weight * max(0.0, min(1.0, o.value(cipro))) for o in guards) / \
        sum(o.weight for o in guards)
    assert s > 0.9, f"ciprofloxacin scored {s:.3f} under the profile built for its own class"

    narrow = [Objective("tpsa", "p25–p75", 1.0, "descriptor", "tpsa",
                        "window", 98.0, 181.6, soft=41.8)]
    assert narrow[0].value(cipro) < 0.5, (
        "the 25-75 window that this test exists to rule out no longer penalises "
        "ciprofloxacin, so the regression it guards has changed shape")


def test_the_property_terms_alone_cannot_rank_and_that_is_the_point():
    """Section 19's sixth verdict-function error, pinned as a property of the profile.

    Saturation was measured on the property terms alone, which read 68.8% and condemned
    the correct window. The property terms are guards; potency ranks, at weight 3.0
    against the guards' 5.0. Two compounds that are both comfortably inside every guard
    must be told apart by potency and by nothing else -- if a future edit lets the guards
    dominate, the measurement that went wrong becomes the right one again.
    """
    guards = [o for o in ANTIBACTERIAL_PROFILE if o.key != "potency"]
    potency = next(o for o in ANTIBACTERIAL_PROFILE if o.key == "potency")
    assert potency.weight / sum(o.weight for o in guards) >= 0.5, (
        "potency no longer carries enough weight to rank compounds the guards all accept")

    weak, strong = _compound(), _compound()
    weak.predicted, strong.predicted = 6.0, 8.0
    assert score(strong, ANTIBACTERIAL_PROFILE) > score(weak, ANTIBACTERIAL_PROFILE)
