"""Can a generator rediscover compounds a real programme went on to make?

The known-answer test for de novo design. A generator is easy to admire and hard to
judge: it emits plausible molecules whether or not it is any good, and nobody can
measure a molecule that was never made. So this asks a question with a known answer.

    Take a real programme. Split it by date. Cut up ONLY the earlier compounds.
    Recombine the pieces. Then ask whether the generated space contains the
    compounds the programme actually made later — and whether it preferentially
    contains the POTENT ones.

That is the 2021 Chem. Commun. claim in retrospective form, and it is falsifiable.

Three controls, because "the generator produced a known drug" is a weak claim on its
own:

  1. **Potency stratification.** Recovery of the later top-decile compounds against
     the rest. A generator that reaches the good ones no more often than the bad ones
     has found chemistry, not activity.
  2. **Foreign fragments.** The same pipeline fed fragments from a different target.
     If that recovers just as much, the generator is exploiting generic medicinal
     chemistry and the programme's own fragments contributed nothing.
  3. **The fixed library.** What `library.py`'s 924 enumerated compounds recover,
     which is the capability this module is meant to replace.

    python3 denovo_validate.py [target] [n_candidates] [foreign_target]
"""
from __future__ import annotations

import sys
import time

import numpy as np

from dmta.chembl import load_compounds
from dmta.generate import fragment_library, generate, recovery


def split_by_date(compounds, quantile: float = 0.75):
    dated = [c for c in compounds if getattr(c, "year", None)]
    if len(dated) < 100:
        return [], [], None
    dated.sort(key=lambda c: c.year)
    cut = dated[int(len(dated) * quantile)].year
    return ([c for c in dated if c.year < cut],
            [c for c in dated if c.year >= cut], cut)


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "CHEMBL203"
    n_candidates = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    foreign_target = sys.argv[3] if len(sys.argv) > 3 else "CHEMBL233"

    compounds = load_compounds(target, "IC50")
    early, late, cut = split_by_date(compounds)
    if not early:
        sys.exit(f"{target} has too few dated compounds for a time split")
    print(f"{target}: {len(compounds)} compounds, time split at {cut}")
    print(f"  {len(early)} earlier (fragmented), {len(late)} later (held out)\n")

    truth = np.array([c.truth for c in late], dtype=float)
    threshold = float(np.percentile(truth, 90))
    potent = truth >= threshold
    print(f'held-out "potent" = its own top decile, pIC50 >= {threshold:.2f} '
          f"({int(potent.sum())} compounds)")

    # Match the generated size range to what the programme actually made, so recovery
    # is not lost to a filter rather than to the generator.
    mws = np.array([c.desc["mw"] for c in late], dtype=float)
    mw_range = (float(np.percentile(mws, 2)), float(np.percentile(mws, 98)))
    print(f"generated MW window matched to held-out compounds: "
          f"{mw_range[0]:.0f}-{mw_range[1]:.0f}\n")

    known = {c.inchikey for c in early if c.inchikey}

    def build(label, source, seed=0):
        t0 = time.time()
        frags = fragment_library(source)
        gen = generate(frags["fragments"], n=n_candidates, seed=seed,
                       mw_range=mw_range, exclude_inchikeys=known, max_depth=5)
        print(f"  {label}: {frags['n_fragments_kept']} fragments -> "
              f"{len(gen)} molecules in {time.time() - t0:.0f}s", flush=True)
        return gen

    print("building:")
    generated = build(f"{target} fragments", early, seed=0)
    foreign_source = load_compounds(foreign_target, "Ki")
    if not foreign_source:
        foreign_source = load_compounds(foreign_target, "IC50")
    foreign = build(f"{foreign_target} fragments (control)",
                    foreign_source[:len(early)], seed=0)

    from dmta.campaign import fresh_library
    fixed = fresh_library()
    print(f"  fixed enumerated library: {len(fixed)} compounds\n")

    rows = []
    for label, gen in ((f"generated from {target} (early)", generated),
                       (f"generated from {foreign_target} (control)", foreign),
                       ("fixed enumerated library", fixed)):
        if not gen:
            continue
        r = recovery(gen, late)
        rows.append((label, r))

    header = (f"{'source of candidates':<36}{'n':>7}{'exact':>8}{'near':>8}"
              f"{'median sim':>12}")
    print(header); print("-" * len(header))
    for label, r in rows:
        print(f"{label:<36}{r['n_generated']:>7}{r['exact_recovered']:>8}"
              f"{r['near_recovered']:>8}{r['median_nearest_similarity']:>12.3f}")

    # --- control 1: does it reach the potent ones preferentially? --------------
    print(f"\npotency stratification, {rows[0][0]}:")
    r = rows[0][1]
    near = np.array(r["_near"])
    nearest = np.array(r["_nearest"])
    print(f"  {'group':<22}{'n':>6}{'near-recovered':>16}{'median sim':>12}")
    for name, mask in (("top decile", potent), ("the rest", ~potent)):
        if mask.sum() == 0:
            continue
        print(f"  {name:<22}{int(mask.sum()):>6}"
              f"{f'{near[mask].sum()} ({near[mask].mean():.0%})':>16}"
              f"{np.median(nearest[mask]):>12.3f}")

    # Is any difference real, or 10 compounds of noise?
    if potent.sum() and (~potent).sum():
        a, b = near[potent].mean(), near[~potent].mean()
        pooled = near.mean()
        rng = np.random.default_rng(0)
        draws = rng.binomial(int(potent.sum()), pooled, 20000) / potent.sum()
        draws2 = rng.binomial(int((~potent).sum()), pooled, 20000) / (~potent).sum()
        p = float((np.abs(draws - draws2) >= abs(a - b)).mean())
        print(f"\n  difference {a - b:+.3f}; a gap that large arises {p:.0%} of the "
              f"time by chance")

    # --- verdict --------------------------------------------------------------
    own, control = rows[0][1], (rows[1][1] if len(rows) > 1 else None)
    print("\n" + "=" * 72)
    if own["near_recovered"] == 0 and own["exact_recovered"] == 0:
        print("  VERDICT: FAILS. The generator reached none of the compounds the")
        print("  programme actually made, even loosely. It is producing molecules,")
        print("  not designs.")
    elif control and own["near_rate"] <= control["near_rate"]:
        print(f"  VERDICT: FAILS the control. Fragments from {foreign_target} recover")
        print(f"  as much ({control['near_rate']:.1%}) as the programme's own")
        print(f"  ({own['near_rate']:.1%}). The programme's chemistry added nothing.")
    else:
        print(f"  VERDICT: the generator reaches {own['near_rate']:.1%} of held-out")
        print(f"  compounds within 0.7 Tanimoto against {control['near_rate']:.1%} for")
        print("  foreign fragments. Read the potency stratification above before")
        print("  calling this useful: reaching the chemistry is not reaching the")
        print("  activity.")


if __name__ == "__main__":
    main()
