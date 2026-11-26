"""Do predicted pharmacophore maps add signal where docking did not?

The known-answer test, run exactly as docking was judged
([RESULTS.md](RESULTS.md) 4-5): real EGFR compounds with real measured potencies,
enrichment measured as top-decile hit rate over the library's own base rate. Docking
scored 0.37x on 2RGP and 0.00x on 4HJO and the ensemble, against 3.51x for the plain
2D model.

Two bars, and the method has to clear both to be worth building on:

  1. **enrichment > 1.5x** -- better than picking at random, by a margin.
  2. **signal on the time split**, where the 2D model collapses to rho 0.11. That
     split is the entire reason to want a structure-based method: it is the case
     where the model has never seen the chemotype. A method that only works where
     the 2D model already works has bought us nothing.

**Three baselines, because the obvious failure mode here is measuring size.** A
bigger molecule has more pharmacophore features, so it can satisfy more of the query
by having more chances -- and bigger molecules are, on average, more potent. If the
pharmacophore score does not beat heavy-atom count, it is an expensive molecular
weight. QED is included as a second cheap descriptor, and the random forest as the
incumbent.

    python3 pharmacophore_validate.py <field.npz> [n_compounds] [target]
"""
from __future__ import annotations

import sys
import time

import numpy as np

from dmta.chembl import load_compounds
from dmta.models import PropertyModel
from dmta.pharmacophore import load_field, score_compound
from dmta.validate import _spearman, scaffold_split


def enrichment(scores: np.ndarray, truth: np.ndarray, threshold: float,
               higher_is_better: bool = True) -> dict:
    """Top-decile hit rate over base rate — the same protocol docking was judged by."""
    active = truth >= threshold
    base = float(active.mean())
    k = max(1, len(scores) // 10)
    order = np.argsort(scores)[::-1] if higher_is_better else np.argsort(scores)
    top = order[:k]
    return {
        "n": len(scores),
        "spearman": _spearman(scores if higher_is_better else -scores, truth),
        "top_decile_hit_rate": round(float(active[top].mean()), 3),
        "base_hit_rate": round(base, 3),
        "enrichment": round(float(active[top].mean()) / base, 2) if base > 0 else None,
    }


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    field_path = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
    target = sys.argv[3] if len(sys.argv) > 3 else "CHEMBL203"

    field = load_field(field_path)
    print("pharmacophore query:", field.summary(), "\n")
    if len(field.spots) < 6:
        sys.exit("fewer than 6 hot spots — the query cannot constrain anything")

    compounds = load_compounds(target, "IC50")
    rng = np.random.default_rng(0)
    if limit and limit < len(compounds):
        idx = rng.choice(len(compounds), size=limit, replace=False)
        compounds = [compounds[i] for i in sorted(idx)]
    print(f"{len(compounds)} real {target} compounds with measured IC50")

    # Dataset-relative, never a hardcoded "active >= 8": on a programme whose best
    # compound is 7.71 a fixed cut reports 0% actives and looks like a failure.
    truth_all = np.array([c.truth for c in compounds], dtype=float)
    threshold = float(np.percentile(truth_all, 90))
    print(f'"active" = top decile of this dataset, pIC50 >= {threshold:.2f}\n')

    t0 = time.time()
    scored = []
    for i, c in enumerate(compounds):
        r = score_compound(field, c.smiles, n_conformers=3)
        if r["score"] is not None:
            scored.append((c, r))
        if (i + 1) % 200 == 0:
            rate = (i + 1) / (time.time() - t0)
            print(f"  {i + 1}/{len(compounds)} at {rate:.1f}/s "
                  f"({len(compounds) - i - 1:.0f} left, "
                  f"~{(len(compounds) - i - 1) / rate / 60:.1f} min)", flush=True)
    print(f"\nscored {len(scored)}/{len(compounds)} in {time.time() - t0:.0f}s "
          f"({len(compounds) - len(scored)} could not be embedded)\n")
    if len(scored) < 50:
        sys.exit("too few compounds scored to conclude anything")

    cs = [c for c, _ in scored]
    truth = np.array([c.truth for c in cs], dtype=float)
    phc4 = np.array([r["score"] for _, r in scored], dtype=float)
    matched = np.array([r["matched"] for _, r in scored], dtype=float)
    heavy = np.array([c.desc["heavy_atoms"] for c in cs], dtype=float)
    qed = np.array([c.desc["qed"] for c in cs], dtype=float)

    rows = [
        ("pharmacophore score", enrichment(phc4, truth, threshold)),
        ("pharmacophore matches", enrichment(matched, truth, threshold)),
        ("heavy atom count", enrichment(heavy, truth, threshold)),
        ("QED", enrichment(qed, truth, threshold)),
    ]

    # The incumbent, on a scaffold split, for context. It is not a like-for-like
    # comparison -- the forest trains on measured values and the pharmacophore
    # score needs none -- and it is reported for scale, not as a fair contest.
    train, test = scaffold_split(cs, test_frac=0.25)
    if len(train) >= 30 and len(test) >= 10:
        for c in train:
            c.measured = c.truth
        try:
            mu, _ = PropertyModel(n_estimators=300).fit(train).predict(test)
        finally:
            for c in train:
                c.measured = None
        t_truth = np.array([c.truth for c in test], dtype=float)
        rows.append(("random forest (scaffold split)",
                     enrichment(np.asarray(mu), t_truth, threshold)))

    header = (f"{'method':<34}{'n':>6}{'rho':>8}{'top decile':>12}"
              f"{'base':>8}{'enrichment':>12}")
    print(header); print("-" * len(header))
    for label, m in rows:
        rho = m["spearman"]
        print(f"{label:<34}{m['n']:>6}{(rho if rho is not None else float('nan')):>8.3f}"
              f"{m['top_decile_hit_rate']:>12.3f}{m['base_hit_rate']:>8.3f}"
              f"{(m['enrichment'] if m['enrichment'] is not None else float('nan')):>12.2f}")

    # Is the score just molecular size wearing a hat?
    size_rho = _spearman(phc4, heavy)
    print(f"\n  pharmacophore score vs heavy atom count: rho {size_rho:.3f}")

    # --- the test that matters: unseen chemotypes over time -------------------
    dated = [c for c in cs if getattr(c, "year", None)]
    phc4_by_id = {c.id: r["score"] for c, r in scored}
    print()
    if len(dated) > 100:
        dated.sort(key=lambda c: c.year)
        cut = dated[int(len(dated) * 0.75)].year
        early = [c for c in dated if c.year < cut]
        late = [c for c in dated if c.year >= cut]
        print(f"time split at {cut}: {len(early)} earlier, {len(late)} later")
        if len(early) >= 30 and len(late) >= 20:
            late_truth = np.array([c.truth for c in late], dtype=float)
            late_phc4 = np.array([phc4_by_id[c.id] for c in late], dtype=float)
            print(f"  {'method':<32}{'rho':>8}{'enrichment':>12}")
            m = enrichment(late_phc4, late_truth, threshold)
            print(f"  {'pharmacophore (no training)':<32}"
                  f"{(m['spearman'] or float('nan')):>8.3f}"
                  f"{(m['enrichment'] or float('nan')):>12.2f}")
            for c in early:
                c.measured = c.truth
            try:
                mu, _ = PropertyModel(n_estimators=300).fit(early).predict(late)
            finally:
                for c in early:
                    c.measured = None
            m2 = enrichment(np.asarray(mu), late_truth, threshold)
            print(f"  {'random forest (trained early)':<32}"
                  f"{(m2['spearman'] or float('nan')):>8.3f}"
                  f"{(m2['enrichment'] or float('nan')):>12.2f}")
    else:
        print("too few dated compounds for a time split")

    # --- verdict --------------------------------------------------------------
    phc4_metrics = rows[0][1]
    best_cheap = max(rows[2][1]["enrichment"] or 0, rows[3][1]["enrichment"] or 0)
    e = phc4_metrics["enrichment"] or 0
    print("\n" + "=" * 72)
    if e < 1.5:
        print(f"  VERDICT: FAILS the first bar. Enrichment {e:.2f}x, needed > 1.5x.")
        print("  Better than docking managed (0.37x) but not enough to build on.")
    elif e <= best_cheap:
        print(f"  VERDICT: FAILS. Enrichment {e:.2f}x does not beat the cheapest")
        print(f"  descriptor baseline ({best_cheap:.2f}x). This is measuring molecular")
        print("  size, not pocket complementarity.")
    else:
        print(f"  VERDICT: PASSES the first bar at {e:.2f}x, beating both the")
        print(f"  descriptor baselines ({best_cheap:.2f}x). Check the time split above")
        print("  before committing to Phase B.")


if __name__ == "__main__":
    main()
