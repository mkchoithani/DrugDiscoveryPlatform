"""Hit finding on LIT-PCBA, the benchmark DUD-E's bias made necessary.

[RESULTS.md](RESULTS.md) 12 ran docking and pharmacophore maps on DUD-E and found the
result could not be trusted as a comparison: a 2D fingerprint model reached AUC 0.996
and the arithmetic ceiling on enrichment, because DUD-E's actives are literature
compounds and its decoys are ZINC catalogue molecules. A fingerprint separates those
two populations without knowing anything about the target.

LIT-PCBA (Tran-Nguyen, Jacquemard and Rognan, *J. Chem. Inf. Model.* 2020, 60, 4263)
exists for that reason:

  - **Actives and inactives are both experimental.** They come from the same PubChem
    confirmatory bioassays, so both populations were actually tested against the
    target. An inactive is a measured inactive, not an untested catalogue molecule.
  - **AVE-unbiased splits.** Asymmetric validation embedding removes the
    train/validation redundancy that lets a model win by recognising analogues.
  - **Realistic class imbalance.** MAPK1: 308 actives against 61,567 inactives — a
    0.5% hit rate, against DUD-E's 10%.

The protocol here is the archive's own: train on the T split, report on the V split.
The pharmacophore and the descriptor baselines need no training, so they see no
advantage from that; the 2D forest does, and the question is whether the AVE unbiasing
has taken its DUD-E advantage away.

    python3 litpcba_benchmark.py <litpcba_target_dir> <field.npz> [max_inactives]
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

from dmta.chem import descriptors, inchikey
from dmta.chembl import murcko
from dmta.hitfinding import hit_finding_report, load_ism
from dmta.library import Compound
from dmta.models import PropertyModel
from dmta.pharmacophore import load_field, score_compound


def as_compound(smiles: str, cid: str, index: int, active: bool) -> Compound | None:
    try:
        desc = descriptors(smiles)
    except Exception:
        return None
    c = Compound(
        id=f"{cid}-{index}", smiles=smiles, inchikey=inchikey(smiles),
        scaffold_id=murcko(smiles), intermediate_key=cid, intermediate_smiles="",
        route=[], desc=desc, alerts=[], overall_yield=1.0, route_confidence=None,
        cost=0.0, lead_days=0,
    )
    c.truth = 1.0 if active else 0.0
    return c


def read_split(directory: Path, name: str, active: bool, limit: int | None,
               seed: int = 0) -> list[Compound]:
    rows = load_ism(directory / name)
    if limit is not None and limit < len(rows):
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(rows), size=limit, replace=False)
        rows = [rows[i] for i in sorted(idx)]
    out = []
    for i, (smiles, cid) in enumerate(rows):
        c = as_compound(smiles, cid, i, active)
        if c is not None:
            out.append(c)
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    directory = Path(sys.argv[1])
    field_path = sys.argv[2]
    max_inactives = (int(sys.argv[3]) if len(sys.argv) > 3
                     and not sys.argv[3].startswith("-") else None)
    want_dock = "--dock" in sys.argv
    dock_pdb = "4QTE"

    target = directory.name
    print(f"LIT-PCBA target {target} (AVE-unbiased split)")
    t0 = time.time()
    active_t = read_split(directory, "active_T.smi", True, None)
    active_v = read_split(directory, "active_V.smi", True, None)
    inactive_t = read_split(directory, "inactive_T.smi", False, max_inactives)
    inactive_v = read_split(directory, "inactive_V.smi", False, max_inactives)
    print(f"  train: {len(active_t)} actives / {len(inactive_t)} inactives")
    print(f"  test:  {len(active_v)} actives / {len(inactive_v)} inactives")
    if max_inactives:
        print(f"  (inactives capped at {max_inactives} per split for runtime; "
              "stated because it raises the apparent hit rate)")
    test = active_v + inactive_v
    labels = np.array([c.truth == 1.0 for c in test])
    print(f"  test hit rate: {labels.mean():.2%}   parsed in {time.time()-t0:.0f}s\n")
    if labels.sum() < 20:
        sys.exit("fewer than 20 actives in the test split — nothing to measure")

    rows = []

    field = load_field(field_path)
    print(f"pharmacophore query: {len(field.spots)} hot spots "
          f"({field.summary()['per_feature']})")
    t0 = time.time()
    scores = []
    for i, c in enumerate(test):
        r = score_compound(field, c.smiles, n_conformers=3)
        scores.append(r["score"] if r["score"] is not None else np.nan)
        if (i + 1) % 2000 == 0:
            rate = (i + 1) / (time.time() - t0)
            print(f"  {i+1}/{len(test)} at {rate:.1f}/s, "
                  f"~{(len(test)-i-1)/rate/60:.0f} min left", flush=True)
    phc4 = np.array(scores, float)
    print(f"  scored in {time.time()-t0:.0f}s "
          f"({int(np.isnan(phc4).sum())} unembeddable)\n")
    rows.append(hit_finding_report(phc4, labels, f"pharmacophore ({target})"))

    rows.append(hit_finding_report(
        np.array([c.desc["heavy_atoms"] for c in test], float), labels, "heavy atom count"))
    rows.append(hit_finding_report(
        np.array([c.desc["qed"] for c in test], float), labels, "QED"))

    train = active_t + inactive_t
    for c in train:
        c.measured = c.truth
    try:
        mu, _ = PropertyModel(n_estimators=300).fit(train).predict(test)
    finally:
        for c in train:
            c.measured = None
    rows.append(hit_finding_report(np.asarray(mu, float), labels, "2D forest"))

    if want_dock:
        from dmta.docking import apply_scores, dock, engine_available
        from dmta.structure import prepare
        if not engine_available():
            print("smina not found — skipping the docking arm\n")
        else:
            print(f"docking {len(test)} compounds against {dock_pdb} "
                  f"(~7s each, so hours)...", flush=True)
            t0 = time.time()
            result = dock(test, prepare(dock_pdb), verbose=True)
            n_applied = apply_scores(test, result)
            print(f"  docked {n_applied}/{len(test)} in {time.time()-t0:.0f}s")
            ds = np.array([c.docking_score if c.docking_score is not None else np.nan
                           for c in test], float)
            if np.isfinite(ds).sum() >= 50:
                # lower kcal/mol is a better fit, so negate for higher = better
                rows.append(hit_finding_report(-ds, labels,
                                               f"docking (smina, {dock_pdb})"))
            else:
                print("  too few compounds docked to report an arm")

    header = f"{'method':<30}{'n':>7}{'AUC':>7}{'EF 1%':>8}{'EF 5%':>8}{'BEDROC':>9}"
    print(header); print("-" * len(header))
    for r in rows:
        print(f"{r['method']:<30}{r['n']:>7}{r['roc_auc']:>7.3f}"
              f"{(r['ef_1pct'] if r['ef_1pct'] is not None else float('nan')):>8.2f}"
              f"{(r['ef_5pct'] if r['ef_5pct'] is not None else float('nan')):>8.2f}"
              f"{r['bedroc_20']:>9.3f}")

    # With few actives, an enrichment factor is a small-count statistic. Say how small.
    n_act = int(labels.sum())
    k1 = max(1, int(round(len(test) * 0.01)))
    print(f"\n  the top 1% is {k1} compounds; at chance it holds "
          f"{k1 * labels.mean():.1f} of the {n_act} actives, so every EF 1% here is "
          f"built from a handful of compounds and is correspondingly noisy")

    phc, size, forest = rows[0], rows[1], rows[3]
    print("\n" + "=" * 69)
    print("  Did the unbiasing work? Compare the 2D forest against DUD-E's "
          "AUC 0.996 / EF1% 10.23:")
    print(f"    2D forest here: AUC {forest['roc_auc']:.3f}, "
          f"EF1% {forest['ef_1pct']}, BEDROC {forest['bedroc_20']:.3f}")
    if forest["roc_auc"] < 0.9:
        print("    -> yes. The fingerprint advantage is gone, so this benchmark can")
        print("       actually compare a structure-based method against a 2D one.")
    else:
        print("    -> no. The 2D model still wins easily, so this benchmark cannot")
        print("       settle the question either.")

    # An enrichment ratio is not evidence until it is converted back into a count of
    # molecules. This verdict declared the pharmacophore a winner on EF1% 2.60 vs
    # 0.00, which was two actives against zero out of 153 compounds -- something that
    # happens 18% of the time with no method at all. Significance is now part of the
    # verdict rather than a caveat printed above it.
    from math import comb

    def tail_probability(ef: float | None) -> tuple[int, float]:
        """Actives implied by an enrichment, and P(at least that many) by chance."""
        if not ef:
            return 0, 1.0
        found = int(round(ef * k1 * labels.mean()))
        if found <= 0:
            return 0, 1.0
        total, act = len(test), int(labels.sum())
        p = sum(comb(act, i) * comb(total - act, k1 - i) / comb(total, k1)
                for i in range(found, min(act, k1) + 1))
        return found, float(p)

    print()
    print(f"  {'method':<24}{'EF 1%':>7}{'actives found':>15}{'P by chance':>14}")
    for r in rows:
        found, p = tail_probability(r["ef_1pct"])
        print(f"  {r['method']:<24}{(r['ef_1pct'] or 0):>7.2f}{found:>15}{p:>14.4f}")

    phc_found, phc_p = tail_probability(phc["ef_1pct"])
    print()
    if phc_p < 0.05 and (phc["ef_1pct"] or 0) > (size["ef_1pct"] or 0):
        print(f"  The pharmacophore finds {phc_found} actives in the top 1% "
              f"(p = {phc_p:.3f}), beating")
        print(f"  molecular size. On an unbiased benchmark that is a real result.")
    else:
        print(f"  The pharmacophore finds {phc_found} actives in the top 1%, which "
              f"chance alone")
        print(f"  produces {phc_p:.0%} of the time. Its EF1% of "
              f"{phc['ef_1pct']} is not evidence of anything.")
        print("  Three benchmarks now agree, and this one is unbiased by construction.")


if __name__ == "__main__":
    main()
