"""Can docking tell a real EGFR binder from a molecule that has no business binding?

This is the job docking is actually known to be good at, and the one that would make
it useful as a cheap first filter over a huge library. It is a different question from
`dock_validate.py`, which asked it to rank potency inside one series — and which it
failed (RESULTS.md section 4).

The test is built so it cannot cheat:

  * Decoys come from unrelated, non-kinase targets (Factor Xa, COX-2) — real drug-like
    molecules with no reason to bind EGFR.
  * Decoys are PROPERTY-MATCHED to the actives on molecular weight and cLogP. Docking
    scores rise with molecular size, so unmatched decoys would let it "win" by simply
    preferring bigger molecules.
  * Molecular weight alone is scored as a control. If docking cannot beat the size
    baseline, it has learned nothing beyond size.

    python3 dock_screen_test.py [pdb_id] [n_actives] [decoys_per_active]
"""
from __future__ import annotations

import random
import statistics
import sys

import numpy as np

from dmta.chembl import load_compounds
from dmta.docking import apply_scores, dock, engine_available
from dmta.structure import prepare

# Non-kinase targets only: a kinase decoy could genuinely bind EGFR and would make
# the test unfairly hard. (target, standard_type, description)
DECOY_TARGETS = [
    ("CHEMBL244", "IC50", "Factor Xa (protease)"),
    ("CHEMBL205", "Ki", "carbonic anhydrase II (lyase)"),
    ("CHEMBL233", "Ki", "mu opioid receptor (GPCR)"),
    ("CHEMBL230", "IC50", "COX-2 (oxidoreductase)"),
]
ACTIVE_THRESHOLD = 8.0


def property_matched_decoys(actives, pool, per_active, mw_tol=35.0, logp_tol=0.8, seed=0):
    """For each active, pick decoys of similar size and greasiness."""
    rng = random.Random(seed)
    available = list(pool)
    rng.shuffle(available)
    chosen, used = [], set()
    unmatched = 0
    for a in actives:
        found = 0
        for d in available:
            if d.inchikey in used:
                continue
            if abs(d.desc["mw"] - a.desc["mw"]) <= mw_tol and \
               abs(d.desc["clogp"] - a.desc["clogp"]) <= logp_tol:
                chosen.append(d)
                used.add(d.inchikey)
                found += 1
                if found >= per_active:
                    break
        if found < per_active:
            unmatched += per_active - found
    return chosen, unmatched


def roc_auc(scores, labels, lower_is_better=True):
    """Probability a random active outranks a random decoy. 0.5 = useless."""
    s = np.array(scores, dtype=float)
    if lower_is_better:
        s = -s
    y = np.array(labels, dtype=bool)
    if y.all() or not y.any():
        return None
    order = s.argsort()
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    n_pos, n_neg = int(y.sum()), int((~y).sum())
    auc = (ranks[y].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return round(float(auc), 3)


def enrichment_at(scores, labels, fraction, lower_is_better=True):
    s = np.array(scores, dtype=float)
    if lower_is_better:
        s = -s
    y = np.array(labels, dtype=bool)
    k = max(1, int(round(len(s) * fraction)))
    top = np.argsort(-s)[:k]
    base = y.mean()
    return round(float(y[top].mean() / base), 2) if base > 0 else None


def main():
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "2RGP"
    n_actives = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    per_active = int(sys.argv[3]) if len(sys.argv) > 3 else 3

    if not engine_available():
        sys.exit("no docking engine found")

    egfr = load_compounds("CHEMBL203")
    actives_pool = [c for c in egfr if c.truth >= ACTIVE_THRESHOLD]
    actives = random.Random(1).sample(actives_pool, min(n_actives, len(actives_pool)))

    decoy_pool, sources = [], []
    for target, std_type, label in DECOY_TARGETS:
        try:
            loaded = load_compounds(target, std_type)
            decoy_pool += loaded
            sources.append(f"{label} n={len(loaded)}")
        except Exception:
            continue          # dataset not downloaded; the others still make a pool
    egfr_keys = {c.inchikey for c in egfr}
    decoy_pool = [c for c in decoy_pool if c.inchikey not in egfr_keys]   # no overlap

    decoys, unmatched = property_matched_decoys(actives, decoy_pool, per_active)
    print(f"{pdb_id}: {len(actives)} EGFR actives (pChEMBL >= {ACTIVE_THRESHOLD}) "
          f"vs {len(decoys)} property-matched decoys")
    for src in sources:
        print(f"    decoy source: {src}")
    print(f"  decoy pool {len(decoy_pool)} compounds, {unmatched} slots unmatched")

    def dist(name, group, key):
        values = [c.desc[key] for c in group]
        return f"  {name:<8} {key:<6} mean {statistics.mean(values):7.1f}  sd {statistics.pstdev(values):6.1f}"

    print("\nproperty match check (these must be close, or the test is rigged):")
    for key in ("mw", "clogp"):
        print(dist("actives", actives, key))
        print(dist("decoys", decoys, key))

    everything = actives + decoys
    labels = [True] * len(actives) + [False] * len(decoys)

    structure = prepare(pdb_id)
    print(f"\ndocking {len(everything)} compounds into {pdb_id}...")
    result = dock(everything, structure, exhaustiveness=8, verbose=True)
    apply_scores(everything, result)

    scored = [(c, lab) for c, lab in zip(everything, labels) if c.docking_score is not None]
    dock_scores = [c.docking_score for c, _ in scored]
    labs = [lab for _, lab in scored]
    mw_scores = [c.desc["mw"] for c, _ in scored]
    rng = random.Random(7)
    rand_scores = [rng.random() for _ in scored]

    print(f"\nscored {len(scored)} of {len(everything)} "
          f"({sum(labs)} actives, {len(labs) - sum(labs)} decoys)\n")
    header = f"{'ranked by':<26}{'ROC AUC':>9}{'EF@1%':>8}{'EF@10%':>9}"
    print(header); print("-" * len(header))
    rows = [
        ("docking score", dock_scores, True),
        ("molecular weight (control)", mw_scores, False),
        ("random (control)", rand_scores, False),
    ]
    results = {}
    for name, values, lower_better in rows:
        auc = roc_auc(values, labs, lower_better)
        ef1 = enrichment_at(values, labs, 0.01, lower_better)
        ef10 = enrichment_at(values, labs, 0.10, lower_better)
        results[name] = auc
        print(f"{name:<26}{auc if auc is not None else '—':>9}"
              f"{ef1 if ef1 is not None else '—':>8}{ef10 if ef10 is not None else '—':>9}")

    auc = results["docking score"]
    mw_auc = results["molecular weight (control)"]
    print("\nVERDICT")
    if auc is None:
        print("  inconclusive")
    elif auc >= 0.7 and auc > (mw_auc or 0) + 0.05:
        print(f"  docking separates binders from non-binders (AUC {auc}, beating the "
              f"size control at {mw_auc}).")
        print("  Usable as a cheap first filter over a large diverse library.")
    elif auc >= 0.7:
        print(f"  AUC {auc} looks good, but the molecular-weight control scores {mw_auc}.")
        print("  Docking is mostly rediscovering molecular size — not worth the compute.")
    else:
        print(f"  docking does not separate binders from non-binders here (AUC {auc}; "
              f"0.5 is chance).")
        print("  Combined with the ranking failure, docking into this structure has no use "
              "in this pipeline.")


if __name__ == "__main__":
    main()
