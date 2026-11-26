"""Were docking and pharmacophore maps rejected unfairly?

Both were judged on potency ranking inside a mature EGFR series, and both failed
([RESULTS.md](RESULTS.md) 4-5 and 10). Neither claims to do that job. This runs them
on the job they do claim: pulling actives out of a library of mostly-inactive
molecules, on DUD-E's EGFR set (542 actives, 35,050 property-matched decoys).

Arms, all scored on the same held-out compounds:

  pharmacophore   predicted field for 2RGP, no training (dmta/pharmacophore.py)
  heavy atoms     the size baseline that beat the pharmacophore last time
  2D forest       Morgan fingerprints, trained on a scaffold-split half of the
                  actives plus a random half of the decoys
  docking         smina against 2RGP, with --dock (slow, so off by default)

**Read the 2D forest result with suspicion, not admiration.** DUD-E's actives come
from the literature in analogue series, and its decoys are property-matched but
structurally distant, so a fingerprint model can score well by recognising the series
rather than the binding. That bias is why LIT-PCBA exists. A big 2D number here is
evidence about the benchmark.

    python3 hitfinding_benchmark.py <dude_dir> <field.npz> [n_decoys] [--dock]
"""
from __future__ import annotations

import sys
import time

import numpy as np

from dmta.chem import descriptors, inchikey
from dmta.chembl import murcko
from dmta.hitfinding import hit_finding_report, load_ism
from dmta.library import Compound
from dmta.models import PropertyModel
from dmta.pharmacophore import load_field, score_compound


def as_compound(smiles: str, cid: str, index: int) -> Compound | None:
    try:
        desc = descriptors(smiles)
    except Exception:
        return None
    return Compound(
        # a real InChIKey, because the docking pose cache is keyed on it
        id=f"{cid}-{index}", smiles=smiles, inchikey=inchikey(smiles),
        scaffold_id=murcko(smiles),
        intermediate_key=cid, intermediate_smiles="", route=[], desc=desc,
        alerts=[], overall_yield=1.0, route_confidence=None, cost=0.0, lead_days=0,
    )


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    dude_dir, field_path = sys.argv[1], sys.argv[2]
    n_decoys = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith("-") else 5000
    want_dock = "--dock" in sys.argv

    actives_raw = load_ism(f"{dude_dir}/actives_final.ism")
    decoys_raw = load_ism(f"{dude_dir}/decoys_final.ism")
    print(f"DUD-E: {len(actives_raw)} actives, {len(decoys_raw)} decoys")

    rng = np.random.default_rng(0)
    if n_decoys < len(decoys_raw):
        idx = rng.choice(len(decoys_raw), size=n_decoys, replace=False)
        decoys_raw = [decoys_raw[i] for i in sorted(idx)]
        print(f"  decoys subsampled to {len(decoys_raw)} to keep the run tractable")

    actives = [c for i, (s, cid) in enumerate(actives_raw)
               if (c := as_compound(s, cid, i)) is not None]
    decoys = [c for i, (s, cid) in enumerate(decoys_raw)
              if (c := as_compound(s, cid, i)) is not None]
    for c in actives:
        c.truth = 1.0
    for c in decoys:
        c.truth = 0.0
    print(f"  parsed {len(actives)} actives, {len(decoys)} decoys\n")

    # The 2D model needs training data, so everything is reported on a held-out
    # portion. Actives are split by scaffold, which is the honest split: a model
    # that has seen a chemotype can recognise it trivially.
    scaffolds = {}
    for c in actives:
        scaffolds.setdefault(c.scaffold_id, []).append(c)
    keys = sorted(scaffolds)
    rng.shuffle(keys)
    train_actives, test_actives = [], []
    for k in keys:
        (train_actives if len(train_actives) < len(actives) * 0.5
         else test_actives).extend(scaffolds[k])
    cut = len(decoys) // 2
    train_decoys, test_decoys = decoys[:cut], decoys[cut:]
    test = test_actives + test_decoys
    labels = np.array([c.truth == 1.0 for c in test])
    print(f"train: {len(train_actives)} actives / {len(train_decoys)} decoys "
          f"({len(scaffolds)} active scaffolds, split by scaffold)")
    print(f"test:  {len(test_actives)} actives / {len(test_decoys)} decoys "
          f"= {labels.mean():.2%} active\n")

    rows = []

    # --- pharmacophore, no training ------------------------------------------
    field = load_field(field_path)
    print(f"pharmacophore query: {len(field.spots)} hot spots")
    t0 = time.time()
    scores = []
    for i, c in enumerate(test):
        r = score_compound(field, c.smiles, n_conformers=3)
        scores.append(r["score"] if r["score"] is not None else np.nan)
        if (i + 1) % 1000 == 0:
            rate = (i + 1) / (time.time() - t0)
            print(f"  {i + 1}/{len(test)} at {rate:.1f}/s, "
                  f"~{(len(test) - i - 1) / rate / 60:.1f} min left", flush=True)
    phc4 = np.array(scores, dtype=float)
    print(f"  scored in {time.time() - t0:.0f}s "
          f"({int(np.isnan(phc4).sum())} unembeddable)\n")
    rows.append(hit_finding_report(phc4, labels, "pharmacophore (2RGP field)"))

    # --- baselines ------------------------------------------------------------
    rows.append(hit_finding_report(
        np.array([c.desc["heavy_atoms"] for c in test], float), labels, "heavy atom count"))
    rows.append(hit_finding_report(
        np.array([c.desc["qed"] for c in test], float), labels, "QED"))

    # --- the 2D incumbent -----------------------------------------------------
    train = train_actives + train_decoys
    for c in train:
        c.measured = c.truth
    try:
        model = PropertyModel(n_estimators=300).fit(train)
        mu, _ = model.predict(test)
    finally:
        for c in train:
            c.measured = None
    rows.append(hit_finding_report(np.asarray(mu, float), labels,
                                  "2D forest (see the caveat)"))

    # --- docking, optional ----------------------------------------------------
    if want_dock:
        from dmta.docking import apply_scores, dock, engine_available
        from dmta.structure import prepare
        if not engine_available():
            print("smina not found — skipping the docking arm\n")
        else:
            print("docking the test set against 2RGP (slow)...", flush=True)
            t0 = time.time()
            structure = prepare("2RGP")
            result = dock(test, structure, verbose=True)
            n_applied = apply_scores(test, result)
            print(f"  docked {n_applied}/{len(test)} in {time.time() - t0:.0f}s")
            if n_applied == 0:
                print("  no compound docked — skipping the arm rather than "
                      "reporting an empty result")
            ds = np.array([c.docking_score if c.docking_score is not None else np.nan
                           for c in test], float)
            if np.isfinite(ds).sum() >= 50:
                # lower kcal/mol is a better fit, so negate to make higher = better
                rows.append(hit_finding_report(-ds, labels, "docking (smina, 2RGP)"))

    header = (f"{'method':<30}{'n':>7}{'AUC':>7}{'EF 1%':>8}{'EF 5%':>8}"
              f"{'BEDROC':>9}")
    print(header); print("-" * len(header))
    for r in rows:
        print(f"{r['method']:<30}{r['n']:>7}{r['roc_auc']:>7.3f}"
              f"{(r['ef_1pct'] if r['ef_1pct'] is not None else float('nan')):>8.2f}"
              f"{(r['ef_5pct'] if r['ef_5pct'] is not None else float('nan')):>8.2f}"
              f"{r['bedroc_20']:>9.3f}")

    print("\n" + "=" * 69)
    # Judge on EARLY enrichment, not AUC. AUC measures the whole ranked list, and
    # nobody screens the whole list — a method can hold a respectable AUC by sorting
    # the bottom 90% well while putting nothing useful in the top 1%. An earlier
    # version of this verdict compared AUC alone and declared the pharmacophore the
    # winner at 0.580 vs 0.516, while it was losing 1.46 to 6.57 on EF at 1%.
    def better(a, b) -> bool:
        # If either method has pinned EF 1% to its ceiling, that column carries no
        # information and BEDROC plus AUC have to decide. Comparing 3.74 against
        # 3.77 when both mean "ten actives out of ten" is not a comparison.
        if a.get("ef_1pct_saturated") or b.get("ef_1pct_saturated"):
            print(f"    (EF 1% is saturated at its ceiling of {a['ef_ceiling']} — "
                  "judging on AUC and BEDROC instead)")
            return a["bedroc_20"] > b["bedroc_20"] and a["roc_auc"] > b["roc_auc"]
        return ((a["ef_1pct"] or 0) > (b["ef_1pct"] or 0)
                and a["bedroc_20"] > b["bedroc_20"])

    phc, size = rows[0], rows[1]
    structure_based = [r for r in rows if "pharmacophore" in r["method"]
                       or "docking" in r["method"]]
    print("  judged on EF at 1% and BEDROC — early enrichment is what a screen uses:")
    for r in structure_based:
        if better(r, size):
            print(f"    {r['method']}: beats the size baseline "
                  f"(EF1% {r['ef_1pct']} vs {size['ef_1pct']}). Worth revisiting as a")
            print("    screening filter.")
        else:
            print(f"    {r['method']}: does NOT beat molecular size "
                  f"(EF1% {r['ef_1pct']} vs {size['ef_1pct']}, "
                  f"BEDROC {r['bedroc_20']} vs {size['bedroc_20']}).")
            print("    The earlier rejection stands; it was not an artefact of the task.")

    forest = next((r for r in rows if "forest" in r["method"]), None)
    if forest and forest["roc_auc"] > 0.95:
        print(f"\n  NOTE ON THE BENCHMARK: the 2D forest reaches AUC "
              f"{forest['roc_auc']:.3f} and")
        print(f"  EF1% {forest['ef_1pct']}, which is the ceiling at this active rate.")
        print("  That is DUD-E's known bias measured, not a model result: its actives")
        print("  are literature compounds and its decoys are ZINC molecules, so a")
        print("  fingerprint separates them almost perfectly without knowing anything")
        print("  about binding. Treat this row as a property of the benchmark.")


if __name__ == "__main__":
    main()
