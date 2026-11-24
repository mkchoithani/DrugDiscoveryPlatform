"""Does docking into this structure actually rank real compounds?

The question that decides whether docking earns a place in the loop. We already know
the answer for every compound in the ChEMBL set, so we dock a random sample and check
whether the score agrees with reality — overall, and split by binding mechanism.

    python3 dock_validate.py [pdb_id] [n_compounds] [target]
"""
from __future__ import annotations

import random
import sys

from dmta.chembl import load_compounds
from dmta.docking import apply_scores, cache_size, dock, engine_available
from dmta.structure import prepare
from dmta.validate import (docking_report_by_mechanism, print_docking_report,
                           random_split, evaluate)


def main():
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "2RGP"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 250
    target = sys.argv[3] if len(sys.argv) > 3 else "CHEMBL203"

    if not engine_available():
        sys.exit("no docking engine found — see dmta/docking.py for the install line")

    structure = prepare(pdb_id)
    print(f"{structure.pdb_id}: {structure.title}")
    print(f"  {structure.resolution} A, chain {structure.chain}, "
          f"{structure.n_protein_atoms} protein atoms")
    print(f"  removed: {structure.removed}")
    print(f"  box from ligand {structure.site.reference_ligand} "
          f"({structure.site.n_ligand_atoms} atoms) at {structure.site.center}, "
          f"size {structure.site.size}")

    library = load_compounds(target)
    sample = random.Random(0).sample(library, min(n, len(library)))
    actives = sum(1 for c in sample if c.truth >= 8.0)
    print(f"\n{target}: docking a random sample of {len(sample)} of {len(library)} "
          f"compounds ({actives} active at pChEMBL >= 8)")
    print(f"  docking cache holds {cache_size()} previous results")

    result = dock(sample, structure, exhaustiveness=8, verbose=True)
    apply_scores(sample, result)
    if result.failed_embed or result.failed_dock:
        print(f"  could not embed: {result.failed_embed[:5]}")
        print(f"  could not dock:  {result.failed_dock[:5]}")

    print("\nDOCKING vs MEASURED POTENCY")
    print("(rho is sign-flipped, so positive means the score agrees with reality)\n")
    report = docking_report_by_mechanism(sample)
    print_docking_report(report)

    # For context: how well does the 2D model do on the same compounds?
    train, _ = random_split([c for c in library if c.id not in {s.id for s in sample}], test_frac=0.01)
    ml = evaluate(train[:3000], sample, n_estimators=150)
    print(f"\nfor comparison, the fingerprint model on these same {ml['n_test']} compounds:")
    print(f"  rho {ml['spearman']}   top-10% hit rate {ml['top_decile_hit_rate']:.0%} "
          f"vs {ml['base_hit_rate']:.0%} base   enrichment {ml['enrichment']}x")

    d = report["all"]
    print("\nVERDICT")
    if "error" in d:
        print(f"  inconclusive: {d['error']}")
    elif d["spearman"] is not None and d["spearman"] >= 0.3 and (d["enrichment"] or 0) >= 1.5:
        print(f"  docking adds signal (rho {d['spearman']}, {d['enrichment']}x enrichment).")
        print("  Worth adding as an objective, especially for new chemotypes.")
    elif d["spearman"] is not None and d["spearman"] > 0.1:
        print(f"  weak signal (rho {d['spearman']}, {d['enrichment']}x enrichment).")
        print("  Not useful alone; may still help on chemotypes the 2D model has never seen.")
    else:
        print(f"  no usable signal (rho {d['spearman']}, {d['enrichment']}x enrichment).")
        print("  Do NOT wire this into the loop. Try another structure before spending more.")


if __name__ == "__main__":
    main()
