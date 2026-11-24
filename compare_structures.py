"""Compare docking structures against the same compounds, plus their ensemble.

Reads scores already in the docking cache, so it costs nothing to re-run. Ensemble
docking takes each compound's best score across structures: the usual remedy when a
single rigid receptor suits only the chemotype it was crystallised with.

    python3 compare_structures.py [target] 2RGP 4HJO ...
"""
from __future__ import annotations

import random
import sys

import numpy as np

from dmta.chembl import load_compounds
from dmta.docking import _load_cache
from dmta.validate import _spearman

THRESHOLD = 8.0


def scores_for(pdb_id: str, cache: dict) -> dict[str, float]:
    """inchikey -> score, for whatever settings this structure was run with."""
    prefix = f"{pdb_id.upper()}|"
    return {k.rsplit("|", 1)[-1]: v for k, v in cache.items() if k.startswith(prefix)}


def metrics(scores: dict[str, float], truth: dict[str, float]) -> dict:
    keys = [k for k in scores if k in truth]
    if len(keys) < 10:
        return {"n": len(keys), "error": "too few compounds"}
    dock = np.array([scores[k] for k in keys])
    real = np.array([truth[k] for k in keys])
    active = real >= THRESHOLD
    k = max(1, len(keys) // 10)
    top = np.argsort(dock)[:k]
    base = float(active.mean())
    return {
        "n": len(keys),
        "spearman": _spearman(-dock, real),
        "top_decile": float(active[top].mean()),
        "base": base,
        "enrichment": round(float(active[top].mean()) / base, 2) if base else None,
        "separation": round(float(dock[~active].mean() - dock[active].mean()), 2)
                      if active.any() and (~active).any() else None,
    }


def main():
    args = sys.argv[1:]
    target = args[0] if args and args[0].startswith("CHEMBL") and len(args[0]) > 7 else "CHEMBL203"
    pdbs = [a for a in args if a != target] or ["2RGP", "4HJO"]

    cache = _load_cache()
    library = load_compounds(target)
    # Exactly the sample dock_validate.py uses, so these numbers are comparable with
    # it — and so hand-picked compounds added later cannot skew the base rate.
    sample = random.Random(0).sample(library, min(250, len(library)))
    truth = {c.inchikey: c.truth for c in sample}

    per_structure = {p: scores_for(p, cache) for p in pdbs}
    for p, s in per_structure.items():
        if not s:
            print(f"warning: no cached scores for {p} — run dock_validate.py {p} first")

    # ensemble: best (most negative) score each compound achieved anywhere
    common = set.intersection(*[set(s) for s in per_structure.values() if s]) if any(per_structure.values()) else set()
    ensemble = {k: min(s[k] for s in per_structure.values() if k in s) for k in common}

    print(f"{target}: docking structures compared on the compounds each one scored\n")
    header = f"{'structure':<14}{'n':>6}{'rho':>8}{'top-10%':>10}{'base':>8}{'enrich':>9}{'separation':>12}"
    print(header); print("-" * len(header))
    rows = {**per_structure, f"ensemble({len(pdbs)})": ensemble}
    for name, s in rows.items():
        if not s:
            print(f"{name:<14}  no data")
            continue
        m = metrics(s, truth)
        if "error" in m:
            print(f"{name:<14}{m['n']:>6}  {m['error']}")
            continue
        sep = "—" if m["separation"] is None else f"{m['separation']:+.2f}"
        rho = "—" if m["spearman"] is None else f"{m['spearman']:.3f}"
        print(f"{name:<14}{m['n']:>6}{rho:>8}{m['top_decile']:>10.0%}"
              f"{m['base']:>8.0%}{str(m['enrichment']):>9}{sep:>12}")

    print("\n(rho is sign-flipped: positive means the score agrees with measured potency)")
    print("A structure beats random only above 1.0x enrichment; the fingerprint model")
    print("reached 3.51x on this same sample, so docking has to clear that to be worth using.")

    best = max((m for m in ((n, metrics(s, truth)) for n, s in rows.items() if s)
                if "error" not in m[1]), key=lambda t: t[1]["enrichment"] or 0, default=None)
    if best:
        name, m = best
        print(f"\nbest of these: {name} at {m['enrichment']}x enrichment, rho {m['spearman']}")
        if (m["enrichment"] or 0) < 1.5:
            print("None of them is usable. Docking stays switched off.")


if __name__ == "__main__":
    main()
