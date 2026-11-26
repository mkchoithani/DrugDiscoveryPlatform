#!/usr/bin/env python3
"""Join LIT-PCBA's screening labels to JUMP Cell Painting's morphological profiles.

This is step one of the F4 test (FEATURES.md). It answers the question that decides
whether the rest of the test can be run at all: **how many compounds have both a
measured screening label and a Cell Painting profile?**

The join is on full InChIKey. Matching on the connectivity block alone would merge
stereoisomers, and stereoisomers routinely differ by orders of magnitude in potency,
so the looser match is counted and reported but never used to build a dataset.

WHY THE FORMULA PRE-FILTER
--------------------------
LIT-PCBA ships ~1.5M inactive SMILES across its targets and JUMP holds 115,796
compounds, so almost every molecule we parse is a miss. An InChIKey costs roughly
2 ms; a molecular formula costs roughly a tenth of that. Computing the formula first
and only paying for the InChIKey when the formula exists in JUMP skips the ~95% of
molecules that cannot match, and changes nothing about which pairs are found: two
molecules with the same InChIKey necessarily have the same formula, so the filter
cannot discard a true match.

Usage:
    python3 jump_match.py                 # every target with a LIT-PCBA directory
    python3 jump_match.py ALDH1 PKM2      # named targets only
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

LITPCBA = Path("/private/tmp/litpcba/ave")
CACHE = Path("data/jump_cache")
COMPOUND_PARQUET = CACHE / "compound.parquet"


def jump_index() -> tuple[dict[str, str], set[str]]:
    """InChIKey -> JCP2022 id, plus the set of molecular formulas JUMP contains."""
    cp = pd.read_parquet(COMPOUND_PARQUET)
    by_key = dict(zip(cp.Metadata_InChIKey, cp.Metadata_JCP2022))
    formulas = set()
    for smi in cp.Metadata_SMILES.dropna():
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            formulas.add(CalcMolFormula(m))
    return by_key, formulas


def scan(path: Path, by_key: dict[str, str], formulas: set[str]) -> list[tuple[str, str]]:
    """Return (inchikey, jcp) for every molecule in `path` that JUMP also holds."""
    found: list[tuple[str, str]] = []
    for line in path.open():
        parts = line.split()
        if not parts:
            continue
        m = Chem.MolFromSmiles(parts[0])
        if m is None:
            continue
        if CalcMolFormula(m) not in formulas:
            continue                      # cannot be a match; skip the InChIKey cost
        key = Chem.MolToInchiKey(m)
        jcp = by_key.get(key)
        if jcp is not None:
            found.append((key, jcp))
    return found


def run(targets: list[str]) -> pd.DataFrame:
    CACHE.mkdir(parents=True, exist_ok=True)
    print("indexing JUMP compounds ...", flush=True)
    t0 = time.time()
    by_key, formulas = jump_index()
    print(f"  {len(by_key)} compounds, {len(formulas)} distinct formulas "
          f"in {time.time() - t0:.0f}s", flush=True)

    rows = []
    for name in targets:
        d = LITPCBA / name
        out = CACHE / f"{name}-matched.parquet"
        if out.exists():
            got = pd.read_parquet(out)
            print(f"{name}: cached, {len(got)} matched", flush=True)
        else:
            t0 = time.time()
            recs = []
            for split, active in (("active_T", 1), ("active_V", 1),
                                  ("inactive_T", 0), ("inactive_V", 0)):
                p = d / f"{split}.smi"
                if not p.exists():
                    continue
                hits = scan(p, by_key, formulas)
                recs += [{"inchikey": k, "jcp": j, "active": active,
                          "split": "T" if split.endswith("_T") else "V"}
                         for k, j in hits]
                print(f"  {name}/{split}: {len(hits)} matched", flush=True)
            got = pd.DataFrame(recs)
            # One compound can appear in both the training and validation split of
            # LIT-PCBA; keep the first and record the duplicate count, because a
            # compound present twice would otherwise sit on both sides of a split.
            dupes = int(got.inchikey.duplicated().sum()) if len(got) else 0
            got = got.drop_duplicates("inchikey")
            got.to_parquet(out)
            print(f"{name}: {len(got)} matched, {dupes} duplicate keys dropped, "
                  f"{time.time() - t0:.0f}s", flush=True)
        if len(got):
            rows.append({"target": name, "matched": len(got),
                         "actives": int(got.active.sum()),
                         "inactives": int((got.active == 0).sum())})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    names = sys.argv[1:] or sorted(
        p.name for p in LITPCBA.iterdir()
        if p.is_dir() and (p / "active_V.smi").exists())
    table = run(names)
    print()
    print(table.to_string(index=False) if len(table) else "nothing matched")
