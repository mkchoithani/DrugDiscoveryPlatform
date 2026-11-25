"""Sweep many targets to find the boundary of the claim.

One result says the loop works. Twenty say it usually works. A hundred, across
different target families, say *when* it works — which is the thing worth knowing
before promising anything to anyone.

Results are appended to a CSV as they finish, and programmes already in that file
are skipped, so an interrupted run resumes instead of starting over.

    python3 mega_sweep.py [--targets A,B,C] [--min-size 40] [--out data/sweep.csv]
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

from dmta.chembl import SUGGESTED_TARGETS, build_dataset, cache_path, programmes
from programme_sweep import MIN_SPREAD, assess

# A spread of target families: kinases, a protease, a channel, an enzyme, a GPCR.
# Breadth matters more than depth — a rule derived from kinases alone would not
# be a rule, it would be a fact about kinases.
DEFAULT_TARGETS = [
    ("CHEMBL203", "IC50"),    # EGFR — kinase
    ("CHEMBL279", "IC50"),    # VEGFR2 — kinase
    ("CHEMBL1824", "IC50"),   # HER2 — kinase
    ("CHEMBL1862", "IC50"),   # ABL1 — kinase
    ("CHEMBL4005", "IC50"),   # PI3K alpha — lipid kinase
    ("CHEMBL244", "IC50"),    # Factor Xa — protease
    ("CHEMBL240", "IC50"),    # hERG — ion channel
    ("CHEMBL230", "IC50"),    # COX-2 — oxidoreductase
    ("CHEMBL205", "Ki"),      # carbonic anhydrase II — lyase
    ("CHEMBL233", "Ki"),      # mu opioid receptor — GPCR
]

FIELDS = ["target", "standard_type", "document", "n", "scaffolds", "spread", "year",
          "rho", "enrichment", "lift", "speedup", "tested", "found", "note"]


def needs_refresh(target: str, standard_type: str) -> bool:
    """A dataset cached before paper ids were kept cannot be split into programmes."""
    path = cache_path(target, standard_type)
    if not path.exists():
        return False
    with open(path, newline="") as fh:
        header = fh.readline()
    return "document" not in header


def done_keys(out: Path) -> set[tuple[str, str]]:
    if not out.exists():
        return set()
    with open(out, newline="") as fh:
        return {(r["target"], r["document"]) for r in csv.DictReader(fh)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", help="comma-separated ChEMBL target ids")
    ap.add_argument("--min-size", type=int, default=40)
    ap.add_argument("--max-per-target", type=int, default=12)
    ap.add_argument("--out", default="data/programme_sweep.csv")
    args = ap.parse_args()

    targets = ([(t.strip(), "IC50") for t in args.targets.split(",")]
               if args.targets else DEFAULT_TARGETS)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    already = done_keys(out)
    if already:
        print(f"resuming — {len(already)} programmes already evaluated\n")

    if not out.exists():
        with open(out, "w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=FIELDS).writeheader()

    t0 = time.time()
    for target, std_type in targets:
        label = SUGGESTED_TARGETS.get(target, target)
        print(f"\n=== {target} ({std_type}) — {label} ===", flush=True)
        try:
            build_dataset(target, std_type, refresh=needs_refresh(target, std_type),
                          verbose=False)
            found = programmes(target, std_type, min_size=args.min_size)
        except Exception as exc:
            print(f"  skipped: {exc}", flush=True)
            continue

        usable = [p for p in found if p["spread"] >= MIN_SPREAD]
        print(f"  {len(found)} programmes of {args.min_size}+, "
              f"{len(usable)} with enough spread", flush=True)

        for prog in usable[:args.max_per_target]:
            if (target, prog["document"]) in already:
                continue
            try:
                row = assess(prog)
            except Exception as exc:
                row = {"document": prog["document"], "n": prog["n"],
                       "scaffolds": prog["scaffolds"], "spread": prog["spread"],
                       "year": prog["year"], "note": f"failed: {exc}"}
            row["target"] = target
            row["standard_type"] = std_type
            with open(out, "a", newline="") as fh:       # append as we go
                csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore").writerow(row)
            print(f"    {row['document']:<16}{row['n']:>5}  "
                  f"rho {str(row.get('rho', '—')):>7}  lift {str(row.get('lift', '—')):>6}  "
                  f"speed {str(row.get('speedup', '—')):>6}"
                  + (f"  [{row['note']}]" if row.get("note") else ""), flush=True)

    print(f"\nswept in {time.time() - t0:.0f}s -> {out}")
    print("now run:  python3 sweep_analysis.py")


if __name__ == "__main__":
    main()
