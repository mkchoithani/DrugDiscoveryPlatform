"""Find the boundary of the claim: when does the loop help, and when does it not?

Takes the swept programmes and looks for a simple rule that predicts success from
facts you know BEFORE running anything — how many compounds, how much the potency
varies, how many chemotypes. A rule needing the answer first would be useless.

The rule is derived on half the programmes and scored on the other half. A rule
fitted and graded on the same data always looks good; that is the error this whole
project exists to avoid, and it applies to this analysis too.

    python3 sweep_analysis.py [path/to/programme_sweep.csv]
"""
from __future__ import annotations

import csv
import random
import sys
from itertools import product
from pathlib import Path


def helped(row: dict) -> bool:
    lift = _f(row.get("lift"))
    speed = _f(row.get("speedup"))
    return (lift or 0) >= 1.3 or (speed or 0) >= 1.5


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if not r.get("note")]
    for r in rows:
        r["_n"] = int(r["n"])
        r["_spread"] = float(r["spread"])
        r["_scaffolds"] = int(r["scaffolds"])
        r["_density"] = r["_n"] / max(1, r["_scaffolds"])   # compounds per chemotype
        r["_rho"] = _f(r.get("rho"))
        r["_helped"] = helped(r)
    return rows


def rate(rows: list[dict]) -> str:
    if not rows:
        return "no programmes"
    k = sum(1 for r in rows if r["_helped"])
    return f"{k}/{len(rows)} ({k / len(rows):.0%})"


def bucket_report(rows: list[dict], key: str, edges: list[float], label: str) -> None:
    print(f"\n  by {label}:")
    for lo, hi in zip([float('-inf')] + edges, edges + [float('inf')]):
        sel = [r for r in rows if lo <= r[key] < hi]
        if not sel:
            continue
        name = (f"< {hi:g}" if lo == float('-inf') else
                f"{lo:g}+" if hi == float('inf') else f"{lo:g}–{hi:g}")
        print(f"    {name:<12} {rate(sel):<16} n={len(sel)}")


def find_rule(train: list[dict]) -> tuple[int, float, float]:
    """The simplest rule that holds: enough compounds, enough potency variation."""
    best, best_score = (100, 2.0, 0.0), -1.0
    for min_n, min_spread in product([40, 60, 80, 100, 120, 150], [1.5, 2.0, 2.5, 3.0]):
        sel = [r for r in train if r["_n"] >= min_n and r["_spread"] >= min_spread]
        if len(sel) < max(5, len(train) * 0.25):      # a rule that accepts almost nobody
            continue                                   # is not a useful rule
        success = sum(1 for r in sel if r["_helped"]) / len(sel)
        coverage = len(sel) / len(train)
        score = success + 0.25 * coverage              # favour accuracy, keep coverage
        if score > best_score:
            best, best_score = (min_n, min_spread, success), score
    return best


def main():
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "data/programme_sweep.csv")
    if not path.exists():
        sys.exit(f"{path} not found — run mega_sweep.py first")
    rows = load(path)
    if len(rows) < 20:
        sys.exit(f"only {len(rows)} programmes evaluated; the sweep is still running")

    targets = sorted({r["target"] for r in rows})
    print(f"{len(rows)} programmes across {len(targets)} targets: {', '.join(targets)}")
    print(f"\noverall, the loop helped in {rate(rows)}")

    bucket_report(rows, "_n", [60, 100, 150, 250], "programme size (compounds)")
    bucket_report(rows, "_spread", [2.0, 3.0, 4.0], "potency spread (log units)")
    bucket_report(rows, "_density", [1.5, 2.5, 4.0], "compounds per chemotype")
    rho_rows = [r for r in rows if r["_rho"] is not None]
    bucket_report(rho_rows, "_rho", [0.2, 0.4, 0.6], "model ranking of unseen chemotypes (rho)")

    # derive on half, grade on the other half
    shuffled = rows[:]
    random.Random(0).shuffle(shuffled)
    cut = len(shuffled) // 2
    train, test = shuffled[:cut], shuffled[cut:]
    min_n, min_spread, train_rate = find_rule(train)

    qualified = [r for r in test if r["_n"] >= min_n and r["_spread"] >= min_spread]
    rejected = [r for r in test if not (r["_n"] >= min_n and r["_spread"] >= min_spread)]

    print(f"\n{'=' * 62}\nQUALIFICATION RULE (derived on {len(train)}, graded on {len(test)})\n{'=' * 62}")
    print(f"\n  Take the programme when it has {min_n}+ compounds "
          f"and {min_spread}+ log units of potency spread.\n")
    print(f"    programmes that qualify:  {rate(qualified):<18} n={len(qualified)}")
    print(f"    programmes it turns away: {rate(rejected):<18} n={len(rejected)}")
    if qualified and rejected:
        q = sum(1 for r in qualified if r["_helped"]) / len(qualified)
        j = sum(1 for r in rejected if r["_helped"]) / len(rejected)
        print(f"\n    the rule separates {q:.0%} success from {j:.0%} — "
              f"a {q - j:+.0%} difference on data it was not fitted to."
              if q != j else "\n    the rule does not separate the two groups.")
        if q - j < 0.15:
            print("    That is too small to qualify customers on. Report the overall")
            print("    rate instead, and say that programme size does not predict success.")

    print("\n  honest limits of this rule:")
    print(f"    - derived from {len(rows)} programmes across {len(targets)} targets, all public")
    print("    - a retrospective replay flatters every method; prospective results will be lower")
    print("    - 'helped' means 1.3x enrichment or reaching the best compound 1.5x sooner")


if __name__ == "__main__":
    main()
