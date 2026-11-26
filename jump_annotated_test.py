#!/usr/bin/env python3
"""F4, the powered version: Cell Painting against fingerprints on JUMP's own labels.

WHY THIS EXISTS ALONGSIDE jump_cellpaint_test.py
------------------------------------------------
The LIT-PCBA version of this test has the better labels and no power. Its inactives
were measured rather than assumed and its split is AVE-unbiased, but only 382 of
LIT-PCBA's actives across fifteen targets are also in JUMP, and exactly one target
(ALDH1) clears a 20-active threshold. One target cannot carry a verdict.

This version has the power and the weaker labels. JUMP ships target annotations for
17,891 of its own compounds across 103 proteins, every one of which has a profile.

THE WEAKNESS IN THESE LABELS, STATED BEFORE THE RESULT
------------------------------------------------------
A 0 here means "not annotated against this protein", which is mostly "never tested"
and only sometimes "tested and inactive". LIT-PCBA exists precisely because that
distinction changes benchmark results -- it is the same defect as DUD-E's untested
decoys, which section 15 showed could reverse a docking verdict. So:

  * this benchmark is easier than a real screen, for every arm equally;
  * the comparison BETWEEN arms is still fair, because every arm sees the same labels;
  * an absolute AUC from this script is not comparable to one from LIT-PCBA, and
    neither file should be quoted as though it were the other.

The question F4 asks is a comparison, not an absolute, so the weak labels damage it
less than having one target would. Both runs are reported.

THE SPLIT
---------
Bemis-Murcko scaffold split, 75/25, scaffolds disjoint across the two sides. JUMP
carries no dates, so a time split is impossible here and the scaffold split is the
harder of the two available. FEATURES.md asked for a time split; this is not one, and
a win here would be necessary rather than sufficient for F4's actual claim.

ARMS AND CONTROLS: identical to jump_cellpaint_test.py, and imported from it, so the
two runs cannot drift apart. `profile-scrambled` is why the first run's pre-specified
comparison was void -- 614 extra columns cost the forest the same AUC whether they
held morphology or noise -- and `fused` is the dilution-immune replacement.

Usage:
    python3 jump_annotated_test.py            # targets with >=MIN_ACTIVES
    python3 jump_annotated_test.py 8          # the 8 largest targets
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).parent))
from dmta.chem import feature_vector                                    # noqa: E402
from jump_cellpaint_test import (ARMS, paired_auc_ci, score_arms,       # noqa: E402
                                 verdict, CACHE, RESAMPLES)

ANNOTATIONS = CACHE / "target_annotations.csv"
PROFILES = CACHE / "pert_median.parquet"
SCORES = CACHE / "annotated_scores.npz"

MIN_ACTIVES = 150          # in the whole set, before the split
MIN_TEST_ACTIVES = 20      # the same floor section 17 used
MAX_INACTIVES = 4000       # keeps one target's fit to seconds; every active is kept
SEED = 0


def scaffold_key(smiles: str) -> str:
    """The Bemis-Murcko scaffold, via the platform's own murcko().

    The first version of this called rdkit.Chem.Scaffolds.MurckoScaffold and wrapped
    it in `except Exception: return ""`. In RDKit 2026.09.2 that call raises "bad any
    cast" on EVERY molecule, not just the acyclic ones, so the guard turned a total
    failure into a silent one: every compound got an empty key, every compound became
    its own group, and the "scaffold split" was a split on row order. The run built on
    it was discarded.

    dmta.chembl.murcko already works around this exact RDKit defect with rdMolHash,
    and has since before this section existed. Reimplementing a thing the codebase
    already does correctly is how the bug got in; delegating is the fix.

    Returns "acyclic" for a molecule with no ring system -- a real key, shared by all
    such molecules, rather than an empty string standing in for an error.
    """
    from dmta.chembl import murcko
    return murcko(smiles)


def scaffold_split(smiles: list[str], frac_train: float = 0.75) -> np.ndarray:
    """True where the compound is in the training side. Scaffolds never cross."""
    groups: dict[str, list[int]] = {}
    for i, s in enumerate(smiles):
        key = scaffold_key(s)
        # "unknown" means the SMILES would not parse, so it is not a shared scaffold
        # and each such compound stands alone. "acyclic" IS a shared key: every
        # ring-free compound belongs to one group and goes to one side of the split.
        groups.setdefault(f"_bad_{i}" if key == "unknown" else key, []).append(i)

    # Largest scaffolds into training first, and no randomness at all: the test side
    # must be the rarer chemistry, deterministically, or the split is a different
    # experiment on every run. (The first version of this line mixed `and` with `or`
    # in one condition, which Python reads as `A or (B and C)` -- it put a coin flip
    # in the middle of a split that is supposed to be reproducible. The same
    # precedence slip already cost this project a bad test assertion once.)
    order = sorted(groups.values(), key=len, reverse=True)
    train = np.zeros(len(smiles), bool)
    n_train_target = int(frac_train * len(smiles))
    n = 0
    for grp in order:
        if n >= n_train_target:
            break
        train[grp] = True
        n += len(grp)
    return train


def load() -> tuple[pd.DataFrame, np.ndarray, list[str], list[str]]:
    ann = pd.read_csv(ANNOTATIONS, low_memory=False)
    prof = pd.read_parquet(PROFILES).rename(columns={"jcp": "JUMP_ID"})
    merged = ann.merge(prof, on="JUMP_ID", how="inner")
    smiles = merged["Standardized_SMILES"].tolist()
    feat_cols = [c for c in prof.columns if c != "JUMP_ID"]
    target_cols = [c for c in ann.columns
                   if c not in ("JUMP_ID", "Standardized_SMILES")]
    P = merged[feat_cols].to_numpy(dtype=float)
    col_mean = np.nanmean(P, axis=0)
    bad = np.isnan(P)
    P[bad] = np.take(col_mean, np.where(bad)[1])
    return merged, P, smiles, target_cols


def main(limit: int | None) -> None:
    if not ANNOTATIONS.exists():
        sys.exit(f"missing {ANNOTATIONS}")
    print("loading annotations and profiles ...", flush=True)
    merged, P, smiles, target_cols = load()
    print(f"  {len(merged)} compounds with both a profile and an annotation row, "
          f"{P.shape[1]} profile features, {len(target_cols)} proteins", flush=True)

    counts = merged[target_cols].apply(pd.to_numeric, errors="coerce").sum()
    chosen = counts[counts >= MIN_ACTIVES].sort_values(ascending=False)
    if limit:
        chosen = chosen.head(limit)
    print(f"  {len(chosen)} proteins with >= {MIN_ACTIVES} actives\n", flush=True)

    print("featurising once for all targets ...", flush=True)
    F = np.array([feature_vector(s) for s in smiles], dtype=float)
    train_mask = scaffold_split(smiles)
    print(f"  fingerprint {F.shape}, scaffold split "
          f"{train_mask.sum()} train / {(~train_mask).sum()} test\n", flush=True)

    rng = np.random.default_rng(SEED)
    results, skipped = {}, []
    for prot in chosen.index:
        y = pd.to_numeric(merged[prot], errors="coerce").fillna(0).to_numpy(int)
        act = np.where(y == 1)[0]
        ina = np.where(y == 0)[0]
        if len(ina) > MAX_INACTIVES:
            ina = rng.choice(ina, MAX_INACTIVES, replace=False)
        keep = np.sort(np.concatenate([act, ina]))

        d = {"target": prot, "F": F[keep], "P": P[keep], "y": y[keep],
             "split": np.where(train_mask[keep], "T", "V"),
             "heavy": np.array([_heavy(smiles[i]) for i in keep], float)}
        n_te_a = int(d["y"][d["split"] == "V"].sum())
        if n_te_a < MIN_TEST_ACTIVES:
            skipped.append((prot, f"{n_te_a} test actives"))
            continue
        s = score_arms(d, seed=SEED)
        if s is None:
            skipped.append((prot, "too few actives on one side"))
            continue
        results[prot] = s
        print(f"  {prot}  n={len(keep):5d}  test actives={n_te_a:4d}  "
              + "  ".join(f"{a}={roc_auc_score(s['_y'], s[a]):.3f}"
                          for a in ("fingerprint", "profile", "fused")), flush=True)

    for p, why in skipped:
        print(f"  skipped {p}: {why} -- not counted as a loss")
    if not results:
        sys.exit("\nno protein cleared the threshold")

    print("\n" + "=" * 78)
    print(f"AUC per protein, scaffold-split held-out set ({len(results)} proteins)\n")
    print(f"{'protein':10s}" + "".join(f"{a:>20s}" for a in ARMS))
    for prot, s in results.items():
        y = s["_y"]
        print(f"{prot:10s}" + "".join(f"{roc_auc_score(y, s[a]):20.3f}" for a in ARMS))

    pairs = {
        "both_minus_fingerprint": ("both", "fingerprint",
                                   "PRE-SPECIFIED PRIMARY (concatenated, against fingerprint)"),
        "both_minus_scrambled": ("both", "profile-scrambled",
                                 "THE CONTROL: real profile against the same columns scrambled"),
        "fused_minus_fingerprint": ("fused", "fingerprint",
                                    "DILUTION-IMMUNE: rank-fused, against fingerprint"),
        "profile_minus_scram": ("profile", "profile-only-scram",
                                "IS THE PROFILE ABOVE CHANCE AT ALL"),
        "profile_minus_size": ("profile", "size",
                               "profile alone, against heavy-atom count"),
    }

    per: dict[str, list[float]] = {}
    for key, (a, b, label) in pairs.items():
        print("\n" + "=" * 78)
        print(label + "\n")
        vals = []
        for prot, s in results.items():
            point, lo, hi = paired_auc_ci(s["_y"], s[a], s[b])
            vals.append(point)
            flag = "  " if lo <= 0 <= hi else ("+ " if lo > 0 else "- ")
            print(f"  {flag}{prot:10s} {point:+.3f}   [{lo:+.3f}, {hi:+.3f}]")
        per[key] = vals

    print("\n" + "=" * 78)
    print(f"ACROSS {len(results)} PROTEINS -- paired over proteins, "
          f"the structure section 17 used\n")
    for key, (a, b, label) in pairs.items():
        v = np.array(per[key])
        boot = np.array([np.mean(np.random.default_rng(i).choice(v, len(v), replace=True))
                         for i in range(RESAMPLES)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        wins = int((v > 0).sum())
        print(f"  {label}")
        print(f"    mean {v.mean():+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  "
              f"wins {wins}/{len(v)}")
        print(f"    -> {verdict(v.mean(), lo, hi, a)}\n")

    np.savez_compressed(SCORES, **{f"{k}": np.array(v) for k, v in per.items()})
    print(f"per-protein differences cached to {SCORES}")


def _heavy(s: str) -> float:
    from rdkit import Chem
    m = Chem.MolFromSmiles(s)
    return float(m.GetNumHeavyAtoms()) if m else 0.0


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
