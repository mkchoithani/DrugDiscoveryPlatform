#!/usr/bin/env python3
"""F4: does a Cell Painting profile beat a structural fingerprint at finding hits?

THE QUESTION, AND WHY IT IS THIS ONE
------------------------------------
Section 1 of RESULTS.md measured the finding that limits this whole platform: sigma
ranks error at 0.44 on a random split, 0.36 on a scaffold split and 0.12 on a time
split, inverting to -0.13 on mu-opioid. A fingerprint model only knows resemblance to
what it has already seen, so it fails exactly where a project needs it.

Two candidate fixes have been measured and neither worked. Docking was rejected four
times (sections 4-6, 15). Chemprop tied (section 16). A morphological profile is the
only remaining candidate that is different in *kind* rather than in degree: a cell's
response to a compound is not a statement about that compound's substructures.

WHAT THIS TEST IS NOT
---------------------
It is not the time-split test FEATURES.md asked for. LIT-PCBA carries no dates, so
this substitutes LIT-PCBA's AVE-unbiased training/validation split, which separates
dissimilar chemotypes. That is the scaffold split's cousin, not the time split's, and
the time split is the harder of the two. A win here would therefore be necessary but
not sufficient evidence for the claim F4 actually makes. This is stated here, before
any result exists, rather than discovered in the discussion.

THE FOUR ARMS, AND THE ONE THAT MATTERS
---------------------------------------
    fingerprint          Morgan + physchem, via dmta.chem.feature_vector -- the incumbent
    profile              the 614-dim Cell Painting median profile, alone
    both                 the two concatenated
    profile-scrambled    THE CONTROL: fingerprint + profiles permuted across compounds
    profile-only-scram   profile alone, permuted -- is `profile` above chance at all?
    fused                rank-average of the fingerprint and profile scores
    size                 heavy-atom count, which has beaten us before (section 10)

`fused` exists because of what the control found on the first run. Concatenating 614
columns onto a 2,060-column fingerprint costs the forest roughly 0.12 of AUC whether
those columns carry morphology or pure noise, so `both` cannot answer the question it
was built for: it measures dilution. Rank-averaging two separately fitted models adds
the profile's information without diluting the fingerprint's, so a difference between
`fused` and `fingerprint` is attributable to the profile's content. `fused` was added
after seeing the control fire, which is a change to the pre-specified plan and is
labelled as one everywhere it appears.

`profile-scrambled` is the arm this test is built around. `both` has 614 more columns
than `fingerprint`, and more columns can win on capacity alone. The scrambled arm has
exactly the same shape and exactly the same marginal distribution per column, and zero
information about which compound is which. If `both` beats `fingerprint` and the
scrambled arm beats it by a similar margin, the gain is dimensionality, not morphology,
and the result is void. Permutation happens once per target under the run's seed and is
applied to train and test together, so the scrambled arm is a real model fitted to real
noise rather than a model handed a broken test set.

PRE-SPECIFIED, BEFORE THE FIRST RUN
-----------------------------------
Primary comparison:  AUC(both) - AUC(fingerprint) on the held-out split.
Pass requires BOTH of:
    (a) the paired interval for `both` - `fingerprint` excludes zero on the positive
        side, per target, and the sign of the effect agrees across a majority of
        targets that clear MIN_TEST_ACTIVES;
    (b) the interval for `profile-scrambled` - `fingerprint` contains zero.
A tie is a tie. Section 16 reported one and kept the incumbent on cost; section 17
reported one and declined to reject a paper over it. The same rule applies here, and a
tie must not be written up as "Cell Painting does not work".

Section 17 cost us a fifth verdict-function error by pre-specifying a modelling choice
on theory and never checking it. So this script reports every arm against every
reference, labelled, rather than only the pre-specified pair.

POWER
-----
JUMP holds 115,796 compounds; LIT-PCBA's actives are mostly not among them. The matched
counts decide whether any of this is answerable, and the run prints them first and
refuses targets that cannot support a verdict instead of quoting an interval that
cannot resolve anything.

Usage:
    python3 jump_cellpaint_test.py            # all matched targets
    python3 jump_cellpaint_test.py ALDH1      # one target
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).parent))
from dmta.chem import feature_vector                                    # noqa: E402

CACHE = Path("data/jump_cache")
LITPCBA = Path("/private/tmp/litpcba/ave")
PROFILES = CACHE / "pert_median.parquet"

RESAMPLES = 4000
SEED = 0
# A verdict needs enough positives in the test split for an AUC to mean anything.
# Section 17 used 20 and that is kept here; lowering a threshold after seeing which
# targets fail it is how a benchmark becomes a selection of its own best cases.
MIN_TEST_ACTIVES = 20
MIN_TRAIN_ACTIVES = 20

ARMS = ("fingerprint", "profile", "both", "profile-scrambled",
        "profile-only-scram", "fused", "size")


# --------------------------------------------------------------- data assembly

def load_profiles() -> pd.DataFrame:
    """JCP2022 id -> median morphological profile, indexed by id."""
    df = pd.read_parquet(PROFILES)
    df = df.rename(columns={"jcp": "Metadata_JCP2022"})
    return df.set_index("Metadata_JCP2022")


def smiles_for(target: str) -> dict[str, str]:
    """InChIKey -> SMILES for every molecule in the target's LIT-PCBA files."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out: dict[str, str] = {}
    for split in ("active_T", "active_V", "inactive_T", "inactive_V"):
        p = LITPCBA / target / f"{split}.smi"
        if not p.exists():
            continue
        for line in p.open():
            parts = line.split()
            if not parts:
                continue
            m = Chem.MolFromSmiles(parts[0])
            if m is None:
                continue
            out.setdefault(Chem.MolToInchiKey(m), parts[0])
    return out


def build(target: str, profiles: pd.DataFrame) -> dict | None:
    matched = CACHE / f"{target}-matched.parquet"
    if not matched.exists():
        return None
    m = pd.read_parquet(matched)
    m = m[m.jcp.isin(profiles.index)]
    if not len(m):
        return None

    smi = smiles_for(target)
    m = m[m.inchikey.isin(smi)]
    if not len(m):
        return None

    P = profiles.loc[m.jcp.values].to_numpy(dtype=float)
    # A profile column that is all-NaN for this subset carries nothing; a few scattered
    # NaNs are filled with the column mean so a single missing well does not delete a
    # compound that is otherwise usable.
    keep = ~np.isnan(P).all(axis=0)
    P = P[:, keep]
    col_mean = np.nanmean(P, axis=0)
    bad = np.isnan(P)
    P[bad] = np.take(col_mean, np.where(bad)[1])

    F = np.array([feature_vector(smi[k]) for k in m.inchikey], dtype=float)
    y = m.active.to_numpy(int)
    split = m.split.to_numpy()
    heavy = np.array([_heavy(smi[k]) for k in m.inchikey], dtype=float)
    return {"target": target, "F": F, "P": P, "y": y, "split": split, "heavy": heavy}


def _heavy(s: str) -> float:
    from rdkit import Chem
    m = Chem.MolFromSmiles(s)
    return float(m.GetNumHeavyAtoms()) if m else 0.0


# ------------------------------------------------------------------- the arms

def score_arms(d: dict, seed: int = SEED) -> dict[str, np.ndarray] | None:
    tr = d["split"] == "T"
    te = ~tr
    y = d["y"]
    if y[tr].sum() < MIN_TRAIN_ACTIVES or y[te].sum() < MIN_TEST_ACTIVES:
        return None
    if len(np.unique(y[te])) < 2:
        return None

    rng = np.random.default_rng(seed)
    # One permutation of the profile rows, applied to train and test alike. The
    # scrambled arm is then an honest model fitted to a real, information-free matrix.
    perm = rng.permutation(len(y))
    P_scram = d["P"][perm]

    feats = {
        "fingerprint": d["F"],
        "profile": d["P"],
        "both": np.hstack([d["F"], d["P"]]),
        "profile-scrambled": np.hstack([d["F"], P_scram]),
        "profile-only-scram": P_scram,
    }
    out: dict[str, np.ndarray] = {}
    for name, X in feats.items():
        clf = RandomForestClassifier(
            n_estimators=300, min_samples_leaf=2, n_jobs=-1,
            class_weight="balanced", random_state=seed)
        clf.fit(X[tr], y[tr])
        out[name] = clf.predict_proba(X[te])[:, 1]
    # Late fusion: rank-average, so neither model's features dilute the other's.
    def ranks(v):
        order = np.empty(len(v), float)
        order[np.argsort(v, kind="stable")] = np.arange(len(v))
        return order / max(len(v) - 1, 1)
    out["fused"] = (ranks(out["fingerprint"]) + ranks(out["profile"])) / 2
    out["size"] = d["heavy"][te]
    out["_y"] = y[te]
    return out


# ------------------------------------------------------------------ statistics

def paired_auc_ci(y: np.ndarray, a: np.ndarray, b: np.ndarray,
                  resamples: int = RESAMPLES, seed: int = SEED) -> tuple[float, float, float]:
    """Stratified paired bootstrap of AUC(a) - AUC(b) on the same compounds."""
    rng = np.random.default_rng(seed)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    diffs = np.empty(resamples)
    for i in range(resamples):
        idx = np.concatenate([rng.choice(pos, len(pos), replace=True),
                              rng.choice(neg, len(neg), replace=True)])
        yy = y[idx]
        diffs[i] = roc_auc_score(yy, a[idx]) - roc_auc_score(yy, b[idx])
    point = roc_auc_score(y, a) - roc_auc_score(y, b)
    return point, float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def verdict(point: float, lo: float, hi: float, what: str) -> str:
    if lo > 0:
        return f"{what} WINS. The interval excludes zero."
    if hi < 0:
        return f"{what} LOSES. The interval excludes zero in the other direction."
    bound = max(abs(lo), abs(hi))
    return (f"A TIE, NOT A REJECTION. The interval contains zero, so this test did not "
            f"separate them.\n     An effect up to {bound:.3f} of AUC, either way, is "
            f"still inside it: a statement about what this could resolve, not a "
            f"measurement of no effect.")


# ------------------------------------------------------------------------ main

def main(targets: list[str]) -> None:
    if not PROFILES.exists():
        sys.exit(f"missing {PROFILES} -- run the profile download first")
    print("loading median profiles ...", flush=True)
    profiles = load_profiles()
    print(f"  {len(profiles)} perturbations x {profiles.shape[1]} features\n", flush=True)

    built, skipped = [], []
    for t in targets:
        d = build(t, profiles)
        if d is None:
            skipped.append((t, "no matched compounds with a profile"))
            continue
        n_tr_a = int(d["y"][d["split"] == "T"].sum())
        n_te_a = int(d["y"][d["split"] == "V"].sum())
        print(f"{d['target']:9s} matched {len(d['y']):5d}  "
              f"train {n_tr_a:4d}a/{(d['split']=='T').sum()-n_tr_a:5d}i  "
              f"test {n_te_a:4d}a/{(d['split']=='V').sum()-n_te_a:5d}i  "
              f"profile dim {d['P'].shape[1]}", flush=True)
        if n_te_a < MIN_TEST_ACTIVES or n_tr_a < MIN_TRAIN_ACTIVES:
            skipped.append((t, f"{n_tr_a} train / {n_te_a} test actives "
                               f"(need {MIN_TRAIN_ACTIVES}/{MIN_TEST_ACTIVES})"))
            continue
        built.append(d)

    print()
    for t, why in skipped:
        print(f"  skipped {t}: {why} -- not counted as a loss")

    if not built:
        print("\n" + "=" * 78)
        print("NO TARGET CAN SUPPORT A VERDICT.")
        print("""
This is the result, not a failure to produce one. The join between a screening
label and a morphological profile does not exist at the size a verdict needs.
That is a measurement of the data, and it is reported as one. It says nothing
about whether Cell Painting predicts activity -- that question is untouched by
this run, and writing it up as a negative would be the error section 17 named.
""")
        return

    print("\n" + "=" * 78)
    print("AUC per target, held-out (AVE-unbiased validation) split\n")
    head = f"{'target':10s}" + "".join(f"{a:>20s}" for a in ARMS)
    print(head)
    results = {}
    for d in built:
        s = score_arms(d)
        if s is None:
            continue
        results[d["target"]] = s
        y = s["_y"]
        row = f"{d['target']:10s}" + "".join(
            f"{roc_auc_score(y, s[a]):20.3f}" for a in ARMS)
        print(row, flush=True)

    pairs = {
        "both_minus_fingerprint": ("both", "fingerprint"),
        "scrambled_minus_fingerprint": ("profile-scrambled", "fingerprint"),
        "both_minus_scrambled": ("both", "profile-scrambled"),
        "fused_minus_fingerprint": ("fused", "fingerprint"),
        "profile_minus_scram": ("profile", "profile-only-scram"),
        "profile_minus_size": ("profile", "size"),
    }
    labels = {
        "both_minus_fingerprint":
            "PRE-SPECIFIED PRIMARY  (fingerprint + profile, against fingerprint alone)",
        "scrambled_minus_fingerprint":
            "THE CONTROL  (fingerprint + SCRAMBLED profile, against fingerprint alone)",
        "both_minus_scrambled":
            "WHAT THE CONTROL LEAVES  (real profile against scrambled, same column count)\n"
            "  -- if this is a tie, the pre-specified primary measured dilution, not cells",
        "fused_minus_fingerprint":
            "ADDED AFTER THE CONTROL FIRED  (rank-fused, against fingerprint alone)\n"
            "  -- late fusion cannot dilute, so a difference here is the profile's content",
        "profile_minus_scram":
            "IS THE PROFILE ABOVE CHANCE AT ALL  (profile against the same profile permuted)",
        "profile_minus_size":
            "profile alone, against heavy-atom count",
    }

    per_target: dict[str, list[float]] = {k: [] for k in pairs}
    for key, (a, b) in pairs.items():
        print("\n" + "=" * 78)
        print(labels[key] + "\n")
        for t, s in results.items():
            point, lo, hi = paired_auc_ci(s["_y"], s[a], s[b])
            per_target[key].append(point)
            print(f"  {t:10s} {point:+.3f}   95% CI [{lo:+.3f}, {hi:+.3f}]")
            print(f"             -> {verdict(point, lo, hi, a)}")

    print("\n" + "=" * 78)
    print("ACROSS TARGETS\n")
    n = len(results)
    for key in pairs:
        v = np.array(per_target[key])
        wins = int((v > 0).sum())
        print(f"  {key:30s} mean {v.mean():+.3f}  wins {wins}/{n}")
    if n < 5:
        print(f"""
  {n} target(s) cleared the threshold. Section 17 used nine and resolved an effect
  of only +/-0.028 AUC; with {n} there is no across-target interval worth quoting,
  so the per-target intervals above are the whole result. Read the control line
  before the primary one: if the scrambled arm moved as much as the real one, the
  primary number is measuring columns rather than cells.""")


if __name__ == "__main__":
    names = sys.argv[1:] or sorted(
        p.stem.replace("-matched", "") for p in CACHE.glob("*-matched.parquet"))
    if not names:
        sys.exit("no matched files -- run jump_match.py first")
    main(names)
