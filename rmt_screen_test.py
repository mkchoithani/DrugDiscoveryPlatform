"""RMT on the task its authors actually validated: virtual screening.

WHY THIS SCRIPT EXISTS
----------------------
[RESULTS.md](RESULTS.md) section 8 rejected `dmta/rmt_model.py` on 24 ChEMBL
programmes: forest median rho 0.562 against RMT's 0.389, a large and properly
measured effect. [PAPERS.md](PAPERS.md) then recorded, against our own result,
that the rejection does not fall where we put it:

    We used it as a predictor. The paper is about cleaning a covariance matrix,
    and never claimed to beat a random forest at QSAR. The effect we measured is
    large, so rejecting it *as our predictor* is safe -- but its actual
    contribution is untested here.

Lee et al., PNAS 2019, 116, 3373 (DOI 10.1073/pnas.1810847116) clean the
correlations between chemical features computed from known actives AND known
inactives, and use the cleaned matrix to pick candidates out of a library. The
reported success is a virtual screen -- prospectively, four confirmed M1
muscarinic agonists -- not a Spearman correlation against continuous potency.
Section 8 therefore graded the method on a task it was never offered for.

This script runs the task it was offered for, and in the data regime the paper
addresses: a handful of known actives, vastly more features, a large library to
rank. LIT-PCBA supplies it with experimentally measured inactives and
AVE-unbiased splits, so a 2D fingerprint cannot win by recognising catalogue
molecules the way it did on DUD-E (section 12).

WHAT IS BEING COMPARED, AND WHY
-------------------------------
The headline is NOT "RMT versus our forest". A method can be worth less than the
incumbent and still have a real contribution; those are different questions, and
section 16 is the record of what conflating them costs. The headline here is:

    rmt-cleaned  versus  rmt-uncleaned

Both arms take the identical path through `RmtModel` -- same features, same
standardisation, same generalised-least-squares solve, same prediction. The only
difference is whether the Marchenko-Pastur clipping happens. Whatever separates
them IS the cleaning, which is the paper's claim and nothing else.

The control needs a ridge to exist at all (p > n makes the raw correlation
matrix singular), and that ridge is the ad-hoc regularisation RMT is proposed as
a principled replacement for. Choosing it badly would hand the comparison to
RMT by default, so it is chosen by cross-validation on the TRAINING split only,
per target, and the chosen value is printed.

Two further arms are reported because they bound the result from each side:

  ols      plain least squares on the same features -- the contrast the paper
           itself draws, and the thing that explodes in this regime
  forest   our incumbent, so the cost decision can be read off the same table
  size     heavy atom count, the null that has already beaten a pharmacophore
           here (section 10). Any arm that loses to it has found nothing.

Mode is "clip", not "project". Section 8 chose "project" because it ranks better,
and said so. But clipping inverts the cleaned matrix literally, which whitens --
and a whitened linear solve against 0/1 labels is a discriminant, which is what
a screen is. The mode section 8 set aside as "right for a discriminant, wrong for
ranking" is the right one here. Both are run anyway, because asserting that
without measuring it would be the same error in a new place.

    python3 rmt_screen_test.py <litpcba_ave_dir> [max_inactives]
"""
from __future__ import annotations

import sys
import time
from math import comb
from pathlib import Path

import numpy as np

from dmta.chem import feature_vector
from dmta.hitfinding import enrichment_factor, roc_auc
from dmta.models import PropertyModel
from dmta.rmt_model import RmtModel
from litpcba_benchmark import read_split
from model_paired_test import sign_test, signed_rank_p

CACHE = Path("data/rmt_screen_cache")
# Every arm's scores are written here per target. The models take ~15 minutes over
# the 15 targets, and the first run of this script threw away all of it to a KeyError
# in a print statement. Scoring and reporting are now separable: fix the report, rerun,
# get the tables back in seconds. Delete the directory to force a real re-score.

RESAMPLES = 4000           # within-target AUC bootstrap
PERMUTATIONS = 20000       # across-target signed-rank permutation
MIN_TEST_ACTIVES = 20      # below this an EF is a coin flip, per section 15
RIDGES = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
# The grid must bracket the optimum at both ends. A control pinned to the largest
# ridge on offer is a control we handicapped, and the benchmark prints the chosen
# value per target so that can be checked rather than assumed.


# ----------------------------------------------------------------- statistics
def paired_auc_ci(a: np.ndarray, b: np.ndarray, labels: np.ndarray,
                  rng: np.random.Generator) -> tuple[float, float, float]:
    """Bootstrap CI for AUC(a) - AUC(b), paired and stratified.

    Paired: both arms are resampled on the SAME compounds every draw. They score
    the same molecules, so their errors are correlated, and an unpaired interval
    would describe the variance of each arm rather than of their difference --
    far too wide, and it would turn every real gap into "no evidence".

    Stratified: actives and inactives are resampled separately, holding the class
    ratio fixed. With a 1% hit rate an unstratified draw varies the number of
    actives from draw to draw, which moves the AUC for a reason that has nothing
    to do with either method.
    """
    pos = np.flatnonzero(labels)
    neg = np.flatnonzero(~labels)
    observed = roc_auc(a, labels) - roc_auc(b, labels)
    draws = np.empty(RESAMPLES)
    for i in range(RESAMPLES):
        idx = np.concatenate([rng.choice(pos, pos.size, replace=True),
                              rng.choice(neg, neg.size, replace=True)])
        lab = labels[idx]
        draws[i] = roc_auc(a[idx], lab) - roc_auc(b[idx], lab)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return observed, float(lo), float(hi)


def actives_in_top(scores: np.ndarray, labels: np.ndarray,
                   fraction: float = 0.01) -> tuple[int, int, float]:
    """Count of actives in the top slice, and P(at least that many) by chance.

    An enrichment factor is a ratio, and a ratio built from three molecules looks
    identical to one built from three hundred. Section 15 was decided on this
    conversion and not on the ratio, so it is computed here rather than left as a
    footnote.
    """
    k = max(1, int(round(len(scores) * fraction)))
    order = np.argsort(-scores, kind="stable")
    found = int(labels[order[:k]].sum())
    total, act = len(scores), int(labels.sum())
    p = sum(comb(act, i) * comb(total - act, k - i) / comb(total, k)
            for i in range(found, min(act, k) + 1)) if found else 1.0
    return found, k, float(p)


# --------------------------------------------------------------------- arms
def choose_ridge(X: np.ndarray, y: np.ndarray, mode: str,
                 seed: int = 0, folds: int = 3) -> float:
    """Pick the control's ridge by AUC on held-out folds of the TRAINING split.

    Never on the test split, and never by hand: the control's whole job is to be
    the best version of "no cleaning", so a ridge chosen to make it look bad
    would make the comparison meaningless.
    """
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(y))
    cuts = np.array_split(order, folds)
    best, best_auc = RIDGES[0], -np.inf
    for ridge in RIDGES:
        aucs = []
        for held in cuts:
            keep = np.setdiff1d(order, held)
            if y[keep].sum() < 2 or y[held].sum() < 1:
                continue
            m = RmtModel(mode=mode, clean=False, shrinkage=ridge)
            m._fit_matrix(X[keep], y[keep])
            mu, _ = m._predict_matrix(X[held])
            aucs.append(roc_auc(np.asarray(mu, float), y[held] == 1.0))
        if aucs and float(np.mean(aucs)) > best_auc:
            best, best_auc = ridge, float(np.mean(aucs))
    return best


def run_target(directory: Path, max_inactives: int, seed: int = 0) -> dict | None:
    name = directory.name
    cached = CACHE / f"{name}-{max_inactives}.npz"
    if cached.exists():
        blob = np.load(cached, allow_pickle=False)
        labels = blob["labels"].astype(bool)
        scores = {k: blob[f"s_{k}"] for k in blob["arms"]}
        timing = {k: float(v) for k, v in zip(blob["arms"], blob["seconds"])}
        print(f"{name}: {len(labels)} compounds, {int(labels.sum())} actives "
              f"(from cache)", flush=True)
        return summarise(name, scores, labels, timing, float(blob["ridge"]), seed)

    t0 = time.time()
    act_t = read_split(directory, "active_T.smi", True, None)
    act_v = read_split(directory, "active_V.smi", True, None)
    ina_t = read_split(directory, "inactive_T.smi", False, max_inactives, seed)
    ina_v = read_split(directory, "inactive_V.smi", False, max_inactives, seed)
    if len(act_v) < MIN_TEST_ACTIVES:
        print(f"{name}: only {len(act_v)} actives in the test split "
              f"(need {MIN_TEST_ACTIVES}) -- skipped, not counted as a loss")
        return None

    train, test = act_t + ina_t, act_v + ina_v
    ytr = np.array([c.truth for c in train], float)
    labels = np.array([c.truth == 1.0 for c in test])
    Xtr = np.array([feature_vector(c.smiles) for c in train], float)
    Xte = np.array([feature_vector(c.smiles) for c in test], float)
    print(f"{name}: train {len(act_t)}a/{len(ina_t)}i  "
          f"test {len(act_v)}a/{len(ina_v)}i  "
          f"hit rate {labels.mean():.2%}  p={Xtr.shape[1]}  "
          f"featurised in {time.time()-t0:.0f}s", flush=True)

    scores: dict[str, np.ndarray] = {}
    timing: dict[str, float] = {}

    for mode in ("clip", "project"):
        t = time.time()
        m = RmtModel(mode=mode, clean=True).fit_matrix_only(Xtr, ytr)
        scores[f"rmt-clean-{mode}"] = np.asarray(m._predict_matrix(Xte)[0], float)
        timing[f"rmt-clean-{mode}"] = time.time() - t

    ridge = choose_ridge(Xtr, ytr, "clip", seed)
    t = time.time()
    m = RmtModel(mode="clip", clean=False, shrinkage=ridge).fit_matrix_only(Xtr, ytr)
    scores["rmt-uncleaned"] = np.asarray(m._predict_matrix(Xte)[0], float)
    timing["rmt-uncleaned"] = time.time() - t
    print(f"  control ridge chosen on the training split: {ridge}", flush=True)

    t = time.time()
    Xd = np.hstack([Xtr, np.ones((len(Xtr), 1))])
    w, *_ = np.linalg.lstsq(Xd, ytr, rcond=None)
    scores["ols"] = np.hstack([Xte, np.ones((len(Xte), 1))]) @ w
    timing["ols"] = time.time() - t

    t = time.time()
    for c in train:
        c.measured = c.truth
    try:
        mu, _ = PropertyModel(n_estimators=300).fit(train).predict(test)
    finally:
        for c in train:
            c.measured = None
    scores["forest"] = np.asarray(mu, float)
    timing["forest"] = time.time() - t

    scores["size"] = np.array([c.desc["heavy_atoms"] for c in test], float)
    timing["size"] = 0.0

    CACHE.mkdir(parents=True, exist_ok=True)
    arms = list(scores)
    np.savez_compressed(
        cached, labels=labels, arms=np.array(arms), ridge=np.array(ridge),
        seconds=np.array([timing[a] for a in arms]),
        **{f"s_{a}": scores[a] for a in arms})

    return summarise(name, scores, labels, timing, ridge, seed)


def summarise(name: str, scores: dict[str, np.ndarray], labels: np.ndarray,
              timing: dict[str, float], ridge: float, seed: int) -> dict | None:
    """Everything downstream of the models: filtering, metrics and the intervals."""
    test_n = len(labels)
    usable = np.ones(test_n, bool)
    for v in scores.values():
        usable &= np.isfinite(v)        # one compound set, or it is not paired
    dropped = int((~usable).sum())
    if dropped:
        print(f"  {dropped} compounds dropped: not every arm scored them, and a "
              f"paired comparison needs one common set", flush=True)
    labels = labels[usable]
    scores = {k: v[usable] for k, v in scores.items()}
    if labels.sum() < MIN_TEST_ACTIVES:
        print(f"  {name}: too few actives survived -- skipped")
        return None

    rng = np.random.default_rng(seed)
    rows = {}
    for key, v in scores.items():
        found, k, p = actives_in_top(v, labels)
        rows[key] = {
            "auc": roc_auc(v, labels),
            "ef1": enrichment_factor(v, labels, 0.01),
            "found": found, "top_k": k, "p_chance": p, "seconds": timing[key],
        }
    # Both modes are carried against both rivals. "clip" was pre-specified as the
    # headline on a theoretical argument -- a screen is a discriminant, and clipping
    # is the discriminant solve. Reporting only the pre-specified arm when the other
    # one is plainly better would overstate the loss; reporting only the better one
    # would hide that the prediction failed. So both are printed, labelled.
    pairs = {
        "clip_minus_uncleaned": ("rmt-clean-clip", "rmt-uncleaned"),
        "project_minus_uncleaned": ("rmt-clean-project", "rmt-uncleaned"),
        "clip_minus_forest": ("rmt-clean-clip", "forest"),
        "project_minus_forest": ("rmt-clean-project", "forest"),
    }
    out = {"target": name, "rows": rows, "n": int(len(labels)),
           "n_actives": int(labels.sum()), "ridge": ridge}
    for key, (a, b) in pairs.items():
        out[key] = paired_auc_ci(scores[a], scores[b], labels, rng)
    return out


# ---------------------------------------------------------------------- main
def read_verdict(diffs: np.ndarray, lo: float, hi: float, wilcoxon: float,
                 sign_p: float, label: str) -> list[str]:
    """Three-valued, as section 16 requires: a win, a loss, or a tie.

    A tie is not a rejection. If the interval contains zero, the only honest
    statement is that this experiment did not separate the two, together with how
    large an effect would still fit inside the interval.
    """
    mean = float(diffs.mean())
    out = [f"  mean difference {mean:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}], "
           f"Wilcoxon P = {wilcoxon:.3f}, sign test P = {sign_p:.3f}"]
    if lo > 0:
        out.append(f"  -> {label} WINS. The interval excludes zero, so the "
                   f"cleaning earned its place on this task.")
    elif hi < 0:
        out.append(f"  -> {label} LOSES. The interval excludes zero in the other "
                   f"direction.")
    else:
        bound = max(abs(lo), abs(hi))
        out.append(f"  -> A TIE, NOT A REJECTION. The interval contains zero, so "
                   f"this test did not separate them.")
        out.append(f"     An effect up to {bound:.3f} of AUC, either way, is still "
                   f"inside the interval: this is a statement about what we "
                   f"could resolve with {len(diffs)} targets, not a measurement "
                   f"of no effect.")
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    root = Path(sys.argv[1])
    max_inactives = int(sys.argv[2]) if len(sys.argv) > 2 else 2000

    dirs = sorted(d for d in root.iterdir()
                  if d.is_dir() and (d / "active_V.smi").exists())
    print(f"RMT as a virtual screen: {len(dirs)} LIT-PCBA targets, inactives "
          f"capped at {max_inactives} per split\n")

    results = []
    for d in dirs:
        try:
            r = run_target(d, max_inactives)
        except Exception as exc:                       # a target must not kill the run
            print(f"{d.name}: FAILED -- {type(exc).__name__}: {exc}", flush=True)
            continue
        if r:
            results.append(r)
            print(f"  done in {sum(v['seconds'] for v in r['rows'].values()):.0f}s "
                  f"of model time\n", flush=True)
    if len(results) < 3:
        sys.exit(f"only {len(results)} targets produced a result -- too few to "
                 f"say anything")

    arms = list(results[0]["rows"])
    col = max(14, max(len(a) for a in arms) + 2)
    print("\n" + "=" * 78)
    print("AUC per target\n")
    print(f"{'target':<12}{'n':>7}{'act':>5}" + "".join(f"{a:>{col}}" for a in arms))
    for r in results:
        print(f"{r['target']:<12}{r['n']:>7}{r['n_actives']:>5}"
              + "".join(f"{r['rows'][a]['auc']:>{col}.3f}" for a in arms))

    print("\nActives in the top 1%, with the chance probability of that many\n")
    print(f"{'target':<12}{'top 1%':>8}" + "".join(f"{a:>{col}}" for a in arms))
    for r in results:
        cells = "".join(
            f"{str(r['rows'][a]['found']) + ' (p=' + format(r['rows'][a]['p_chance'], '.3f') + ')':>{col}}"
            for a in arms)
        print(f"{r['target']:<12}{results[0]['rows'][arms[0]]['top_k']:>8}" + cells)

    print("\n  The inactive cap keeps every active but samples the inactives, so the")
    print("  hit rates above are inflated -- ALDH1 reaches 40%, where a \"top 1%\" is")
    print("  barely a screen at all. AUC is a rank statistic over active/inactive")
    print("  pairs and is unaffected by that ratio, which is why the verdict below")
    print("  is read off AUC and the top-1% counts are context only.")
    print("\nWithin-target paired AUC intervals, cleaned (clip) minus uncleaned\n")
    print(f"{'target':<12}{'ridge':>8}{'difference':>14}{'95% CI':>24}")
    for r in results:
        obs, lo, hi = r["clip_minus_uncleaned"]
        print(f"{r['target']:<12}{r['ridge']:>8}{obs:>+14.3f}"
              f"{f'[{lo:+.3f}, {hi:+.3f}]':>24}")
    if any(r["ridge"] == max(RIDGES) for r in results):
        pinned = [r["target"] for r in results if r["ridge"] == max(RIDGES)]
        print(f"\n  WARNING: the control's ridge ran to the top of the grid on "
              f"{', '.join(pinned)}.")
        print(f"  On those targets the control may be handicapped by the grid "
              f"rather than beaten.")

    rng = np.random.default_rng(0)
    print("\n" + "=" * 78)
    for label, key, lhs, rival in (
            ("THE PAPER'S CLAIM, pre-specified arm", "clip_minus_uncleaned",
             "rmt-clean-clip", "rmt-uncleaned"),
            ("THE PAPER'S CLAIM, better arm", "project_minus_uncleaned",
             "rmt-clean-project", "rmt-uncleaned"),
            ("pre-specified arm against our incumbent", "clip_minus_forest",
             "rmt-clean-clip", "forest"),
            ("better arm against our incumbent", "project_minus_forest",
             "rmt-clean-project", "forest")):
        diffs = np.array([r[key][0] for r in results])
        # The unit of resampling is the TARGET, not the compound. A method that wins
        # on one huge target and loses on eight small ones has not won; bootstrapping
        # compounds would let the huge one carry the result.
        draws = np.array([rng.choice(diffs, diffs.size, replace=True).mean()
                          for _ in range(RESAMPLES)])
        lo, hi = (float(x) for x in np.percentile(draws, [2.5, 97.5]))
        wins = int((diffs > 0).sum())
        print(f"\n{label}  ({lhs} minus {rival}, over {len(diffs)} targets)")
        print(f"  wins {wins}/{len(diffs)}")
        for line in read_verdict(diffs, lo, hi,
                                 signed_rank_p(diffs, rng, PERMUTATIONS),
                                 sign_test(wins, len(diffs)),
                                 "cleaning" if "uncleaned" in key else "RMT"):
            print(line)

    print("\nFit time per target, median seconds")
    for a in arms:
        print(f"  {a:<22}{np.median([r['rows'][a]['seconds'] for r in results]):>8.1f}")


if __name__ == "__main__":
    main()
