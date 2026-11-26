#!/usr/bin/env python3
"""Does the platform's model predict an ADME endpoint, or is it reading off cLogP?

Every number this platform reports is a potency number, and RESULTS.md says in its own
words that potency is not what programmes die on. This runs the same model, the same
splits and the same calibration check against two public ADME endpoints instead:

  logd         lipophilicity (logD 7.4), 4,200 compounds from ChEMBL via MoleculeNet
  solubility   aqueous solubility (log mol/L), 1,128 compounds, the Delaney/ESOL set

THE TRAP THIS TEST IS BUILT AROUND
----------------------------------
`dmta.chem.feature_vector` concatenates a Morgan fingerprint with a twelve-element
physchem block, and the second element of that block is Crippen cLogP. For a potency
endpoint that is one weak feature among thousands. For *these two endpoints it is very
close to the answer*: logD is lipophilicity measured, and cLogP is lipophilicity
computed; solubility is dominated by the same quantity with the sign flipped.

So "our model beats cLogP" would be a claim about a model that was handed cLogP. The
honest comparison needs a fourth arm -- the identical forest with the physchem block
removed -- to separate having learned something from having read the baseline off its
own input. Both arms are reported. If the fingerprint-only forest collapses to the
constant predictor, what we have is a cLogP calculator with extra steps, and that is a
finding worth more than a win.

The arms:

  forest            the platform's PropertyModel, exactly as a campaign uses it
  forest_fp_only    the same forest on the fingerprint alone (physchem block removed)
  descriptor_ridge  ridge on the twelve physchem descriptors, no fingerprint at all
  clogp             Crippen cLogP alone, signed for the endpoint -- the trivial answer
  heavy_atoms       counting heavy atoms, the baseline that beat our pharmacophore screen
  mean              the constant predictor, the floor any claim has to clear

Splits are the platform's own `random_split` and `scaffold_split`. There is no time
split: neither public set carries dates. That is worth stating loudly, because the time
split is the column where this platform's uncertainty estimate collapsed (0.44 -> 0.36
-> 0.12, inverting to -0.13 on mu-opioid), and it is the column a production model
lives in. These results therefore cover the two easier columns only.

Usage:  python3 adme_test.py [logd|solubility|both]
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV

from dmta.chem import descriptors, feature_vector
from dmta.ingest import load_csv
from dmta.models import DescriptorRidge, PropertyModel, choose_model, regression_metrics
from dmta.validate import random_split, scaffold_split

RAW = Path("data/adme_raw")
SEEDS = (0, 1, 2, 3, 4)
RESAMPLES = 4000
# The physchem block is the last 12 columns of feature_vector; cLogP is its second
# element. Both facts are asserted against dmta.chem at startup rather than trusted.
N_PHYS = 12
CLOGP_AT = 1      # within the physchem block
HEAVY_AT = -3     # within the whole vector

DATASETS = {
    "logd": {
        "file": "Lipophilicity.csv",
        "smiles": "smiles", "value": "exp", "id": "CMPD_CHEMBLID",
        "label": "logD 7.4 (lipophilicity)",
        # higher cLogP means higher logD: the trivial predictor points the same way
        "clogp_sign": +1.0,
        "source": "MoleculeNet Lipophilicity, curated from ChEMBL",
    },
    "solubility": {
        "file": "delaney-processed.csv",
        "smiles": "smiles", "value": "measured log solubility in mols per litre",
        "id": "Compound ID",
        "label": "aqueous solubility (log mol/L)",
        # more lipophilic means less soluble: the trivial predictor is NEGATED cLogP,
        # and getting this sign wrong would hand the baseline a guaranteed loss
        "clogp_sign": -1.0,
        "source": "MoleculeNet ESOL / Delaney",
    },
}


# --------------------------------------------------------------------- features
def check_feature_layout() -> None:
    """Fail loudly if the physchem block is not where this script thinks it is.

    Everything below depends on being able to strip cLogP from the feature vector.
    If dmta.chem changes shape, a silent mismatch would make the fingerprint-only
    arm a lie -- it would still contain the baseline it is supposed to exclude.
    """
    probe = "CCOc1ccccc1C(=O)O"
    vec = feature_vector(probe)
    d = descriptors(probe)
    phys = vec[-N_PHYS:]
    if not np.isclose(phys[CLOGP_AT], d["clogp"] / 5.0, atol=1e-5):
        raise SystemExit(
            f"feature layout changed: expected cLogP/5 = {d['clogp'] / 5.0:.4f} at "
            f"position -{N_PHYS - CLOGP_AT}, found {phys[CLOGP_AT]:.4f}. "
            "Fix the constants in this file before trusting any arm."
        )
    if not np.isclose(vec[HEAVY_AT], d["heavy_atoms"] / 40.0, atol=1e-5):
        raise SystemExit(
            f"heavy-atom column moved: expected {d['heavy_atoms'] / 40.0:.4f} at "
            f"{HEAVY_AT}, found {vec[HEAVY_AT]:.4f}."
        )


def build_features(compounds) -> tuple[np.ndarray, np.ndarray]:
    X = np.array([feature_vector(c.smiles) for c in compounds], dtype=np.float32)
    return X, X[:, :-N_PHYS]          # full, fingerprint-only


# ------------------------------------------------------------------------ arms
def predict_arms(train, test, Xtr, Xte, fp_tr, fp_te, ytr, sign) -> dict:
    """Every arm's predictions on the same test set, plus the model's sigma."""
    out: dict[str, np.ndarray] = {}

    model = PropertyModel(n_estimators=300, seed=0).fit(train, target="measured")
    mu, sigma = model.predict(test)
    out["forest"] = mu

    # The model actually shipped for this endpoint, not a bare RidgeCV: it carries
    # the bootstrap sigma the acquisition functions need, and sigma is the half of
    # it that the rho and MAE tables below cannot see.
    shipped = DescriptorRidge(seed=0).fit(train, target="measured")
    mu_r, sigma_r = shipped.predict(test)
    out["platform_ridge"] = mu_r

    fp = RandomForestRegressor(n_estimators=300, max_features="sqrt", n_jobs=-1,
                               random_state=0).fit(fp_tr, ytr)
    out["forest_fp_only"] = fp.predict(fp_te)

    phys_tr, phys_te = Xtr[:, -N_PHYS:], Xte[:, -N_PHYS:]
    ridge = RidgeCV(alphas=np.logspace(-3, 3, 13)).fit(phys_tr, ytr)
    out["descriptor_ridge"] = ridge.predict(phys_te)

    # The trivial arms predict a *ranking*, not a value, so they are put on the
    # endpoint's scale by a least-squares fit on the training split only. Without
    # that they would lose on MAE by construction and the comparison would be rigged.
    for name, raw_tr, raw_te in (
        ("clogp", sign * phys_tr[:, CLOGP_AT], sign * phys_te[:, CLOGP_AT]),
        ("heavy_atoms", Xtr[:, HEAVY_AT], Xte[:, HEAVY_AT]),
    ):
        a, b = np.polyfit(raw_tr, ytr, 1)
        out[name] = a * raw_te + b

    out["mean"] = np.full(len(test), float(ytr.mean()))
    return out, sigma, sigma_r, model, shipped


# ------------------------------------------------------------------- statistics
def spearman(p: np.ndarray, o: np.ndarray) -> float | None:
    """None when the prediction is constant: no ranking, not a zero one.

    Running this benchmark is what found the tie-handling bug in dmta.models._rank
    that made exactly this case come back as a small random number. The constant
    arm is kept -- it is the floor any MAE claim has to clear -- but it is scored
    on MAE only, and excluded from every rank comparison below.
    """
    return regression_metrics(list(p), list(o))["spearman"]


def paired_ci(pred_a: np.ndarray, pred_b: np.ndarray, obs: np.ndarray,
              seed: int = 0) -> tuple[float, float, float]:
    """Bootstrap the difference in Spearman rho over TEST COMPOUNDS.

    The platform's other benchmarks resample whole targets, because the question
    there is whether a method generalises across programmes. Here there is one
    dataset, so the resampling unit is the compound and the interval describes
    uncertainty *within this endpoint only*. It says nothing about the next
    endpoint, and the seed-to-seed table below is the other half of the picture.
    """
    rng = np.random.default_rng(seed)
    n = len(obs)
    diffs = []
    for _ in range(RESAMPLES):
        idx = rng.integers(0, n, n)
        a, b = spearman(pred_a[idx], obs[idx]), spearman(pred_b[idx], obs[idx])
        if a is not None and b is not None:
            diffs.append(a - b)
    a, b = spearman(pred_a, obs), spearman(pred_b, obs)
    if a is None or b is None or len(diffs) < RESAMPLES // 2:
        raise ValueError("one arm has no ranking on this split — not comparable by rho")
    lo, hi = np.percentile(np.array(diffs), [2.5, 97.5])
    return a - b, float(lo), float(hi)


def read_verdict(point: float, lo: float, hi: float, what: str) -> str:
    if lo > 0:
        return f"{what} WINS — the interval excludes zero."
    if hi < 0:
        return f"{what} LOSES — the interval excludes zero in the other direction."
    bound = max(abs(lo), abs(hi))
    return (f"A TIE, NOT A WIN. The interval contains zero, so this test did not "
            f"separate them; an effect up to {bound:.3f} of rho, either way, is "
            f"still inside it.")


# ------------------------------------------------------------------------- run
def run(key: str) -> None:
    cfg = DATASETS[key]
    path = RAW / cfg["file"]
    if not path.exists():
        raise SystemExit(f"missing {path} — see the download step in RESULTS.md section 18")

    # units="log" makes ingest pass the values through untouched. Without it the
    # potency path would treat a logD of 3.54 as a concentration and invert it.
    #
    # mw_range is widened, but that does NOT keep every small solute: dmta.chembl
    # .standardize drops anything under 8 heavy atoms, and on ESOL that removes
    # about a fifth of the file (benzene, chloroform and the like). It is left in
    # place because it is what the platform really does to an uploaded CSV, and
    # because those solutes are predictable from size alone -- keeping them would
    # inflate every arm here, the heavy-atom baseline most of all. The count is
    # printed below so the reader can see the size of the exclusion.
    compounds, report = load_csv(
        path, smiles_column=cfg["smiles"], value_column=cfg["value"],
        id_column=cfg["id"], units="log", mw_range=(0.0, 2000.0),
    )
    for c in compounds:
        c.measured = c.truth

    y_all = np.array([c.measured for c in compounds])
    print(f"\n{'=' * 78}\n{key}: {cfg['label']}\n{'=' * 78}")
    print(f"  source     {cfg['source']}")
    print(f"  rows read  {report.rows_read}, usable compounds {len(compounds)}")
    if report.dropped:
        for reason, n in sorted(report.dropped.items(), key=lambda kv: -kv[1]):
            print(f"             dropped {n}: {reason}")
    print(f"  range      {y_all.min():.2f} to {y_all.max():.2f}, sd {y_all.std():.2f}")
    print(f"  trivial arm is {'+' if cfg['clogp_sign'] > 0 else '-'}cLogP")

    t0 = time.time()
    X_all, fp_all = build_features(compounds)
    print(f"  featurised in {time.time() - t0:.0f}s, p = {X_all.shape[1]}")
    index = {id(c): i for i, c in enumerate(compounds)}

    arms = ["forest", "platform_ridge", "forest_fp_only", "descriptor_ridge",
            "clogp", "heavy_atoms", "mean"]
    results: dict[str, dict[str, list[float]]] = {
        split: {a: [] for a in arms} for split in ("random", "scaffold")
    }
    sigma_rank: dict[str, list[float]] = {"random": [], "scaffold": []}
    sigma_rank_ridge: dict[str, list[float]] = {"random": [], "scaffold": []}
    chosen: dict[str, list[str]] = {"random": [], "scaffold": []}
    keep: dict[str, dict] = {}

    for split_name, splitter in (("random", random_split), ("scaffold", scaffold_split)):
        for seed in SEEDS:
            train, test = splitter(compounds, seed=seed)
            if len(test) < 30 or len(train) < 50:
                print(f"  {split_name} seed {seed}: split too small — skipped")
                continue
            tr_i = [index[id(c)] for c in train]
            te_i = [index[id(c)] for c in test]
            ytr = np.array([c.measured for c in train])
            yte = np.array([c.measured for c in test])

            preds, sigma, sigma_r, model, shipped = predict_arms(
                train, test, X_all[tr_i], X_all[te_i], fp_all[tr_i], fp_all[te_i],
                ytr, cfg["clogp_sign"],
            )
            # what the selector would pick, seeing only the training split
            _, pick = choose_model(train, target="measured", seed=0)
            chosen[split_name].append(pick["chosen"])
            for a in arms:
                results[split_name][a].append(spearman(preds[a], yte))

            # the platform's own diagnostic: does sigma rank the error?
            err = np.abs(preds["forest"] - yte)
            sr = spearman(sigma, err)
            if sr is not None:
                sigma_rank[split_name].append(sr)
            err_r = np.abs(preds["platform_ridge"] - yte)
            sr_r = spearman(sigma_r, err_r)
            if sr_r is not None:
                sigma_rank_ridge[split_name].append(sr_r)

            if seed == SEEDS[0]:
                keep[split_name] = {
                    "preds": preds, "obs": yte, "n_train": len(train),
                    "n_test": len(test), "model": model.report(),
                }

    # ------------------------------------------------------------------ tables
    for split_name in ("random", "scaffold"):
        rows = results[split_name]
        if not rows["forest"]:
            continue
        print(f"\nSpearman rho, {split_name} split, {len(rows['forest'])} seeds")
        print(f"  {'arm':<18}{'mean':>8}{'min':>8}{'max':>8}")
        for a in arms:
            v = [x for x in rows[a] if x is not None]
            if not v:
                print(f"  {a:<18}{'—':>8}{'—':>8}{'—':>8}   (constant: no ranking)")
                continue
            v = np.array(v)
            print(f"  {a:<18}{v.mean():>8.3f}{v.min():>8.3f}{v.max():>8.3f}")
        for label, store in (("forest", sigma_rank), ("ridge ", sigma_rank_ridge)):
            v = np.array(store[split_name])
            if not len(v):
                continue
            print(f"  sigma ranks error, {label}: {v.mean():+.3f} mean over seeds "
                  f"({'usable' if v.mean() > 0.2 else 'weak or uninformative'})")
        picks = chosen[split_name]
        print(f"  choose_model picked: "
              f"{', '.join(f'{p} x{picks.count(p)}' for p in sorted(set(picks)))}")

    # -------------------------------------------------------------- the verdicts
    for split_name in ("random", "scaffold"):
        if split_name not in keep:
            continue
        k = keep[split_name]
        obs, preds = k["obs"], k["preds"]
        print(f"\n{'-' * 78}\nVERDICTS on the {split_name} split "
              f"(seed {SEEDS[0]}: {k['n_train']} train / {k['n_test']} test)")
        print(f"  platform model: CV R2 {k['model']['cv_r2']}, "
              f"sigma calibration {k['model']['sigma_calibration']}")
        # "mean" is absent here on purpose: a constant predictor has no ranking,
        # so a rho comparison against it is not defined. It appears in the MAE
        # table below, which is where a constant floor actually means something.
        for rival in ("platform_ridge", "clogp", "descriptor_ridge", "heavy_atoms"):
            p, lo, hi = paired_ci(preds["forest"], preds[rival], obs)
            print(f"\n  forest vs {rival}")
            print(f"    d-rho {p:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}]")
            print(f"    -> {read_verdict(p, lo, hi, 'the forest')}")
        p, lo, hi = paired_ci(preds["forest_fp_only"], preds["clogp"], obs)
        print(f"\n  THE ONE THAT MATTERS: forest WITHOUT the physchem block, vs cLogP")
        print(f"    d-rho {p:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}]")
        print(f"    -> {read_verdict(p, lo, hi, 'the fingerprint-only forest')}")
        p, lo, hi = paired_ci(preds["forest"], preds["forest_fp_only"], obs)
        print(f"\n  how much the physchem block is worth: forest minus forest_fp_only")
        print(f"    d-rho {p:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}]")

        print("\n  MAE, same split (lower is better)")
        for a in arms:
            m = regression_metrics(list(preds[a]), list(obs))
            print(f"    {a:<18}{m['mae']:>8}   R2 {m['r2']}")


def main() -> None:
    check_feature_layout()
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    keys = list(DATASETS) if which == "both" else [which]
    for k in keys:
        if k not in DATASETS:
            raise SystemExit(f"unknown dataset {k!r}; choose from {list(DATASETS)} or 'both'")
        run(k)


if __name__ == "__main__":
    main()
