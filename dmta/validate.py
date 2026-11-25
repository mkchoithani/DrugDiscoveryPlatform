"""Honest validation on real data.

Three splits, in increasing order of difficulty and realism:

  random    — the number everyone quotes. Optimistic, because analogues of a test
              compound sit in the training set.
  scaffold  — whole Bemis–Murcko scaffolds held out. Asks: can it generalise to a
              chemotype it has not seen?
  time      — train on everything published before a cutoff, test on what came
              after. This is the only split that matches how the model is used:
              predicting compounds nobody has made yet.

A platform that reports only the random split is lying to its users by omission.
"""
from __future__ import annotations

import random
from collections import defaultdict

import numpy as np

from .library import Compound
from .models import PropertyModel, regression_metrics


def random_split(compounds: list[Compound], test_frac: float = 0.2, seed: int = 0):
    idx = list(range(len(compounds)))
    random.Random(seed).shuffle(idx)
    cut = int(len(idx) * (1 - test_frac))
    return [compounds[i] for i in idx[:cut]], [compounds[i] for i in idx[cut:]]


def scaffold_split(compounds: list[Compound], test_frac: float = 0.2, seed: int = 0):
    """Hold out whole scaffolds, largest groups first into train (standard practice)."""
    groups: dict[str, list[Compound]] = defaultdict(list)
    for c in compounds:
        groups[c.scaffold_id].append(c)
    ordered = sorted(groups.values(), key=len, reverse=True)
    rng = random.Random(seed)
    rng.shuffle(ordered := ordered[:])          # shuffle ties, keep rough size order
    ordered.sort(key=len, reverse=True)
    n_test_target = int(len(compounds) * test_frac)
    train, test = [], []
    for group in ordered:
        (test if len(test) < n_test_target and len(group) < 50 else train).extend(group)
    return train, test


def time_split(compounds: list[Compound], cutoff_year: int | None = None, test_frac: float = 0.2):
    dated = [c for c in compounds if getattr(c, "year", None)]
    if len(dated) < 50:
        raise ValueError("not enough dated compounds for a time split")
    if cutoff_year is None:
        years = sorted(c.year for c in dated)
        cutoff_year = years[int(len(years) * (1 - test_frac))]
    train = [c for c in dated if c.year < cutoff_year]
    test = [c for c in dated if c.year >= cutoff_year]
    return train, test, cutoff_year


def evaluate(train: list[Compound], test: list[Compound], n_estimators: int = 300,
             threshold: float | None = None) -> dict:
    """threshold defines "active". Left as None it is the 90th percentile of the
    TRAINING values — dataset-relative, and set without looking at the test set.
    A hardcoded cut (say pIC50 >= 8) silently reports 0% on any programme whose
    compounds never reach it, which is most of them."""
    for c in train:
        c.measured = c.truth
    model = PropertyModel(n_estimators=n_estimators).fit(train)
    mu, sigma = model.predict(test)
    truth = [c.truth for c in test]
    metrics = regression_metrics([float(x) for x in mu], truth)
    metrics["n_train"] = len(train)
    metrics["n_test"] = len(test)
    metrics["mean_sigma"] = round(float(sigma.mean()), 3)
    # does the uncertainty actually rank errors? (the property that makes it useful)
    errors = np.abs(np.array(mu) - np.array(truth))
    metrics["sigma_error_spearman"] = _spearman(sigma, errors)
    # top-decile enrichment: of the 10% the model likes most, how many are truly active?
    if threshold is None:
        threshold = float(np.percentile([c.truth for c in train], 90))
    metrics["active_threshold"] = round(threshold, 3)
    k = max(1, len(test) // 10)
    top_idx = np.argsort(mu)[::-1][:k]
    active = np.array(truth) >= threshold
    metrics["top_decile_hit_rate"] = round(float(active[top_idx].mean()), 3)
    metrics["base_hit_rate"] = round(float(active.mean()), 3)
    metrics["enrichment"] = (round(metrics["top_decile_hit_rate"] / metrics["base_hit_rate"], 2)
                             if metrics["base_hit_rate"] > 0 else None)
    for c in train:
        c.measured = None
    return metrics


def _spearman(a, b) -> float | None:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 3:
        return None
    ra, rb = _rank(a), _rank(b)
    sd = ra.std() * rb.std()
    return round(float(((ra - ra.mean()) * (rb - rb.mean())).mean() / sd), 3) if sd > 1e-9 else None


def _rank(a: np.ndarray) -> np.ndarray:
    order = a.argsort()
    r = np.empty_like(order, dtype=float)
    r[order] = np.arange(len(a), dtype=float)
    return r


def full_report(compounds: list[Compound], n_estimators: int = 300, verbose: bool = True) -> dict:
    out = {}
    tr, te = random_split(compounds)
    out["random"] = evaluate(tr, te, n_estimators)
    tr, te = scaffold_split(compounds)
    out["scaffold"] = evaluate(tr, te, n_estimators)
    try:
        tr, te, cutoff = time_split(compounds)
        out["time"] = evaluate(tr, te, n_estimators)
        out["time"]["cutoff_year"] = cutoff
    except ValueError as e:
        out["time"] = {"error": str(e)}

    if verbose:
        print(f"{'split':<10}{'n_train':>9}{'n_test':>8}{'R2':>7}{'MAE':>7}{'rho':>7}"
              f"{'sigma~err':>11}{'top-10% hits':>14}{'enrich':>8}")
        for name, m in out.items():
            if "error" in m:
                print(f"{name:<10}  {m['error']}")
                continue
            hits = "{:.0%} vs {:.0%}".format(m["top_decile_hit_rate"], m["base_hit_rate"])
            print(f"{name:<10}{m['n_train']:>9}{m['n_test']:>8}"
                  f"{_f(m['r2']):>7}{_f(m['mae']):>7}{_f(m['spearman']):>7}"
                  f"{_f(m['sigma_error_spearman']):>11}{hits:>14}{_f(m['enrichment']):>8}")
        for name, m in out.items():
            if "error" not in m:
                label, note = sigma_verdict(m["sigma_error_spearman"])
                m["sigma_verdict"] = label
                print(f"  σ on the {name} split: {label} — {note}")
        beta, why = recommend_beta(out)
        out["recommended_beta"] = beta
        print(f"\n  recommended exploration β for this target: {beta} ({why})")
    return out


SIGMA_VERDICTS = (
    (0.30, "usable", "σ ranks error well enough to explore on"),
    (0.10, "weak", "σ carries little information about error"),
    (-9.9, "unreliable", "σ is uninformative or inverted — do not explore on it"),
)


def sigma_verdict(rho: float | None) -> tuple[str, str]:
    """How much the uncertainty estimate can be trusted to rank errors.

    This matters because every acquisition function except `greedy` spends part of
    the batch on σ. If σ does not rank error, that spend buys nothing.
    """
    if rho is None:
        return "unknown", "not enough data"
    for threshold, label, note in SIGMA_VERDICTS:
        if rho >= threshold:
            return label, note
    return "unreliable", "σ is inverted"


def recommend_beta(report: dict) -> tuple[float, str]:
    """A starting exploration weight, grounded only in what has been measured.

    Two separate measurements inform this, and they are not known to be causally
    linked — do not let the tidy story substitute for your own sweep:

    1. σ ranks error at ρ≈0.44 on a random split, ≈0.36 on a scaffold split and
       ≈0.12 on a time split (inverting on one of four targets). A production model
       sits in the third regime.
    2. On EGFR at a 2.7% budget, exploration hurt monotonically: β=0 found 117
       actives, β=1 found 75. That is a budget effect as much as a σ effect.

    So the ceiling here is 0.25, not 1.0, and the honest advice is always to sweep.
    """
    rhos = [(report.get(k) or {}).get("sigma_error_spearman")
            for k in ("scaffold", "time")]
    rhos = [r for r in rhos if r is not None]
    worst = min(rhos) if rhos else None
    label, _ = sigma_verdict(worst)
    if label == "usable":
        return 0.25, ("σ still ranks error under shift, so a little exploration may pay — "
                      "but sweep β, it was worse than 0 on EGFR at a 2.7% budget")
    if label == "weak":
        return 0.1, "σ is weakly informative under shift — explore barely, and sweep β"
    return 0.0, "σ does not rank error under shift — exploring on it buys nothing; use greedy or EI"


def _f(v) -> str:
    return "—" if v is None else f"{v:.3f}" if isinstance(v, float) and abs(v) < 10 else str(v)


# --------------------------------------------------------- docking validation

def docking_report(compounds: list[Compound], threshold: float = 8.0) -> dict:
    """Does the docking score actually track measured potency?

    Docking returns kcal/mol where LOWER is a better predicted fit, while higher
    pChEMBL means more potent — so a working scorer gives a NEGATIVE correlation.
    It is reported sign-flipped, so positive always means 'agrees with reality'.
    """
    scored = [c for c in compounds if c.docking_score is not None and c.truth is not None]
    if len(scored) < 10:
        return {"error": f"only {len(scored)} compounds have both a score and a measurement"}

    dock = np.array([c.docking_score for c in scored], dtype=float)
    truth = np.array([c.truth for c in scored], dtype=float)
    active = truth >= threshold

    rho = _spearman(-dock, truth)                     # sign-flipped: positive = good
    k = max(1, len(scored) // 10)
    best_by_dock = np.argsort(dock)[:k]               # most negative = best predicted
    top_rate = float(active[best_by_dock].mean())
    base_rate = float(active.mean())

    return {
        "n": len(scored),
        "spearman": rho,
        "top_decile_hit_rate": round(top_rate, 3),
        "base_hit_rate": round(base_rate, 3),
        "enrichment": round(top_rate / base_rate, 2) if base_rate > 0 else None,
        "mean_score_actives": round(float(dock[active].mean()), 2) if active.any() else None,
        "mean_score_inactives": round(float(dock[~active].mean()), 2) if (~active).any() else None,
        "separation": (round(float(dock[~active].mean() - dock[active].mean()), 2)
                       if active.any() and (~active).any() else None),
        "score_range": (round(float(dock.min()), 2), round(float(dock.max()), 2)),
    }


def docking_report_by_mechanism(compounds: list[Compound], threshold: float = 8.0) -> dict:
    """The same report, split by whether a compound binds covalently.

    Standard docking models reversible binding. Covalent binders form a permanent
    bond, so their measured potency reflects chemistry the score cannot see. If
    docking works at all, it should work better on the reversible subset — and if
    it does not, the whole result is suspect.
    """
    from rdkit import Chem
    warheads = [Chem.MolFromSmarts(s) for s in
                ("[NX3][CX3](=O)[CX3]=[CX3]", "[NX3][CX3](=O)C#C", "[NX3][CX3](=O)[CH2]Cl")]

    covalent, reversible = [], []
    for c in compounds:
        mol = Chem.MolFromSmiles(c.smiles)
        if mol is None:
            continue
        (covalent if any(mol.HasSubstructMatch(w) for w in warheads if w)
         else reversible).append(c)

    return {
        "all": docking_report(compounds, threshold),
        "reversible": docking_report(reversible, threshold),
        "covalent": docking_report(covalent, threshold),
        "n_covalent": len(covalent),
        "n_reversible": len(reversible),
    }


def print_docking_report(report: dict) -> None:
    print(f"{'subset':<14}{'n':>6}{'rho':>8}{'top-10%':>10}{'base':>8}{'enrich':>9}{'separation':>12}")
    print("-" * 67)
    for name in ("all", "reversible", "covalent"):
        m = report.get(name, {})
        if "error" in m or not m:
            print(f"{name:<14}  {m.get('error', 'no data')}")
            continue
        sep = "—" if m["separation"] is None else f"{m['separation']:+.2f}"
        print(f"{name:<14}{m['n']:>6}{_f(m['spearman']):>8}"
              f"{m['top_decile_hit_rate']:>10.0%}{m['base_hit_rate']:>8.0%}"
              f"{_f(m['enrichment']):>9}{sep:>12}")
