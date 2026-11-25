"""Compare models on identical data: does anything beat the random forest?

Paired comparison across many real programmes. Every model sees exactly the same
scaffold split of exactly the same compounds, so the difference between them is the
model and nothing else. One dataset could not settle this; dozens can.

Contenders: the random forest (incumbent), Chemprop's D-MPNN, and the RMT-denoised
linear model of Lee et al., PNAS 2019. Chemprop is skipped automatically when it is
not installed, so the RMT comparison runs without a torch install.

The 2 s.e. rule below is a screen, not the verdict. `--save results.json` writes the
per-programme scores so `model_paired_test.py` can bootstrap the paired difference
properly, without refitting anything.

    python3 model_benchmark.py [n_programmes] [model,model,...] [--save results.json]
"""
from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path

from dmta.chembl import programmes
from dmta.models import PropertyModel, regression_metrics
from dmta.validate import _spearman, scaffold_split

TARGETS = [("CHEMBL203", "IC50"), ("CHEMBL279", "IC50"), ("CHEMBL1824", "IC50"),
           ("CHEMBL244", "IC50"), ("CHEMBL205", "Ki"), ("CHEMBL233", "Ki")]
MIN_SPREAD = 1.5


def evaluate_model(model, train, test) -> dict:
    """Train on `train`, score `test`. Identical protocol for every model."""
    for c in train:
        c.measured = c.truth
    try:
        model.fit(train)
        mu, sigma = model.predict(test)
    finally:
        for c in train:
            c.measured = None
    truth = [c.truth for c in test]
    m = regression_metrics([float(x) for x in mu], truth)

    import numpy as np
    k = max(1, len(test) // 10)
    top = np.argsort(mu)[::-1][:k]
    threshold = float(np.percentile([c.truth for c in train], 90))
    active = np.array(truth) >= threshold
    base = float(active.mean())
    m["enrichment"] = round(float(active[top].mean()) / base, 2) if base else None
    m["sigma_ranks_error"] = _spearman(sigma, np.abs(np.array(mu) - np.array(truth)))
    return m


def contenders(wanted: list[str] | None) -> list[tuple[str, callable]]:
    """(name, factory) for every model that can run here. Order sets the columns."""
    from dmta.rmt_model import RmtModel
    entries = [("forest", lambda: PropertyModel(n_estimators=300)),
               ("rmt-project", lambda: RmtModel(mode="project")),
               ("rmt-clip", lambda: RmtModel(mode="clip"))]
    try:
        from dmta.chemprop_model import ChempropModel, chemprop_available
        if chemprop_available():
            entries.append(("chemprop", lambda: ChempropModel(epochs=40)))
        else:
            print("chemprop not installed — skipping it\n")
    except ImportError:
        print("chemprop_model not importable — skipping it\n")
    if wanted:
        known = {n for n, _ in entries}
        missing = [w for w in wanted if w not in known]
        if missing:
            sys.exit(f"unknown model(s): {', '.join(missing)}; have {sorted(known)}")
        entries = [e for e in entries if e[0] in wanted]
    return entries


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    limit = int(args[0]) if args else 24
    wanted = args[1].split(",") if len(args) > 1 else None
    save_to = None
    if "--save" in sys.argv:
        i = sys.argv.index("--save")
        if i + 1 >= len(sys.argv):
            sys.exit("--save needs a filename")
        save_to = Path(sys.argv[i + 1])
    models = contenders(wanted)
    names = [n for n, _ in models]

    pool = []
    for target, std in TARGETS:
        try:
            pool += [(target, p) for p in programmes(target, std, min_size=60)
                     if p["spread"] >= MIN_SPREAD]
        except Exception as exc:
            print(f"  skipping {target}: {exc}")
    pool.sort(key=lambda tp: -tp[1]["n"])
    pool = pool[:limit]
    print(f"comparing on {len(pool)} programmes "
          f"from {len({t for t, _ in pool})} targets\n")

    col = max(16, max(len(n) for n in names) + 5)   # "rmt-project rho" is 15 wide
    header = (f"{'paper':<16}{'n':>5}"
              + "".join(f"{name + ' rho':>{col}}" for name in names)
              + f"{'winner':>10}")
    print(header); print("-" * len(header))

    rows, t0 = [], time.time()
    timing = {n: 0.0 for n in names}
    for target, prog in pool:
        train, test = scaffold_split(prog["compounds"], test_frac=0.25)
        if len(train) < 30 or len(test) < 10:
            continue
        scores = {}
        for name, factory in models:
            t1 = time.time()
            scores[name] = evaluate_model(factory(), train, test)
            timing[name] += time.time() - t1
        winner = max(names, key=lambda n: scores[n]["spearman"] if scores[n]["spearman"] is not None else -9)
        rows.append({"document": prog["document"], "n": prog["n"],
                     "scores": scores, "winner": winner})
        print(f"{prog['document']:<16}{prog['n']:>5}"
              + "".join(f"{_f(scores[n]['spearman']):>{col}}" for n in names)
              + f"{winner:>10}", flush=True)

    if not rows:
        sys.exit("no programme produced a usable split")

    def column(name, metric):
        return [r["scores"][name][metric] for r in rows
                if r["scores"][name][metric] is not None]

    print(f"\n{len(rows)} programmes in {time.time() - t0:.0f}s\n")
    print(f"  {'model':<12}{'median rho':>12}{'median MAE':>12}{'wins':>8}{'fit seconds':>14}")
    for name in names:
        wins = sum(1 for r in rows if r["winner"] == name)
        print(f"  {name:<12}{statistics.median(column(name, 'spearman')):>12.3f}"
              f"{statistics.median(column(name, 'mae')):>12.3f}"
              f"{wins:>8}{timing[name]:>14.1f}")

    # Paired differences say more than two medians: is a gap consistent, or noise?
    print("\n  paired against the forest, programme by programme:")
    for name in names:
        if name == "forest":
            continue
        paired = [r["scores"][name]["spearman"] - r["scores"]["forest"]["spearman"]
                  for r in rows
                  if r["scores"][name]["spearman"] is not None
                  and r["scores"]["forest"]["spearman"] is not None]
        if len(paired) < 2:
            continue
        med = statistics.median(paired)
        spread = statistics.pstdev(paired)
        stderr = spread / (len(paired) ** 0.5)
        verdict = ("no reliable difference" if abs(med) < 2 * stderr
                   else ("BETTER than the forest" if med > 0 else "WORSE than the forest"))
        print(f"    {name:<10} median {med:+.3f}   spread {spread:.3f}"
              f"   2 s.e. {2 * stderr:.3f}   -> {verdict}")

    if save_to is not None:
        # Every programme's scores, so the paired test never has to refit a model.
        payload = {
            "models": names,
            "programmes": [
                {"document": r["document"], "n": r["n"],
                 "scores": {n: {k: v for k, v in r["scores"][n].items()
                                if isinstance(v, (int, float, type(None)))}
                            for n in names}}
                for r in rows],
        }
        save_to.write_text(json.dumps(payload, indent=1))
        print(f"\n  per-programme scores -> {save_to}")
        print(f"  now run: python3 model_paired_test.py {save_to}")


def _f(v, dp=3) -> str:
    return "—" if v is None else f"{v:.{dp}f}"


if __name__ == "__main__":
    main()
