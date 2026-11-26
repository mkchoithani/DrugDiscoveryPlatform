"""PostEra's DEL Simulator as a second ground-truth landscape.

Menzeleev, Chitturi, Davies, Schroeder and Lee, *Digital Discovery* 2026, 5, 4976
(`github.com/postera-ai/del-simulator`). Their simulator exists for the same reason
`dmta/assay.py`'s SimulatedAssay does: you cannot honestly benchmark a model or a
selection strategy on noisy real data, because you never know the right answer for
the compounds you did not make.

Why we load it. The headline number in RESULTS.md 2 is that UCB selection finds
2.75x as many of the top 5% as random screening. That was measured on a landscape
**we wrote**. If the 2.75x is a property of our landscape rather than of the method,
nothing in this repository would reveal it. Theirs is an independently authored
landscape, and its affinities are anchored on 2,660 real soluble-epoxide-hydrolase
binders rather than on a random weighting we chose. Re-running our own acquisition
comparison inside it is the cleanest available audit of our own result.

What this module does NOT do: any of the DEL-specific machinery. Their simulator also
models synthesis yield, multi-round selection, PCR amplification and NGS readout. We
use none of it — only the library and its true affinities, because the thing being
tested is our selection loop, not DEL sequencing.

Generating the inputs (their repository, their Python 3.11 environment):

    export DEL_SIMULATOR_EXPERIMENT_PATH=$(pwd)/data/experiments/example
    python scripts/generate_library.py   data/experiments/example/config.yaml
    python scripts/generate_affinities.py data/experiments/example/config.yaml

That writes `library.csv` (nsynthon_id, smiles) and `affinity_seh_512.csv`. Point
`load_del_landscape` at the pair.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .chem import descriptors
from .chembl import murcko
from .library import Compound


def _read_affinity_column(path: Path) -> list[float]:
    """Their affinity CSV, in order.

    It carries a single column (`pKd`) and **no identifier**: the script builds a
    one-column frame from the affinity vector and writes it with `index=False`, so
    row i of this file belongs to row i of the library CSV it was generated from.
    That makes the join positional, which is also the one way to get a silently
    wrong answer here — hence the length check in `load_del_landscape`.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        col = next((f for f in fields
                    if f.lower() in ("pkd", "p_activity", "pactivity", "affinity", "pki")),
                   fields[0] if fields else None)
        if col is None:
            raise ValueError(f"no columns in {path}")
        values = []
        for row in reader:
            try:
                values.append(float(row[col]))
            except (TypeError, ValueError):
                values.append(float("nan"))
    if not values:
        raise ValueError(f"no affinities in {path}")
    return values


def load_del_landscape(library_csv: str | Path, affinity_csv: str | Path,
                       n: int = 1000, seed: int = 0,
                       sampling: str = "uniform") -> list[Compound]:
    """A size-matched sample of their library, as platform Compounds carrying truth.

    `n` defaults to 1000 deliberately. Their example library has 105,197 distinct
    molecules and ours has 924; the benchmark metric is "how many of the top 5% did
    you find", and the top 5% of 105k compounds is a far harder target than the top
    5% of 924. A size-matched sample keeps the difficulty comparable, which is the
    whole point of the comparison. Pass a larger `n` to make it harder on purpose.

    `sampling` decides the *shape* of the sample, and it matters more than it looks:

      "uniform"  draw at random from the whole library. Over a 105k-member DEL this
                 gives a sample where almost every molecule is its own chemotype —
                 903 distinct Murcko scaffolds in 924 compounds, measured.
      "dense"    take whole groups of products that share their first two reaction
                 cycles, until `n` is reached. Products within a group differ only in
                 the third cycle, so they are close neighbours — the same shape as our
                 own enumerated library, which has 924 compounds over 5 scaffolds.

    Neither is biased towards actives; both draw blind to affinity. The distinction
    exists because active learning needs near-neighbours to generalise from, and a
    uniform draw over a huge combinatorial space has almost none. Comparing the two
    separates "their landscape is different" from "their *sample* is sparse", which
    are very different conclusions about our own results.
    """
    library_csv, affinity_csv = Path(library_csv), Path(affinity_csv)
    affinities = _read_affinity_column(affinity_csv)

    with open(library_csv, newline="") as fh:
        library_rows = list(csv.DictReader(fh))
    if len(library_rows) != len(affinities):
        raise ValueError(
            f"{library_csv.name} has {len(library_rows)} rows but "
            f"{affinity_csv.name} has {len(affinities)}. The affinity file carries no "
            "identifier, so the join is positional and these must match exactly — "
            "they were generated from different libraries.")

    # With yield noise one n-synthon yields several product rows (major product plus
    # side products), so neither nsynthon_id nor smiles is unique. Keep the highest
    # relative_fraction per molecule: the intended product, not a trace impurity.
    best: dict[str, tuple[float, str, float]] = {}
    for row, affinity in zip(library_rows, affinities):
        smiles = row.get("smiles")
        if not smiles or affinity != affinity:          # NaN affinity
            continue
        fraction = float(row.get("relative_fraction") or 1.0)
        key = row.get("nsynthon_id") or smiles
        prev = best.get(smiles)
        if prev is None or fraction > prev[0]:
            best[smiles] = (fraction, key, affinity)
    if not best:
        raise ValueError(f"no usable molecule in {library_csv.name}")

    items = sorted(best.items())                        # deterministic before sampling
    rng = np.random.default_rng(seed)
    if sampling not in ("uniform", "dense"):
        raise ValueError("sampling must be 'uniform' or 'dense'")
    if sampling == "dense":
        groups: dict[str, list] = {}
        for smiles, payload in items:
            key = payload[1]
            groups.setdefault("_".join(key.split("_")[:2]), []).append((smiles, payload))
        order = rng.permutation(sorted(groups))
        items, taken = [], 0
        for key in order:
            if n and taken >= n:
                break
            items.extend(groups[key])
            taken += len(groups[key])
        items = items[:n] if n else items
    elif n and n < len(items):
        items = [items[i] for i in rng.choice(len(items), size=n, replace=False)]

    out: list[Compound] = []
    for smiles, (_fraction, key, affinity) in items:
        try:
            desc = descriptors(smiles)
        except Exception:
            continue                      # a molecule RDKit will not parse is not usable
        # nsynthon_id is "c1_c2_c3": one index per reaction cycle. The first two
        # cycles are shared between every product that differs only in the third,
        # which is exactly our notion of a shared intermediate — so the
        # intermediate-sharing bonus stays meaningful in their library too.
        parts = key.split("_")
        intermediate = "_".join(parts[:2]) if len(parts) >= 3 else key
        comp = Compound(
            id=f"DEL-{key}-{len(out)}", smiles=smiles, inchikey="",
            scaffold_id=murcko(smiles), intermediate_key=intermediate,
            intermediate_smiles="", route=[], desc=desc, alerts=[],
            overall_yield=1.0, route_confidence=None, cost=0.0, lead_days=0,
        )
        comp.truth = float(affinity)                   # revealed only when "made"
        out.append(comp)
    if not out:
        raise ValueError("every sampled molecule failed to parse")
    return out


class NoisyReplayAssay:
    """Replay their affinities with assay noise added, so only the landscape differs.

    RetrospectiveAssay replays a stored value exactly. The 2.75x was measured with
    SimulatedAssay at noise=0.2, and a noiseless oracle makes selection easier, so
    comparing against a noiseless replay would flatter the loop for the wrong reason.
    This adds matching noise. `truth` stays clean — the benchmark scores against the
    real landscape, not against the noisy reading, which is how a real programme is
    judged too.
    """
    name = "del-simulator"

    def __init__(self, library: list[Compound], noise: float = 0.2, seed: int = 42):
        self._truth = {c.id: float(c.truth) for c in library}
        self.noise = noise
        self._rng = np.random.default_rng(seed)

    def truth(self, compound_id: str) -> float:
        return self._truth[compound_id]

    def measure(self, batch: list[Compound]) -> dict[str, float]:
        return {c.id: round(self._truth[c.id] + float(self._rng.normal(0, self.noise)), 3)
                for c in batch}

    def ceiling(self) -> float:
        return max(self._truth.values())

    def actives(self, threshold: float) -> set[str]:
        return {k for k, v in self._truth.items() if v >= threshold}

    def top_fraction(self, frac: float = 0.05) -> set[str]:
        k = max(1, int(len(self._truth) * frac))
        ranked = sorted(self._truth.items(), key=lambda kv: kv[1], reverse=True)
        return {cid for cid, _ in ranked[:k]}

    def threshold_at(self, frac: float = 0.05) -> float:
        k = max(1, int(len(self._truth) * frac))
        return round(sorted(self._truth.values(), reverse=True)[k - 1], 3)


def landscape_summary(compounds: list[Compound], assay=None) -> dict:
    """Enough to say whether the two landscapes are even comparable.

    Truth lives on the compound for a replayed landscape and inside the oracle for a
    simulated one, so either may be the source; pass the assay for the latter.
    """
    if assay is not None:
        v = np.array([assay.truth(c.id) for c in compounds], dtype=float)
    else:
        v = np.array([c.truth for c in compounds], dtype=float)
    return {
        "n": len(v), "min": round(float(v.min()), 2), "max": round(float(v.max()), 2),
        "median": round(float(np.median(v)), 2), "spread": round(float(v.max() - v.min()), 2),
        "sd": round(float(v.std()), 2),
        "top5pct_threshold": round(float(np.percentile(v, 95)), 2),
        "scaffolds": len({c.scaffold_id for c in compounds}),
        "intermediates": len({c.intermediate_key for c in compounds}),
    }
