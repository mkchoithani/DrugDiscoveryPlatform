"""Score compounds against a predicted pharmacophore field.

Boe, Robinson, Lee et al., *Digital Discovery* 2026, 5, 5001
(DOI 10.1039/d6dd00226a), "Learning probabilistic pharmacophore maps enables
efficient hit finding". Their model takes a rigid protein pocket and predicts, over
a 3D grid, the probability that each of six interaction feature types belongs at
each point. Their argument is that docking is expensive and pose-dependent and that
co-folding models carry training-set bias, so a coarser, cheaper representation is
the better instrument. That matters here because docking failed outright on this
platform: 0.0x to 0.37x enrichment across three EGFR structures
([RESULTS.md](../RESULTS.md) 4-5).

**What they release and what they do not.** `postera-ai/phc4-field-generator` gives
featurization, model and inference — the *field*. It does not include the screening
side: taking that field and ranking a library against it. That is this module.

## How a compound is scored

Not by rigid-body optimisation. Searching six pose parameters per conformer against a
grid is docking again, with a different scoring function and the same pose-dependence
the paper is trying to escape. Instead this does classical pharmacophore matching,
which is what the field is shaped for:

1. Reduce the field to discrete typed **hot spots** — local maxima per feature
   channel, above a probability floor, thinned so two spots cannot sit on top of
   each other.
2. Annotate each ligand conformer with RDKit's `BaseFeatures` factory. Its families
   are Donor, Acceptor, NegIonizable, PosIonizable, Aromatic and Hydrophobe — the
   same six channels the model predicts, in the same vocabulary, so no mapping
   between schemes is needed.
3. Find the largest geometrically consistent correspondence between ligand features
   and hot spots: same feature type, and every pairwise distance agreeing within a
   tolerance. This is a maximum-weight clique in a correspondence graph, and it is
   translation- and rotation-invariant, so no pose search happens at all.
4. Score the compound by the best such match over its conformers.

A match of size k means the molecule can place k of its interaction features on k
predicted hot spots simultaneously, at the right distances from each other. The
weight is the sum of those spots' predicted probabilities, so a match onto confident
spots counts for more than one onto marginal ones.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field as dc_field
from pathlib import Path

import numpy as np

# The six channels the model predicts, and the RDKit feature family that means the
# same thing. Identical names, but stated explicitly rather than assumed: a silent
# mismatch here would score every compound against the wrong channel.
FEATURE_CHANNELS: tuple[str, ...] = (
    "Donor", "Acceptor", "NegIonizable", "PosIonizable", "Aromatic", "Hydrophobe",
)


@dataclass
class HotSpot:
    feature: str
    position: np.ndarray          # (3,) in the PDB's coordinate frame
    probability: float


@dataclass
class PharmacophoreField:
    """A predicted field, reduced to typed hot spots."""
    spots: list[HotSpot]
    spacing: float
    origin: np.ndarray
    shape: tuple[int, int, int]
    features: list[str]
    threshold: float
    meta: dict = dc_field(default_factory=dict)

    @property
    def positions(self) -> np.ndarray:
        return (np.stack([s.position for s in self.spots]) if self.spots
                else np.zeros((0, 3)))

    def summary(self) -> dict:
        per_feature = {f: sum(1 for s in self.spots if s.feature == f)
                       for f in self.features}
        return {"spots": len(self.spots), "threshold": self.threshold,
                "spacing": self.spacing, "per_feature": per_feature,
                **self.meta}


def load_field(npz_path: str | Path, threshold: float = 0.9,
               min_separation: float = 3.0, max_spots_per_feature: int = 4
               ) -> PharmacophoreField:
    """Reduce `phc4_field.npz` to a sparse typed pharmacophore query.

    The grid is reduced by non-maximum suppression: take the highest-probability
    point, discard everything within `min_separation` of it, repeat.

    **The defaults are strict for a measured reason.** The field is not sparse. On
    the paper's own sample pocket, `Aromatic` exceeds probability 0.5 at 23,595 of
    55,986 grid points — 42% of the box — because the quantity predicted is "could a
    feature of this type sit here", and inside a pocket the answer is usually yes. A
    0.5 threshold therefore yields a query so dense that any geometry satisfies it:
    benzene, octane and adenosine all scored within 0.07 of each other on the first
    attempt. Keeping only confidently predicted points, well separated, a few per
    channel, is what makes the query actually constrain anything.
    """
    data = np.load(npz_path, allow_pickle=True)
    probs = np.asarray(data["probs"], dtype=np.float32)
    origin = np.asarray(data["origin"], dtype=float).reshape(3)
    shape = tuple(int(v) for v in np.asarray(data["shape"]).reshape(3))
    spacing = float(np.asarray(data["spacing"]).reshape(-1)[0])
    features = [str(f) for f in np.asarray(data["features"]).reshape(-1)]

    if probs.shape[0] != int(np.prod(shape)):
        raise ValueError(f"field has {probs.shape[0]} rows but shape {shape} implies "
                         f"{int(np.prod(shape))}")
    if probs.shape[1] != len(features):
        raise ValueError(f"field has {probs.shape[1]} channels but {len(features)} names")
    unknown = set(features) - set(FEATURE_CHANNELS)
    if unknown:
        raise ValueError(f"unexpected feature channels {sorted(unknown)}; this module "
                         f"knows {FEATURE_CHANNELS}")

    # Row order matches the lattice, so an index unpacks straight into coordinates.
    idx = np.arange(probs.shape[0])
    ijk = np.stack(np.unravel_index(idx, shape), axis=1)
    coords = origin + ijk * spacing

    spots: list[HotSpot] = []
    for channel, name in enumerate(features):
        column = probs[:, channel]
        keep = np.flatnonzero(column >= threshold)
        if keep.size == 0:
            continue
        order = keep[np.argsort(column[keep])[::-1]]
        chosen: list[int] = []
        for i in order:
            if len(chosen) >= max_spots_per_feature:
                break
            if all(np.linalg.norm(coords[i] - coords[j]) >= min_separation
                   for j in chosen):
                chosen.append(int(i))
        spots.extend(HotSpot(name, coords[i], float(column[i])) for i in chosen)

    return PharmacophoreField(
        spots=spots, spacing=spacing, origin=origin, shape=shape, features=features,
        threshold=threshold,
        meta={"grid_points": int(probs.shape[0]),
              "max_probability": round(float(probs.max()), 3)},
    )


# ------------------------------------------------------------------ ligand side

_FACTORY = None


def _feature_factory():
    """RDKit's BaseFeatures definitions — the same vocabulary the model predicts."""
    global _FACTORY
    if _FACTORY is None:
        from rdkit import RDConfig
        from rdkit.Chem import ChemicalFeatures
        _FACTORY = ChemicalFeatures.BuildFeatureFactory(
            str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    return _FACTORY


def conformer_features(smiles: str, n_conformers: int = 5, seed: int = 0,
                       max_features: int = 24
                       ) -> list[tuple[list[str], np.ndarray]]:
    """Per conformer: the feature types present and where they sit.

    A feature's position is the centroid of the atoms that define it, which is what
    makes an aromatic ring one point rather than six.

    Returns [] when RDKit cannot embed the molecule — a compound with no 3D structure
    cannot be scored, and silently returning a zero would look like a weak binder.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.5
    if AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params) == 0:
        return []
    try:
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200)
    except Exception:
        pass                      # an unoptimised conformer is still a conformer

    raw = [f for f in _feature_factory().GetFeaturesForMol(mol)
           if f.GetFamily() in FEATURE_CHANNELS]
    if not raw:
        return []

    out = []
    for conf in mol.GetConformers():
        positions = np.array(conf.GetPositions())
        types, points = [], []
        for f in raw:
            ids = list(f.GetAtomIds())
            types.append(f.GetFamily())
            points.append(positions[ids].mean(axis=0))
        if len(types) > max_features:          # keep the molecule tractable
            types, points = types[:max_features], points[:max_features]
        out.append((types, np.array(points)))
    return out


# ------------------------------------------------------------------- matching

def _max_weight_clique(compatible: np.ndarray, weights: np.ndarray,
                       max_size: int) -> tuple[float, int]:
    """Largest geometrically consistent correspondence, by total spot probability.

    Depth-first over the compatibility graph with an optimistic bound: if the best
    conceivable completion cannot beat the incumbent, stop descending. Node weights
    are the hot spots' predicted probabilities, so size is rewarded but confidence
    decides between matches of equal size.
    """
    n = len(weights)
    if n == 0:
        return 0.0, 0
    order = np.argsort(weights)[::-1]
    best_score, best_size = 0.0, 0
    top = float(weights[order[0]]) if n else 0.0

    def descend(current: list[int], candidates: np.ndarray, score: float):
        nonlocal best_score, best_size
        if score > best_score:
            best_score, best_size = score, len(current)
        if len(current) >= max_size or candidates.size == 0:
            return
        headroom = min(max_size - len(current), candidates.size) * top
        if score + headroom <= best_score:
            return
        for position, v in enumerate(candidates):
            remaining = candidates[position + 1:]
            headroom = min(max_size - len(current) - 1, remaining.size + 1) * top
            if score + headroom <= best_score:
                return
            descend(current + [int(v)],
                    remaining[compatible[v, remaining]],
                    score + float(weights[v]))

    descend([], order, 0.0)
    return best_score, best_size


def match_score(field: PharmacophoreField, types: list[str], points: np.ndarray,
                tolerance: float = 1.0, max_size: int = 8) -> tuple[float, int]:
    """Score one conformer against the field."""
    if not field.spots or len(types) == 0:
        return 0.0, 0
    spot_types = [s.feature for s in field.spots]
    spot_pos = field.positions
    spot_prob = np.array([s.probability for s in field.spots])

    pairs = [(li, si) for li, lt in enumerate(types)
             for si, st in enumerate(spot_types) if lt == st]
    if not pairs:
        return 0.0, 0

    ligand_idx = np.array([p[0] for p in pairs])
    spot_idx = np.array([p[1] for p in pairs])
    weights = spot_prob[spot_idx]

    # Two correspondences can coexist if they use different ligand features and
    # different spots, and the distance between the ligand pair matches the distance
    # between the spot pair. That constraint is what makes this a 3D match rather
    # than a bag-of-features count, and it needs no alignment to evaluate.
    ligand_d = np.linalg.norm(points[ligand_idx][:, None] - points[ligand_idx][None], axis=-1)
    spot_d = np.linalg.norm(spot_pos[spot_idx][:, None] - spot_pos[spot_idx][None], axis=-1)
    compatible = (np.abs(ligand_d - spot_d) <= tolerance)
    compatible &= ligand_idx[:, None] != ligand_idx[None, :]
    compatible &= spot_idx[:, None] != spot_idx[None, :]
    np.fill_diagonal(compatible, False)

    return _max_weight_clique(compatible, weights, max_size)


def score_compound(field: PharmacophoreField, smiles: str, n_conformers: int = 5,
                   tolerance: float = 1.0, max_size: int = 8, seed: int = 0
                   ) -> dict:
    """Best pharmacophore match for a molecule, over its conformers.

    `coverage` is the fraction of the query satisfied, which is the number to compare
    across pockets; `score` weights each satisfied point by the model's confidence in
    it, which is the number to rank within one pocket.
    """
    conformers = conformer_features(smiles, n_conformers=n_conformers, seed=seed)
    if not conformers:
        return {"score": None, "matched": 0, "coverage": None, "conformers": 0}
    best = (0.0, 0)
    for types, points in conformers:
        result = match_score(field, types, points, tolerance=tolerance,
                             max_size=max_size)
        if result > best:
            best = result
    n_spots = max(len(field.spots), 1)
    return {"score": round(best[0], 4), "matched": best[1],
            "coverage": round(best[1] / n_spots, 4),
            "conformers": len(conformers)}


def score_library(field: PharmacophoreField, compounds, n_conformers: int = 5,
                  tolerance: float = 1.0, max_size: int = 8,
                  progress_every: int = 250) -> list[dict]:
    """Score many compounds, keeping count of the ones that could not be embedded."""
    out = []
    for i, c in enumerate(compounds):
        result = score_compound(field, c.smiles, n_conformers=n_conformers,
                                tolerance=tolerance, max_size=max_size)
        result["id"] = c.id
        out.append(result)
        if progress_every and (i + 1) % progress_every == 0:
            failed = sum(1 for r in out if r["score"] is None)
            print(f"  scored {i + 1}/{len(compounds)} ({failed} unembeddable)",
                  flush=True)
    return out
