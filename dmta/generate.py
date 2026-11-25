"""Design new molecules from the fragments of measured ones.

Follows the approach of Yang, Lee et al., *Chem. Commun.* 2021, 57, 5909 — the
closest published work to this platform's loop. From 42 compounds with measurable
IC50 and 515 inactives their model designed 5 new compounds predicted to beat the
incumbent, *together with predicted synthesis routes*; all 5 were made and 3 had
measurable activity.

The gap this closes. `dmta/library.py` enumerates a fixed 924-compound space from a
building-block catalogue: every compound is synthesisable by construction, but the
loop can only ever propose what was enumerated up front. This module lets it propose
molecules that were not in any catalogue, by taking apart the compounds a programme
has already measured and putting the pieces back together differently.

## Why BRICS rather than a generative network

BRICS (Degen et al. 2008) cuts molecules at bonds that correspond to known
retrosynthetic disconnections, and labels each cut end with its environment type. A
fragment may only be rejoined to a compatible partner, so recombination stays inside
chemistry that is plausibly makeable. That is the same principle as the paper's —
synthesis is part of the design rather than a filter afterwards — and it needs no
training, no GPU, and no second failure surface.

**What this does NOT give you, stated plainly.** The paper designs compounds *with
predicted routes*. BRICS bond types are retrosynthetically motivated cuts; they are
not a verified route, a reagent list, or a yield. So generated compounds carry
`route=[]` and `route_confidence=None`, exactly as literature compounds pulled from
ChEMBL do, and the scoring profile must be one that does not score route objectives
(`CHEMBL_PROFILE` or `POTENCY_LED_PROFILE`). Turning a generated structure into a
real route needs retrosynthesis — AiZynthFinder or equivalent — and is not done here.
Calling these "synthesisable" would be a claim this code does not support.

## STATUS: incomplete. The known-answer test has not been run.

`denovo_validate.py` is written but has not produced a result, because plain
`BRICSBuild` recombination is the wrong engine for lead-like molecules and the
measurement says so: asking for molecules in the 300-600 Da range that real EGFR
compounds occupy yielded **46 molecules in 210 seconds**, because the builder starts
from a single fragment and grows, so the overwhelming majority of what it emits is
far too small and is thrown away by the size filter.

The fix is not more compute. It is to **seed the build with the programme's own
cores** — `BRICSBuild(fragments, seeds=[...])` starts from a given scaffold and
decorates it, which is both faster and closer to what a medicinal chemist does when
elaborating a series. Until that is done, nothing in this module has been shown to
work, and the docstring above describes an intention rather than a measured result.
"""
from __future__ import annotations

import itertools
import random
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, Crippen, Descriptors

from .chem import alerts, descriptors, inchikey
from .chembl import murcko
from .library import Compound

RDLogger.DisableLog("rdApp.*")        # BRICS products are noisy mid-construction


def fragment_library(compounds: list[Compound], min_occurrences: int = 1,
                     target: str = "truth") -> dict:
    """Cut measured compounds into BRICS fragments.

    `min_occurrences` keeps only fragments seen in at least that many compounds,
    which is a crude way to prefer motifs the programme actually explored over
    one-off decorations.

    Returns the fragments plus enough bookkeeping to say where they came from, since
    "the generator rediscovered a compound" means nothing if the fragments came from
    that compound.
    """
    counts: Counter[str] = Counter()
    for c in compounds:
        mol = Chem.MolFromSmiles(c.smiles)
        if mol is None:
            continue
        for frag in BRICS.BRICSDecompose(mol):
            counts[frag] += 1

    kept = sorted(f for f, n in counts.items() if n >= min_occurrences)
    values = [getattr(c, target, None) for c in compounds]
    values = [v for v in values if v is not None]
    return {
        "fragments": kept,
        "counts": dict(counts),
        "n_source_compounds": len(compounds),
        "n_fragments_total": len(counts),
        "n_fragments_kept": len(kept),
        "source_value_range": ((round(min(values), 2), round(max(values), 2))
                               if values else None),
    }


def _acceptable(mol, smiles: str, mw_range: tuple[float, float],
                max_sa: float | None, drop_alerts: bool
                ) -> tuple[bool, dict | None, list[str]]:
    """Cheap rejections first.

    `descriptors()` computes QED and the Ertl synthetic-accessibility score, which
    together cost far more than everything else in the generate loop. Most candidates
    are rejected on size alone, so weight is checked from the molecule directly before
    any of that runs — this is ~5x the whole generator's throughput, measured.
    """
    weight = Descriptors.MolWt(mol)
    if not (mw_range[0] <= weight <= mw_range[1]):
        return False, None, []
    try:
        desc = descriptors(smiles)
    except Exception:
        return False, None, []
    if max_sa is not None and desc["sa_score"] is not None and desc["sa_score"] > max_sa:
        return False, None, []
    flags = alerts(smiles)
    if drop_alerts and flags:
        return False, None, []
    return True, desc, flags


def generate(fragments: list[str], n: int = 2000, seed: int = 0, max_depth: int = 4,
             mw_range: tuple[float, float] = (200.0, 650.0), max_sa: float | None = 6.0,
             drop_alerts: bool = True, exclude_inchikeys: set[str] | None = None,
             scan_limit: int | None = None) -> list[Compound]:
    """Lazily recombine fragments into novel Compounds.

    BRICS.BRICSBuild is itself a generator, so the space is never enumerated — this
    walks it until `n` molecules pass the filters or the walk is exhausted. That is
    the point: the space is far too large to enumerate, which is precisely why the
    fixed 924-compound library had to be replaced rather than extended.

    `exclude_inchikeys` drops molecules already known, so a "designed" compound is
    never one that was handed in as input.
    """
    mols = [m for m in (Chem.MolFromSmiles(f) for f in fragments) if m is not None]
    if len(mols) < 2:
        return []
    exclude = exclude_inchikeys or set()
    scan_limit = scan_limit if scan_limit is not None else max(n * 60, 20_000)

    # BRICSBuild takes no RNG seed — `scrambleReagents` shuffles via the `random`
    # module's global state, so reproducibility has to be bought by seeding it here.
    # Its `seeds` argument means starting fragments, not randomness.
    random.seed(seed)
    built = BRICS.BRICSBuild(mols, onlyCompleteMols=True, maxDepth=max_depth,
                             scrambleReagents=True)
    out: list[Compound] = []
    seen: set[str] = set()
    for scanned, mol in enumerate(itertools.islice(built, scan_limit)):
        if len(out) >= n:
            break
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            continue
        if not smiles or "*" in smiles:        # an unfinished build, not a molecule
            continue
        key = inchikey(smiles)
        if key in seen or key in exclude:
            continue
        ok, desc, flags = _acceptable(mol, smiles, mw_range, max_sa, drop_alerts)
        if not ok:
            continue
        seen.add(key)
        scaffold = murcko(smiles)
        out.append(Compound(
            id=f"GEN-{len(out):05d}", smiles=smiles, inchikey=key,
            scaffold_id=scaffold,
            # No route, so no shared intermediate either; grouping falls back to the
            # scaffold, which is what the intermediate bonus can honestly use here.
            intermediate_key=scaffold, intermediate_smiles="",
            route=[], desc=desc, alerts=flags, overall_yield=1.0,
            route_confidence=None, cost=0.0, lead_days=0,
        ))
    return out


def design_batch(measured: list[Compound], model, k: int = 5, n_candidates: int = 2000,
                 seed: int = 0, min_occurrences: int = 1, **kwargs) -> dict:
    """The paper's loop in one call: measure, fragment, generate, rank, propose.

    `model` is anything with `.fit(compounds)` and `.apply(compounds)` — the random
    forest by default. Returns the top `k` designs and the bookkeeping needed to
    judge them.
    """
    frags = fragment_library(measured, min_occurrences=min_occurrences)
    known = {c.inchikey for c in measured if c.inchikey}
    known |= {inchikey(c.smiles) for c in measured}
    candidates = generate(frags["fragments"], n=n_candidates, seed=seed,
                          exclude_inchikeys=known, **kwargs)
    if not candidates:
        return {"designs": [], "fragments": frags, "n_candidates": 0}

    for c in measured:
        if c.measured is None:
            c.measured = getattr(c, "truth", None)
    model.fit([c for c in measured if c.measured is not None])
    model.apply(candidates)
    ranked = sorted(candidates, key=lambda c: c.predicted or -99, reverse=True)
    best_measured = max((c.measured for c in measured if c.measured is not None),
                        default=None)
    return {
        "designs": ranked[:k],
        "fragments": frags,
        "n_candidates": len(candidates),
        "incumbent": round(best_measured, 2) if best_measured is not None else None,
        "predicted": [round(c.predicted, 2) for c in ranked[:k]],
        "beats_incumbent": (sum(1 for c in ranked[:k]
                                if best_measured is not None
                                and (c.predicted or -99) > best_measured)),
    }


def recovery(generated: list[Compound], held_out: list[Compound],
             similarity_cutoff: float = 0.7) -> dict:
    """Did the generator rediscover compounds the programme actually made?

    The honest known-answer test. Exact recovery compares InChIKey **connectivity
    blocks** (the first 14 characters), not full keys: BRICS recombination does not
    preserve stereochemistry, so demanding a full-key match would fail a molecule
    that is right in every respect this platform models.

    Near recovery counts a held-out compound as reached if any generated molecule is
    within `similarity_cutoff` Tanimoto of it, which is the weaker claim that the
    generator got into the right region of chemistry.
    """
    import numpy as np
    from .chem import bulk_tanimoto, fingerprint

    if not generated or not held_out:
        return {"error": "need both generated and held-out compounds"}

    def skeleton(c) -> str:
        key = c.inchikey or inchikey(c.smiles)
        return key.split("-")[0]

    gen_keys = {skeleton(c) for c in generated}
    gen_fps = np.stack([fingerprint(c.smiles) for c in generated])

    exact, near, nearest = [], [], []
    for c in held_out:
        is_exact = skeleton(c) in gen_keys
        sims = bulk_tanimoto(fingerprint(c.smiles), gen_fps)
        best = float(sims.max())
        exact.append(is_exact)
        near.append(best >= similarity_cutoff)
        nearest.append(best)

    return {
        "n_generated": len(generated),
        "n_held_out": len(held_out),
        "exact_recovered": int(sum(exact)),
        "exact_rate": round(float(np.mean(exact)), 4),
        "near_recovered": int(sum(near)),
        "near_rate": round(float(np.mean(near)), 4),
        "median_nearest_similarity": round(float(np.median(nearest)), 3),
        "_exact": exact, "_near": near, "_nearest": nearest,
    }
