"""Real bioactivity data from ChEMBL.

Fetches measured potencies for a target, cleans them the way a computational chemist
would, and presents them as Compounds the platform can run its loop over.

Cleaning is the part that matters. Public bioactivity data is noisy: mixed assay
formats, salts, replicate measurements that disagree, censored values ("> 10 µM").
Everything dropped here is dropped for a stated reason, and the counts are reported.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from statistics import median

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolHash

from .chem import descriptors, inchikey
from .library import Compound

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
USER_AGENT = "dmta-platform/1.0 (research prototype)"

# Well-populated targets, for convenience. Any ChEMBL target id works.
SUGGESTED_TARGETS = {
    "CHEMBL203": "EGFR (epidermal growth factor receptor)",
    "CHEMBL279": "VEGFR2 / KDR",
    "CHEMBL1824": "HER2 / ERBB2",
    "CHEMBL4005": "PI3K alpha",
    "CHEMBL244": "Factor Xa",
    "CHEMBL240": "hERG",
    "CHEMBL230": "COX-2",
    "CHEMBL1862": "ABL1",
}


# ------------------------------------------------------------------- fetching

def fetch_activities(target: str, standard_type: str = "IC50", max_records: int = 25_000,
                     page_size: int = 1000, verbose: bool = True) -> list[dict]:
    """Page through the ChEMBL activity endpoint. Public API, no credentials."""
    rows, offset = [], 0
    while len(rows) < max_records:
        query = urllib.parse.urlencode({
            "target_chembl_id": target,
            "standard_type": standard_type,
            "pchembl_value__isnull": "false",
            "limit": page_size,
            "offset": offset,
        })
        req = urllib.request.Request(f"{API}?{query}", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as resp:
            payload = json.load(resp)
        batch = payload.get("activities", [])
        rows.extend(batch)
        total = payload["page_meta"]["total_count"]
        if verbose:
            print(f"  fetched {len(rows):>6} / {min(total, max_records)}")
        offset += page_size
        if not batch or offset >= total:
            break
        time.sleep(0.2)          # be polite to a public service
    return rows[:max_records]


# ------------------------------------------------------------------ cleaning

_LARGEST = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()


def standardize(smiles: str) -> str | None:
    """Largest fragment, neutralised, canonical. Returns None if unusable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _UNCHARGER.uncharge(_LARGEST.choose(mol))
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    if mol.GetNumHeavyAtoms() < 8:
        return None
    return Chem.MolToSmiles(mol)


def clean(rows: list[dict], max_replicate_spread: float = 1.0,
          mw_range: tuple[float, float] = (150.0, 750.0), verbose: bool = True) -> list[dict]:
    reasons: dict[str, int] = defaultdict(int)
    kept: dict[str, dict] = {}
    measurements: dict[str, list[float]] = defaultdict(list)
    years: dict[str, list[int]] = defaultdict(list)

    for r in rows:
        if r.get("standard_relation") != "=":
            reasons["censored value (>, <, ~)"] += 1; continue
        if r.get("data_validity_comment"):
            reasons["flagged by ChEMBL"] += 1; continue
        if r.get("assay_type") != "B":
            reasons["not a binding assay"] += 1; continue
        if r.get("target_organism") != "Homo sapiens":
            reasons["non-human target"] += 1; continue
        if r.get("pchembl_value") in (None, ""):
            reasons["no pChEMBL value"] += 1; continue
        smiles = standardize(r.get("canonical_smiles") or "")
        if smiles is None:
            reasons["unparsable / too small"] += 1; continue
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        if not (mw_range[0] <= mw <= mw_range[1]):
            reasons[f"MW outside {mw_range[0]:.0f}–{mw_range[1]:.0f}"] += 1; continue

        key = inchikey(smiles)
        measurements[key].append(float(r["pchembl_value"]))
        if r.get("document_year"):
            years[key].append(int(r["document_year"]))
        kept.setdefault(key, {
            "inchikey": key, "smiles": smiles,
            "chembl_id": r.get("molecule_chembl_id"),
            "document": r.get("document_chembl_id") or "",
            "assay": r.get("assay_chembl_id") or "",
        })

    out = []
    for key, rec in kept.items():
        values = measurements[key]
        if len(values) > 1 and (max(values) - min(values)) > max_replicate_spread:
            reasons[f"replicates disagree by >{max_replicate_spread} log"] += 1
            continue
        rec["pchembl"] = round(median(values), 3)
        rec["n_measurements"] = len(values)
        rec["year"] = min(years[key]) if years[key] else None
        out.append(rec)

    if verbose:
        print(f"  {len(rows)} records → {len(out)} unique compounds")
        for reason, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
            print(f"    dropped {n:>6}  {reason}")
    return sorted(out, key=lambda r: (r["year"] or 0, r["chembl_id"] or ""))


# -------------------------------------------------------------------- caching

def cache_path(target: str, standard_type: str = "IC50") -> Path:
    return DATA_DIR / f"{target}-{standard_type}.csv"


def build_dataset(target: str, standard_type: str = "IC50", max_records: int = 25_000,
                  refresh: bool = False, verbose: bool = True) -> Path:
    path = cache_path(target, standard_type)
    if path.exists() and not refresh:
        if verbose:
            print(f"using cached {path}")
        return path
    DATA_DIR.mkdir(exist_ok=True)
    if verbose:
        print(f"fetching {standard_type} data for {target} "
              f"({SUGGESTED_TARGETS.get(target, 'target')}) from ChEMBL")
    rows = clean(fetch_activities(target, standard_type, max_records, verbose=verbose), verbose=verbose)
    if not rows:
        raise RuntimeError(f"no usable data for {target}/{standard_type}")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["chembl_id", "inchikey", "smiles", "pchembl",
                                           "n_measurements", "year", "document", "assay"])
        w.writeheader()
        w.writerows(rows)
    if verbose:
        print(f"wrote {len(rows)} compounds to {path}")
    return path


# ----------------------------------------------------- platform integration

def load_compounds(target: str, standard_type: str = "IC50", limit: int | None = None,
                   refresh: bool = False, verbose: bool = False) -> list[Compound]:
    """Real ChEMBL compounds as platform Compounds.

    These have no enumerated route — they came from the literature, not from our
    building blocks — so route fields are None and the ChEMBL profile skips them,
    leaning on the Ertl SA score for synthesisability instead.
    """
    path = build_dataset(target, standard_type, refresh=refresh, verbose=verbose)
    out: list[Compound] = []
    with open(path, newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            if limit and i >= limit:
                break
            smiles = row["smiles"]
            comp = Compound(
                id=row["chembl_id"], smiles=smiles, inchikey=row["inchikey"],
                scaffold_id=murcko(smiles), intermediate_key=row["chembl_id"],
                intermediate_smiles=smiles, route=[], desc=descriptors(smiles), alerts=[],
                overall_yield=1.0, route_confidence=None, cost=0.0,
                lead_days=0,
            )
            comp.truth = float(row["pchembl"])             # revealed only when "made"
            comp.year = int(row["year"]) if row["year"] else None
            comp.document = row.get("document") or ""
            comp.assay = row.get("assay") or ""
            out.append(comp)
    return out


def murcko(smiles: str) -> str:
    """Bemis–Murcko scaffold — the grouping that makes a split honest.

    Uses rdMolHash rather than rdkit.Chem.Scaffolds.MurckoScaffold: in RDKit
    2026.09.2 the latter raises "bad any cast" on ordinary molecules.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "unknown"
        core = rdMolHash.MolHash(mol, rdMolHash.HashFunction.MurckoScaffold)
        return core or "acyclic"
    except Exception:
        return "unknown"


def dataset_summary(compounds: list[Compound]) -> dict:
    values = [c.truth for c in compounds]
    years = [c.year for c in compounds if c.year]
    return {
        "compounds": len(compounds),
        "scaffolds": len({c.scaffold_id for c in compounds}),
        "pchembl_min": round(min(values), 2),
        "pchembl_median": round(median(values), 2),
        "pchembl_max": round(max(values), 2),
        "actives_ge_8": sum(1 for v in values if v >= 8),
        "year_range": (min(years), max(years)) if years else None,
        "mean_mw": round(sum(c.desc["mw"] for c in compounds) / len(compounds), 1),
    }


# ------------------------------------------------- programmes within a target

def programmes(target: str, standard_type: str = "IC50", min_size: int = 40,
               max_programmes: int | None = None) -> list[dict]:
    """Split a target's data into single-paper campaigns.

    One ChEMBL document is usually one medicinal-chemistry paper: a single team
    optimising a single series. That is far closer to a live programme than the
    whole target, which is thousands of unrelated compounds from hundreds of groups
    and reads as a scatter rather than a campaign.
    """
    from collections import defaultdict
    compounds = load_compounds(target, standard_type)
    if compounds and not any(getattr(c, "document", "") for c in compounds):
        # A cache written before paper ids were kept yields zero programmes and no
        # error — exactly the silent emptiness worth refusing to return.
        raise RuntimeError(
            f"the cached dataset for {target} has no paper ids, so it cannot be split "
            f"into programmes. Re-fetch it:\n"
            f"    python3 -m dmta.cli chembl-fetch --target {target} --refresh")
    by_doc: dict[str, list] = defaultdict(list)
    for c in compounds:
        if getattr(c, "document", ""):
            by_doc[c.document].append(c)

    out = []
    for doc, members in by_doc.items():
        if len(members) < min_size:
            continue
        values = [m.truth for m in members]
        years = [m.year for m in members if m.year]
        out.append({
            "document": doc,
            "target": target,
            "n": len(members),
            "scaffolds": len({m.scaffold_id for m in members}),
            "spread": round(max(values) - min(values), 2),
            "median": round(sorted(values)[len(values) // 2], 2),
            "year": min(years) if years else None,
            "compounds": members,
        })
    out.sort(key=lambda d: -d["n"])
    return out[:max_programmes] if max_programmes else out
