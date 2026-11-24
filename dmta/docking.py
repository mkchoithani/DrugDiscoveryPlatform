"""Molecular docking: score a molecule by physically fitting it into the protein.

Uses smina (a Vina fork) because it reads PDB/SDF directly and can take the binding
box straight from the crystal ligand — fewer conversion steps, fewer silent errors.

Docking is the one scorer here that does NOT depend on having seen similar molecules
before, which is why it is the candidate fix for the time-split failure in RESULTS.md.
It is also the one that fails quietly, so nothing in this module is trusted until
`validate.py` has shown it ranks known compounds correctly.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .library import Compound
from .structure import PreparedStructure

CACHE = Path(__file__).resolve().parent.parent / "data" / "docking_cache.json"
POSES = Path(__file__).resolve().parent.parent / "data" / "poses"
ENGINE = "smina"


class DockingUnavailable(RuntimeError):
    """Raised when no docking engine is installed — never silently skipped."""


def _find_engine() -> str | None:
    """PATH first, then conda envs — smina's arm64 build pins an old Python, so it
    usually lives in an env of its own rather than alongside this code."""
    found = shutil.which(ENGINE)
    if found:
        return found
    import os
    roots = [os.environ.get("CONDA_PREFIX", ""), os.path.expanduser("~/anaconda3"),
             os.path.expanduser("~/miniconda3"), "/opt/anaconda3", "/opt/miniconda3"]
    for root in filter(None, roots):
        for candidate in (Path(root) / "envs" / "dmta-docking" / "bin" / ENGINE,
                          Path(root) / "bin" / ENGINE):
            if candidate.exists():
                return str(candidate)
    return None


def engine_available() -> bool:
    return _find_engine() is not None


def require_engine() -> str:
    path = _find_engine()
    if path is None:
        raise DockingUnavailable(
            f"{ENGINE} not found. Install it in its own environment:\n"
            f"    conda create -y -n dmta-docking -c conda-forge {ENGINE}\n"
            "or run docking in a Linux container (see ROADMAP.md)."
        )
    return path


# ---------------------------------------------------------------- 3D ligands

def embed_3d(smiles: str, name: str, seed: int = 0xF00D) -> Chem.Mol | None:
    """SMILES -> a single reasonable 3D conformer with hydrogens."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True            # second chance for awkward ring systems
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=400)
    except Exception:
        pass                                      # an unoptimised conformer still docks
    mol.SetProp("_Name", name)
    return mol


def write_sdf(compounds: list[Compound], path: Path) -> tuple[Path, list[str]]:
    """Returns the SDF path and the ids that failed to embed — never silently dropped."""
    failed = []
    writer = Chem.SDWriter(str(path))
    for c in compounds:
        mol = embed_3d(c.smiles, c.id)
        if mol is None:
            failed.append(c.id)
            continue
        writer.write(mol)
    writer.close()
    return path, failed


# -------------------------------------------------------------------- docking

@dataclass
class DockingResult:
    scores: dict[str, float]          # compound id -> affinity (kcal/mol, lower = better)
    failed_embed: list[str]
    failed_dock: list[str]
    seconds: float
    engine: str
    settings: dict

    def summary(self) -> str:
        return (f"{len(self.scores)} docked in {self.seconds:.0f}s "
                f"({self.seconds / max(1, len(self.scores)):.1f}s each), "
                f"{len(self.failed_embed)} embed failures, {len(self.failed_dock)} dock failures")


def dock(compounds: list[Compound], structure: PreparedStructure,
         exhaustiveness: int = 8, autobox_add: float = 5.0, cpu: int = 0,
         use_cache: bool = True, verbose: bool = False, chunk_size: int = 25,
         save_poses: bool = True) -> DockingResult:
    """Dock a batch. Scores are kcal/mol — more negative means a better predicted fit.

    Work is done in chunks and the cache is written after each one, so an
    interrupted run keeps everything finished so far. Docking a few hundred
    compounds takes long enough that losing all of it to one crash is a real cost.
    """
    import time
    engine = require_engine()
    cache = _load_cache()
    key_prefix = f"{structure.pdb_id}|{exhaustiveness}|{autobox_add}"

    def cache_key(c: Compound) -> str:
        return f"{key_prefix}|{c.inchikey}"

    todo = [c for c in compounds if not (use_cache and cache_key(c) in cache)]
    scores = {c.id: cache[cache_key(c)] for c in compounds
              if use_cache and cache_key(c) in cache}

    failed_embed: list[str] = []
    failed_dock: list[str] = []
    t0 = time.time()
    pose_dir = POSES / structure.pdb_id
    if save_poses:
        pose_dir.mkdir(parents=True, exist_ok=True)

    for start in range(0, len(todo), chunk_size):
        batch = todo[start:start + chunk_size]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            ligands, embed_failures = write_sdf(batch, tmp / "ligands.sdf")
            failed_embed.extend(embed_failures)
            out = tmp / "docked.sdf"
            cmd = [
                engine,
                "-r", structure.receptor_pdb,
                "-l", str(ligands),
                "--autobox_ligand", structure.reference_ligand_pdb,
                "--autobox_add", str(autobox_add),
                "--exhaustiveness", str(exhaustiveness),
                "--num_modes", "1",
                "--seed", "42",
                "-o", str(out),
                "--quiet",
            ]
            if cpu:
                cmd += ["--cpu", str(cpu)]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"{ENGINE} failed: {proc.stderr[-800:]}")

            by_id = {c.id: c for c in batch}
            docked_ids = set()
            if out.exists():
                for mol in Chem.SDMolSupplier(str(out), removeHs=False):
                    if mol is None or not mol.HasProp("_Name"):
                        continue
                    name = mol.GetProp("_Name")
                    if not mol.HasProp("minimizedAffinity"):
                        continue
                    value = float(mol.GetProp("minimizedAffinity"))
                    scores[name] = value
                    docked_ids.add(name)
                    compound = by_id.get(name)
                    if compound is not None:
                        cache[cache_key(compound)] = value
                        if save_poses:
                            writer = Chem.SDWriter(str(pose_dir / f"{compound.inchikey}.sdf"))
                            writer.write(mol)
                            writer.close()

            failed_dock.extend(c.id for c in batch
                               if c.id not in docked_ids and c.id not in embed_failures)
        if use_cache:
            _save_cache(cache)        # progress survives an interruption
        if verbose:
            done = min(start + chunk_size, len(todo))
            print(f"    docked {done}/{len(todo)}  ({time.time() - t0:.0f}s elapsed)", flush=True)

    result = DockingResult(
        scores=scores, failed_embed=failed_embed, failed_dock=failed_dock,
        seconds=time.time() - t0, engine=ENGINE,
        settings={"pdb_id": structure.pdb_id, "exhaustiveness": exhaustiveness,
                  "autobox_add": autobox_add, "site": structure.site.as_dict()},
    )
    if verbose:
        print("  " + result.summary())
    return result


def pose_path(pdb_id: str, inchikey: str) -> Path | None:
    """Where a docked pose was saved, if we have one."""
    path = POSES / pdb_id.upper() / f"{inchikey}.sdf"
    return path if path.exists() else None


def apply_scores(compounds: list[Compound], result: DockingResult) -> int:
    """Attach docking scores to compounds; returns how many were set."""
    n = 0
    for c in compounds:
        if c.id in result.scores:
            c.docking_score = result.scores[c.id]
            n += 1
        else:
            c.docking_score = None
    return n


# ---------------------------------------------------------------------- cache

def _load_cache() -> dict:
    if CACHE.exists():
        try:
            return json.loads(CACHE.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def _save_cache(cache: dict) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache))


def cache_size() -> int:
    return len(_load_cache())
