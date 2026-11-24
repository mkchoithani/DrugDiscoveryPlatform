"""Runpod Serverless handler: dock a batch of SMILES against a receptor.

Request:
    {
      "input": {
        "receptor_pdb":  "<contents of the prepared receptor .pdb>",
        "reference_pdb": "<contents of the crystal ligand .pdb, defines the box>",
        "compounds":     [{"id": "CHEMBL123", "smiles": "..."}, ...],
        "exhaustiveness": 8,
        "autobox_add": 5.0
      }
    }

Response:
    {
      "scores": {"CHEMBL123": -9.42, ...},
      "failed_embed": [...], "failed_dock": [...],
      "engine": "unidock" | "smina",
      "seconds": 41.2
    }

Design notes that matter in production:
  * Failures are returned per compound, never swallowed. A missing score must be
    distinguishable from a bad score, or the loop will treat "we could not run this"
    as "this molecule is poor".
  * The receptor travels with the request, so the worker holds no state and any
    structure change is visible in the request that caused it.
  * The engine actually used is reported, because a silent CPU fallback would make
    a GPU bill look inexplicable and change throughput by orders of magnitude.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import runpod
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

UNIDOCK = shutil.which("unidock")
SMINA = shutil.which("smina")


def embed_3d(smiles: str, name: str, seed: int = 0xF00D):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=400)
    except Exception:
        pass
    mol.SetProp("_Name", name)
    return mol


def box_from_reference(pdb_text: str, padding: float):
    coords = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise ValueError("reference ligand contained no atoms")
    xs, ys, zs = zip(*coords)
    centre = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    size = (max(xs) - min(xs) + 2 * padding,
            max(ys) - min(ys) + 2 * padding,
            max(zs) - min(zs) + 2 * padding)
    return centre, size


def run_smina(work: Path, receptor: Path, reference: Path, ligands: Path,
              exhaustiveness: int, autobox_add: float) -> dict:
    out = work / "docked.sdf"
    cmd = [SMINA, "-r", str(receptor), "-l", str(ligands),
           "--autobox_ligand", str(reference), "--autobox_add", str(autobox_add),
           "--exhaustiveness", str(exhaustiveness), "--num_modes", "1",
           "--seed", "42", "-o", str(out), "--quiet",
           "--cpu", str(os.cpu_count() or 4)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"smina failed: {proc.stderr[-500:]}")
    return read_scores(out)


def run_unidock(work: Path, receptor: Path, reference: Path, ligands: Path,
                exhaustiveness: int, autobox_add: float) -> dict:
    """Uni-Dock takes one file per ligand, so the SDF is split first."""
    centre, size = box_from_reference(reference.read_text(), autobox_add)
    lig_dir = work / "ligands"
    lig_dir.mkdir(exist_ok=True)
    paths = []
    for mol in Chem.SDMolSupplier(str(ligands), removeHs=False):
        if mol is None:
            continue
        path = lig_dir / f"{mol.GetProp('_Name')}.sdf"
        writer = Chem.SDWriter(str(path))
        writer.write(mol)
        writer.close()
        paths.append(path)

    index = work / "ligands.txt"
    index.write_text("\n".join(str(p) for p in paths))
    out_dir = work / "out"
    out_dir.mkdir(exist_ok=True)
    cmd = [UNIDOCK, "--receptor", str(receptor), "--ligand_index", str(index),
           "--center_x", str(centre[0]), "--center_y", str(centre[1]), "--center_z", str(centre[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exhaustiveness), "--num_modes", "1",
           "--dir", str(out_dir), "--search_mode", "balance"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"unidock failed: {proc.stderr[-500:]}")

    scores: dict[str, float] = {}
    for sdf in out_dir.glob("*.sdf"):
        scores.update(read_scores(sdf))
    return scores


def read_scores(sdf: Path) -> dict:
    scores: dict[str, float] = {}
    if not sdf.exists():
        return scores
    for mol in Chem.SDMolSupplier(str(sdf), removeHs=False):
        if mol is None or not mol.HasProp("_Name"):
            continue
        for prop in ("minimizedAffinity", "Uni-Dock RESULT", "docking_score"):
            if mol.HasProp(prop):
                raw = mol.GetProp(prop)
                try:
                    scores[mol.GetProp("_Name")] = float(str(raw).split()[0])
                except ValueError:
                    continue
                break
    return scores


def handler(event):
    payload = event.get("input", {})
    compounds = payload.get("compounds") or []
    if not compounds:
        return {"error": "no compounds supplied"}
    if not payload.get("receptor_pdb") or not payload.get("reference_pdb"):
        return {"error": "receptor_pdb and reference_pdb are both required"}

    exhaustiveness = int(payload.get("exhaustiveness", 8))
    autobox_add = float(payload.get("autobox_add", 5.0))
    started = time.time()

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        receptor = work / "receptor.pdb"
        reference = work / "reference.pdb"
        receptor.write_text(payload["receptor_pdb"])
        reference.write_text(payload["reference_pdb"])

        ligands = work / "ligands.sdf"
        writer = Chem.SDWriter(str(ligands))
        failed_embed = []
        for c in compounds:
            mol = embed_3d(c["smiles"], c["id"])
            if mol is None:
                failed_embed.append(c["id"])
                continue
            writer.write(mol)
        writer.close()

        engine = "unidock" if UNIDOCK else "smina"
        try:
            if UNIDOCK:
                scores = run_unidock(work, receptor, reference, ligands, exhaustiveness, autobox_add)
            elif SMINA:
                scores = run_smina(work, receptor, reference, ligands, exhaustiveness, autobox_add)
            else:
                return {"error": "no docking engine present in the image"}
        except RuntimeError as exc:
            if UNIDOCK and SMINA:          # GPU engine failed: say so, then fall back
                engine = "smina (unidock failed)"
                scores = run_smina(work, receptor, reference, ligands, exhaustiveness, autobox_add)
            else:
                return {"error": str(exc)}

    failed_dock = [c["id"] for c in compounds
                   if c["id"] not in scores and c["id"] not in failed_embed]
    return {
        "scores": scores,
        "failed_embed": failed_embed,
        "failed_dock": failed_dock,
        "engine": engine,
        "seconds": round(time.time() - started, 1),
        "n_requested": len(compounds),
        "n_scored": len(scores),
    }


runpod.serverless.start({"handler": handler})
