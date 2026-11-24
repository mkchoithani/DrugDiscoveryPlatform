"""Protein structure preparation for docking.

Takes a PDB entry and produces: a clean receptor, and a binding box derived from
the ligand that was co-crystallised with it. Using the crystal ligand to define
the site removes the single most error-prone manual step — you are not guessing
where the pocket is, you are reading it off the experiment.

Every choice made here is recorded in the returned manifest, because a silently
different preparation is the classic way docking results stop being reproducible.
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "structures"
USER_AGENT = "dmta-platform/1.0 (research prototype)"

# Things routinely present in a PDB file that are not the protein.
CRYSTALLISATION_JUNK = {
    "HOH", "DOD",                      # water
    "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "DMS", "ACT", "FMT", "TRS",
    "CL", "NA", "K", "MG", "CA", "ZN",  # ions (keep deliberately if catalytic)
}


@dataclass
class BindingSite:
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    reference_ligand: str
    n_ligand_atoms: int

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class PreparedStructure:
    pdb_id: str
    receptor_pdb: str
    reference_ligand_pdb: str
    site: BindingSite
    chain: str
    resolution: float | None
    title: str
    n_protein_atoms: int
    removed: dict
    kept_waters: list[str]

    def as_dict(self) -> dict:
        d = asdict(self)
        d["site"] = self.site.as_dict()
        return d


def fetch_pdb(pdb_id: str, refresh: bool = False) -> Path:
    pdb_id = pdb_id.upper()
    DATA.mkdir(parents=True, exist_ok=True)
    path = DATA / f"{pdb_id}.pdb"
    if path.exists() and not refresh:
        return path
    req = urllib.request.Request(f"https://files.rcsb.org/download/{pdb_id}.pdb",
                                 headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as resp:
        path.write_bytes(resp.read())
    return path


def _atom_lines(text: str):
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            yield line


def _field(line: str):
    return {
        "record": line[:6].strip(),
        "name": line[12:16].strip(),
        "altloc": line[16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21].strip(),
        "resseq": line[22:26].strip(),
        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
        "element": line[76:78].strip(),
    }


def list_ligands(pdb_id: str) -> list[dict]:
    """What non-protein molecules are in this entry, and how big."""
    text = fetch_pdb(pdb_id).read_text()
    counts: dict[tuple, int] = {}
    for line in _atom_lines(text):
        if not line.startswith("HETATM"):
            continue
        f = _field(line)
        key = (f["resname"], f["chain"], f["resseq"])
        counts[key] = counts.get(key, 0) + 1
    out = [{"resname": k[0], "chain": k[1], "resseq": k[2], "atoms": v}
           for k, v in counts.items()]
    return sorted(out, key=lambda d: -d["atoms"])


def pick_reference_ligand(pdb_id: str, min_atoms: int = 12) -> dict:
    """The biggest non-junk heteroatom group — almost always the bound inhibitor."""
    candidates = [l for l in list_ligands(pdb_id)
                  if l["resname"] not in CRYSTALLISATION_JUNK and l["atoms"] >= min_atoms]
    if not candidates:
        raise ValueError(f"{pdb_id}: no co-crystallised ligand big enough to define a site. "
                         "Supply the site centre manually.")
    return candidates[0]


def prepare(pdb_id: str, chain: str | None = None, ligand_resname: str | None = None,
            padding: float = 5.0, keep_waters: list[str] | None = None,
            refresh: bool = False) -> PreparedStructure:
    """Clean the structure and define the box from the crystal ligand."""
    pdb_id = pdb_id.upper()
    path = fetch_pdb(pdb_id, refresh=refresh)
    text = path.read_text()

    title = ""
    resolution = None
    for line in text.splitlines():
        if line.startswith("TITLE"):
            title += " " + line[10:].strip()
        elif line.startswith("REMARK   2 RESOLUTION."):
            try:
                resolution = float(line.split("RESOLUTION.")[1].split("ANGSTROM")[0])
            except (IndexError, ValueError):
                pass
        elif line.startswith("ATOM"):
            break

    ref = (pick_reference_ligand(pdb_id) if ligand_resname is None
           else next(l for l in list_ligands(pdb_id) if l["resname"] == ligand_resname.upper()))
    chain = chain or ref["chain"]
    keep_waters = [w.upper() for w in (keep_waters or [])]

    receptor_lines, ligand_lines = [], []
    removed: dict[str, int] = {}
    for line in _atom_lines(text):
        f = _field(line)
        if f["altloc"] not in ("", "A"):          # keep one conformer only
            removed["alternate conformers"] = removed.get("alternate conformers", 0) + 1
            continue
        is_ref = (f["record"] == "HETATM" and f["resname"] == ref["resname"]
                  and f["chain"] == ref["chain"] and f["resseq"] == ref["resseq"])
        if is_ref:
            ligand_lines.append(line)
            continue
        if f["chain"] != chain:
            removed[f"other chains ({f['chain']})"] = removed.get(f"other chains ({f['chain']})", 0) + 1
            continue
        if f["record"] == "HETATM":
            rn = f["resname"]
            if rn in ("HOH", "DOD") and f["resseq"] in keep_waters:
                receptor_lines.append(line)
                continue
            removed[rn] = removed.get(rn, 0) + 1
            continue
        receptor_lines.append(line)

    if not ligand_lines:
        raise ValueError(f"{pdb_id}: reference ligand {ref['resname']} not found after filtering")

    coords = [( _field(l)["x"], _field(l)["y"], _field(l)["z"]) for l in ligand_lines]
    xs, ys, zs = zip(*coords)
    center = (round(sum(xs) / len(xs), 3), round(sum(ys) / len(ys), 3), round(sum(zs) / len(zs), 3))
    size = (round(max(xs) - min(xs) + 2 * padding, 2),
            round(max(ys) - min(ys) + 2 * padding, 2),
            round(max(zs) - min(zs) + 2 * padding, 2))

    rec_path = DATA / f"{pdb_id}_receptor.pdb"
    lig_path = DATA / f"{pdb_id}_ref_ligand.pdb"
    rec_path.write_text("\n".join(receptor_lines) + "\nEND\n")
    lig_path.write_text("\n".join(ligand_lines) + "\nEND\n")

    prepared = PreparedStructure(
        pdb_id=pdb_id,
        receptor_pdb=str(rec_path),
        reference_ligand_pdb=str(lig_path),
        site=BindingSite(center=center, size=size,
                         reference_ligand=ref["resname"], n_ligand_atoms=len(ligand_lines)),
        chain=chain, resolution=resolution, title=title.strip(),
        n_protein_atoms=len(receptor_lines),
        removed=removed, kept_waters=keep_waters,
    )
    (DATA / f"{pdb_id}_manifest.json").write_text(json.dumps(prepared.as_dict(), indent=2))
    return prepared
