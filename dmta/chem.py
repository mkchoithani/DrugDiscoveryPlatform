"""RDKit-backed chemistry: descriptors, fingerprints, synthetic accessibility."""
from __future__ import annotations

import os
import sys
from functools import lru_cache

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, QED, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem import RDConfig

RDLogger.DisableLog("rdApp.*")

# RDKit ships Ertl's synthetic-accessibility score in Contrib, not the main API.
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
try:
    import sascorer  # type: ignore
except Exception:  # pragma: no cover
    sascorer = None

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=False)

DESCRIPTOR_KEYS = (
    "mw", "clogp", "tpsa", "hbd", "hba", "rotb", "rings", "aromatic_rings",
    "fsp3", "heavy_atoms", "formal_charge", "qed", "sa_score", "stereocenters",
)


def mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparsable SMILES: {smiles!r}")
    return mol


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(mol_from_smiles(smiles))


def inchikey(smiles: str) -> str:
    return Chem.MolToInchiKey(mol_from_smiles(smiles))


@lru_cache(maxsize=100_000)
def descriptors(smiles: str) -> dict:
    m = mol_from_smiles(smiles)
    return {
        "mw": round(Descriptors.MolWt(m), 2),
        "clogp": round(Crippen.MolLogP(m), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(m), 1),
        "hbd": rdMolDescriptors.CalcNumHBD(m),
        "hba": rdMolDescriptors.CalcNumHBA(m),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(m),
        "rings": rdMolDescriptors.CalcNumRings(m),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "fsp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
        "heavy_atoms": m.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(m),
        "qed": round(QED.qed(m), 3),
        "sa_score": round(sascorer.calculateScore(m), 2) if sascorer else None,
        "stereocenters": len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)),
    }


@lru_cache(maxsize=100_000)
def fingerprint(smiles: str) -> np.ndarray:
    fp = _MORGAN.GetFingerprintAsNumPy(mol_from_smiles(smiles))
    return fp.astype(np.float32)


def feature_vector(smiles: str) -> np.ndarray:
    """Fingerprint + scaled physchem block — what the models actually see."""
    d = descriptors(smiles)
    phys = np.array([
        d["mw"] / 500.0, d["clogp"] / 5.0, d["tpsa"] / 120.0, d["hbd"] / 5.0,
        d["hba"] / 10.0, d["rotb"] / 10.0, d["rings"] / 5.0, d["aromatic_rings"] / 4.0,
        d["fsp3"], d["heavy_atoms"] / 40.0, d["qed"], (d["sa_score"] or 3.0) / 10.0,
    ], dtype=np.float32)
    return np.concatenate([fingerprint(smiles), phys])


def tanimoto(a: np.ndarray, b: np.ndarray) -> float:
    both = float(np.dot(a, b))
    either = float(a.sum() + b.sum() - both)
    return both / either if either else 0.0


def bulk_tanimoto(query: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    both = matrix @ query
    either = matrix.sum(axis=1) + query.sum() - both
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(either > 0, both / either, 0.0)


# --- structural alerts: a real design filter, not a nicety ---
_ALERTS = {
    "nitro_aromatic": "[c][N+](=O)[O-]",
    "michael_acceptor": "[CX3]=[CX3][CX3]=[OX1]",
    "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "thiol": "[SX2H]",
    "hydrazine": "[NX3][NX3]",
    "epoxide": "C1OC1",
    "quaternary_N": "[N+;!$([N+][O-])]",
    "long_alkyl": "CCCCCCC",
}
_ALERT_PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in _ALERTS.items()}


def alerts(smiles: str) -> list[str]:
    m = mol_from_smiles(smiles)
    return [name for name, patt in _ALERT_PATTERNS.items() if patt is not None and m.HasSubstructMatch(patt)]
