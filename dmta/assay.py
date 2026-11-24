"""Assay providers — where measurements come from.

Swap SimulatedAssay for CsvAssay (results handed back by a real lab) or write your
own against an ELN/LIMS API. The rest of the platform only needs `measure(batch)`.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Protocol

import numpy as np
from rdkit import Chem

from .chem import fingerprint
from .library import Compound


class AssayProvider(Protocol):
    name: str

    def measure(self, batch: list[Compound]) -> dict[str, float]: ...


class SimulatedAssay:
    """A hidden structure–activity landscape, for benchmarking selection strategies.

    Potency is driven by: a random weighting over fingerprint bits (the unknown
    pharmacophore), a lipophilicity term that helps then hurts, a size penalty, and
    bonuses for two substructures the 'target' happens to like. Plus assay noise.
    """
    name = "simulated"

    MOTIFS = {
        "sulfone": "S(=O)(=O)",
        "aminopyrimidine": "Nc1ncccn1",
        "pyrazole": "c1cn[nH,n]c1",
        "fluorine": "[F]",
    }
    MOTIF_EFFECT = {"sulfone": 0.9, "aminopyrimidine": 1.3, "pyrazole": 0.5, "fluorine": 0.25}

    def __init__(self, library: list[Compound], seed: int = 42, noise: float = 0.2):
        self.noise = noise
        self.seed = seed
        rng = np.random.default_rng(seed)
        n_bits = len(fingerprint(library[0].smiles))
        w = rng.normal(0, 0.35, n_bits)
        hot = rng.choice(n_bits, 15, replace=False)
        w[hot] += rng.normal(0, 2.0, 15)

        patterns = {k: Chem.MolFromSmarts(v) for k, v in self.MOTIFS.items()}
        self._truth: dict[str, float] = {}
        for c in library:
            fp = fingerprint(c.smiles)
            z = float(fp @ w) / max(1.0, float(np.sqrt(fp.sum())))
            lipo = 0.5 * min(c.desc["clogp"], 4.5) - 0.12 * max(0.0, c.desc["clogp"] - 4.5) ** 2
            size = -0.004 * max(0.0, c.desc["mw"] - 430)
            mol = Chem.MolFromSmiles(c.smiles)
            motif = sum(eff for k, eff in self.MOTIF_EFFECT.items()
                        if patterns[k] is not None and mol.HasSubstructMatch(patterns[k]))
            self._truth[c.id] = float(np.clip(4.6 + 1.9 * z + lipo + size + motif, 4.0, 11.0))
        self._rng = np.random.default_rng(seed + 1)

    def truth(self, compound_id: str) -> float:
        return self._truth[compound_id]

    def measure(self, batch: list[Compound]) -> dict[str, float]:
        return {
            c.id: round(self._truth[c.id] + float(self._rng.normal(0, self.noise)), 3)
            for c in batch
        }

    # benchmarking helpers
    def ceiling(self) -> float:
        return max(self._truth.values())

    def actives(self, threshold: float) -> set[str]:
        return {k for k, v in self._truth.items() if v >= threshold}

    def top_fraction(self, frac: float = 0.05) -> set[str]:
        """The objectively best `frac` of the library — the benchmark target set."""
        k = max(1, int(len(self._truth) * frac))
        ranked = sorted(self._truth.items(), key=lambda kv: kv[1], reverse=True)
        return {cid for cid, _ in ranked[:k]}

    def threshold_at(self, frac: float = 0.05) -> float:
        k = max(1, int(len(self._truth) * frac))
        return round(sorted(self._truth.values(), reverse=True)[k - 1], 3)


class RetrospectiveAssay:
    """Replay real measurements that already exist, revealed only when selected.

    This is how active learning is honestly benchmarked on public data: the loop
    may not look at a compound's value until it decides to make it.
    """
    name = "retrospective"

    def __init__(self, library: list[Compound]):
        self._truth = {c.id: float(c.truth) for c in library}

    def truth(self, compound_id: str) -> float:
        return self._truth[compound_id]

    def measure(self, batch: list[Compound]) -> dict[str, float]:
        return {c.id: self._truth[c.id] for c in batch}

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


class CsvAssay:
    """Results returned by a real lab: a CSV of compound_id,value."""
    name = "csv"

    def __init__(self, path: str | Path, id_column: str = "compound_id", value_column: str = "pIC50"):
        self.path, self.id_column, self.value_column = Path(path), id_column, value_column

    def measure(self, batch: list[Compound]) -> dict[str, float]:
        wanted = {c.id for c in batch}
        out: dict[str, float] = {}
        with open(self.path, newline="") as fh:
            for row in csv.DictReader(fh):
                cid = row[self.id_column].strip()
                if cid in wanted and row[self.value_column].strip():
                    out[cid] = float(row[self.value_column])
        missing = wanted - out.keys()
        if missing:
            raise ValueError(f"{len(missing)} compounds have no result in {self.path}: {sorted(missing)[:5]}…")
        return out


def write_batch_csv(batch: list[Compound], path: str | Path) -> Path:
    """The worklist a chemistry team actually receives."""
    path = Path(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["compound_id", "smiles", "intermediate", "step1", "reagent1",
                    "step2", "reagent2", "est_yield", "cost_usd", "lead_days", "pIC50"])
        for c in batch:
            w.writerow([c.id, c.smiles, c.intermediate_smiles,
                        c.route[0].reaction, c.route[0].reagent_name,
                        c.route[1].reaction, c.route[1].reagent_name,
                        c.overall_yield, c.cost, c.lead_days, ""])
    return path
