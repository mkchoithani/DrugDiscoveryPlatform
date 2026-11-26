"""Synthesisable library enumeration: real reaction SMARTS over a building-block catalogue.

Every product exists because two known reactions joined catalogue reagents, so each
compound carries its route, yield estimate, reagent cost and lead time. This is the
make-on-demand model (Enamine REAL and friends) rather than retrosynthesis of an
arbitrary target.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field, asdict

from rdkit import Chem
from rdkit.Chem import AllChem

from .chem import alerts, canonical, descriptors, inchikey


# ---------------------------------------------------------------- reactions

@dataclass(frozen=True)
class Reaction:
    name: str
    smarts: str
    typical_yield: float      # median literature/plate yield
    reliability: float        # P(route works first time) — feeds route confidence
    cost_per_mmol: float      # reagent + workup, USD

    def run(self, a: Chem.Mol, b: Chem.Mol) -> Chem.Mol | None:
        rxn = _compiled(self.smarts)
        for reactants in ((a, b), (b, a)):
            products = rxn.RunReactants(reactants)
            if products:
                p = products[0][0]
                try:
                    Chem.SanitizeMol(p)
                except Exception:
                    continue
                return p
        return None


_RXN_CACHE: dict[str, AllChem.ChemicalReaction] = {}


def _compiled(smarts: str) -> AllChem.ChemicalReaction:
    if smarts not in _RXN_CACHE:
        _RXN_CACHE[smarts] = AllChem.ReactionFromSmarts(smarts)
    return _RXN_CACHE[smarts]


REACTIONS = {
    "amide_coupling": Reaction(
        "amide coupling",
        "[C:1](=[O:2])[OX2H1].[N;H2,H1;!$(N-[C,S]=[O,S,N]);!$(N-c):3]>>[C:1](=[O:2])[N:3]",
        0.86, 0.93, 11.0),
    "suzuki": Reaction(
        "Suzuki coupling",
        "[c:1][Br,I].[#6:2][BX3]([OX2H1])[OX2H1]>>[c:1][#6:2]",
        0.74, 0.82, 26.0),
    "snar": Reaction(
        "SNAr",
        "[c:1][Cl,F].[N;H2,H1;!$(N-[C,S]=[O,S,N]):2]>>[c:1][N:2]",
        0.63, 0.71, 14.0),
    "reductive_amination": Reaction(
        "reductive amination",
        "[C:1]=[OX1].[N;H2,H1;!$(N-[C,S]=[O,S,N]);!$(N-c):2]>>[C:1][N:2]",
        0.78, 0.88, 9.0),
    "sulfonamide": Reaction(
        "sulfonamide formation",
        "[S:1](=[O:2])(=[O:3])[Cl].[N;H2,H1;!$(N-[C,S]=[O,S,N]):4]>>[S:1](=[O:2])(=[O:3])[N:4]",
        0.81, 0.89, 15.0),
}


# ------------------------------------------------------- building blocks

@dataclass(frozen=True)
class Reagent:
    id: str
    name: str
    smiles: str
    role: str                 # amine | boronic_acid | sulfonyl_chloride
    price_per_g: float
    lead_days: int
    in_stock: bool = True


@dataclass(frozen=True)
class Scaffold:
    id: str
    name: str
    smiles: str
    step1: str                # reaction key for the first (diversity) point
    step2: str                # reaction key for the late-stage point
    step1_role: str
    step2_role: str


# Bifunctional scaffolds: one handle for each step, orthogonal so order is controllable.
SCAFFOLDS = [
    Scaffold("S1", "4-bromobenzoic acid", "OC(=O)c1ccc(Br)cc1",
             "amide_coupling", "suzuki", "amine", "boronic_acid"),
    Scaffold("S2", "6-bromonicotinic acid", "OC(=O)c1ccc(Br)nc1",
             "amide_coupling", "suzuki", "amine", "boronic_acid"),
    Scaffold("S3", "4-chloro-6-bromoquinazoline", "Clc1ncnc2ccc(Br)cc12",
             "snar", "suzuki", "amine", "boronic_acid"),
    Scaffold("S4", "4-bromobenzaldehyde", "O=Cc1ccc(Br)cc1",
             "reductive_amination", "suzuki", "amine", "boronic_acid"),
    Scaffold("S5", "5-bromothiophene-2-carboxylic acid", "OC(=O)c1ccc(Br)s1",
             "amide_coupling", "suzuki", "amine", "boronic_acid"),
]

AMINES = [
    Reagent("N01", "morpholine", "C1COCCN1", "amine", 0.8, 1),
    Reagent("N02", "piperidine", "C1CCNCC1", "amine", 0.9, 1),
    Reagent("N03", "1-methylpiperazine", "CN1CCNCC1", "amine", 1.6, 1),
    Reagent("N04", "cyclopropylamine", "NC1CC1", "amine", 2.4, 2),
    Reagent("N05", "benzylamine", "NCc1ccccc1", "amine", 0.7, 1),
    Reagent("N06", "2-methoxyethylamine", "NCCOC", "amine", 1.1, 1),
    Reagent("N07", "4-aminotetrahydropyran", "NC1CCOCC1", "amine", 12.0, 4),
    Reagent("N08", "trans-4-aminocyclohexanol", "NC1CCC(O)CC1", "amine", 18.0, 6),
    Reagent("N09", "azetidine", "C1CNC1", "amine", 22.0, 5),
    Reagent("N10", "3-aminopyridine", "Nc1cccnc1", "amine", 3.2, 2),
    Reagent("N11", "2-aminoethanesulfonamide", "NCCS(N)(=O)=O", "amine", 34.0, 9),
    Reagent("N12", "(S)-3-aminotetrahydrofuran", "N[C@H]1CCOC1", "amine", 41.0, 7),
    Reagent("N13", "N,N-dimethylethylenediamine", "CN(C)CCN", "amine", 1.9, 1),
    Reagent("N14", "4-fluorobenzylamine", "NCc1ccc(F)cc1", "amine", 4.5, 2),
]

BORONIC_ACIDS = [
    Reagent("B01", "phenylboronic acid", "OB(O)c1ccccc1", "boronic_acid", 1.2, 1),
    Reagent("B02", "4-fluorophenylboronic acid", "OB(O)c1ccc(F)cc1", "boronic_acid", 2.1, 1),
    Reagent("B03", "4-methoxyphenylboronic acid", "OB(O)c1ccc(OC)cc1", "boronic_acid", 2.4, 1),
    Reagent("B04", "4-(methylsulfonyl)phenylboronic acid", "OB(O)c1ccc(cc1)S(C)(=O)=O", "boronic_acid", 16.0, 4),
    Reagent("B05", "pyridin-3-ylboronic acid", "OB(O)c1cccnc1", "boronic_acid", 5.5, 2),
    Reagent("B06", "pyrimidin-5-ylboronic acid", "OB(O)c1cncnc1", "boronic_acid", 28.0, 7),
    Reagent("B07", "4-carbamoylphenylboronic acid", "OB(O)c1ccc(cc1)C(N)=O", "boronic_acid", 12.0, 3),
    Reagent("B08", "3,5-difluorophenylboronic acid", "OB(O)c1cc(F)cc(F)c1", "boronic_acid", 6.8, 2),
    Reagent("B09", "1-methylpyrazol-4-ylboronic acid", "OB(O)c1cnn(C)c1", "boronic_acid", 9.4, 3),
    Reagent("B10", "4-(trifluoromethyl)phenylboronic acid", "OB(O)c1ccc(cc1)C(F)(F)F", "boronic_acid", 7.2, 2),
    Reagent("B11", "4-hydroxyphenylboronic acid", "OB(O)c1ccc(O)cc1", "boronic_acid", 8.1, 3),
    Reagent("B12", "quinolin-6-ylboronic acid", "OB(O)c1ccc2ncccc2c1", "boronic_acid", 44.0, 10),
    Reagent("B13", "4-(dimethylamino)phenylboronic acid", "OB(O)c1ccc(cc1)N(C)C", "boronic_acid", 5.9, 2),
    Reagent("B14", "cyclopropylboronic acid", "OB(O)C1CC1", "boronic_acid", 13.5, 4),
]

CATALOGUE: dict[str, list[Reagent]] = {"amine": AMINES, "boronic_acid": BORONIC_ACIDS}


# ---------------------------------------------------------------- products

@dataclass
class RouteStep:
    n: int
    reaction: str
    reagent_id: str
    reagent_name: str
    product_smiles: str
    yield_: float
    cost: float
    lead_days: int


@dataclass
class Compound:
    id: str
    smiles: str
    inchikey: str
    scaffold_id: str
    intermediate_key: str
    intermediate_smiles: str
    route: list[RouteStep]
    desc: dict
    alerts: list[str]
    overall_yield: float
    route_confidence: float
    cost: float
    lead_days: int
    # filled in by the platform
    score: float | None = None
    objectives: dict = field(default_factory=dict)
    predicted: float | None = None
    uncertainty: float | None = None
    docking_score: float | None = None   # kcal/mol, lower = better fit
    measured: float | None = None
    cycle_made: int | None = None

    def as_dict(self) -> dict:
        d = asdict(self)
        d["route"] = [asdict(s) for s in self.route]
        return d


def enumerate_library(max_compounds: int | None = None, drop_alerts: bool = True) -> list[Compound]:
    """Two-step enumeration: scaffold -> intermediate (step 1) -> product (step 2)."""
    out: list[Compound] = []
    seen: set[str] = set()

    for scaffold in SCAFFOLDS:
        core = Chem.MolFromSmiles(scaffold.smiles)
        if core is None:
            continue
        rxn1, rxn2 = REACTIONS[scaffold.step1], REACTIONS[scaffold.step2]

        for amine in CATALOGUE[scaffold.step1_role]:
            inter = rxn1.run(core, Chem.MolFromSmiles(amine.smiles))
            if inter is None:
                continue
            inter_smi = Chem.MolToSmiles(inter)

            for boron in CATALOGUE[scaffold.step2_role]:
                prod = rxn2.run(inter, Chem.MolFromSmiles(boron.smiles))
                if prod is None:
                    continue
                smi = Chem.MolToSmiles(prod)
                key = inchikey(smi)
                if key in seen:
                    continue
                flags = alerts(smi)
                if drop_alerts and flags:
                    continue
                seen.add(key)

                steps = [
                    RouteStep(1, rxn1.name, amine.id, amine.name, inter_smi,
                              rxn1.typical_yield, rxn1.cost_per_mmol + amine.price_per_g, amine.lead_days),
                    RouteStep(2, rxn2.name, boron.id, boron.name, smi,
                              rxn2.typical_yield, rxn2.cost_per_mmol + boron.price_per_g, boron.lead_days),
                ]
                overall = steps[0].yield_ * steps[1].yield_
                confidence = rxn1.reliability * rxn2.reliability
                out.append(Compound(
                    id=f"{scaffold.id}-{amine.id}-{boron.id}",
                    smiles=smi, inchikey=key,
                    scaffold_id=scaffold.id,
                    intermediate_key=f"{scaffold.id}-{amine.id}",
                    intermediate_smiles=inter_smi,
                    route=steps,
                    desc=descriptors(smi),
                    alerts=flags,
                    overall_yield=round(overall, 3),
                    route_confidence=round(confidence, 3),
                    cost=round(sum(s.cost for s in steps), 2),
                    lead_days=max(s.lead_days for s in steps) + 2 * len(steps),
                ))
                if max_compounds and len(out) >= max_compounds:
                    return out
    return out


def group_by_intermediate(batch: list[Compound]) -> dict:
    """Parallel-synthesis plan: make each intermediate once, then split it."""
    groups: dict[str, dict] = {}
    for c in batch:
        g = groups.setdefault(c.intermediate_key, {
            "key": c.intermediate_key,
            "smiles": c.intermediate_smiles,
            # compounds from the literature carry no enumerated route
            "step1": c.route[0].reaction if c.route else "no enumerated route",
            "reagent": c.route[0].reagent_name if c.route else "—",
            "members": [],
        })
        g["members"].append(c.id)
    steps_each = max((len(c.route) for c in batch), default=1) or 1
    naive = steps_each * len(batch)
    shared = len(groups) + len(batch) if steps_each > 1 else len(batch)
    # Compounds that came from a measurement rather than from our enumeration carry
    # no route, so each one keys its own group and nothing is actually shared. The
    # caller needs to be able to tell that apart from a real parallel plan, and the
    # compound count has to come from the batch rather than be inferred from the
    # step arithmetic -- len(groups) equals len(batch) in exactly that case, so any
    # difference of the two counts collapses to zero and reads as "0 compounds".
    return {
        "groups": list(groups.values()),
        "n_compounds": len(batch),
        "n_shared_intermediates": sum(1 for g in groups.values() if len(g["members"]) > 1),
        "routed": steps_each > 1,
        "naive_steps": naive,
        "shared_steps": shared,
        "steps_saved": naive - shared,
        "reagent_cost": round(sum(c.cost for c in batch), 2),
        "longest_lead_days": max((c.lead_days for c in batch), default=0),
    }
