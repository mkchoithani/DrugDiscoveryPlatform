// design.js — DESIGN + MAKE: enumerate a synthesisable library from real reactions,
// score it multi-objectively, and plan routes with shared intermediates.

import { parseSmiles, descriptors, fingerprint, desire } from './chem.js';

// --- Building blocks. Each carries a catalogue cost/lead time so "Make" is real. ---

export const CORES = [
  { id: 'C1', name: 'quinazoline-4-amine', smiles: 'c1ccc2c(c1)ncnc2N{R1}', note: 'kinase hinge binder' },
  { id: 'C2', name: 'pyrazolo[3,4-d]pyrimidine', smiles: 'c1nc2c(c(n1)N{R1})cn[nH]2', note: 'ATP-site mimetic' },
  { id: 'C3', name: 'benzamide', smiles: 'O=C(N{R1})c1ccccc1', note: 'flat, cheap, permeable' }
];

// R1 goes on by SNAr or reductive amination; the product is the shared intermediate.
export const R1 = [
  { id: 'A1', smiles: 'C1CCN(CC1)C(=O)', name: 'piperidine-1-carbonyl', rxn: 'amide coupling', cost: 42, lead: 3 },
  { id: 'A2', smiles: 'c1ccc(cc1)C', name: 'benzyl', rxn: 'reductive amination', cost: 18, lead: 2 },
  { id: 'A3', smiles: 'C1COCC1', name: 'tetrahydrofuran-3-yl', rxn: 'SNAr', cost: 65, lead: 5 },
  { id: 'A4', smiles: 'c1ccncc1', name: 'pyridin-4-yl', rxn: 'SNAr', cost: 30, lead: 3 },
  { id: 'A5', smiles: 'C1CC1', name: 'cyclopropyl', rxn: 'reductive amination', cost: 22, lead: 2 },
  { id: 'A6', smiles: 'CC(C)C', name: 'isobutyl', rxn: 'reductive amination', cost: 12, lead: 1 },
  { id: 'A7', smiles: 'C1CN(CC1)S(=O)(=O)C', name: '1-(methanesulfonyl)pyrrolidin-3-yl', rxn: 'SNAr', cost: 88, lead: 8 },
  { id: 'A8', smiles: 'c1cc(sc1)C', name: '5-methylthiophen-2-yl', rxn: 'Suzuki', cost: 47, lead: 4 },
  { id: 'A9', smiles: 'CCOC', name: '2-methoxyethyl', rxn: 'reductive amination', cost: 10, lead: 1 },
  { id: 'A10', smiles: 'C1CCC(CC1)O', name: '4-hydroxycyclohexyl', rxn: 'reductive amination', cost: 58, lead: 5 }
];

// R2 is the late-stage diversification point (Suzuki / amide), done in parallel.
export const R2 = [
  { id: 'B1', smiles: 'c1ccc(cc1)F', name: '4-fluorophenyl', rxn: 'Suzuki', cost: 28, lead: 2 },
  { id: 'B2', smiles: 'c1ccc(cc1)S(=O)(=O)N', name: '4-sulfamoylphenyl', rxn: 'Suzuki', cost: 95, lead: 7 },
  { id: 'B3', smiles: 'c1ccc(cc1)O', name: '4-hydroxyphenyl', rxn: 'Suzuki', cost: 34, lead: 3 },
  { id: 'B4', smiles: 'c1cc(ccc1)C(=O)N', name: 'benzamide-4-yl', rxn: 'amide coupling', cost: 40, lead: 3 },
  { id: 'B5', smiles: 'c1ccc(cc1)Cl', name: '4-chlorophenyl', rxn: 'Suzuki', cost: 25, lead: 2 },
  { id: 'B6', smiles: 'C1CCNCC1', name: 'piperidin-4-yl', rxn: 'SNAr', cost: 55, lead: 4 },
  { id: 'B7', smiles: 'c1cnc(nc1)N', name: '2-aminopyrimidin-5-yl', rxn: 'Suzuki', cost: 120, lead: 10 },
  { id: 'B8', smiles: 'CC(=O)N', name: 'acetamide', rxn: 'amide coupling', cost: 9, lead: 1 },
  { id: 'B9', smiles: 'c1ccc(cc1)C(F)(F)F', name: '4-(trifluoromethyl)phenyl', rxn: 'Suzuki', cost: 62, lead: 4 },
  { id: 'B10', smiles: 'c1ccc(cc1)OC', name: '4-methoxyphenyl', rxn: 'Suzuki', cost: 26, lead: 2 },
  { id: 'B11', smiles: 'c1cc2c(cc1)nccc2', name: 'quinolin-6-yl', rxn: 'Suzuki', cost: 110, lead: 9 },
  { id: 'B12', smiles: 'C1CN(CCO1)C', name: '4-methylmorpholin-2-yl', rxn: 'SNAr', cost: 44, lead: 3 },
  { id: 'B13', smiles: 'c1cc(cnc1)C(=O)N', name: 'nicotinamide-5-yl', rxn: 'amide coupling', cost: 71, lead: 6 },
  { id: 'B14', smiles: 'c1ccc(cc1)N(C)C', name: '4-(dimethylamino)phenyl', rxn: 'Suzuki', cost: 38, lead: 3 }
];

// Join the pieces into one SMILES. Shared intermediate = core + R1 (step 1),
// then R2 is appended at the diversification vector (step 2).
function assemble(core, r1, r2) {
  const intermediate = core.smiles.replace('{R1}', `(${r1.smiles})`);
  return { intermediate, product: `${intermediate}${r2.smiles}` };
}

export function enumerateLibrary() {
  const out = [];
  for (const core of CORES)
    for (const a of R1)
      for (const b of R2) {
        const { intermediate, product } = assemble(core, a, b);
        let mol, d;
        try { mol = parseSmiles(product); d = descriptors(mol); } catch { continue; }
        out.push({
          id: `${core.id}-${a.id}-${b.id}`,
          smiles: product,
          core, r1: a, r2: b,
          intermediateSmiles: intermediate,
          intermediateKey: `${core.id}-${a.id}`,
          desc: d,
          fp: fingerprint(mol),
          route: route(core, a, b),
          observed: null,          // assay result, once "made and tested"
          pred: null, unc: null    // model output
        });
      }
  return out;
}

function route(core, a, b) {
  const steps = [
    { n: 1, rxn: a.rxn, text: `${core.name} + ${a.name}`, yield: yieldOf(a.rxn), bb: [core.name, a.name] },
    { n: 2, rxn: b.rxn, text: `intermediate + ${b.name}`, yield: yieldOf(b.rxn), bb: [b.name] }
  ];
  const overall = steps.reduce((p, s) => p * s.yield, 1);
  return {
    steps,
    overallYield: Math.round(overall * 100) / 100,
    cost: a.cost + b.cost,
    leadDays: Math.max(a.lead, b.lead) + 2 * steps.length,
    // route confidence: known reactions with decent yields score high
    confidence: Math.round(Math.min(1, overall * 1.35) * 100) / 100
  };
}

const YIELDS = { 'amide coupling': 0.86, 'Suzuki': 0.74, 'SNAr': 0.63, 'reductive amination': 0.78 };
const yieldOf = r => YIELDS[r] ?? 0.6;

// --- Multi-parameter optimisation score (the "Design" objective) ---

export const OBJECTIVES = [
  { key: 'potency', label: 'Potency (pIC50 ≥ 8)', weight: 3, f: c => desire.above(c.pred ?? c.observed ?? 6, 8, 2.5) },
  { key: 'mw', label: 'MW ≤ 500', weight: 1, f: c => desire.below(c.desc.mw, 500, 150) },
  { key: 'clogp', label: 'cLogP 1–4', weight: 1.5, f: c => desire.window(c.desc.clogp, 1, 4, 2) },
  { key: 'tpsa', label: 'TPSA 40–120', weight: 1, f: c => desire.window(c.desc.tpsa, 40, 120, 40) },
  { key: 'rotb', label: 'Rot. bonds ≤ 8', weight: 0.5, f: c => desire.below(c.desc.rotb, 8, 4) },
  { key: 'hbd', label: 'HBD ≤ 3', weight: 0.5, f: c => desire.below(c.desc.hbd, 3, 2) },
  { key: 'route', label: 'Synthesisable', weight: 2, f: c => c.route.confidence },
  { key: 'cost', label: 'BB cost', weight: 0.5, f: c => desire.below(c.route.cost, 80, 120) }
];

export function scoreCompound(c) {
  let num = 0, den = 0;
  const parts = {};
  for (const o of OBJECTIVES) {
    const v = Math.max(0, Math.min(1, o.f(c)));
    parts[o.key] = v;
    num += o.weight * v;
    den += o.weight;
  }
  c.parts = parts;
  c.score = Math.round((num / den) * 1000) / 1000;
  return c.score;
}

// Group a chosen batch by shared intermediate — the parallel-synthesis win.
export function batchByIntermediate(batch) {
  const groups = new Map();
  for (const c of batch) {
    if (!groups.has(c.intermediateKey)) {
      groups.set(c.intermediateKey, {
        key: c.intermediateKey, smiles: c.intermediateSmiles,
        label: `${c.core.name} + ${c.r1.name}`, members: []
      });
    }
    groups.get(c.intermediateKey).members.push(c);
  }
  const list = [...groups.values()];
  const naive = batch.length * 2;                         // 2 steps each, one at a time
  const shared = list.length + batch.length;              // make each intermediate once, then diversify
  return { groups: list, naiveSteps: naive, sharedSteps: shared, saved: naive - shared };
}
