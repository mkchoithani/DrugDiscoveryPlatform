// run.mjs — headless Design-Make-Test loop. `node run.mjs [cycles] [batch]`
import { enumerateLibrary, scoreCompound, batchByIntermediate } from './design.js';
import { makeOracle, trainEnsemble, applyModel, selectBatch, rng, r2, mae } from './ml.js';
import { parseSmiles, descriptors } from './chem.js';

const CYCLES = Number(process.argv[2] ?? 6);
const BATCH = Number(process.argv[3] ?? 12);

// sanity check the parser on known molecules
const checks = [
  ['CCO', 'ethanol', 46.07],
  ['c1ccccc1', 'benzene', 78.11],
  ['CC(=O)Oc1ccccc1C(=O)O', 'aspirin', 180.16],
  ['CN1C=NC2=C1C(=O)N(C)C(=O)N2C', 'caffeine', 194.19]
];
console.log('parser check (MW):');
for (const [smi, name, ref] of checks) {
  const d = descriptors(parseSmiles(smi));
  const err = Math.abs(d.mw - ref);
  console.log(`  ${name.padEnd(9)} calc ${String(d.mw).padStart(6)}  ref ${ref}  ${err < 0.5 ? 'ok' : 'MISMATCH'}   cLogP ${d.clogp}  TPSA ${d.tpsa}  HBD ${d.hbd}`);
}

const lib = enumerateLibrary();
console.log(`\nlibrary: ${lib.length} virtual compounds from ${new Set(lib.map(c => c.intermediateKey)).size} intermediates`);

const oracle = makeOracle(lib, 42, 0.25);
const r = rng(5);
const ceiling = Math.max(...lib.map(c => oracle.true(c.id)));

// cycle 0: a small diverse seed set, no model yet
lib.forEach(c => scoreCompound(c));
let batch = selectBatch(lib, { size: BATCH, beta: 0, diversity: 1.2, rank: () => r() });
let history = [];

for (let cycle = 0; cycle <= CYCLES; cycle++) {
  // MAKE
  const plan = batchByIntermediate(batch);
  // TEST
  const pre = batch.filter(c => c.pred !== null).map(c => [c.pred, oracle.true(c.id)]);
  batch.forEach(c => (c.observed = oracle.assay(c.id, r)));
  const tested = lib.filter(c => c.observed !== null);
  const best = Math.max(...tested.map(c => c.observed));
  history.push({ cycle, tested: tested.length, best, r2: r2(pre), mae: mae(pre), steps: plan.sharedSteps, saved: plan.saved });
  console.log(
    `cycle ${cycle}: made ${batch.length} (${plan.sharedSteps} steps vs ${plan.naiveSteps} naive) ` +
    `| tested ${String(tested.length).padStart(3)} | best pIC50 ${best.toFixed(2)} / ceiling ${ceiling.toFixed(2)}` +
    (pre.length ? ` | prospective R² ${r2(pre)} MAE ${mae(pre)}` : ' | (seed batch, no model)')
  );
  if (cycle === CYCLES) break;
  // DESIGN next cycle
  const model = trainEnsemble(tested, { members: 8, lambda: 1.0, seed: 11 + cycle });
  applyModel(model, lib);
  lib.forEach(c => scoreCompound(c));
  batch = selectBatch(lib, { size: BATCH, beta: 1.0, diversity: 0.6, rank: c => c.score * 10 });
}

const tested = lib.filter(c => c.observed !== null).sort((a, b) => b.observed - a.observed);
console.log(`\ntop hits found (${tested.length} compounds tested of ${lib.length}):`);
for (const c of tested.slice(0, 5)) {
  console.log(`  ${c.id}  pIC50 ${c.observed.toFixed(2)}  MW ${c.desc.mw}  cLogP ${c.desc.clogp}  TPSA ${c.desc.tpsa}  route ${c.route.steps.map(s => s.rxn).join(' → ')} (${(c.route.overallYield * 100).toFixed(0)}% yield)`);
}

// Baseline: random screening of the same number of compounds, averaged over 200 draws.
const rr = rng(99);
let sumBest = 0, sumHits = 0;
const DRAWS = 200;
const goodIds = new Set(lib.filter(c => oracle.true(c.id) >= 8).map(c => c.id));
for (let d = 0; d < DRAWS; d++) {
  const shuffled = [...lib];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(rr() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  const draw = shuffled.slice(0, tested.length);
  sumBest += Math.max(...draw.map(c => oracle.true(c.id)));
  sumHits += draw.filter(c => goodIds.has(c.id)).length;
}
const alHits = tested.filter(c => goodIds.has(c.id)).length;
console.log(`\nbenchmark over ${tested.length} compounds made (library ${lib.length}, ${goodIds.size} true actives pIC50>=8):`);
console.log(`  random screen : best ${(sumBest / DRAWS).toFixed(2)}  actives found ${(sumHits / DRAWS).toFixed(1)}  (mean of ${DRAWS} draws)`);
console.log(`  active learn  : best ${tested[0].observed.toFixed(2)}  actives found ${alHits}`);
console.log(`  enrichment    : ${(alHits / (sumHits / DRAWS)).toFixed(2)}x   (ceiling pIC50 ${ceiling.toFixed(2)})`);
