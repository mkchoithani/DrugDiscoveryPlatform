// ml.js — TEST: a simulated assay (ground truth), a bagged ridge model with
// uncertainty, and active-learning batch selection.

import { tanimoto } from './chem.js';

// ---------- deterministic RNG so runs are reproducible ----------
export function rng(seed = 7) {
  let s = seed >>> 0;
  return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
}
function gauss(r) {
  const u = Math.max(1e-9, r()), v = r();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// ---------- the hidden truth: a structure–activity landscape ----------
// Unknown to the model. Built from a random weight vector over fingerprint bits
// plus physchem terms, so the model has something learnable but non-trivial.
export function makeOracle(library, seed = 42, noise = 0.25) {
  const r = rng(seed);
  const nbits = library[0].fp.length;
  const w = Float64Array.from({ length: nbits }, () => gauss(r) * 0.42);
  // a few bits matter a lot (a real pharmacophore)
  for (let k = 0; k < 12; k++) w[Math.floor(r() * nbits)] += gauss(r) * 2.2;

  const truth = new Map();
  for (const c of library) {
    let z = 0;
    for (let i = 0; i < nbits; i++) if (c.fp[i] > 0) z += w[i];
    z /= Math.sqrt(c.fp.reduce((s, x) => s + x, 0) || 1);
    // lipophilicity helps potency up to a point (a real and annoying correlation)
    const lipo = 0.55 * Math.min(c.desc.clogp, 5) - 0.09 * Math.max(0, c.desc.clogp - 5) ** 2;
    const size = -0.0025 * Math.max(0, c.desc.mw - 420);
    const pic50 = 6.0 + 1.35 * z + lipo + size;
    truth.set(c.id, Math.max(4, Math.min(11, pic50)));
  }
  return {
    noise,
    true: id => truth.get(id),
    // one assay measurement = truth + experimental error
    assay: (id, r2) => Math.round((truth.get(id) + gauss(r2) * noise) * 100) / 100
  };
}

// ---------- features ----------
function featurize(c) {
  const d = c.desc;
  const extra = [
    1,
    d.mw / 500, d.clogp / 5, d.tpsa / 120, d.hbd / 5, d.hba / 8,
    d.rotb / 10, d.rings / 5, d.fsp3, d.halo / 3, d.aromAtoms / 20
  ];
  const x = new Float64Array(c.fp.length + extra.length);
  x.set(c.fp, 0);
  x.set(extra, c.fp.length);
  return x;
}

// ---------- ridge regression (normal equations, Cholesky) ----------
function ridgeFit(X, y, lambda) {
  const n = X.length, p = X[0].length;
  const A = Array.from({ length: p }, () => new Float64Array(p));
  const b = new Float64Array(p);
  for (let k = 0; k < n; k++) {
    const xk = X[k], yk = y[k];
    for (let i = 0; i < p; i++) {
      if (xk[i] === 0) continue;
      b[i] += xk[i] * yk;
      const Ai = A[i];
      for (let j = i; j < p; j++) if (xk[j] !== 0) Ai[j] += xk[i] * xk[j];
    }
  }
  for (let i = 0; i < p; i++) {
    A[i][i] += lambda;
    for (let j = 0; j < i; j++) A[i][j] = A[j][i];
  }
  return solveSPD(A, b);
}

function solveSPD(A, b) {
  const p = b.length;
  const L = Array.from({ length: p }, () => new Float64Array(p));
  for (let i = 0; i < p; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) L[i][j] = Math.sqrt(Math.max(s, 1e-10));
      else L[i][j] = s / L[j][j];
    }
  }
  const z = new Float64Array(p), x = new Float64Array(p);
  for (let i = 0; i < p; i++) {
    let s = b[i];
    for (let k = 0; k < i; k++) s -= L[i][k] * z[k];
    z[i] = s / L[i][i];
  }
  for (let i = p - 1; i >= 0; i--) {
    let s = z[i];
    for (let k = i + 1; k < p; k++) s -= L[k][i] * x[k];
    x[i] = s / L[i][i];
  }
  return x;
}

// ---------- bagged ensemble → prediction + uncertainty ----------
export function trainEnsemble(trainSet, { members = 8, lambda = 1.0, seed = 11 } = {}) {
  const r = rng(seed);
  const X = trainSet.map(featurize);
  const y = trainSet.map(c => c.observed);
  const models = [];
  for (let m = 0; m < members; m++) {
    const idx = Array.from({ length: X.length }, () => Math.floor(r() * X.length));
    models.push(ridgeFit(idx.map(i => X[i]), idx.map(i => y[i]), lambda));
  }
  const mean = y.reduce((a, b) => a + b, 0) / (y.length || 1);
  return {
    n: trainSet.length,
    predict(c) {
      const x = featurize(c);
      const ps = models.map(w => {
        let s = 0;
        for (let i = 0; i < x.length; i++) if (x[i] !== 0) s += w[i] * x[i];
        return s;
      });
      const mu = ps.reduce((a, b) => a + b, 0) / ps.length;
      const sd = Math.sqrt(ps.reduce((a, b) => a + (b - mu) ** 2, 0) / ps.length);
      return { mu: Number.isFinite(mu) ? mu : mean, sd: Number.isFinite(sd) ? sd : 1 };
    }
  };
}

export function applyModel(model, library) {
  for (const c of library) {
    const { mu, sd } = model.predict(c);
    c.pred = Math.round(mu * 100) / 100;
    c.unc = Math.round(sd * 100) / 100;
  }
}

// ---------- active learning ----------
// Upper-confidence-bound on the design score, penalised for similarity to what
// is already selected — informative *and* diverse, which is the point.
export function selectBatch(library, { size = 12, beta = 1.0, diversity = 0.6, rank }) {
  const pool = library.filter(c => c.observed === null);
  const picked = [];
  const scored = pool.map(c => ({ c, base: rank(c) + beta * (c.unc ?? 1) }));
  scored.sort((a, b) => b.base - a.base);

  while (picked.length < size && scored.length) {
    let bestI = -1, bestV = -Infinity;
    for (let i = 0; i < Math.min(scored.length, 120); i++) {
      const { c, base } = scored[i];
      let sim = 0;
      for (const p of picked) sim = Math.max(sim, tanimoto(c.fp, p.fp));
      const v = base - diversity * sim;
      if (v > bestV) { bestV = v; bestI = i; }
    }
    picked.push(scored.splice(bestI, 1)[0].c);
  }
  return picked;
}

export function metrics(library) {
  const tested = library.filter(c => c.observed !== null);
  const untested = library.filter(c => c.observed === null && c.pred !== null);
  const ys = tested.map(c => c.observed);
  const best = ys.length ? Math.max(...ys) : null;
  // held-out-ish check: model error on the most recent batch before it was used
  return { nTested: tested.length, best, nUntested: untested.length };
}

export function r2(pairs) {
  if (pairs.length < 3) return null;
  const ys = pairs.map(p => p[1]);
  const m = ys.reduce((a, b) => a + b, 0) / ys.length;
  const ssTot = ys.reduce((s, y) => s + (y - m) ** 2, 0);
  const ssRes = pairs.reduce((s, [p, y]) => s + (y - p) ** 2, 0);
  return ssTot ? Math.round((1 - ssRes / ssTot) * 100) / 100 : null;
}

export function mae(pairs) {
  if (!pairs.length) return null;
  return Math.round((pairs.reduce((s, [p, y]) => s + Math.abs(p - y), 0) / pairs.length) * 100) / 100;
}
