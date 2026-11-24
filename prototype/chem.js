// chem.js — minimal cheminformatics: SMILES parser, descriptors, fingerprints.
// Organic subset only (enough for the enumerated library in this prototype).

const MASS = { C: 12.011, N: 14.007, O: 15.999, S: 32.06, F: 18.998, Cl: 35.45, Br: 79.904, I: 126.90, H: 1.008 };
const VAL = { C: 4, N: 3, O: 2, S: 2, F: 1, Cl: 1, Br: 1, I: 1 };

// Crippen-style coarse atomic logP contributions, keyed by element + context.
function atomLogP(a, mol) {
  const nbrs = mol.adj[a.i];
  const el = a.el;
  if (el === 'C') {
    const hetero = nbrs.some(n => mol.atoms[n.j].el !== 'C');
    if (a.arom) return hetero ? 0.10 : 0.16;
    return hetero ? -0.07 : 0.32;
  }
  if (el === 'N') return a.arom ? -0.45 : (a.h > 0 ? -0.95 : -0.55);
  if (el === 'O') return a.h > 0 ? -0.45 : -0.30;
  if (el === 'S') return 0.65;
  if (el === 'F') return 0.40;
  if (el === 'Cl') return 0.75;
  if (el === 'Br') return 0.95;
  if (el === 'I') return 1.10;
  return 0;
}

// Ertl-style TPSA contributions (coarse).
function atomTPSA(a, mol) {
  const deg = mol.adj[a.i].length;
  if (a.el === 'N') {
    if (a.arom) return a.h > 0 ? 15.79 : 12.89;
    if (a.h === 2) return 26.02;
    if (a.h === 1) return 12.03;
    return deg >= 3 ? 3.24 : 12.36;
  }
  if (a.el === 'O') {
    if (a.arom) return 13.14;
    if (a.h > 0) return 20.23;
    return deg >= 2 ? 9.23 : 17.07;
  }
  return 0;
}

export function parseSmiles(smi) {
  const atoms = [], bonds = [];
  const ringOpen = {};
  const stack = [];
  let prev = -1, pendingBond = 0, i = 0;

  const push = (el, arom, hExplicit) => {
    const a = { i: atoms.length, el, arom, h: hExplicit, hExplicit: hExplicit !== null };
    atoms.push(a);
    return a.i;
  };
  const bond = (x, y, order) => {
    if (x < 0 || y < 0 || x === y) return;
    bonds.push({ a: x, b: y, order });
  };

  while (i < smi.length) {
    const c = smi[i];
    if (c === '(') { stack.push(prev); i++; continue; }
    if (c === ')') { prev = stack.pop(); i++; continue; }
    if (c === '-') { pendingBond = 1; i++; continue; }
    if (c === '=') { pendingBond = 2; i++; continue; }
    if (c === '#') { pendingBond = 3; i++; continue; }
    if (c === ':') { pendingBond = 1.5; i++; continue; }
    if (c === '/' || c === '\\' || c === '.') { i++; continue; }
    if (/[0-9]/.test(c)) {
      const n = c;
      if (ringOpen[n] === undefined) { ringOpen[n] = prev; }
      else {
        const other = ringOpen[n];
        const o = pendingBond || (atoms[prev].arom && atoms[other].arom ? 1.5 : 1);
        bond(other, prev, o);
        delete ringOpen[n];
      }
      pendingBond = 0; i++; continue;
    }
    if (c === '[') {
      const end = smi.indexOf(']', i);
      const inner = smi.slice(i + 1, end);
      const m = inner.match(/^([A-Za-z][a-z]?)(?:H(\d?))?/);
      let el = m[1], arom = false;
      if (el[0] === el[0].toLowerCase()) { arom = true; el = el[0].toUpperCase() + el.slice(1); }
      let h = null;
      if (/H/.test(inner)) h = m[2] === '' || m[2] === undefined ? 1 : parseInt(m[2], 10);
      const idx = push(el, arom, h === null ? null : h);
      const o = pendingBond || (prev >= 0 && atoms[prev].arom && arom ? 1.5 : 1);
      bond(prev, idx, o);
      prev = idx; pendingBond = 0; i = end + 1; continue;
    }
    const two = smi.slice(i, i + 2);
    let el = null, arom = false, len = 1;
    if (two === 'Cl' || two === 'Br') { el = two; len = 2; }
    else if (/[A-Z]/.test(c)) { el = c; }
    else if (/[a-z]/.test(c)) { el = c.toUpperCase(); arom = true; }
    else { i++; continue; }
    const idx = push(el, arom, null);
    const o = pendingBond || (prev >= 0 && atoms[prev].arom && arom ? 1.5 : 1);
    bond(prev, idx, o);
    prev = idx; pendingBond = 0; i += len;
  }

  const adj = atoms.map(() => []);
  bonds.forEach((b, k) => {
    adj[b.a].push({ j: b.b, order: b.order, k });
    adj[b.b].push({ j: b.a, order: b.order, k });
  });
  const mol = { atoms, bonds, adj, smiles: smi };

  // implicit hydrogens
  for (const a of atoms) {
    if (a.hExplicit) continue;
    const sum = adj[a.i].reduce((s, n) => s + n.order, 0);
    const target = VAL[a.el] ?? 0;
    a.h = Math.max(0, Math.round(target - sum));
  }

  // ring membership (bonds in cycles) via bridge detection
  markRings(mol);
  return mol;
}

function markRings(mol) {
  const n = mol.atoms.length;
  const disc = new Array(n).fill(-1), low = new Array(n).fill(0);
  let timer = 0;
  mol.bonds.forEach(b => (b.ring = false));
  const dfs = (u, parentEdge) => {
    disc[u] = low[u] = timer++;
    for (const nb of mol.adj[u]) {
      if (nb.k === parentEdge) continue;
      if (disc[nb.j] === -1) {
        dfs(nb.j, nb.k);
        low[u] = Math.min(low[u], low[nb.j]);
        if (low[nb.j] > disc[u]) mol.bonds[nb.k].ring = false;
        else mol.bonds[nb.k].ring = true;
      } else {
        low[u] = Math.min(low[u], disc[nb.j]);
        mol.bonds[nb.k].ring = true;
      }
    }
  };
  for (let i = 0; i < n; i++) if (disc[i] === -1) dfs(i, -1);
  mol.atoms.forEach(a => (a.ring = mol.adj[a.i].some(nb => mol.bonds[nb.k].ring)));
}

export function descriptors(mol) {
  let mw = 0, logp = 0, tpsa = 0, hbd = 0, hba = 0, arom = 0, heavy = 0, halo = 0, sp3 = 0, carbons = 0;
  for (const a of mol.atoms) {
    mw += (MASS[a.el] ?? 0) + a.h * MASS.H;
    logp += atomLogP(a, mol) + a.h * 0.12;
    tpsa += atomTPSA(a, mol);
    if ((a.el === 'N' || a.el === 'O') && a.h > 0) hbd += a.h;
    if (a.el === 'N' || a.el === 'O') hba += 1;
    if (a.arom) arom += 1;
    if (['F', 'Cl', 'Br', 'I'].includes(a.el)) halo += 1;
    if (a.el === 'C') {
      carbons++;
      const maxOrder = Math.max(0, ...mol.adj[a.i].map(nb => nb.order));
      if (!a.arom && maxOrder <= 1) sp3++;
    }
    heavy++;
  }
  let rot = 0;
  for (const b of mol.bonds) {
    if (b.ring || b.order !== 1) continue;
    const da = mol.adj[b.a].length, db = mol.adj[b.b].length;
    const ta = mol.atoms[b.a], tb = mol.atoms[b.b];
    if (da > 1 && db > 1 && !(ta.el === 'C' && ta.h === 3) && !(tb.el === 'C' && tb.h === 3)) rot++;
  }
  const rings = mol.bonds.length - mol.atoms.length + 1;
  return {
    mw: round(mw, 1), clogp: round(logp, 2), tpsa: round(tpsa, 1),
    hbd, hba, rotb: rot, rings: Math.max(0, rings),
    aromAtoms: arom, heavy, halo, fsp3: carbons ? round(sp3 / carbons, 2) : 0
  };
}

const round = (x, d) => Math.round(x * 10 ** d) / 10 ** d;

// Morgan-style circular fingerprint, folded to `bits`.
export function fingerprint(mol, radius = 2, bits = 512) {
  let ids = mol.atoms.map(a => hash(`${a.el}|${a.arom ? 1 : 0}|${mol.adj[a.i].length}|${a.h}|${a.ring ? 1 : 0}`));
  const on = new Set(ids.map(x => x % bits));
  for (let r = 0; r < radius; r++) {
    const next = ids.map((id, i) => {
      const env = mol.adj[i].map(nb => `${nb.order}:${ids[nb.j]}`).sort().join(',');
      return hash(`${r}|${id}|${env}`);
    });
    next.forEach(x => on.add(x % bits));
    ids = next;
  }
  const v = new Float64Array(bits);
  on.forEach(b => (v[b] = 1));
  return v;
}

function hash(s) {
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
  return h >>> 0;
}

export function tanimoto(a, b) {
  let both = 0, either = 0;
  for (let i = 0; i < a.length; i++) {
    const x = a[i] > 0, y = b[i] > 0;
    if (x && y) both++;
    if (x || y) either++;
  }
  return either ? both / either : 0;
}

// Desirability curves — each returns 0..1.
export const desire = {
  window: (x, lo, hi, soft) => {
    if (x >= lo && x <= hi) return 1;
    const d = x < lo ? lo - x : x - hi;
    return Math.max(0, 1 - d / soft);
  },
  below: (x, limit, soft) => (x <= limit ? 1 : Math.max(0, 1 - (x - limit) / soft)),
  above: (x, limit, soft) => (x >= limit ? 1 : Math.max(0, 1 - (limit - x) / soft))
};
