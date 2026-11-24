// app.js — UI driver for the Design-Make-Test loop.
import { enumerateLibrary, scoreCompound, batchByIntermediate, OBJECTIVES } from './design.js';
import { makeOracle, trainEnsemble, applyModel, selectBatch, rng, r2, mae } from './ml.js';

const $ = id => document.getElementById(id);
const SV = 'http://www.w3.org/2000/svg';

let lib, oracle, rand, batch, cycle, history, parity, baseline, ceiling, tab = 'batch';

function cfg() {
  return {
    batch: +$('batch').value, beta: +$('beta').value,
    diversity: +$('div').value, noise: +$('noise').value, seed: +$('seed').value
  };
}

function reset() {
  const c = cfg();
  lib = enumerateLibrary();
  oracle = makeOracle(lib, c.seed, c.noise);
  rand = rng(5);
  ceiling = Math.max(...lib.map(x => oracle.true(x.id)));
  cycle = 0; history = []; parity = [];
  lib.forEach(scoreCompound);
  // seed batch: diverse, no model yet
  const r0 = rng(c.seed + 1);
  batch = selectBatch(lib, { size: c.batch, beta: 0, diversity: 1.2, rank: () => r0() });
  baseline = randomBaseline(lib.length);
  render();
}

// Expected best-of-n under random screening, averaged over draws — the honest comparison.
function randomBaseline(n) {
  const r = rng(999), draws = 120, pts = [];
  const truths = lib.map(x => oracle.true(x.id));
  for (let k = 5; k <= n; k += Math.max(1, Math.round(n / 40))) {
    let s = 0;
    for (let d = 0; d < draws; d++) {
      let best = -Infinity;
      for (let i = 0; i < k; i++) best = Math.max(best, truths[Math.floor(r() * n)]);
      s += best;
    }
    pts.push([k, s / draws]);
  }
  return pts;
}

function step() {
  const c = cfg();
  oracle.noise = c.noise;
  // MAKE
  const plan = batchByIntermediate(batch);
  // TEST — record prospective predictions before revealing results
  for (const m of batch) if (m.pred !== null) parity.push([m.pred, null, m.id]);
  batch.forEach(m => (m.observed = oracle.assay(m.id, rand)));
  parity.forEach(p => { if (p[1] === null) p[1] = lib.find(x => x.id === p[2]).observed; });

  const tested = lib.filter(x => x.observed !== null);
  const recent = parity.slice(-batch.length).map(p => [p[0], p[1]]);
  history.push({
    cycle, made: batch.length, tested: tested.length,
    best: Math.max(...tested.map(x => x.observed)),
    steps: plan.sharedSteps, naive: plan.naiveSteps, saved: plan.saved,
    cost: batch.reduce((s, m) => s + m.route.cost, 0),
    r2: recent.length >= 3 ? r2(recent) : null, mae: recent.length ? mae(recent) : null,
    hits: batch.filter(m => m.observed >= 8).length, plan
  });

  // DESIGN next
  const model = trainEnsemble(tested, { members: 8, lambda: 1.0, seed: 11 + cycle });
  applyModel(model, lib);
  lib.forEach(scoreCompound);
  batch = selectBatch(lib, { size: c.batch, beta: c.beta, diversity: c.diversity, rank: m => m.score * 10 });
  cycle++;
  render();
}

// ---------------- rendering ----------------

function render() {
  const h = history.at(-1);
  const tested = lib.filter(x => x.observed !== null);

  // DESIGN panel
  const top = [...lib].sort((a, b) => b.score - a.score)[0];
  $('dBig').innerHTML = `${top.score.toFixed(3)}<small>top multi-parameter score &middot; ${lib.length} virtual compounds ranked</small>`;
  $('dObj').innerHTML = OBJECTIVES.map(o => {
    const v = top.parts[o.key];
    return `<div class="kv"><span>${o.label}</span><span><span class="bar"><i style="width:${(v * 100).toFixed(0)}%;background:${v > .75 ? 'var(--test)' : v > .4 ? 'var(--make)' : 'var(--hot)'}"></i></span></span></div>`;
  }).join('');

  // MAKE panel
  if (h) {
    const g = h.plan;
    $('mBig').innerHTML = `${h.steps}<small>steps for ${h.made} compounds &middot; ${h.naive} one-at-a-time (&minus;${h.saved})</small>`;
    $('mBody').innerHTML = [
      ['Shared intermediates', `${g.groups.length}`],
      ['Largest parallel set', `${Math.max(...g.groups.map(x => x.members.length))} analogues`],
      ['Building-block spend', `$${h.cost}`],
      ['Mean overall yield', `${(batchMean(h, m => m.route.overallYield) * 100).toFixed(0)}%`],
      ['Longest lead time', `${Math.max(...g.groups.flatMap(x => x.members.map(m => m.route.leadDays)))} d`]
    ].map(kv).join('');
  } else {
    const g = batchByIntermediate(batch);
    $('mBig').innerHTML = `${g.sharedSteps}<small>steps planned for the seed batch</small>`;
    $('mBody').innerHTML = [['Shared intermediates', `${g.groups.length}`], ['Status', 'not yet made']].map(kv).join('');
  }

  // TEST panel
  if (h) {
    $('tBig').innerHTML = `${h.best.toFixed(2)}<small>best pIC50 measured &middot; best in library ${ceiling.toFixed(2)}</small>`;
    $('tBody').innerHTML = [
      ['Cycle', `${h.cycle}`],
      ['Compounds made', `${h.tested} / ${lib.length}`],
      ['Actives this batch', `${h.hits} / ${h.made} (pIC50 ≥ 8)`],
      ['Model R² (prospective)', h.r2 === null ? '—' : h.r2.toFixed(2)],
      ['Model MAE', h.mae === null ? '—' : `${h.mae.toFixed(2)} log units`],
      ['Library explored', `${((h.tested / lib.length) * 100).toFixed(0)}%`]
    ].map(kv).join('');
  } else {
    $('tBody').innerHTML = [['Status', 'no data yet'], ['Next', 'run a cycle']].map(kv).join('');
  }

  drawProgress();
  drawParity();
  drawTable(tested);
  $('run').disabled = $('run5').disabled = lib.every(x => x.observed !== null);
}

const kv = ([k, v]) => `<div class="kv"><span>${k}</span><span>${v}</span></div>`;
const batchMean = (h, f) => h.plan.groups.flatMap(g => g.members).reduce((s, m) => s + f(m), 0) /
  h.plan.groups.flatMap(g => g.members).length;

function axes(svg, xlab, ylab, xr, yr) {
  svg.innerHTML = '';
  const W = 520, H = 240, L = 46, R = 14, T = 14, B = 34;
  const sx = v => L + ((v - xr[0]) / (xr[1] - xr[0])) * (W - L - R);
  const sy = v => H - B - ((v - yr[0]) / (yr[1] - yr[0])) * (H - T - B);
  const el = (n, a) => { const e = document.createElementNS(SV, n); for (const k in a) e.setAttribute(k, a[k]); svg.appendChild(e); return e; };
  el('line', { x1: L, y1: H - B, x2: W - R, y2: H - B, class: 'ax' });
  el('line', { x1: L, y1: T, x2: L, y2: H - B, class: 'ax' });
  for (let i = 0; i <= 4; i++) {
    const vx = xr[0] + (i / 4) * (xr[1] - xr[0]), vy = yr[0] + (i / 4) * (yr[1] - yr[0]);
    el('text', { x: sx(vx), y: H - B + 15, class: 'axlbl', 'text-anchor': 'middle' }).textContent = fmt(vx);
    el('text', { x: L - 7, y: sy(vy) + 3, class: 'axlbl', 'text-anchor': 'end' }).textContent = fmt(vy);
    el('line', { x1: L, y1: sy(vy), x2: W - R, y2: sy(vy), class: 'ax', opacity: i ? .35 : 0 });
  }
  el('text', { x: (W + L) / 2, y: H - 2, class: 'axlbl', 'text-anchor': 'middle' }).textContent = xlab;
  el('text', { x: 10, y: (H - B + T) / 2, class: 'axlbl', 'text-anchor': 'middle', transform: `rotate(-90 10 ${(H - B + T) / 2})` }).textContent = ylab;
  return { sx, sy, el, W, H, L, R, T, B };
}
const fmt = v => (Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(1));

function drawProgress() {
  const yr = [5.5, Math.ceil(ceiling * 2) / 2 + 0.3];
  const xr = [0, lib.length];
  const a = axes($('progress'), 'compounds made', 'best pIC50 found', xr, yr);
  // ceiling
  a.el('line', { x1: a.sx(0), y1: a.sy(ceiling), x2: a.sx(xr[1]), y2: a.sy(ceiling), stroke: 'var(--faint)', 'stroke-dasharray': '4 4' });
  a.el('text', { x: a.sx(xr[1]), y: a.sy(ceiling) - 6, class: 'axlbl', 'text-anchor': 'end' }).textContent = `best in library ${ceiling.toFixed(2)}`;
  // random baseline
  a.el('path', {
    d: baseline.map((p, i) => `${i ? 'L' : 'M'}${a.sx(p[0])},${a.sy(p[1])}`).join(' '),
    fill: 'none', stroke: 'var(--faint)', 'stroke-width': 1.5, opacity: .7
  });
  a.el('text', { x: a.sx(baseline.at(-1)[0]) - 4, y: a.sy(baseline.at(-1)[1]) + 14, class: 'axlbl', 'text-anchor': 'end' }).textContent = 'random screening';
  if (!history.length) return;
  const pts = history.map(h => [h.tested, h.best]);
  a.el('path', {
    d: pts.map((p, i) => `${i ? 'L' : 'M'}${a.sx(p[0])},${a.sy(p[1])}`).join(' '),
    fill: 'none', stroke: 'var(--test)', 'stroke-width': 2.5
  });
  pts.forEach(p => a.el('circle', { cx: a.sx(p[0]), cy: a.sy(p[1]), r: 3.5, fill: 'var(--test)' }));
  a.el('text', { x: a.sx(pts.at(-1)[0]) + 6, y: a.sy(pts.at(-1)[1]) - 6, class: 'axlbl', fill: 'var(--test)' }).textContent = 'this loop';
}

function drawParity() {
  const svg = $('parity');
  if (parity.length < 2) {
    const a = axes(svg, 'predicted pIC50', 'measured pIC50', [5, 10], [5, 10]);
    a.el('text', { x: 260, y: 120, class: 'axlbl', 'text-anchor': 'middle' }).textContent = 'run two cycles to see prospective predictions';
    return;
  }
  const all = parity.flatMap(p => [p[0], p[1]]);
  const lo = Math.floor(Math.min(...all) * 2) / 2 - .2, hi = Math.ceil(Math.max(...all) * 2) / 2 + .2;
  const a = axes(svg, 'predicted pIC50', 'measured pIC50', [lo, hi], [lo, hi]);
  a.el('line', { x1: a.sx(lo), y1: a.sy(lo), x2: a.sx(hi), y2: a.sy(hi), stroke: 'var(--faint)', 'stroke-dasharray': '4 4' });
  parity.forEach((p, i) => {
    const age = i / Math.max(1, parity.length - 1);
    a.el('circle', { cx: a.sx(p[0]), cy: a.sy(p[1]), r: 4, fill: 'var(--design)', opacity: 0.25 + 0.65 * age });
  });
  const pairs = parity.map(p => [p[0], p[1]]);
  a.el('text', { x: a.W - 20, y: 26, class: 'axlbl', 'text-anchor': 'end' })
    .textContent = `all cycles: R² ${r2(pairs)}  MAE ${mae(pairs)}  n=${pairs.length}`;
  a.el('text', { x: a.W - 20, y: 40, class: 'axlbl', 'text-anchor': 'end' })
    .textContent = 'darker = later cycle';
}

function drawTable(tested) {
  const wrap = $('tableWrap');
  if (tab === 'routes') {
    const g = (history.at(-1)?.plan) ?? batchByIntermediate(batch);
    wrap.innerHTML = g.groups.map(grp => `
      <div style="margin:0 0 16px">
        <div style="font-weight:600;margin-bottom:2px">${grp.label}
          <span class="pill">${grp.members.length} analogues off one intermediate</span></div>
        <div class="smiles" style="max-width:none">${grp.smiles}</div>
        <table><thead><tr><th>Compound</th><th>Step 1</th><th>Step 2</th><th class="n">Yield</th><th class="n">BB cost</th><th class="n">Lead</th><th class="n">pIC50</th></tr></thead><tbody>
        ${grp.members.map(m => `<tr${m.observed >= 8 ? ' class="hit"' : ''}>
          <td>${m.id}</td>
          <td class="route">${m.route.steps[0].rxn}</td>
          <td class="route">${m.route.steps[1].rxn} &middot; ${m.r2.name}</td>
          <td class="n">${(m.route.overallYield * 100).toFixed(0)}%</td>
          <td class="n">$${m.route.cost}</td>
          <td class="n">${m.route.leadDays}d</td>
          <td class="n">${m.observed === null ? '—' : m.observed.toFixed(2)}</td></tr>`).join('')}
        </tbody></table>
      </div>`).join('');
    $('note').textContent = `${g.sharedSteps} reactions instead of ${g.naiveSteps}: each intermediate is made once, then split across its analogues.`;
    return;
  }

  let rows, caption;
  if (tab === 'batch') {
    rows = (history.length ? history.at(-1).plan.groups.flatMap(g => g.members) : batch)
      .sort((a, b) => (b.observed ?? b.score) - (a.observed ?? a.score));
    caption = history.length ? `Batch from cycle ${history.at(-1).cycle}, now measured.` : 'Seed batch — diverse, chosen without a model.';
  } else if (tab === 'tested') {
    rows = [...tested].sort((a, b) => b.observed - a.observed);
    caption = `${tested.length} compounds made and tested of ${lib.length} in the library.`;
  } else {
    rows = lib.filter(x => x.observed === null).sort((a, b) => b.score - a.score).slice(0, 60);
    caption = 'Top of the untested ranking — what the loop would make next.';
  }

  wrap.innerHTML = `<table><thead><tr>
    <th>ID</th><th>SMILES</th><th class="n">Score</th><th class="n">Pred</th><th class="n">±σ</th>
    <th class="n">Measured</th><th class="n">MW</th><th class="n">cLogP</th><th class="n">TPSA</th>
    <th class="n">HBD</th><th class="n">RotB</th><th>Route</th><th class="n">Yield</th>
    </tr></thead><tbody>${rows.map(m => `<tr${m.observed >= 8 ? ' class="hit"' : ''}>
      <td>${m.id}</td>
      <td class="smiles" title="${m.smiles}">${m.smiles}</td>
      <td class="n">${m.score.toFixed(3)}</td>
      <td class="n">${m.pred === null ? '—' : m.pred.toFixed(2)}</td>
      <td class="n">${m.unc === null ? '—' : m.unc.toFixed(2)}</td>
      <td class="n">${m.observed === null ? '—' : m.observed.toFixed(2)}</td>
      <td class="n">${m.desc.mw}</td><td class="n">${m.desc.clogp}</td><td class="n">${m.desc.tpsa}</td>
      <td class="n">${m.desc.hbd}</td><td class="n">${m.desc.rotb}</td>
      <td class="route">${m.route.steps.map(s => s.rxn).join(' → ')}</td>
      <td class="n">${(m.route.overallYield * 100).toFixed(0)}%</td>
    </tr>`).join('')}</tbody></table>`;
  $('note').textContent = caption;
}

// ---------------- wiring ----------------
$('run').onclick = step;
$('run5').onclick = () => { for (let i = 0; i < 5 && !lib.every(x => x.observed !== null); i++) step(); };
$('reset').onclick = reset;
for (const [id, out] of [['beta', 'betaV'], ['div', 'divV'], ['noise', 'noiseV']]) {
  $(id).oninput = () => { $(out).textContent = (+$(id).value).toFixed(2).replace(/0$/, ''); };
}
['seed'].forEach(id => ($(id).onchange = reset));
$('tabs').onclick = e => {
  const b = e.target.closest('button[data-tab]');
  if (!b) return;
  tab = b.dataset.tab;
  [...$('tabs').children].forEach(x => x.setAttribute('aria-selected', x === b));
  drawTable(lib.filter(x => x.observed !== null));
};

reset();
