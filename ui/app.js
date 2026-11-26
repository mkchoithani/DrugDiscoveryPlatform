// UI for the DMTA platform. Talks only to /api — no chemistry in the browser.
const $ = id => document.getElementById(id);
const SV = 'http://www.w3.org/2000/svg';

const state = { campaign: null, info: null, batch: null, chemistry: null, tab: 'batch', busy: false };

async function api(path, opts = {}) {
  const res = await fetch(`/api${path}`, {
    headers: { 'Content-Type': 'application/json' }, ...opts,
  });
  if (!res.ok) throw new Error((await res.json().catch(() => ({}))).detail || res.statusText);
  return res.status === 204 ? null : res.json();
}

function toast(msg, isError = false) {
  const t = document.createElement('div');
  t.className = 'toast' + (isError ? ' err' : '');
  t.textContent = msg;
  document.body.appendChild(t);
  setTimeout(() => t.remove(), isError ? 6000 : 2600);
}

function busy(on, label) {
  state.busy = on;
  for (const b of ['run', 'run3', 'propose']) $(b).disabled = on;
  if (on && label) $('run').textContent = label;
  if (!on) $('run').textContent = 'Run cycle';
}

// ---------------------------------------------------------------- bootstrap

async function boot() {
  state.chemistry = await api('/chemistry');
  const list = await api('/campaigns');
  if (!list.length) {
    await api('/campaigns', { method: 'POST', body: JSON.stringify({ name: 'demo', target: 'simulated kinase' }) });
    return boot();
  }
  $('campaign').innerHTML = list.map(c =>
    `<option value="${c.name}">${c.name} — ${c.measured} tested, ${c.cycles} cycles</option>`).join('');
  state.campaign = state.campaign && list.some(c => c.name === state.campaign) ? state.campaign : list[0].name;
  $('campaign').value = state.campaign;
  await refresh();
}

async function refresh() {
  $('tBig').innerHTML = '<span class="spin">◴</span><small>loading campaign — a real ChEMBL library takes a moment to featurise</small>';
  state.info = await api(`/campaigns/${state.campaign}/state`);
  const s = state.info.settings;
  $('batch').value = s.batch_size;
  $('acq').value = s.acquisition;
  $('beta').value = s.beta; $('betaV').textContent = (+s.beta).toFixed(1);
  $('div').value = s.diversity; $('divV').textContent = (+s.diversity).toFixed(2);
  const opt = $('campaign').querySelector(`option[value="${state.campaign}"]`);
  if (opt) opt.textContent = `${state.campaign} — ${state.info.tested} tested, ${state.info.cycle} cycles`;
  $('worklist').href = `/api/campaigns/${state.campaign}/worklist.csv`;
  $('export').href = `/api/campaigns/${state.campaign}/export.csv`;
  render();
  if (state.tab === 'batch' && !state.batch) await loadBatch();
  else drawTable();
}

async function loadBatch() {
  state.batch = await api(`/campaigns/${state.campaign}/propose?size=${+$('batch').value}`, { method: 'POST' });
  render();
  drawTable();
}

// ------------------------------------------------------------------ actions

$('run').onclick = async () => {
  busy(true, 'Running…');
  try {
    const r = await api(`/campaigns/${state.campaign}/cycle?size=${+$('batch').value}`, { method: 'POST' });
    toast(`Cycle ${r.cycle}: made ${r.made}, best ${r.best_measured}, ${r.hits_in_batch} hits`);
    state.batch = null;
    await refresh(); await loadBatch();
  } catch (e) { toast(e.message, true); } finally { busy(false); }
};

$('run3').onclick = async () => {
  busy(true, 'Running…');
  try {
    for (let i = 0; i < 3; i++)
      await api(`/campaigns/${state.campaign}/cycle?size=${+$('batch').value}`, { method: 'POST' });
    state.batch = null;
    await refresh(); await loadBatch();
    toast('3 cycles complete');
  } catch (e) { toast(e.message, true); } finally { busy(false); }
};

$('propose').onclick = async () => {
  busy(true);
  try { await loadBatch(); state.tab = 'batch'; syncTabs(); toast('Batch proposed — nothing made yet'); }
  catch (e) { toast(e.message, true); } finally { busy(false); }
};

$('campaign').onchange = async e => { state.campaign = e.target.value; state.batch = null; await refresh(); await loadBatch(); };

$('newCampaign').onclick = async () => {
  const sel = $('targetSelect');
  if (!sel.options.length) {
    const targets = await api('/targets');
    sel.innerHTML = '<option value="simulated kinase">Simulated landscape — enumerated library</option>'
      + targets.map(t => `<option value="${t.id}">${t.id} — ${t.name}${t.cached ? '' : ' (will download)'}</option>`).join('');
    sel.onchange = () => {
      $('targetNote').textContent = sel.value.startsWith('CHEMBL')
        ? 'Real measured IC50 data, replayed: a value is revealed only when the loop makes that compound.'
        : 'Enumerated make-on-demand library against a simulated assay — routes and costs are modelled.';
    };
    sel.onchange();
  }
  $('dlg').showModal();
};
$('dlg').addEventListener('close', async e => {
  const dlg = $('dlg');
  if (dlg.returnValue !== 'ok') return;
  const f = new FormData(dlg.querySelector('form'));
  const body = {
    name: f.get('name'), target: f.get('target'),
    batch_size: +f.get('batch_size'), oracle_seed: +f.get('oracle_seed'),
    assay_noise: +f.get('assay_noise'),
    acquisition: $('acq').value, beta: +$('beta').value, diversity: +$('div').value,
  };
  try {
    await api('/campaigns', { method: 'POST', body: JSON.stringify(body) });
    state.campaign = body.name; state.batch = null;
    await boot(); await loadBatch();
    toast(`Campaign ${body.name} created`);
  } catch (err) { toast(err.message, true); }
});

for (const [id, out, dp] of [['beta', 'betaV', 1], ['div', 'divV', 2]])
  $(id).oninput = () => ($(out).textContent = (+$(id).value).toFixed(dp));

$('tabs').onclick = e => {
  const b = e.target.closest('button[data-tab]');
  if (!b) return;
  state.tab = b.dataset.tab;
  syncTabs();
  drawTable();
};
const syncTabs = () => [...$('tabs').children].forEach(x => x.setAttribute('aria-selected', x.dataset.tab === state.tab));

// ------------------------------------------------------------------ render

function render() {
  const info = state.info, last = info.history.at(-1);

  // Design
  const top = state.batch?.compounds?.[0];
  $('dBig').innerHTML = top
    ? `${top.score.toFixed(3)}<small>top score of ${info.library_size} enumerated compounds · ${info.tested} already made</small>`
    : `—<small>${info.library_size} enumerated compounds</small>`;
  const profile = info.profile ?? state.chemistry?.profile ?? [];
  $('dObj').innerHTML = top && profile.length
    ? profile.map(o => {
        const v = top.objectives[o.key] ?? 0;
        const col = v > 0.75 ? 'var(--test)' : v > 0.4 ? 'var(--make)' : 'var(--hot)';
        return `<div class="kv"><span>${o.label}</span><span>
          <span class="bar"><i style="width:${(v * 100).toFixed(0)}%;background:${col}"></i></span>
          <b style="font:11px var(--mono);color:var(--faint)">${v.toFixed(2)}</b></span></div>`;
      }).join('')
    : '<p class="muted">Propose a batch to see the profile.</p>';

  // Make
  const plan = state.batch?.plan ?? last?.plan;
  const div = state.batch?.diversity ?? last?.diversity;
  if (plan) {
    // The compound count comes from the plan, not from step arithmetic: when the batch
    // is measured compounds with no enumerated route, every compound keys its own group
    // and shared_steps - groups.length is exactly zero, which read as "0 compounds".
    const shared = plan.n_shared_intermediates > 0
      ? `${plan.naive_steps} one-at-a-time (−${plan.steps_saved})`
      : 'no intermediate shared by two compounds';
    $('mBig').innerHTML = `${plan.shared_steps}<small>reactions for ${plan.n_compounds} compound${plan.n_compounds === 1 ? '' : 's'} · ${shared}</small>`;
    $('mBody').innerHTML = kvs([
      ['Shared intermediates', plan.n_shared_intermediates],
      ['Largest parallel set', (() => {
        const n = Math.max(...plan.groups.map(g => g.members.length));
        return n > 1 ? `${n} analogues` : 'none — one compound each';
      })()],
      ['Reagent cost', `$${plan.reagent_cost.toFixed(0)}`],
      ['Longest lead time', `${plan.longest_lead_days} d`],
      ['Mean pairwise similarity', div?.mean_pairwise_tanimoto ?? '—'],
      ['Scaffolds represented', div?.n_scaffolds ?? '—'],
    ]);
  }

  // Test
  const oracle = info.oracle;
  $('tBig').innerHTML = info.best_measured
    ? `${info.best_measured}<small>best measured pIC50${oracle ? ` · best in library ${oracle.ceiling}` : ''}</small>`
    : `—<small>nothing measured yet</small>`;
  $('tBody').innerHTML = kvs([
    ['Data source', info.data_source === 'retrospective' ? 'real ChEMBL' : 'simulated'],
    ['Cycle', info.cycle],
    ['Compounds made', `${info.tested} / ${info.library_size}`],
    ['Hits (≥ ' + info.settings.hit_threshold + ')', info.hits],
    ['Model CV R²', info.model?.cv_r2 ?? '—'],
    ['Model CV MAE', info.model?.cv_mae ?? '—'],
    ['σ calibration', info.model?.sigma_calibration ?? '—'],
    ['Last cycle prospective ρ', last?.prospective?.spearman ?? '—'],
  ]);

  drawProgress();
  drawAccuracy();
}

const kvs = rows => rows.map(([k, v]) => `<div class="kv"><span>${k}</span><span>${v}</span></div>`).join('');

function axes(svg, xlab, ylab, xr, yr) {
  svg.innerHTML = '';
  const W = 520, H = 230, L = 48, R = 16, T = 16, B = 34;
  const sx = v => L + ((v - xr[0]) / (xr[1] - xr[0] || 1)) * (W - L - R);
  const sy = v => H - B - ((v - yr[0]) / (yr[1] - yr[0] || 1)) * (H - T - B);
  const el = (n, a, text) => {
    const e = document.createElementNS(SV, n);
    for (const k in a) e.setAttribute(k, a[k]);
    if (text !== undefined) e.textContent = text;
    svg.appendChild(e); return e;
  };
  el('line', { x1: L, y1: H - B, x2: W - R, y2: H - B, class: 'ax' });
  el('line', { x1: L, y1: T, x2: L, y2: H - B, class: 'ax' });
  for (let i = 0; i <= 4; i++) {
    const vx = xr[0] + (i / 4) * (xr[1] - xr[0]), vy = yr[0] + (i / 4) * (yr[1] - yr[0]);
    el('text', { x: sx(vx), y: H - B + 15, class: 'axlbl', 'text-anchor': 'middle' }, fmt(vx));
    el('text', { x: L - 7, y: sy(vy) + 3, class: 'axlbl', 'text-anchor': 'end' }, fmt(vy));
    if (i) el('line', { x1: L, y1: sy(vy), x2: W - R, y2: sy(vy), class: 'ax', opacity: .3 });
  }
  el('text', { x: (W + L) / 2, y: H - 2, class: 'axlbl', 'text-anchor': 'middle' }, xlab);
  el('text', { x: 11, y: (H - B + T) / 2, class: 'axlbl', 'text-anchor': 'middle',
    transform: `rotate(-90 11 ${(H - B + T) / 2})` }, ylab);
  return { sx, sy, el, W, H, L, R, T, B };
}
const fmt = v => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2));

function drawProgress() {
  const h = state.info.history;
  const made = h.map(x => x.tested_total);
  const xr = [0, Math.max(state.info.library_size * 0.25, ...made, 10)];
  const yr = [0, Math.max(...h.map(x => x.best_measured ?? 0), state.info.oracle?.ceiling ?? 9)];
  const a = axes($('progress'), 'compounds made', 'best pIC50 found', xr, [yr[1] - 3.2, yr[1] + 0.2]);
  if (state.info.oracle) {
    a.el('line', { x1: a.sx(xr[0]), y1: a.sy(state.info.oracle.ceiling), x2: a.sx(xr[1]),
      y2: a.sy(state.info.oracle.ceiling), stroke: 'var(--faint)', 'stroke-dasharray': '4 4' });
    a.el('text', { x: a.sx(xr[1]), y: a.sy(state.info.oracle.ceiling) - 6, class: 'axlbl', 'text-anchor': 'end' },
      `best in library ${state.info.oracle.ceiling}`);
  }
  if (!h.length) return a.el('text', { x: 260, y: 110, class: 'axlbl', 'text-anchor': 'middle' }, 'run a cycle');
  const pts = h.map(x => [x.tested_total, x.best_measured]);
  a.el('path', { d: pts.map((p, i) => `${i ? 'L' : 'M'}${a.sx(p[0])},${a.sy(p[1])}`).join(' '),
    fill: 'none', stroke: 'var(--test)', 'stroke-width': 2.5 });
  pts.forEach(p => a.el('circle', { cx: a.sx(p[0]), cy: a.sy(p[1]), r: 3.5, fill: 'var(--test)' }));
}

function drawAccuracy() {
  const h = state.info.history.filter(x => x.prospective && x.prospective.spearman !== null);
  const a = axes($('accuracy'), 'cycle', 'prospective agreement', [0, Math.max(3, state.info.cycle)], [-0.4, 1]);
  a.el('line', { x1: a.sx(0), y1: a.sy(0), x2: a.sx(Math.max(3, state.info.cycle)), y2: a.sy(0),
    stroke: 'var(--faint)', opacity: .5 });
  if (!h.length) return a.el('text', { x: 260, y: 110, class: 'axlbl', 'text-anchor': 'middle' },
    'needs two cycles — the first batch has no model behind it');
  for (const [key, colour, label] of [['spearman', 'var(--design)', 'Spearman ρ (ranking)'], ['r2', 'var(--make)', 'R² (absolute)']]) {
    const pts = h.map(x => [x.cycle, Math.max(-0.4, x.prospective[key] ?? 0)]);
    a.el('path', { d: pts.map((p, i) => `${i ? 'L' : 'M'}${a.sx(p[0])},${a.sy(p[1])}`).join(' '),
      fill: 'none', stroke: colour, 'stroke-width': 2 });
    pts.forEach(p => a.el('circle', { cx: a.sx(p[0]), cy: a.sy(p[1]), r: 3, fill: colour }));
    a.el('text', { x: a.sx(pts.at(-1)[0]) + 6, y: a.sy(pts.at(-1)[1]) + 3, class: 'axlbl', fill: colour }, label.split(' ')[0]);
  }
  a.el('text', { x: a.W - 16, y: 24, class: 'axlbl', 'text-anchor': 'end' }, 'ρ = ranking, R² = absolute');
}

// ------------------------------------------------------------------- tables

async function drawTable() {
  const wrap = $('tableWrap');
  if (state.tab === 'chemistry') return drawChemistry(wrap);
  if (state.tab === 'tested') {
    const rows = await api(`/campaigns/${state.campaign}/compounds?status=tested&limit=200`);
    $('note').textContent = `${state.info.tested} compounds made and tested of ${state.info.library_size}.`;
    return wrap.innerHTML = rows.length ? compoundTable(rows, true) : '<p class="muted">Nothing tested yet.</p>';
  }
  if (!state.batch) return wrap.innerHTML = '<p class="muted">Press “Propose batch”.</p>';
  if (state.tab === 'routes') return drawRoutes(wrap);
  $('note').textContent = `Proposed for cycle ${state.batch.cycle} — not yet made. Download the worklist to send it to the lab.`;
  wrap.innerHTML = compoundTable(state.batch.compounds, false);
}

function compoundTable(rows, tested) {
  const th = state.info.settings.hit_threshold;
  return `<table><thead><tr>
    <th>ID</th><th>SMILES</th><th class="n">Score</th><th class="n">Pred</th><th class="n">±σ</th>
    ${tested ? '<th class="n">Measured</th>' : ''}
    <th class="n">MW</th><th class="n">cLogP</th><th class="n">TPSA</th><th class="n">QED</th>
    <th class="n">SA</th><th>Route</th><th class="n">Yield</th><th class="n">$</th>
  </tr></thead><tbody>${rows.map(c => `<tr${tested && c.measured >= th ? ' class="hit"' : ''}>
    <td>${c.id}</td><td class="smiles" title="${c.smiles}">${c.smiles}</td>
    <td class="n">${(c.score ?? 0).toFixed(3)}</td>
    <td class="n">${c.predicted ?? '—'}</td><td class="n">${c.uncertainty ?? '—'}</td>
    ${tested ? `<td class="n">${c.measured}</td>` : ''}
    <td class="n">${c.desc.mw}</td><td class="n">${c.desc.clogp}</td><td class="n">${c.desc.tpsa}</td>
    <td class="n">${c.desc.qed}</td><td class="n">${c.desc.sa_score}</td>
    <td class="route">${c.route.map(s => s.reaction).join(' → ')}</td>
    <td class="n">${(c.overall_yield * 100).toFixed(0)}%</td><td class="n">${c.cost.toFixed(0)}</td>
  </tr>`).join('')}</tbody></table>`;
}

function drawRoutes(wrap) {
  const { plan, compounds } = state.batch;
  const byId = Object.fromEntries(compounds.map(c => [c.id, c]));
  // A batch of measured compounds carries no enumerated route, so there is no step 2
  // to tabulate. Say that, instead of reading route[1] off an empty array.
  if (!plan.routed) {
    wrap.innerHTML = '<p class="muted">These compounds come from measured data rather than '
      + 'our enumeration, so they carry no synthetic route. A parallel plan appears for '
      + 'batches drawn from the enumerated library.</p>';
    $('note').textContent = `${plan.n_compounds} compounds, each made on its own. `
      + `Reagent spend $${plan.reagent_cost.toFixed(0)}, longest lead ${plan.longest_lead_days} days.`;
    return;
  }
  wrap.innerHTML = plan.groups.map(g => `
    <div class="group">
      <h4>${g.step1} with ${g.reagent}<span class="pill">${g.members.length} analogue${g.members.length > 1 ? 's' : ''} from one intermediate</span></h4>
      <div class="smiles" style="max-width:none">${g.smiles}</div>
      <table><thead><tr><th>Compound</th><th>Step 2</th><th>Reagent</th>
        <th class="n">Yield</th><th class="n">$</th><th class="n">Lead</th></tr></thead><tbody>
      ${g.members.map(id => { const c = byId[id]; return `<tr>
        <td>${c.id}</td><td class="route">${c.route[1].reaction}</td><td>${c.route[1].reagent_name}</td>
        <td class="n">${(c.overall_yield * 100).toFixed(0)}%</td>
        <td class="n">${c.cost.toFixed(0)}</td><td class="n">${c.lead_days}d</td></tr>`; }).join('')}
      </tbody></table>
    </div>`).join('');
  $('note').textContent = `${plan.shared_steps} reactions instead of ${plan.naive_steps}: each intermediate is made once, then split across its analogues. Reagent spend $${plan.reagent_cost.toFixed(0)}, longest lead ${plan.longest_lead_days} days.`;
}

function drawChemistry(wrap) {
  const c = state.chemistry;
  wrap.innerHTML = `
    <div class="group"><h4>Reactions</h4>
      <table><thead><tr><th>Reaction</th><th>SMARTS</th><th class="n">Typical yield</th><th class="n">Reliability</th></tr></thead>
      <tbody>${c.reactions.map(r => `<tr><td>${r.name}</td><td class="route">${r.smarts}</td>
        <td class="n">${(r.typical_yield * 100).toFixed(0)}%</td><td class="n">${(r.reliability * 100).toFixed(0)}%</td></tr>`).join('')}
      </tbody></table></div>
    <div class="group"><h4>Scaffolds</h4>
      <table><thead><tr><th>ID</th><th>Name</th><th>SMILES</th><th>Step 1</th><th>Step 2</th></tr></thead>
      <tbody>${c.scaffolds.map(s => `<tr><td>${s.id}</td><td>${s.name}</td><td class="smiles">${s.smiles}</td>
        <td class="route">${s.step1}</td><td class="route">${s.step2}</td></tr>`).join('')}</tbody></table></div>
    ${Object.entries(c.reagents).map(([role, rs]) => `
    <div class="group"><h4>${role.replace('_', ' ')}s <span class="pill">${rs.length}</span></h4>
      <table><thead><tr><th>ID</th><th>Name</th><th>SMILES</th><th class="n">$/g</th><th class="n">Lead</th></tr></thead>
      <tbody>${rs.map(r => `<tr><td>${r.id}</td><td>${r.name}</td><td class="smiles">${r.smiles}</td>
        <td class="n">${r.price_per_g}</td><td class="n">${r.lead_days}d</td></tr>`).join('')}</tbody></table></div>`).join('')}`;
  $('note').textContent = 'Every product in the library comes from one of these reactions applied to catalogue reagents — which is why each has a route.';
}

boot().catch(e => toast('Startup failed: ' + e.message, true));
