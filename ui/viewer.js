// 3D structure and pose viewer. Renders the real coordinates the docking used.
const $ = id => document.getElementById(id);

const state = {
  pdb: '2RGP',
  viewer: null,
  structure: null,
  models: {},           // name -> 3Dmol model or shape
  shown: { cartoon: true, pocket: true, crystal: true, pose: true, box: false, spin: false },
  pocketResi: [],
  selected: null,
};

const COLOURS = {
  protein: '#9AA69F',   // muted grey-green ribbons
  crystal: '#E0A32E',   // the ligand from the crystal structure
  pose: '#58D68D',      // the docked pose
  surface: '#5FA8A0',   // pocket surface
};

// --- geometry helpers: 3Dmol's `within` selector did not match across models,
// so the pocket residues are computed here instead. Deterministic and checkable.
function parsePdbAtoms(text) {
  const atoms = [];
  for (const line of text.split('\n')) {
    if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
    atoms.push({
      resi: parseInt(line.slice(22, 26), 10),
      x: parseFloat(line.slice(30, 38)),
      y: parseFloat(line.slice(38, 46)),
      z: parseFloat(line.slice(46, 54)),
    });
  }
  return atoms;
}

function pocketResidues(receptorText, ligandText, cutoff = 8) {
  const prot = parsePdbAtoms(receptorText);
  const lig = parsePdbAtoms(ligandText);
  const c2 = cutoff * cutoff;
  const near = new Set();
  for (const a of prot) {
    if (near.has(a.resi)) continue;
    for (const b of lig) {
      const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
      if (dx * dx + dy * dy + dz * dz <= c2) { near.add(a.resi); break; }
    }
  }
  return [...near];
}

async function api(path) {
  const res = await fetch(`/api${path}`);
  if (!res.ok) throw new Error((await res.json().catch(() => ({}))).detail || res.statusText);
  return res.json();
}

// ---------------------------------------------------------------- bootstrap

async function boot() {
  state.viewer = $3Dmol.createViewer($('gl'), { backgroundColor: '#0E1512', antialias: true });
  await loadStructure(state.pdb);
  await loadDockedList();
}

async function loadStructure(pdbId) {
  $('glLoading').style.display = 'grid';
  const s = await api(`/structures/${pdbId}`);
  state.structure = s;
  const v = state.viewer;
  v.clear();
  state.models = {};

  // the protein
  state.models.protein = v.addModel(s.receptor_data, 'pdb');
  // the ligand that was in the crystal — this is what defined the pocket
  state.models.crystal = v.addModel(s.reference_ligand_data, 'pdb');
  state.pocketResi = pocketResidues(s.receptor_data, s.reference_ligand_data, 8);

  renderStyles();
  v.zoomTo({ model: state.models.crystal });
  v.zoom(0.45);
  v.render();
  $('glLoading').style.display = 'none';

  $('meta').innerHTML = [
    ['PDB entry', s.pdb_id],
    ['Resolution', `${s.resolution} Å`],
    ['Chain', s.chain],
    ['Protein atoms', s.n_protein_atoms],
    ['Site from', `${s.site.reference_ligand} (${s.site.n_ligand_atoms} atoms)`],
    ['Box centre', s.site.center.map(x => x.toFixed(1)).join(', ')],
    ['Box size', s.site.size.map(x => x.toFixed(1)).join(' × ') + ' Å'],
    ['Removed', Object.entries(s.removed).map(([k, n]) => `${k}×${n}`).join(', ') || 'nothing'],
    ['Pocket residues', state.pocketResi.length],
  ].map(([k, val]) => `<div><span>${k}</span><span>${val}</span></div>`).join('');
  $('meta').innerHTML += `<p class="hint" style="margin-top:10px">${s.title}</p>`;
}

function renderStyles() {
  const v = state.viewer;
  const { protein, crystal, pose } = state.models;
  v.removeAllSurfaces();
  v.removeAllShapes();

  if (protein) {
    v.setStyle({ model: protein }, state.shown.cartoon
      ? { cartoon: { color: COLOURS.protein, opacity: 0.85 } }
      : {});
    if (state.shown.pocket && state.pocketResi && state.pocketResi.length) {
      // surface only on the residues that line the pocket, so you can see inside
      v.addSurface($3Dmol.SurfaceType.VDW,
        { opacity: 0.7, color: COLOURS.surface },
        { model: protein, resi: state.pocketResi });
    }
  }
  if (crystal) {
    v.setStyle({ model: crystal }, state.shown.crystal
      ? { stick: { colorscheme: 'yellowCarbon', radius: 0.18 } }
      : {});
  }
  if (pose) {
    v.setStyle({ model: pose }, state.shown.pose
      ? { stick: { colorscheme: 'greenCarbon', radius: 0.21 } }
      : {});
  }
  if (state.shown.box && state.structure) {
    const { center, size } = state.structure.site;
    v.addBox({
      center: { x: center[0], y: center[1], z: center[2] },
      dimensions: { w: size[0], h: size[1], d: size[2] },
      color: '#7aa2f7', opacity: 0.12, wireframe: false,
    });
  }
  v.spin(state.shown.spin ? 'y' : false);
  v.render();
}

// ------------------------------------------------------------------- poses

async function loadDockedList() {
  let rows = [];
  try {
    rows = await api(`/structures/${state.pdb}/docked?limit=200`);
  } catch (e) { /* none yet */ }

  if (!rows.length) {
    $('pick').innerHTML = `<p class="muted">No docked poses yet for ${state.pdb}.<br><br>
      Run <code>python3 dock_validate.py ${state.pdb} 250</code> — poses are saved
      as it goes, so they appear here while it is still running.</p>`;
    return;
  }
  $('pick').innerHTML = `<table><thead><tr>
      <th>Compound</th><th class="n">Dock</th><th class="n">Real</th>
    </tr></thead><tbody>${rows.map(r => `
      <tr data-id="${r.id}">
        <td>${r.id}</td>
        <td class="n">${r.docking_score === null || r.docking_score === undefined ? '—' : r.docking_score.toFixed(2)}</td>
        <td class="n">${r.measured === null || r.measured === undefined ? '—' : r.measured.toFixed(2)}</td>
      </tr>`).join('')}</tbody></table>`;

  $('pick').onclick = async e => {
    const tr = e.target.closest('tr[data-id]');
    if (!tr) return;
    [...$('pick').querySelectorAll('tr[data-id]')].forEach(x =>
      x.setAttribute('aria-selected', x === tr));
    await showPose(tr.dataset.id);
  };
  // show the best-scoring one straight away
  await showPose(rows[0].id);
  const first = $('pick').querySelector('tr[data-id]');
  if (first) first.setAttribute('aria-selected', 'true');
}

async function showPose(compoundId) {
  let data;
  try {
    data = await api(`/structures/${state.pdb}/pose/${encodeURIComponent(compoundId)}`);
  } catch (e) {
    $('selected').innerHTML = `<p class="muted">${e.message}</p>`;
    return;
  }
  const v = state.viewer;
  if (state.models.pose) v.removeModel(state.models.pose);
  state.models.pose = v.addModel(data.pose_sdf, 'sdf');
  state.selected = data;
  renderStyles();

  $('selected').innerHTML = [
    ['Compound', data.compound_id],
    ['Docking score', data.docking_score === null ? '—' : `${data.docking_score} kcal/mol`],
    ['Measured pIC50', data.measured === null ? '—' : data.measured],
    ['MW', data.desc.mw],
    ['cLogP', data.desc.clogp],
    ['QED', data.desc.qed],
  ].map(([k, val]) => `<div><span>${k}</span><span>${val}</span></div>`).join('')
    + `<p class="hint" style="margin-top:10px;word-break:break-all">${data.smiles}</p>`;
}

// ------------------------------------------------------------------ controls

$('toolbar').onclick = e => {
  const b = e.target.closest('button[data-toggle]');
  if (!b) return;
  const key = b.dataset.toggle;
  state.shown[key] = !state.shown[key];
  b.setAttribute('aria-pressed', String(state.shown[key]));
  renderStyles();
};

$('reset').onclick = () => {
  state.viewer.zoomTo({ model: state.models.crystal });
  state.viewer.zoom(0.45);
  state.viewer.render();
};

$('pdb').onchange = async e => {
  state.pdb = e.target.value;
  await loadStructure(state.pdb);
  await loadDockedList();
};

boot().catch(e => {
  $('glLoading').textContent = 'Failed to load: ' + e.message;
});
