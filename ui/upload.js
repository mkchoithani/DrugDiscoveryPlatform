// Upload a programme CSV, confirm the columns, watch it run, read the report.
const $ = id => document.getElementById(id);
const state = { upload: null, job: null, poll: null };

function show(step) {
  for (const s of ['step1', 'step2', 'step3', 'step4']) $(s).classList.toggle('on', s === step);
  if (step !== 'step4') $('error').innerHTML = '';
}

function fail(message, detail) {
  if (state.poll) { clearInterval(state.poll); state.poll = null; }
  show(state.upload ? 'step2' : 'step1');
  $('error').innerHTML = `<div class="err"><b>${escapeHtml(message)}</b>
    ${detail ? `<code>${escapeHtml(detail)}</code>` : ''}</div>`;
}

const escapeHtml = s => String(s).replace(/[&<>"]/g, c =>
  ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// ------------------------------------------------------------------ step 1

const drop = $('drop');
drop.onclick = () => $('file').click();
drop.ondragover = e => { e.preventDefault(); drop.classList.add('over'); };
drop.ondragleave = () => drop.classList.remove('over');
drop.ondrop = e => {
  e.preventDefault();
  drop.classList.remove('over');
  if (e.dataTransfer.files[0]) inspect(e.dataTransfer.files[0]);
};
$('file').onchange = e => { if (e.target.files[0]) inspect(e.target.files[0]); };

// Sample datasets, so the page can be tried without hunting for a file.
document.addEventListener('click', async e => {
  const link = e.target.closest('a[data-sample]');
  if (!link) return;
  e.preventDefault();
  const name = link.dataset.sample;
  try {
    const res = await fetch(`./samples/${name}`);
    if (!res.ok) throw new Error(`sample not found (${res.status})`);
    inspect(new File([await res.blob()], name, { type: 'text/csv' }));
  } catch (err) {
    fail('Could not load the sample file.', err.message);
  }
});

async function inspect(file) {
  if (!/\.(csv|tsv|txt)$/i.test(file.name)) {
    return fail('That does not look like a CSV file.', file.name);
  }
  const body = new FormData();
  body.append('file', file);
  try {
    const res = await fetch('/api/evaluate/inspect', { method: 'POST', body });
    const data = await res.json();
    if (!res.ok) return fail('Could not read the file.', data.detail || res.statusText);
    state.upload = data;
    fillMapping(data);
    show('step2');
  } catch (e) {
    fail('Could not reach the server.', e.message);
  }
}

function fillMapping(d) {
  const option = (col, selected, blank) =>
    `<option value="${escapeHtml(col)}"${col === selected ? ' selected' : ''}>${escapeHtml(col || blank)}</option>`;
  const fill = (el, selected, optional) => {
    el.innerHTML = (optional ? option('', null, '— none —') : '')
      + d.columns.map(c => option(c, selected)).join('');
  };
  fill($('cSmiles'), d.detected.smiles, false);
  fill($('cValue'), d.detected.value, false);
  fill($('cId'), d.detected.id, true);
  fill($('cDate'), d.detected.date, true);
  $('cUnits').value = d.units;

  $('fileInfo').textContent =
    `${d.filename} · ${d.columns.length} columns` +
    (d.missing.length ? ` · could not detect: ${d.missing.join(', ')} — please pick them` : '');

  $('peek').innerHTML = d.sample.length
    ? `<table><thead><tr>${d.columns.map(c => `<th>${escapeHtml(c)}</th>`).join('')}</tr></thead>
       <tbody>${d.sample.map(r => `<tr>${d.columns.map(c =>
         `<td>${escapeHtml(r[c] ?? '')}</td>`).join('')}</tr>`).join('')}</tbody></table>`
    : '';
  $('timeHint').textContent = 'usually under a minute';
}

$('back').onclick = () => { state.upload = null; show('step1'); };

// ------------------------------------------------------------------ step 3

$('run').onclick = async () => {
  if (!state.upload) return;
  const body = new FormData();
  body.append('upload_id', state.upload.upload_id);
  body.append('smiles', $('cSmiles').value);
  body.append('value', $('cValue').value);
  body.append('id_column', $('cId').value);
  body.append('date', $('cDate').value);
  body.append('units', $('cUnits').value);
  body.append('customer', $('cCustomer').value);
  body.append('filename', state.upload.filename || '');

  show('step3');
  $('stage').textContent = 'starting…';
  try {
    const res = await fetch('/api/evaluate/run', { method: 'POST', body });
    const data = await res.json();
    if (!res.ok) return fail('Could not start the evaluation.', data.detail || res.statusText);
    state.job = data.job_id;
    state.poll = setInterval(checkJob, 1200);
  } catch (e) {
    fail('Could not reach the server.', e.message);
  }
};

async function checkJob() {
  try {
    const res = await fetch(`/api/evaluate/${state.job}/status`);
    const job = await res.json();
    if (!res.ok) return fail('Lost track of the evaluation.', job.detail);
    $('stage').textContent = job.stage || job.status;
    if (job.status === 'error') {
      return fail(job.error || 'The evaluation failed.', job.trace);
    }
    if (job.status === 'done') {
      clearInterval(state.poll); state.poll = null;
      const url = `/api/evaluate/${state.job}/report`;
      $('report').src = url;
      $('openReport').href = url;
      show('step4');
    }
  } catch (e) {
    fail('Could not reach the server.', e.message);
  }
}

$('again').onclick = () => { state.upload = null; state.job = null; $('file').value = ''; show('step1'); };
