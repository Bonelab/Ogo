/* The exported CSV is the study record. Browser storage is only a convenience. */
function decisionKey(row) { return JSON.stringify([row.site, row.model_id]); }

function inclusionRecord(row, decision = {}) {
  const manual = decision.manual_decision || 'automatic';
  if (!['automatic', 'include', 'exclude'].includes(manual)) throw new Error('Invalid manual decision');
  return {
    model_id: row.model_id, site: row.site, model_path: row.model_path || '',
    automatic_qc_status: row.qc_status, automatic_qc_reasons: row.qc_reasons || '',
    manual_decision: manual,
    manual_reason: manual === 'automatic' ? '' : (decision.manual_reason || 'manual_review'),
    manual_note: manual === 'automatic' ? '' : (decision.manual_note || ''),
    reviewer: decision.reviewer || '', review_timestamp: decision.review_timestamp || '',
    final_inclusion: manual === 'automatic' ? ({pass: 'include', fail: 'exclude'}[row.qc_status] || 'pending') : manual,
    solver_eligible: row.generation_status === 'too_short' ? 'False' : '',
    qc_policy_version: row.qc_policy_version || '', measurement_source: row.measurement_source || ''
  };
}

function encodeCSV(records, fields) {
  const quote = value => '"' + String(value ?? '').replaceAll('"', '""') + '"';
  return [fields.map(quote).join(','), ...records.map(row => fields.map(key => quote(row[key])).join(','))].join('\r\n') + '\r\n';
}

function parseCSV(text) {
  text = text.replace(/^\uFEFF/, '');
  const rows = []; let row = [], field = '', quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (c === '"') {
      if (quoted && text[i + 1] === '"') { field += '"'; i++; }
      else quoted = !quoted;
    } else if (!quoted && c === ',') { row.push(field); field = ''; }
    else if (!quoted && (c === '\r' || c === '\n')) {
      if (c === '\r' && text[i + 1] === '\n') i++;
      row.push(field); if (row.some(value => value !== '')) rows.push(row);
      row = []; field = '';
    } else field += c;
  }
  if (quoted) throw new Error('Unterminated CSV quotation');
  if (row.length || field) { row.push(field); rows.push(row); }
  const headers = rows.shift() || [];
  if (new Set(headers).size !== headers.length) throw new Error('Duplicate CSV headers');
  if (!['model_id', 'site', 'manual_decision'].every(key => headers.includes(key))) throw new Error('Missing review CSV columns');
  return rows.map(values => {
    if (values.length !== headers.length) throw new Error('Invalid CSV row length');
    return Object.fromEntries(headers.map((key, index) => [key, values[index]]));
  });
}

function validateReviews(records, rows) {
  const known = new Map(rows.map(row => [decisionKey(row), row]));
  const decisions = new Map();
  for (const record of records) {
    const key = decisionKey(record);
    if (!known.has(key)) throw new Error('Unknown model in review CSV: ' + record.model_id);
    if (decisions.has(key)) throw new Error('Duplicate model in review CSV: ' + record.model_id);
    inclusionRecord(known.get(key), record);
    decisions.set(key, record);
  }
  return decisions;
}

function imageVisible(image, stage, camera) {
  return image.stage === stage && (image.view === camera || image.fallback === 'true');
}

function pageIndices(matches, page, size) { return matches.slice(page * size, (page + 1) * size); }

function excludedByFlags(row, decision, flags) {
  if (flags.has('manual_review_pending') && (!decision.manual_decision || decision.manual_decision === 'automatic')) return true;
  return (row.qc_reasons || '').split(';').some(reason => flags.has(reason));
}

function sortIndices(indices, rows, mode) {
  const patient = (a, b) => rows[a].model_id.localeCompare(rows[b].model_id, undefined, {numeric: true});
  const count = row => new Set((row.qc_reasons || '').split(';').filter(Boolean)).size;
  return [...indices].sort((a, b) => {
    if (mode === 'patient') return patient(a, b);
    if (mode === 'flags') return count(rows[b]) - count(rows[a]) || patient(a, b);
    const [key, order] = mode.split(':');
    const value = row => {
      const raw = row[key];
      return raw === '' || raw == null || !Number.isFinite(Number(raw)) ? null : Math.abs(Number(raw));
    };
    const x = value(rows[a]), y = value(rows[b]);
    if (x === null || y === null) return x === y ? patient(a, b) : x === null ? 1 : -1;
    return (order === 'desc' ? y - x : x - y) || patient(a, b);
  });
}

function preserveScroll(action, viewport = window) {
  const x = viewport.scrollX, y = viewport.scrollY;
  action();
  viewport.scrollTo({left: x, top: y, behavior: 'instant'});
}

if (typeof module !== 'undefined') module.exports = {inclusionRecord, encodeCSV, parseCSV, validateReviews, imageVisible, pageIndices, excludedByFlags, sortIndices, preserveScroll};

if (typeof document !== 'undefined') {
  const data = JSON.parse(document.getElementById('review-data').textContent);
  const rows = data.rowValues.map(([group, values]) => Object.fromEntries(data.rowColumns[group].map((key, i) => [key, values[i]])));
  for (const row of rows) Object.assign(row, (data.outcomes || {})[row.model_id] || {});
  const $ = id => document.getElementById(id);
  let decisions = validateReviews(data.seeds, rows), selected = null, autosave = true;
  let page = 0;
  const previewObserver = typeof IntersectionObserver === 'undefined' ? null : new IntersectionObserver(entries => {
    for (const entry of entries) {
      if (!entry.isIntersecting || entry.target.hidden) continue;
      entry.target.src = entry.target.dataset.thumbnail || entry.target.dataset.src;
      previewObserver.unobserve(entry.target);
    }
  }, {rootMargin: '1200px 0px'});
  const message = text => { $('message').textContent = text; };
  try {
    const saved = localStorage.getItem(data.storageKey);
    if (saved && !data.imported) decisions = validateReviews(JSON.parse(saved), rows);
  } catch (error) { autosave = false; }
  function persist() {
    try { localStorage.setItem(data.storageKey, JSON.stringify([...decisions.values()])); }
    catch (error) { autosave = false; }
  }
  function updateImages(container, full = false) {
    for (const image of container.querySelectorAll('img')) {
      image.hidden = !imageVisible(image.dataset, $('stage').value, $('camera').value);
      image.style.height = $('camera').value === 'all' || image.dataset.fallback === 'true' ? '560px' : '280px';
      if (!image.hidden) {
        if (full || !previewObserver) image.src = full ? image.dataset.src : (image.dataset.thumbnail || image.dataset.src);
        else previewObserver.observe(image);
      } else {
        if (previewObserver) previewObserver.unobserve(image);
        image.removeAttribute('src');
      }
    }
  }
  function measurementTable(row) {
    const table = document.createElement('table');
    for (const [key, value] of Object.entries(row).sort(([a], [b]) => a.localeCompare(b))) {
      if (key.endsWith('_image')) continue;
      const tr = table.insertRow(), th = document.createElement('th');
      th.textContent = key; tr.append(th); tr.insertCell().textContent = value;
    }
    return table;
  }
  function imagesFor(index) {
    return (data.images[index] || []).map(asset => {
      const image = document.createElement('img');
      Object.assign(image.dataset, asset);
      image.loading = 'eager'; image.decoding = 'async';
      image.alt = asset.stage + ': ' + asset.view;
      image.hidden = true;
      return image;
    });
  }
  function makeCard(index) {
    const card = document.createElement('article');
    card.dataset.index = index;
    card.innerHTML = '<button class="review-open"></button><p class="automatic"></p><p class="decision"></p><button class="quick-decision" data-decision="include" aria-pressed="false">Include</button> <button class="quick-decision" data-decision="exclude" aria-pressed="false">Exclude</button> <button class="review-open">Review</button><details><summary>Measurements</summary></details>';
    card.querySelector('.review-open').textContent = rows[index].model_id;
    card.append(...imagesFor(index));
    card.querySelector('details').addEventListener('toggle', event => {
      if (event.target.open && !event.target.querySelector('table')) event.target.append(measurementTable(rows[index]));
    });
    return card;
  }
  function refresh() {
    const pageSize = Number($('page-size').value);
    const search = $('search').value.toLowerCase();
    const excluded = new Set([...$('reason-picker').querySelectorAll('input:checked')].map(input => input.dataset.qcReason));
    $('reason-summary').textContent = 'Exclude flags' + (excluded.size ? ' (' + excluded.size + ')' : '');
    let matches = [];
    rows.forEach((row, index) => {
      if (!row.model_id.toLowerCase().includes(search) ||
        ($('site').value && row.site !== $('site').value) ||
        ($('status').value && row.qc_status !== $('status').value) ||
        excludedByFlags(row, decisions.get(decisionKey(row)) || {}, excluded) ||
        ($('inclusion').value && inclusionRecord(row, decisions.get(decisionKey(row))).final_inclusion !== $('inclusion').value)) return;
      matches.push(index);
    });
    matches = sortIndices(matches, rows, $('sort').value);
    const pages = Math.max(1, Math.ceil(matches.length / pageSize));
    page = Math.min(page, pages - 1);
    if (previewObserver) previewObserver.disconnect();
    const cards = pageIndices(matches, page, pageSize).map(makeCard);
    document.querySelector('main').replaceChildren(...cards);
    for (const card of cards) {
      const row = rows[Number(card.dataset.index)], record = inclusionRecord(row, decisions.get(decisionKey(row)));
      card.querySelector('.automatic').textContent = 'Automatic: ' + row.qc_status + (row.qc_reasons ? ': ' + row.qc_reasons : '');
      card.querySelector('.decision').textContent = record.manual_decision === 'exclude' ? 'Manually excluded' : 'Study: ' + record.final_inclusion + (record.manual_decision !== 'automatic' ? ' (manual)' : '');
      for (const button of card.querySelectorAll('.quick-decision')) {
        button.setAttribute('aria-pressed', String(record.manual_decision === button.dataset.decision));
        button.style.fontWeight = record.manual_decision === button.dataset.decision ? '700' : '400';
      }
      updateImages(card);
    }
    $('previous-page').disabled = page === 0;
    $('next-page').disabled = page === pages - 1;
    $('page-status').textContent = 'Page ' + (page + 1) + ' / ' + pages;
    $('review-stage').value = $('stage').value;
    $('review-camera').value = $('camera').value;
    updateImages($('review-images'), true);
    $('message').textContent = matches.length + ' / ' + rows.length + ' models. ' + (autosave ? '' : 'Browser autosave unavailable. ') + 'Export CSV to preserve the study record.';
  }
  function openReview(index) {
    selected = rows[index];
    const record = inclusionRecord(selected, decisions.get(decisionKey(selected)));
    $('review-title').textContent = selected.model_id + ': ' + selected.site;
    $('review-auto').textContent = 'Automatic: ' + selected.qc_status + '; ' + selected.qc_reasons +
      (selected.generation_status === 'too_short' ? '. Too short: manual inclusion does not make this model solver-eligible.' : '');
    $('manual-reason').value = data.reasons.includes(record.manual_reason) ? record.manual_reason : 'other';
    if (!record.manual_reason) $('manual-reason').value = 'manual_review';
    $('manual-note').value = record.manual_note;
    $('review-images').replaceChildren(...imagesFor(index));
    updateImages($('review-images'), true);
    $('review-measurements').replaceChildren(measurementTable(selected));
    $('review-dialog').showModal();
  }
  function decide(manual) {
    const decision = {manual_decision: manual, manual_reason: $('manual-reason').value,
      manual_note: $('manual-note').value, reviewer: $('reviewer').value.trim() || 'unspecified',
      review_timestamp: new Date().toISOString()};
    decisions.set(decisionKey(selected), inclusionRecord(selected, decision));
    persist();
    preserveScroll(() => { $('review-dialog').close(); refresh(); });
  }
  document.querySelector('main').addEventListener('click', event => {
    const card = event.target.closest('article');
    if (!card) return;
    if (event.target.closest('.review-open,img')) openReview(Number(card.dataset.index));
    const button = event.target.closest('.quick-decision');
    if (button) {
      const row = rows[Number(card.dataset.index)];
      const decision = {manual_decision: button.dataset.decision, manual_reason: 'manual_review',
        manual_note: '', reviewer: $('reviewer').value.trim() || 'unspecified',
        review_timestamp: new Date().toISOString()};
      decisions.set(decisionKey(row), inclusionRecord(row, decision));
      persist(); preserveScroll(refresh);
    }
  });
  $('previous-page').onclick = () => { page--; refresh(); window.scrollTo(0, 0); };
  $('next-page').onclick = () => { page++; refresh(); window.scrollTo(0, 0); };
  for (const id of ['search', 'site', 'status', 'inclusion', 'stage', 'camera', 'page-size', 'sort']) $(id).addEventListener('input', () => { page = 0; refresh(); });
  $('reason-picker').addEventListener('change', () => { page = 0; refresh(); });
  $('clear-flags').onclick = () => {
    for (const input of $('reason-picker').querySelectorAll('input')) input.checked = false;
    page = 0; refresh();
  };
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape') $('reason-picker').open = false;
  });
  for (const id of ['stage', 'camera']) $('review-' + id).addEventListener('change', () => {
    $(id).value = $('review-' + id).value;
    refresh();
  });
  $('close').onclick = () => $('review-dialog').close();
  for (const id of ['include', 'exclude']) $(id).onclick = () => decide(id);
  $('reset').onclick = () => decide('automatic');
  $('export').onclick = () => {
    const csv = encodeCSV(rows.map(row => inclusionRecord(row, decisions.get(decisionKey(row)))), data.fields);
    const url = URL.createObjectURL(new Blob([csv], {type: 'text/csv;charset=utf-8'}));
    const link = document.createElement('a'); link.href = url; link.download = 'study_inclusion.csv'; link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    message('CSV export started: check your browser downloads for study_inclusion.csv. The project CSV is unchanged.');
    window.alert('CSV export started.\n\nCheck your browser downloads for study_inclusion.csv.\n\nThis is a separate file: it does not update the project CSV or the CSV inside the gallery ZIP.');
  };
  $('import').onchange = async event => {
    try {
      const records = parseCSV(await event.target.files[0].text());
      const imported = validateReviews(records, rows);
      for (const [key, record] of imported) decisions.set(key, record);
      persist(); refresh(); message('Imported ' + imported.size + ' decisions. Automatic QC remains unchanged.');
    } catch (error) { message('Import rejected: ' + error.message); }
    event.target.value = '';
  };
  refresh();
}
