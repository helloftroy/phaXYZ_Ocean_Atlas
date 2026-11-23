/* PHA Ocean Atlas explorer: views, filters, table and export.
 *
 * Three views over the same data at different granularity (genomes, phaC
 * proteins, phaC clusters). Each declares its columns, its facets and how to
 * reach genome coordinates for the map; everything else -- filtering, faceted
 * counts, sorting, the virtualised table, CSV export -- is shared.
 *
 * Facet counts are computed with that facet's own selection removed, so
 * choosing one phylum does not zero out the others. That costs one extra pass
 * over the rows per facet, which at this size is a few milliseconds and worth
 * it: static counts make a filter UI feel broken.
 */
(() => {
  'use strict';

  const BLANK = -1;
  const ROW_HEIGHT = 27;
  const WINDOW_ROWS = 60;

  const el = (id) => document.getElementById(id);
  const loading = el('loading');
  const loadingText = el('loading-text');

  let atlas = null;
  let view = null;
  const views = {};
  const state = { view: 'genomes', search: '', sort: null, scaleByCount: true };
  let rows = [];
  let detailRows = null;      // rows behind the open map-cell panel, if any
  let map = null;

  /* ------------------------------------------------------------- facets */

  /* Each facet owns its widget, its current selection and a predicate over row
   * indices. `counts` is optional: facets that list coded values provide a
   * column so the shared code can tally them. */

  function checksFacet({ id, title, column, table, searchable = false, openByDefault = false }) {
    const selected = new Set();
    return {
      id, title, kind: 'checks', column, table, searchable, selected, openByDefault,
      active: () => selected.size > 0,
      clear: () => selected.clear(),
      predicate() {
        if (!selected.size) return null;
        const values = column;
        return (row) => selected.has(values[row]);
      },
      summary: () => (selected.size ? `${selected.size}` : ''),
    };
  }

  function chipsFacet({ id, title, options, test, mode = 'any', note, openByDefault = false }) {
    const selected = new Set();
    return {
      id, title, kind: 'chips', options, selected, mode, note, openByDefault,
      active: () => selected.size > 0,
      clear: () => selected.clear(),
      predicate() {
        if (!selected.size) return null;
        const picked = [...selected];
        return mode === 'all'
          ? (row) => picked.every((value) => test(row, value))
          : (row) => picked.some((value) => test(row, value));
      },
      summary: () => (selected.size ? `${selected.size}` : ''),
    };
  }

  function rangeFacet({ id, title, column, unit = '', step = 'any', places = 0, openByDefault = false }) {
    const bounds = { min: null, max: null };
    return {
      id, title, kind: 'range', column, unit, step, places, bounds, openByDefault,
      active: () => bounds.min !== null || bounds.max !== null,
      clear: () => { bounds.min = null; bounds.max = null; },
      predicate() {
        if (bounds.min === null && bounds.max === null) return null;
        const values = column;
        const lo = bounds.min === null ? -Infinity : bounds.min;
        const hi = bounds.max === null ? Infinity : bounds.max;
        return (row) => {
          const v = values[row];
          if (v === BLANK || Number.isNaN(v)) return false;
          return v >= lo && v <= hi;
        };
      },
      summary: () => {
        if (bounds.min !== null && bounds.max !== null) return `${bounds.min}–${bounds.max}`;
        if (bounds.min !== null) return `≥ ${bounds.min}`;
        if (bounds.max !== null) return `≤ ${bounds.max}`;
        return '';
      },
    };
  }

  function boolFacet({ id, title, column, yes = 'Yes', no = 'No', test, openByDefault = false }) {
    const choice = { value: 'any' };
    const check = test || ((row) => column[row] === 1);
    return {
      id, title, kind: 'bool', yes, no, choice, openByDefault,
      active: () => choice.value !== 'any',
      clear: () => { choice.value = 'any'; },
      predicate() {
        if (choice.value === 'any') return null;
        const want = choice.value === 'yes';
        return (row) => check(row) === want;
      },
      summary: () => (choice.value === 'any' ? '' : choice.value === 'yes' ? yes : no),
    };
  }

  /* ------------------------------------------------------------- filtering */

  function select(excludeFacetId) {
    const predicates = [];
    for (const facet of view.facets) {
      if (facet.id === excludeFacetId) continue;
      const predicate = facet.predicate();
      if (predicate) predicates.push(predicate);
    }
    const needle = state.search.trim().toLowerCase();
    const matchText = needle && view.searchText ? view.searchText(needle) : null;

    const out = [];
    const n = view.count();
    for (let row = 0; row < n; row++) {
      let ok = true;
      for (let p = 0; p < predicates.length; p++) {
        if (!predicates[p](row)) { ok = false; break; }
      }
      if (ok && matchText && !matchText(row)) ok = false;
      if (ok) out.push(row);
    }
    return out;
  }

  function countsFor(facet) {
    if (!facet.column) return null;
    const visible = select(facet.id);
    const tally = new Map();
    const values = facet.column;
    for (let i = 0; i < visible.length; i++) {
      const code = values[visible[i]];
      if (code === BLANK) continue;
      tally.set(code, (tally.get(code) || 0) + 1);
    }
    return tally;
  }

  /* ------------------------------------------------------------- rendering */

  function renderFacets() {
    const host = el('facet-list');
    host.textContent = '';
    for (const facet of view.facets) {
      const details = document.createElement('details');
      details.className = 'facet';
      // Open state lives on the facet, not the element: the panel is rebuilt on
      // every refresh and a section the user opened must stay open.
      if (facet.openState === undefined) facet.openState = !!facet.openByDefault;
      details.open = facet.openState || facet.active();
      details.addEventListener('toggle', () => { facet.openState = details.open; });

      const summary = document.createElement('summary');
      summary.append(document.createTextNode(facet.title));
      const badge = document.createElement('span');
      badge.className = 'count';
      badge.textContent = facet.summary();
      summary.append(badge);
      details.append(summary);

      const body = document.createElement('div');
      body.className = 'facet-body';
      details.append(body);
      buildFacetBody(facet, body);
      host.append(details);
    }
  }

  function buildFacetBody(facet, body) {
    if (facet.kind === 'checks') {
      const tally = countsFor(facet) || new Map();
      const list = document.createElement('div');
      list.className = 'checks';

      const paint = () => {
        const query = facet.query || '';
        list.textContent = '';
        const entries = [...tally.entries()]
          .map(([code, n]) => ({ code, n, name: facet.table[code] || '' }))
          .filter((e) => e.name && (!query || e.name.toLowerCase().includes(query)))
          .sort((a, b) => b.n - a.n || a.name.localeCompare(b.name));
        // Selected values must stay visible even when the search or the other
        // filters would otherwise hide them, or they cannot be unticked.
        for (const code of facet.selected) {
          if (!entries.some((e) => e.code === code)) {
            entries.unshift({ code, n: tally.get(code) || 0, name: facet.table[code] || String(code) });
          }
        }
        const shown = entries.slice(0, 220);
        for (const entry of shown) {
          const label = document.createElement('label');
          const box = document.createElement('input');
          box.type = 'checkbox';
          box.checked = facet.selected.has(entry.code);
          box.addEventListener('change', () => {
            box.checked ? facet.selected.add(entry.code) : facet.selected.delete(entry.code);
            refresh();
          });
          const text = document.createElement('span');
          text.textContent = entry.name;
          const n = document.createElement('span');
          n.className = 'n';
          n.textContent = Data.int(entry.n);
          label.append(box, text, n);
          list.append(label);
        }
        if (!shown.length) {
          const empty = document.createElement('p');
          empty.className = 'hint';
          empty.textContent = 'Nothing matches.';
          list.append(empty);
        } else if (entries.length > shown.length) {
          const more = document.createElement('p');
          more.className = 'hint';
          more.textContent = `${Data.int(entries.length - shown.length)} more — use the search box.`;
          list.append(more);
        }
      };

      if (facet.searchable) {
        const search = document.createElement('input');
        search.type = 'search';
        search.placeholder = 'Filter this list…';
        search.value = facet.query || '';
        search.addEventListener('input', () => {
          facet.query = search.value.trim().toLowerCase();
          paint();
        });
        body.append(search);
      }
      body.append(list);
      paint();
    } else if (facet.kind === 'chips') {
      const wrap = document.createElement('div');
      wrap.className = 'chips';
      for (const option of facet.options) {
        const chip = document.createElement('button');
        chip.type = 'button';
        chip.className = 'chip';
        chip.textContent = option.label;
        chip.setAttribute('aria-pressed', facet.selected.has(option.value) ? 'true' : 'false');
        if (option.title) chip.title = option.title;
        chip.addEventListener('click', () => {
          facet.selected.has(option.value)
            ? facet.selected.delete(option.value)
            : facet.selected.add(option.value);
          refresh();
        });
        wrap.append(chip);
      }
      body.append(wrap);
      if (facet.note) {
        const note = document.createElement('p');
        note.className = 'hint';
        note.textContent = facet.note;
        body.append(note);
      }
    } else if (facet.kind === 'range') {
      const wrap = document.createElement('div');
      wrap.className = 'range';
      const make = (which, placeholder) => {
        const input = document.createElement('input');
        input.type = 'number';
        input.step = facet.step;
        input.placeholder = placeholder;
        if (facet.bounds[which] !== null) input.value = facet.bounds[which];
        input.addEventListener('change', () => {
          const raw = input.value.trim();
          facet.bounds[which] = raw === '' ? null : Number(raw);
          refresh();
        });
        return input;
      };
      const dash = document.createElement('span');
      dash.textContent = '–';
      wrap.append(make('min', 'min'), dash, make('max', 'max'));
      body.append(wrap);
      if (facet.unit) {
        const note = document.createElement('p');
        note.className = 'hint';
        note.textContent = facet.unit;
        body.append(note);
      }
    } else if (facet.kind === 'bool') {
      for (const [value, text] of [['any', 'Any'], ['yes', facet.yes], ['no', facet.no]]) {
        const label = document.createElement('label');
        label.className = 'inline';
        const radio = document.createElement('input');
        radio.type = 'radio';
        radio.name = `bool-${facet.id}`;
        radio.checked = facet.choice.value === value;
        radio.addEventListener('change', () => { facet.choice.value = value; refresh(); });
        const span = document.createElement('span');
        span.textContent = text;
        label.append(radio, span);
        body.append(label);
      }
    }
  }

  function renderSummary() {
    const host = el('summary');
    host.textContent = '';
    for (const stat of view.stats(rows)) {
      const box = document.createElement('div');
      box.className = 'stat';
      const v = document.createElement('div');
      v.className = 'v';
      v.textContent = stat.v;
      const k = document.createElement('div');
      k.className = 'k';
      k.textContent = stat.k;
      if (stat.title) box.title = stat.title;
      box.append(v, k);
      host.append(box);
    }
  }

  function renderTable() {
    const head = el('thead-row');
    if (head.dataset.view !== state.view) {
      head.textContent = '';
      view.columns.forEach((column, index) => {
        const th = document.createElement('th');
        th.textContent = column.label;
        if (column.title) th.title = column.title;
        th.addEventListener('click', () => {
          const current = state.sort;
          state.sort = current && current.index === index
            ? { index, dir: -current.dir }
            : { index, dir: column.type === 'text' ? 1 : -1 };
          sortRows();
          renderTable();
        });
        head.append(th);
      });
      head.dataset.view = state.view;
    }
    [...head.children].forEach((th, index) => {
      if (state.sort && state.sort.index === index) {
        th.setAttribute('aria-sort', state.sort.dir > 0 ? 'ascending' : 'descending');
      } else {
        th.removeAttribute('aria-sort');
      }
    });
    paintWindow();

    el('table-note').textContent = view.tableNote
      ? view.tableNote(rows)
      : `${Data.int(rows.length)} rows match.`;
  }

  /* Only the visible slice is in the DOM; spacer rows carry the rest of the
   * scroll height, so 241,864 rows scroll without 241,864 elements. */
  function paintWindow() {
    const scroller = el('table-scroll');
    const body = el('tbody');
    const total = rows.length;
    const first = Math.max(0, Math.floor(scroller.scrollTop / ROW_HEIGHT) - 10);
    const last = Math.min(total, first + WINDOW_ROWS + 20);

    body.textContent = '';
    const columns = view.columns;

    if (first > 0) {
      const spacer = document.createElement('tr');
      const cell = document.createElement('td');
      cell.colSpan = columns.length;
      cell.style.height = `${first * ROW_HEIGHT}px`;
      cell.style.padding = '0';
      cell.style.borderBottom = 'none';
      spacer.append(cell);
      body.append(spacer);
    }

    for (let i = first; i < last; i++) {
      const row = rows[i];
      const tr = document.createElement('tr');
      for (const column of columns) {
        const td = document.createElement('td');
        if (column.cls) td.className = column.cls;
        const value = column.cell ? column.cell(row) : column.value(row);
        if (value instanceof Node) td.append(value);
        else td.textContent = value === null || value === undefined ? '' : value;
        tr.append(td);
      }
      body.append(tr);
    }

    if (last < total) {
      const spacer = document.createElement('tr');
      const cell = document.createElement('td');
      cell.colSpan = columns.length;
      cell.style.height = `${(total - last) * ROW_HEIGHT}px`;
      cell.style.padding = '0';
      cell.style.borderBottom = 'none';
      spacer.append(cell);
      body.append(spacer);
    }
  }

  function sortRows() {
    if (!state.sort) return;
    const column = view.columns[state.sort.index];
    const dir = state.sort.dir;
    const key = column.sortKey || column.value;
    const numeric = column.type !== 'text';
    rows.sort((a, b) => {
      let va = key(a);
      let vb = key(b);
      if (numeric) {
        va = va === '' || va === null || va === undefined ? -Infinity : Number(va);
        vb = vb === '' || vb === null || vb === undefined ? -Infinity : Number(vb);
        if (Number.isNaN(va)) va = -Infinity;
        if (Number.isNaN(vb)) vb = -Infinity;
        return (va - vb) * dir;
      }
      return String(va).localeCompare(String(vb)) * dir;
    });
  }

  function renderMap() {
    el('map-title').textContent = view.mapTitle;
    const points = view.mapPoints(rows);
    map.setPoints(points, {
      scaleByCount: state.scaleByCount,
      describe: view.describePoint,
      onSelect: (bin) => (bin ? renderDetail(bin) : hideDetail()),
    });
    hideDetail();
    updateMapNote();
  }

  function updateMapNote() {
    const plotted = map.plotted || 0;
    el('map-note').textContent = plotted
      ? `${Data.int(plotted)} occurrences in ${Data.int(map.bins.length)} cells of ` +
        `${map.resolution}°. Scroll to zoom, drag to pan, click a cell for what is in it.`
      : 'No records in this selection carry coordinates.';
    el('zoom-level').textContent = `${map.zoom.toFixed(1)}×`;
  }

  /* ------------------------------------------------------------- detail */

  /* A clicked cell holds genome rows. In the genomes view those are the records
   * themselves; in the other two they are the genomes an occurrence links to,
   * so the proteins or clusters present there have to be looked up and then
   * narrowed to the ones the current filter kept. */
  function rowsInCell(bin) {
    if (state.view === 'genomes') return bin.rows;
    const allowed = new Set(rows);
    const group = state.view === 'proteins' ? atlas.proteinsIn : atlas.clustersIn;
    return Data.keysIn(group, bin.rows, allowed);
  }

  const DETAIL_LIMIT = 200;

  function renderDetail(bin) {
    const found = rowsInCell(bin);
    detailRows = found;

    const north = bin.lat < 0 ? 'S' : 'N';
    const east = bin.lon < 0 ? 'W' : 'E';
    el('detail-title').textContent =
      `${Data.int(found.length)} ${found.length === 1 ? view.unit.slice(0, -1) : view.unit} ` +
      `near ${Math.abs(bin.lat).toFixed(1)}°${north}, ${Math.abs(bin.lon).toFixed(1)}°${east}` +
      (state.view === 'genomes' ? '' : ` (in ${Data.int(bin.count)} genome records)`);

    const body = el('detail-body');
    body.textContent = '';
    if (!found.length) {
      const empty = document.createElement('p');
      empty.className = 'hint';
      empty.textContent = 'Nothing here matches the current filters.';
      body.append(empty);
    } else {
      const columns = view.columns.slice(0, 6);
      const table = document.createElement('table');
      const head = document.createElement('thead');
      const headRow = document.createElement('tr');
      for (const column of columns) {
        const th = document.createElement('th');
        th.textContent = column.label;
        headRow.append(th);
      }
      head.append(headRow);
      table.append(head);

      const tbody = document.createElement('tbody');
      for (const row of found.slice(0, DETAIL_LIMIT)) {
        const tr = document.createElement('tr');
        for (const column of columns) {
          const td = document.createElement('td');
          if (column.cls) td.className = column.cls;
          const value = column.cell ? column.cell(row) : column.value(row);
          if (value instanceof Node) td.append(value);
          else td.textContent = value === null || value === undefined ? '' : value;
          tr.append(td);
        }
        tbody.append(tr);
      }
      table.append(tbody);
      body.append(table);

      if (found.length > DETAIL_LIMIT) {
        const more = document.createElement('p');
        more.className = 'hint';
        more.textContent =
          `Showing the first ${DETAIL_LIMIT}. Download for all ${Data.int(found.length)}.`;
        body.append(more);
      }
    }
    el('detail').hidden = false;
  }

  function hideDetail() {
    detailRows = null;
    el('detail').hidden = true;
  }

  function refresh() {
    rows = select(null);
    sortRows();
    renderFacets();
    renderSummary();
    renderTable();
    renderMap();
  }

  /* ------------------------------------------------------------- views */

  function buildViews() {
    const { lookups, genomes, proteins, clusters } = atlas;
    const L = lookups;
    const g = genomes.cols;
    const p = proteins.cols;
    const c = clusters.cols;
    const familyIndex = new Map(L.families.map((name, i) => [name, i]));

    const famCount = (row, family) => g.fam[row * 16 + familyIndex.get(family)];

    /* ---- genomes ---- */
    views.genomes = {
      title: 'Genomes', unit: 'genomes',
      mapTitle: 'Where these genomes were sampled',
      count: () => genomes.n,
      facets: [
        chipsFacet({
          id: 'fam', openByDefault: true,
          title: 'Must carry these genes',
          mode: 'all',
          options: L.families.map((name) => ({
            value: name,
            label: name.replace('_synthase', ' (syn)').replace('_regulator', ' (reg)'),
            title: `Genomes with at least one ${name} gene`,
          })),
          test: (row, family) => famCount(row, family) > 0,
        }),
        checksFacet({ id: 'ph', title: 'Phylum', column: g.ph, table: L.phylum, searchable: true, openByDefault: true }),
        checksFacet({ id: 'gn', title: 'Genus', column: g.gn, table: L.genus, searchable: true }),
        checksFacet({ id: 'hb', title: 'Habitat', column: g.hb, table: L.habitat, searchable: true }),
        rangeFacet({ id: 'dep', title: 'Depth', column: g.dep, unit: 'metres below surface', step: '1' }),
        checksFacet({ id: 'arch', title: 'Pathway architecture', column: g.arch, table: L.architecture }),
        checksFacet({ id: 'st', title: 'Study', column: g.st, table: L.study, searchable: true }),
        boolFacet({
          id: 'coords', title: 'Has coordinates', column: null, yes: 'Mapped', no: 'No coordinates',
          test: (row) => !Number.isNaN(g.lat[row]) && !Number.isNaN(g.lon[row]),
        }),
      ],
      columns: [
        {
          label: 'Genome', type: 'text', cls: 'mono',
          value: (row) => (atlas.names ? atlas.names[row] : `row ${row}`),
        },
        { label: 'Phylum', type: 'text', value: (row) => Data.label(L.phylum, g.ph[row]) },
        { label: 'Genus', type: 'text', value: (row) => Data.label(L.genus, g.gn[row]) },
        { label: 'Habitat', type: 'text', cls: 'wrap', value: (row) => Data.label(L.habitat, g.hb[row]) },
        {
          label: 'Depth (m)', type: 'number', cls: 'num',
          value: (row) => Data.dec(g.dep[row], 0), sortKey: (row) => g.dep[row],
        },
        {
          label: 'PHA genes', type: 'text', cls: 'wrap',
          title: 'Families present; phaC is shown dark',
          cell: (row) => {
            const frag = document.createDocumentFragment();
            L.families.forEach((family, i) => {
              const n = g.fam[row * 16 + i];
              if (!n) return;
              const chip = document.createElement('span');
              chip.className = family === 'phaC' ? 'fam c' : 'fam';
              chip.textContent = n > 1 ? `${shortFamily(family)}×${n}` : shortFamily(family);
              frag.append(chip);
            });
            return frag;
          },
          sortKey: (row) => {
            let total = 0;
            for (let i = 0; i < 16; i++) if (g.fam[row * 16 + i]) total++;
            return total;
          },
          exportValue: (row) => L.families
            .map((family, i) => {
              const n = g.fam[row * 16 + i];
              return n ? (n > 1 ? `${family}x${n}` : family) : null;
            })
            .filter(Boolean)
            .join(' '),
        },
        { label: 'Architecture', type: 'text', cls: 'wrap', value: (row) => Data.label(L.architecture, g.arch[row]) },
        { label: 'Study', type: 'text', value: (row) => Data.label(L.study, g.st[row]) },
      ],
      stats: (matched) => {
        let phac = 0;
        let mapped = 0;
        const phyla = new Set();
        for (let i = 0; i < matched.length; i++) {
          const row = matched[i];
          if (g.phac[row] === 1) phac++;
          if (!Number.isNaN(g.lat[row])) mapped++;
          if (g.ph[row] !== BLANK) phyla.add(g.ph[row]);
        }
        return [
          { v: Data.int(matched.length), k: 'genomes' },
          { v: Data.int(phac), k: 'carry phaC' },
          { v: Data.int(phyla.size), k: 'phyla' },
          { v: Data.int(mapped), k: 'with coordinates' },
        ];
      },
      mapPoints: (matched) => ({ lat: g.lat, lon: g.lon, rows: matched }),
      describePoint: (sample) => {
        const names = new Set(sample.map((row) => Data.label(L.phylum, g.ph[row])).filter(Boolean));
        return names.size ? [...names].slice(0, 3).join(', ') : '';
      },
      searchText: (needle) => (row) => {
        const name = atlas.names ? atlas.names[row] : '';
        if (name && name.toLowerCase().includes(needle)) return true;
        const genus = Data.label(L.genus, g.gn[row]);
        return genus ? genus.toLowerCase().includes(needle) : false;
      },
      exportName: 'pha_atlas_genomes',
    };

    /* ---- phaC proteins ---- */
    views.proteins = {
      title: 'phaC proteins', unit: 'proteins',
      mapTitle: 'Where these phaC proteins occur',
      count: () => proteins.n,
      facets: [
        chipsFacet({
          id: 'cls', title: 'Synthase class', openByDefault: true,
          options: [
            { value: 'I', label: 'Class I' },
            { value: 'II', label: 'Class II' },
            { value: 'III', label: 'Class III' },
          ],
          note: 'Best-scoring class HMM. Class III and IV cannot be separated reliably and are reported together.',
          test: (row, value) => p.cls[row] === value,
        }),
        boolFacet({ id: 'clsok', title: 'Clears the class trusted cutoff', column: p.clsok, yes: 'Yes', no: 'No' }),
        boolFacet({
          id: 'tr', title: 'Catalytic triad complete', column: p.tr, openByDefault: true,
          yes: 'Complete', no: 'Not resolved by alignment',
        }),
        rangeFacet({ id: 'pid', title: 'Identity to nearest reference', column: p.pid, unit: 'per cent', step: '1' }),
        rangeFacet({ id: 'tcov', title: 'Reference coverage', column: p.tcov, unit: 'fraction, 0–1', step: '0.05' }),
        rangeFacet({ id: 'qtm', title: 'Structural match (qTM)', column: p.qtm, unit: 'fraction, 0–1; measured for representatives only', step: '0.05' }),
        chipsFacet({
          id: 'th', title: 'Predicted thermostable at',
          options: L.thresholds.map((t) => ({ value: t, label: `≥${t}°C` })),
          note: 'TemStaPro prediction, not a measured melting temperature.',
          test: (row, value) => p.th[row] !== BLANK && p.th[row] >= value,
        }),
        rangeFacet({ id: 'ng', title: 'Genomes carrying it', column: p.ng, unit: 'uncapped genome count', step: '1' }),
        rangeFacet({ id: 'len', title: 'Protein length', column: p.len, unit: 'amino acids', step: '10' }),
        checksFacet({ id: 'ph', title: 'Phylum', column: p.ph, table: L.phylum, searchable: true }),
        checksFacet({ id: 'gn', title: 'Genus', column: p.gn, table: L.genus, searchable: true }),
      ],
      columns: [
        { label: 'Protein', type: 'text', cls: 'mono', value: (row) => p.id[row] },
        { label: 'Length', type: 'number', cls: 'num', value: (row) => Data.int(p.len[row]), sortKey: (row) => p.len[row] },
        { label: 'Class', type: 'text', value: (row) => p.cls[row] || '' },
        {
          label: 'Cutoff', type: 'number', title: 'Clears the class HMM trusted cutoff',
          cell: (row) => tick(p.clsok[row] === 1), sortKey: (row) => p.clsok[row],
          exportValue: (row) => (p.clsok[row] === 1 ? 'yes' : 'no'),
        },
        {
          label: 'Triad', type: 'number', title: 'Catalytic triad complete by HMM alignment columns',
          cell: (row) => tick(p.tr[row] === 1), sortKey: (row) => p.tr[row],
          exportValue: (row) => (p.tr[row] === 1 ? 'complete' : 'not resolved'),
        },
        { label: '% id', type: 'number', cls: 'num', value: (row) => Data.dec(p.pid[row], 1), sortKey: (row) => p.pid[row] },
        { label: 'Ref cov', type: 'number', cls: 'num', value: (row) => Data.dec(p.tcov[row], 2), sortKey: (row) => p.tcov[row] },
        {
          label: 'qTM', type: 'number', cls: 'num', title: 'Structural similarity to its best reference',
          value: (row) => Data.dec(p.qtm[row], 2), sortKey: (row) => p.qtm[row],
        },
        {
          label: 'Thermo', type: 'number', cls: 'num', title: 'Highest TemStaPro threshold passed',
          value: (row) => (p.th[row] > 0 ? `≥${p.th[row]}°C` : ''), sortKey: (row) => p.th[row],
        },
        { label: 'Genomes', type: 'number', cls: 'num', value: (row) => Data.int(p.ng[row]), sortKey: (row) => p.ng[row] },
        { label: 'Cluster', type: 'text', cls: 'mono', value: (row) => p.cl[row] || '' },
        { label: 'Phylum', type: 'text', value: (row) => Data.label(L.phylum, p.ph[row]) },
        { label: 'Genus', type: 'text', value: (row) => Data.label(L.genus, p.gn[row]) },
      ],
      stats: (matched) => {
        let triad = 0;
        let thermo = 0;
        let genomes = 0;
        const clusterSet = new Set();
        for (let i = 0; i < matched.length; i++) {
          const row = matched[i];
          if (p.tr[row] === 1) triad++;
          if (p.th[row] >= 60) thermo++;
          genomes += Math.max(p.ng[row], 0);
          if (p.cl[row]) clusterSet.add(p.cl[row]);
        }
        return [
          { v: Data.int(matched.length), k: 'proteins' },
          { v: Data.int(clusterSet.size), k: 'clusters' },
          { v: Data.int(triad), k: 'triad complete' },
          { v: Data.int(thermo), k: 'stable ≥60°C' },
          { v: Data.int(genomes), k: 'genome occurrences', title: 'Sum over proteins; a genome with two paralogs counts twice' },
        ];
      },
      mapPoints: (matched) => ({ lat: g.lat, lon: g.lon, rows: Data.genomesFor(atlas.byProtein, matched) }),
      describePoint: (sample) => {
        const names = new Set(sample.map((row) => Data.label(L.genus, g.gn[row])).filter(Boolean));
        return names.size ? [...names].slice(0, 3).join(', ') : '';
      },
      searchText: (needle) => (row) =>
        p.id[row].includes(needle) ||
        (p.cl[row] || '').includes(needle) ||
        (Data.label(L.genus, p.gn[row]) || '').toLowerCase().includes(needle),
      tableNote: (matched) =>
        `${Data.int(matched.length)} of ${Data.int(proteins.n)} phaC proteins match. ` +
        'A blank qTM means the protein was not among the dereplication representatives, not that it scored low.',
      exportName: 'pha_atlas_phac_proteins',
    };

    /* ---- phaC clusters ---- */
    views.clusters = {
      title: 'phaC clusters', unit: 'clusters',
      mapTitle: 'Where these phaC clusters occur',
      count: () => clusters.n,
      facets: [
        rangeFacet({ id: 'ng', title: 'Genomes in the cluster', column: c.ng, unit: 'uncapped genome count', step: '1', openByDefault: true }),
        rangeFacet({ id: 'nloc', title: 'Distinct sampling locations', column: c.nloc, unit: 'count', step: '1' }),
        rangeFacet({ id: 'nph', title: 'Phyla spanned', column: c.nph, unit: 'count', step: '1' }),
        rangeFacet({ id: 'np', title: 'Distinct protein sequences', column: c.np, unit: 'count', step: '1' }),
        checksFacet({ id: 'ph', title: 'Dominant phylum', column: c.ph, table: L.phylum, searchable: true }),
        checksFacet({ id: 'gen', title: 'Top genus', column: c.gen, table: L.genus, searchable: true }),
      ],
      columns: [
        { label: 'Cluster', type: 'text', cls: 'mono', value: (row) => c.id[row] },
        { label: 'Genomes', type: 'number', cls: 'num', value: (row) => Data.int(c.ng[row]), sortKey: (row) => c.ng[row] },
        { label: 'Proteins', type: 'number', cls: 'num', value: (row) => Data.int(c.np[row]), sortKey: (row) => c.np[row] },
        { label: 'Locations', type: 'number', cls: 'num', value: (row) => Data.int(c.nloc[row]), sortKey: (row) => c.nloc[row] },
        { label: 'Classes', type: 'number', cls: 'num', value: (row) => Data.int(c.ncls[row]), sortKey: (row) => c.ncls[row] },
        { label: 'Phyla', type: 'number', cls: 'num', value: (row) => Data.int(c.nph[row]), sortKey: (row) => c.nph[row] },
        { label: 'Dominant phylum', type: 'text', value: (row) => Data.label(L.phylum, c.ph[row]) },
        { label: 'Top genus', type: 'text', value: (row) => Data.label(L.genus, c.gen[row]) },
      ],
      stats: (matched) => {
        let singletons = 0;
        let oneLocation = 0;
        let genomes = 0;
        for (let i = 0; i < matched.length; i++) {
          const row = matched[i];
          if (c.ng[row] === 1) singletons++;
          if (c.nloc[row] === 1) oneLocation++;
          genomes += Math.max(c.ng[row], 0);
        }
        return [
          { v: Data.int(matched.length), k: 'clusters' },
          { v: Data.int(singletons), k: 'in one genome' },
          { v: Data.int(oneLocation), k: 'at one location' },
          { v: Data.int(genomes), k: 'genome memberships' },
        ];
      },
      mapPoints: (matched) => ({ lat: g.lat, lon: g.lon, rows: Data.genomesFor(atlas.byCluster, matched) }),
      describePoint: (sample) => {
        const names = new Set(sample.map((row) => Data.label(L.genus, g.gn[row])).filter(Boolean));
        return names.size ? [...names].slice(0, 3).join(', ') : '';
      },
      searchText: (needle) => (row) =>
        c.id[row].includes(needle) ||
        (Data.label(L.genus, c.gen[row]) || '').toLowerCase().includes(needle),
      exportName: 'pha_atlas_phac_clusters',
    };
  }

  function shortFamily(name) {
    return name.replace('pha', '').replace('_synthase', 'sy').replace('_regulator', 'rg');
  }

  function tick(on) {
    const span = document.createElement('span');
    span.className = on ? 'yes' : 'no';
    span.textContent = on ? '✓' : '—';
    return span;
  }

  /* ------------------------------------------------------------- export */

  function exportCSV(which, tag) {
    const target = which || rows;
    const columns = view.columns.map((column) => column.label);
    const text = Data.toCSV(columns, target, (row, index) => {
      const column = view.columns[index];
      if (column.exportValue) return column.exportValue(row);
      if (column.cell) return column.sortKey ? column.sortKey(row) : '';
      const value = column.value ? column.value(row) : '';
      // Strip the grouping separators the table adds, so the CSV stays numeric.
      return typeof value === 'string' && column.type === 'number'
        ? value.replace(/,/g, '') : value;
    });
    const stamp = new Date().toISOString().slice(0, 10);
    const suffix = tag ? `_${tag}` : '';
    Data.download(`${view.exportName}${suffix}_${target.length}rows_${stamp}.csv`, text);
  }

  /* ------------------------------------------------------------- about */

  function renderAbout() {
    const host = el('about-body');
    const counts = atlas.meta.counts;
    host.textContent = '';

    const intro = document.createElement('p');
    intro.textContent =
      `This site indexes ${Data.int(counts.genomes_in_matrix)} ocean genomes carrying at least one ` +
      `polyhydroxyalkanoate biosynthesis gene, drawn from a survey of 274,282 genomes in the Ocean ` +
      `Microbiomics Database. ${Data.int(counts.phac_positive_genomes)} carry a PHA synthase (phaC), ` +
      `represented by ${Data.int(counts.phac_proteins)} distinct protein sequences in ` +
      `${Data.int(counts.phac_clusters)} clusters at 70% identity. Everything shown is filtered and ` +
      `plotted in your browser; nothing is sent anywhere.`;
    host.append(intro);

    const caveatHead = document.createElement('h3');
    caveatHead.textContent = 'Read before you use these numbers';
    host.append(caveatHead);
    const list = document.createElement('ol');
    for (const caveat of atlas.meta.caveats) {
      const item = document.createElement('li');
      item.textContent = caveat;
      list.append(item);
    }
    host.append(list);

    const dataHead = document.createElement('h3');
    dataHead.textContent = 'Downloads and source';
    host.append(dataHead);
    const links = document.createElement('p');
    links.innerHTML =
      'Any filtered selection exports as CSV from the button above the table. ' +
      'The full archive — per-family annotation tables, cluster assignments, catalytic-site ' +
      'scoring, predicted structures, thermostability predictions and the figure source tables ' +
      '— is deposited on Zenodo (DOI pending). Analysis code and the scripts that generate ' +
      'this site are at <a href="https://github.com/helloftroy/phaXYZ_Ocean_Atlas">github.com/helloftroy/phaXYZ_Ocean_Atlas</a>.';
    host.append(links);

    const built = document.createElement('p');
    built.className = 'hint';
    built.textContent = `Data generated ${atlas.meta.built}. Coastline: Natural Earth 110 m (public domain).`;
    host.append(built);
  }

  /* ------------------------------------------------------------- startup */

  function switchView(name) {
    state.view = name;
    view = views[name];
    state.sort = null;
    state.search = '';
    el('search').value = '';
    for (const facet of view.facets) facet.clear();
    document.querySelectorAll('.tabs button').forEach((button) => {
      button.setAttribute('aria-selected', button.dataset.view === name ? 'true' : 'false');
    });
    el('table-title').textContent = view.title;
    el('thead-row').dataset.view = '';
    el('table-scroll').scrollTop = 0;
    refresh();
  }

  function wire() {
    document.querySelectorAll('.tabs button').forEach((button) => {
      button.addEventListener('click', () => switchView(button.dataset.view));
    });
    el('reset').addEventListener('click', () => {
      for (const facet of view.facets) { facet.clear(); facet.query = ''; }
      state.search = '';
      el('search').value = '';
      refresh();
    });

    let searchTimer = null;
    el('search').addEventListener('input', (event) => {
      clearTimeout(searchTimer);
      const value = event.target.value;
      searchTimer = setTimeout(() => { state.search = value; refresh(); }, 180);
    });

    el('download').addEventListener('click', () => exportCSV());
    el('map-scale').addEventListener('change', (event) => {
      state.scaleByCount = event.target.checked;
      map.setScaleByCount(state.scaleByCount);
    });

    el('zoom-in').addEventListener('click', () => { map.zoomBy(1.6); updateMapNote(); });
    el('zoom-out').addEventListener('click', () => { map.zoomBy(1 / 1.6); updateMapNote(); });
    el('zoom-reset').addEventListener('click', () => { map.resetView(); hideDetail(); updateMapNote(); });
    el('map').addEventListener('wheel', () => requestAnimationFrame(updateMapNote), { passive: true });
    el('map').addEventListener('dblclick', () => requestAnimationFrame(updateMapNote));

    el('detail-close').addEventListener('click', () => { hideDetail(); map.clearSelection(); });
    el('detail-download').addEventListener('click', () => {
      if (detailRows && detailRows.length) exportCSV(detailRows, 'map_cell');
    });

    let scrollFrame = null;
    el('table-scroll').addEventListener('scroll', () => {
      if (scrollFrame) return;
      scrollFrame = requestAnimationFrame(() => { scrollFrame = null; paintWindow(); });
    });
  }

  async function main() {
    try {
      let done = 0;
      atlas = await Data.load(() => {
        done++;
        loadingText.textContent = `Loading the atlas… (${done}/7)`;
      });
    } catch (error) {
      loading.innerHTML = '';
      const message = document.createElement('p');
      message.textContent = `Could not load the data: ${error.message}`;
      loading.append(message);
      return;
    }

    buildViews();
    map = new AtlasMap.Map2D(el('map'), el('tip'));
    map.setLand(atlas.land);

    document.querySelectorAll('[data-stat="genomes"]').forEach((node) => {
      node.textContent = Data.int(atlas.meta.counts.genomes_in_matrix);
    });
    document.querySelectorAll('[data-stat="built"]').forEach((node) => {
      node.textContent = atlas.meta.built;
    });

    renderAbout();
    wire();
    loading.hidden = true;
    el('app').hidden = false;
    el('about').hidden = false;
    switchView('genomes');

    // Genome names arrive after first paint; repaint the table once they land.
    atlas.namesPromise.then((names) => {
      if (!names) return;
      atlas.names = names;
      if (state.view === 'genomes') paintWindow();
    });
  }

  main();
})();
