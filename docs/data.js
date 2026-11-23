/* Loading and parsing for the PHA Ocean Atlas explorer.
 *
 * The data files are plain CSV with integer codes in place of repeated strings
 * (taxonomy, habitat, study), and the code tables live in lookups.json. Parsing
 * goes straight into typed arrays, one per column, because every filter is a
 * linear scan over 241,864 genome rows and a scan over typed arrays is fast
 * enough that no index is needed -- which keeps the whole site static.
 *
 * Splitting on "," without quote handling is safe here and only here: the
 * generator (site_build/build_site_data.py) emits no string column into these
 * files. Genome names live one-per-line in their own file, and every other
 * label is an integer code resolved through lookups.json.
 */
const Data = (() => {
  'use strict';

  const BLANK_INT = -1;

  async function fetchText(path, onProgress) {
    const response = await fetch(path);
    if (!response.ok) throw new Error(`${path}: ${response.status} ${response.statusText}`);
    if (onProgress) onProgress(path);
    return response.text();
  }

  async function fetchJSON(path, onProgress) {
    return JSON.parse(await fetchText(path, onProgress));
  }

  /* Parse a generated CSV into one typed array per column.
   * types maps a header name to 'f32' | 'i32' | 'str' | 'fam'. A column absent
   * from types is skipped, so adding a column to the generator cannot break
   * an older copy of this file. */
  function parseColumnar(text, types) {
    const lines = text.split('\n');
    const header = lines[0].split(',');
    let last = lines.length;
    while (last > 1 && lines[last - 1] === '') last--;   // tolerate trailing newline
    const n = last - 1;

    const cols = {};
    const order = [];
    header.forEach((name, position) => {
      const type = types[name];
      if (!type) return;
      order.push({ name, position, type });
      if (type === 'f32') cols[name] = new Float32Array(n).fill(NaN);
      else if (type === 'i32') cols[name] = new Int32Array(n).fill(BLANK_INT);
      else if (type === 'fam') cols[name] = new Uint8Array(n * 16);
      else cols[name] = new Array(n);
    });

    for (let row = 0; row < n; row++) {
      const parts = lines[row + 1].split(',');
      for (let f = 0; f < order.length; f++) {
        const { name, position, type } = order[f];
        const raw = parts[position];
        if (type === 'fam') {
          // 16 single-digit copy numbers, one per PHA family, capped at 9
          const base = row * 16;
          if (raw) for (let k = 0; k < 16 && k < raw.length; k++) cols[name][base + k] = raw.charCodeAt(k) - 48;
        } else if (raw === undefined || raw === '') {
          continue;                                       // keep the fill value
        } else if (type === 'f32') {
          cols[name][row] = +raw;
        } else if (type === 'i32') {
          cols[name][row] = raw.charCodeAt(0) === 45 ? BLANK_INT : +raw;
        } else {
          cols[name][row] = raw;
        }
      }
    }
    return { n, cols };
  }

  const GENOME_TYPES = {
    lat: 'f32', lon: 'f32', ph: 'i32', gn: 'i32', hb: 'i32', eco: 'i32',
    st: 'i32', dep: 'f32', dz: 'i32', fam: 'fam', arch: 'i32', phac: 'i32',
  };

  const PROTEIN_TYPES = {
    id: 'str', len: 'i32', cl: 'str', cls: 'str', clsok: 'i32', tr: 'i32',
    pid: 'f32', tcov: 'f32', qtm: 'f32', th: 'i32', ng: 'i32',
    ph: 'i32', gn: 'i32', ref: 'i32',
  };

  const CLUSTER_TYPES = {
    id: 'str', ng: 'i32', np: 'i32', nloc: 'i32', ncls: 'i32', nph: 'i32',
    ph: 'i32', gen: 'i32',
  };

  const OCCURRENCE_TYPES = { p: 'str', c: 'str', g: 'i32' };

  async function load(onProgress) {
    const [meta, lookups, land, genomeText, proteinText, clusterText, occurrenceText] =
      await Promise.all([
        fetchJSON('data/meta.json', onProgress),
        fetchJSON('data/lookups.json', onProgress),
        fetchJSON('data/land.json', onProgress),
        fetchText('data/genomes.csv', onProgress),
        fetchText('data/proteins.csv', onProgress),
        fetchText('data/clusters.csv', onProgress),
        fetchText('data/occurrences.csv', onProgress),
      ]);

    const genomes = parseColumnar(genomeText, GENOME_TYPES);
    const proteins = parseColumnar(proteinText, PROTEIN_TYPES);
    const clusters = parseColumnar(clusterText, CLUSTER_TYPES);
    const occurrences = parseColumnar(occurrenceText, OCCURRENCE_TYPES);

    // Protein and cluster ids to row numbers, so an occurrence resolves without scanning.
    const proteinRow = new Map();
    for (let i = 0; i < proteins.n; i++) proteinRow.set(proteins.cols.id[i], i);
    const clusterRow = new Map();
    for (let i = 0; i < clusters.n; i++) clusterRow.set(clusters.cols.id[i], i);

    // Occurrences grouped by protein row and by cluster row, as flat index lists.
    const byProtein = groupOccurrences(occurrences, proteinRow, 'p');
    const byCluster = groupOccurrences(occurrences, clusterRow, 'c');

    // Genome names are 8 MB and only the table and CSV export need them, so they
    // arrive separately and the first paint does not wait on them.
    const namesPromise = fetchText('data/genome_names.csv')
      .then((text) => {
        const names = text.split('\n');
        if (names.length && names[names.length - 1] === '') names.pop();
        return names;
      })
      .catch(() => null);

    return {
      meta, lookups, land, genomes, proteins, clusters, occurrences,
      proteinRow, clusterRow, byProtein, byCluster, namesPromise,
      names: null,
    };
  }

  /* Build a CSR-style grouping: offsets[k]..offsets[k+1] indexes items for key k. */
  function groupOccurrences(occurrences, keyRow, keyColumn) {
    const keys = occurrences.cols[keyColumn];
    const genomeRows = occurrences.cols.g;
    const counts = new Int32Array(keyRow.size + 1);
    const resolved = new Int32Array(occurrences.n).fill(BLANK_INT);

    for (let i = 0; i < occurrences.n; i++) {
      const row = keyRow.get(keys[i]);
      if (row === undefined) continue;
      resolved[i] = row;
      counts[row + 1]++;
    }
    for (let k = 1; k < counts.length; k++) counts[k] += counts[k - 1];

    const offsets = counts;
    const cursor = offsets.slice();
    const items = new Int32Array(offsets[offsets.length - 1]);
    for (let i = 0; i < occurrences.n; i++) {
      const row = resolved[i];
      if (row === BLANK_INT) continue;
      items[cursor[row]++] = genomeRows[i];
    }
    return { offsets, items };
  }

  /* Every genome row reachable from a filtered set of protein or cluster rows. */
  function genomesFor(group, rows) {
    const seen = new Set();
    const out = [];
    for (let r = 0; r < rows.length; r++) {
      const key = rows[r];
      for (let i = group.offsets[key]; i < group.offsets[key + 1]; i++) {
        const genome = group.items[i];
        if (!seen.has(genome)) { seen.add(genome); out.push(genome); }
      }
    }
    return out;
  }

  function countFor(group, key) {
    return group.offsets[key + 1] - group.offsets[key];
  }

  /* ----------------------------------------------------------- formatting */

  const groupedNumber = new Intl.NumberFormat('en-US');

  function int(value) {
    return value === null || value === undefined || value === BLANK_INT || Number.isNaN(value)
      ? '' : groupedNumber.format(value);
  }

  function dec(value, places) {
    return value === null || value === undefined || Number.isNaN(value)
      ? '' : value.toFixed(places);
  }

  function label(table, code) {
    return code === BLANK_INT || code === undefined || code === null ? '' : (table[code] || '');
  }

  /* ----------------------------------------------------------- CSV export */

  function toCSV(columns, rows, valueAt) {
    const quote = (value) => {
      const text = value === null || value === undefined ? '' : String(value);
      return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
    };
    const out = [columns.join(',')];
    for (let r = 0; r < rows.length; r++) {
      const row = rows[r];
      out.push(columns.map((_, c) => quote(valueAt(row, c))).join(','));
    }
    return out.join('\n') + '\n';
  }

  function download(filename, text) {
    const blob = new Blob([text], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  return { load, BLANK_INT, genomesFor, countFor, int, dec, label, toCSV, download };
})();
