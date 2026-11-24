/* Checks for the explorer site in docs/, to run after every rebuild.
 *
 *   node site_build/check_site.mjs
 *
 * There is no browser here, so this drives the real docs/*.js in a Node vm
 * against the real generated files, with stubs for fetch, canvas and pointer
 * events. It covers the three things that have actually broken: generated data
 * disagreeing with its own manifest, page wiring drifting from the scripts, and
 * map geometry. Exit status is non-zero on any failure.
 */
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import vm from 'node:vm';

const DOCS = join(dirname(fileURLToPath(import.meta.url)), '..', 'docs');
const read = (name) => readFileSync(join(DOCS, name), 'utf8');

let failures = 0;
function fail(message) { failures++; console.log(`FAIL  ${message}`); }
function check(label, got, want, tolerance = 0) {
  const ok = tolerance
    ? Math.abs(got - want) <= tolerance
    : String(got) === String(want);
  if (ok) console.log(`PASS  ${label}: ${got}`);
  else fail(`${label}: ${got} (expected ${want}${tolerance ? ` +-${tolerance}` : ''})`);
}
function section(title) { console.log(`\n--- ${title} ---`); }

/* ------------------------------------------------------------ load the data */

function makeSandbox(extra = {}) {
  const sandbox = {
    console, Intl, Math, JSON, Map, Set, Number, Array, String, Object, Promise,
    Float32Array, Int32Array, Uint8Array, setTimeout, clearTimeout, Blob, URL,
    requestAnimationFrame: (fn) => { fn(); return 1; },
    document: {
      documentElement: {},
      createElement: () => ({ style: {}, click() {}, remove() {} }),
      body: { appendChild() {} },
    },
    getComputedStyle: () => ({ getPropertyValue: () => '#123456' }),
    // The loader appends ?v=<build stamp> for cache busting; strip it to a path.
    fetch: async (path) => ({
      ok: true, status: 200, statusText: 'OK',
      text: async () => read(String(path).split('?')[0]),
    }),
    ...extra,
  };
  sandbox.globalThis = sandbox;
  vm.createContext(sandbox);
  return sandbox;
}

const sandbox = makeSandbox({ window: { devicePixelRatio: 1, addEventListener() {} } });
vm.runInContext(read('map.js'), sandbox, { filename: 'map.js' });
vm.runInContext(read('data.js'), sandbox, { filename: 'data.js' });

const atlas = await vm.runInContext('Data.load()', sandbox);
const names = await atlas.namesPromise;
const { meta, genomes, proteins, clusters, occurrences, byProtein, byCluster, proteinsIn, clustersIn } = atlas;
const g = genomes.cols, p = proteins.cols, c = clusters.cols;
const keysIn = vm.runInContext('Data.keysIn', sandbox);
const genomesFor = vm.runInContext('Data.genomesFor', sandbox);
const counts = meta.counts;

section('the data matches its own manifest');
check('genome rows', genomes.n, counts.genomes_total);
check('protein rows', proteins.n, counts.phac_proteins);
check('cluster rows', clusters.n, counts.phac_clusters);
check('occurrence rows', occurrences.n, counts.mappable_occurrences);
check('genome names', names.length, genomes.n);
check('proteins with qTM', [...p.qtm].filter((v) => !Number.isNaN(v)).length, counts.proteins_with_qtm);

let withCoords = 0, withHabitat = 0, withDepth = 0, phacFlag = 0, phacFamily = 0, unprofiled = 0;
for (let i = 0; i < genomes.n; i++) {
  if (!Number.isNaN(g.lat[i]) && !Number.isNaN(g.lon[i])) withCoords++;
  if (g.hb[i] !== -1) withHabitat++;
  if (!Number.isNaN(g.dep[i])) withDepth++;
  if (g.phac[i] === 1) phacFlag++;
  if (g.fam[i * 16 + 2] > 0) phacFamily++;          // index 2 is phaC
  if (g.prof[i] === 0) unprofiled++;
}
check('with coordinates', withCoords, counts.genomes_with_coordinates);
check('with a habitat', withHabitat, counts.genomes_with_habitat);
check('with a depth', withDepth, counts.genomes_with_depth);
check('phaC carriers', phacFlag, counts.phac_positive_genomes);
check('phaC flag agrees with the family column', phacFamily, phacFlag);
check('unprofiled genomes', unprofiled, counts.genomes_unprofiled);
check('nothing unplaceable', counts.unmappable_occurrences, 0);

section('no placeholder leaked through as a value');
check('blank protein ids', [...p.id].filter((v) => !v).length, 0);
for (const junk of ['nan', 'NA', 'DEPRECATED', 'Unknown', 'None']) {
  const inLabels = ['phylum', 'genus', 'habitat', 'ecosystem', 'study', 'depth_zone']
    .flatMap((key) => atlas.lookups[key])
    .filter((label) => label === junk).length;
  if (inLabels) fail(`"${junk}" is being shown as a real label`);
}
console.log('PASS  lookup tables carry no placeholder labels');

section('index integrity');
const sumGroup = (group) => group.offsets[group.offsets.length - 1];
check('byProtein items', sumGroup(byProtein), occurrences.n);
check('byCluster items', sumGroup(byCluster), occurrences.n);
check('proteinsIn items', sumGroup(proteinsIn), occurrences.n);
let outOfRange = 0;
for (const group of [byProtein, byCluster]) {
  for (let i = 0; i < group.items.length; i++) {
    if (group.items[i] < 0 || group.items[i] >= genomes.n) outOfRange++;
  }
}
check('genome rows in range', outOfRange, 0);

// proteinsIn must be the exact inverse of byProtein: the popup depends on it.
const forward = new Set();
for (let i = 0; i < proteins.n; i++) {
  for (let k = byProtein.offsets[i]; k < byProtein.offsets[i + 1]; k++) forward.add(`${i}:${byProtein.items[k]}`);
}
let missingEdges = 0;
for (let row = 0; row < genomes.n; row++) {
  for (let k = proteinsIn.offsets[row]; k < proteinsIn.offsets[row + 1]; k++) {
    if (!forward.has(`${proteinsIn.items[k]}:${row}`)) missingEdges++;
  }
}
check('reverse index is an exact inverse', missingEdges, 0);
check('every phaC genome reachable from proteins',
  genomesFor(byProtein, [...Array(proteins.n).keys()]).length, counts.phac_positive_genomes);

section('known facts from the project');
const deepDive = names.indexOf('CARD22-1_SAMN24292811_MAG_00000010');
check('Desulfoluna deep-dive genome present', deepDive >= 0, true);
check('  phylum', atlas.lookups.phylum[g.ph[deepDive]], 'Desulfobacterota');
check('  genus', atlas.lookups.genus[g.gn[deepDive]], 'Desulfoluna');
check('  habitat', atlas.lookups.habitat[g.hb[deepDive]], 'Coral tissue');
check('  latitude', g.lat[deepDive], 22.18, 0.01);
check('  longitude', g.lon[deepDive], 39.02, 0.01);
const copies = keysIn(proteinsIn, [deepDive], new Set([...Array(proteins.n).keys()]));
check('  phaC copies', copies.length, 7);
check('  triad complete by alignment', copies.filter((i) => p.tr[i] === 1).length, 5);

let singletons = 0, oneLocation = 0, largest = 0;
for (let i = 0; i < clusters.n; i++) {
  if (c.ng[i] === 1) singletons++;
  if (c.nloc[i] === 1) oneLocation++;
  largest = Math.max(largest, c.ng[i]);
}
check('singleton clusters', singletons, 1752);
check('single-location clusters', oneLocation, 1882);
check('largest cluster', largest, 1201);

section('page wiring');
const html = read('index.html');
const app = read('app.js');
const mapSrc = read('map.js');
const css = read('style.css');

const htmlIds = new Set([...html.matchAll(/\bid="([^"]+)"/g)].map((m) => m[1]));
const usedIds = new Set([...app.matchAll(/\bel\('([^']+)'\)/g)].map((m) => m[1]));
for (const id of usedIds) if (!htmlIds.has(id)) fail(`app.js uses #${id}, absent from the page`);
console.log(`PASS  ${usedIds.size} element ids all present`);

const scripts = [...html.matchAll(/<script src="([^"?]+)/g)].map((m) => m[1]);
if (scripts.indexOf('app.js') < Math.max(scripts.indexOf('data.js'), scripts.indexOf('map.js'))) {
  fail('app.js is loaded before its dependencies');
} else console.log(`PASS  script order: ${scripts.join(' -> ')}`);

for (const name of new Set([...mapSrc.matchAll(/style\('(--[a-z0-9-]+)'/g)].map((m) => m[1]))) {
  if (!css.includes(`${name}:`)) fail(`map.js reads ${name}, undefined in the stylesheet`);
}
console.log('PASS  css variables used by the map are defined');

// The virtualised table sizes its spacer rows from ROW_HEIGHT.
check('row height agrees with the stylesheet',
  Number(app.match(/ROW_HEIGHT\s*=\s*(\d+)/)[1]),
  Number(css.match(/tbody td \{[^}]*height:\s*(\d+)px/)[1]));

// hidden only gets display:none from the UA sheet, so author rules beat it.
const override = css.match(/\[hidden\]\s*\{[^}]*display:\s*none\s*!important/);
if (!override) fail('no "[hidden] { display: none !important }" rule');
else {
  const at = css.indexOf(override[0]);
  for (const selector of ['\\.loading', 'main']) {
    const rule = css.match(new RegExp(`(^|\\n)${selector}\\s*\\{[^}]*display:`));
    if (rule && css.indexOf(rule[0]) < at) fail(`${selector} sets display before the [hidden] override`);
  }
  console.log('PASS  hidden elements really hide');
}

const views = [...app.matchAll(/views\.(\w+) = \{/g)].map((m) => m[1]);
for (const tab of [...html.matchAll(/data-view="([^"]+)"/g)].map((m) => m[1])) {
  if (!views.includes(tab)) fail(`tab "${tab}" has no view`);
}
console.log(`PASS  views: ${views.join(', ')}`);

section('map geometry');
// The projection must agree with figures/scripts/_simple_world_map.py, so that
// a point lands where the manuscript figures would put it.
const EXPECTED = [
  [0, 0, 0, 0], [180, 0, 2.666270, 0], [0, 90, 0, 1.352300],
  [39.02, 22.18, 0.565381, 0.371926], [-70.5, -45.25, -0.934414, -0.757457],
];
for (const [lon, lat, wantX, wantY] of EXPECTED) {
  const [x, y] = vm.runInContext(`AtlasMap.project(${lon}, ${lat})`, sandbox);
  check(`project(${lon}, ${lat})`, `${x.toFixed(6)},${y.toFixed(6)}`,
    `${wantX.toFixed(6)},${wantY.toFixed(6)}`);
}

const resolutions = [1, 2, 4, 8, 20].map((z) => vm.runInContext(`AtlasMap.resolutionForZoom(${z})`, sandbox));
for (let i = 1; i < resolutions.length; i++) {
  if (resolutions[i] > resolutions[i - 1]) fail('grid coarsens as zoom increases');
}
console.log(`PASS  grid follows zoom: ${resolutions.join(', ')} degrees`);

section('map interaction');
const handlers = { canvas: {}, window: {} };
const listen = (bag) => (type, fn) => { (bag[type] ||= []).push(fn); };
const fire = (bag, type, event) => (handlers[bag][type] || []).forEach((fn) => fn(event));
const WIDTH = 900, HEIGHT = 420;
const canvas = {
  width: 0, height: 0, clientWidth: WIDTH, clientHeight: HEIGHT, style: {},
  parentElement: { clientWidth: WIDTH },
  getContext: () => new Proxy({}, { get: () => () => {}, set: () => true }),
  addEventListener: listen(handlers.canvas),
  getBoundingClientRect: () => ({ left: 0, top: 0, width: WIDTH, height: HEIGHT }),
};
const mapBox = makeSandbox({
  window: { devicePixelRatio: 1, addEventListener: listen(handlers.window) },
  __canvas: canvas,
  __tip: { hidden: true, style: {}, getBoundingClientRect: () => ({ width: 120, height: 40 }) },
  __selected: [],
});
vm.runInContext(read('map.js'), mapBox, { filename: 'map.js' });

const lat = [], lon = [], rows = [];
for (const [la, lo] of [[22.18, 39.02], [-45.25, -70.5], [60.1, -20.4], [0.5, 120.9]]) {
  for (let k = 0; k < 25; k++) { lat.push(la + k * 0.004); lon.push(lo + k * 0.004); rows.push(lat.length - 1); }
}
mapBox.__pts = { lat: Float32Array.from(lat), lon: Float32Array.from(lon), rows };
const run = (expr) => vm.runInContext(expr, mapBox);
run(`globalThis.m = new AtlasMap.Map2D(__canvas, __tip);
     m.setLand([]); m.setPoints(__pts, { scaleByCount: true, onSelect: (b) => __selected.push(b) });`);

check('all points binned', run('m.plotted'), 100);
check('cells at full extent', run('m.bins.length'), 4);

// Zooming must hold the geography under the pointer. Track a fixed coordinate,
// not a cell centre: cell centres legitimately move as the grid gets finer.
const screenOf = (ln, lt) => run(`(() => {
  const [px, py] = AtlasMap.project(${ln}, ${lt}); const t = m.transform;
  return [t.cx + m.panX + px * t.scale * m.zoom, t.cy + m.panY - py * t.scale * m.zoom];
})()`);
for (const [ln, lt, factor] of [[39.02, 22.18, 4], [-150, -60, 1.8]]) {
  const before = screenOf(ln, lt);
  run(`m._zoomAt(${before[0]}, ${before[1]}, ${factor})`);
  const after = screenOf(ln, lt);
  check(`zoom holds ${ln},${lt} in x`, after[0], before[0], 0.01);
  check(`zoom holds ${ln},${lt} in y`, after[1], before[1], 0.01);
}
check('grid finer when zoomed', run('m.resolution') < 2, true);
run('m.resetView()');
check('reset restores zoom', run('m.zoom'), 1);
check('reset restores grid', run('m.resolution'), 2);

// At zoom 1 the whole map fits, so panning is correctly clamped to nothing.
run('m.panX = 400; m.redraw();');
check('no panning at full extent', run('m.panX'), 0);

const spot = run('(() => { const b = m.bins[m.bins.length - 1]; return [b.x, b.y, b.count]; })()');
fire('canvas', 'mousedown', { button: 0, clientX: spot[0], clientY: spot[1] });
fire('window', 'mouseup', { clientX: spot[0], clientY: spot[1] });
check('click selects a cell', run('__selected.length'), 1);
check('cell keeps every row for the panel', run('m.selected.rows.length'), spot[2]);

run('m.zoomBy(4);');
const panBefore = run('m.panX');
fire('canvas', 'mousedown', { button: 0, clientX: spot[0], clientY: spot[1] });
fire('window', 'mousemove', { clientX: spot[0] + 60, clientY: spot[1] + 20 });
fire('window', 'mouseup', { clientX: spot[0] + 60, clientY: spot[1] + 20 });
check('drag pans', run('m.panX') !== panBefore, true);
check('drag does not select', run('__selected.length'), 1);

fire('canvas', 'mousedown', { button: 0, clientX: spot[0], clientY: spot[1] });
fire('window', 'mousemove', { clientX: spot[0] + 2, clientY: spot[1] + 1 });
fire('window', 'mouseup', { clientX: spot[0] + 2, clientY: spot[1] + 1 });
check('a small twitch is still a click', run('__selected.length'), 2);

console.log(`\n${failures === 0 ? 'all checks passed' : `${failures} FAILURE(S)`}`);
process.exit(failures === 0 ? 0 : 1);
