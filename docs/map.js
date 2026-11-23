/* Canvas occurrence map for the PHA Ocean Atlas explorer.
 *
 * Robinson projection, the same one the manuscript figures use, with the same
 * interpolation table -- so a point sits where the paper would put it. Drawing
 * is plain canvas with a locally generated coastline: no tile service, nothing
 * to rate-limit, and the map works offline.
 *
 * Occurrences are binned onto a degree grid before drawing. A filtered set can
 * hold 227,000 points, most of them at repeated station coordinates, so drawing
 * them individually would be both slow and misleading -- a hundred genomes from
 * one cruise station would look like one dot. Binning makes the count visible
 * and the hit-test cheap.
 *
 * The grid follows the zoom: at 2 degrees a cell is a sensible unit of the
 * whole ocean, but zoomed in it would merge stations that are plainly separate,
 * so the cell shrinks as you magnify and the points redistribute.
 */
const AtlasMap = (() => {
  'use strict';

  const X_SCALE = 0.8487;
  const Y_SCALE = 1.3523;
  const X_MAX = X_SCALE * Math.PI;
  const Y_MAX = Y_SCALE;

  const ROBINSON_X = [
    1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216,
    0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322,
  ];
  const ROBINSON_Y = [
    0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958,
    0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000,
  ];

  /* Linear interpolation between the 5-degree table entries, as in the figures. */
  function project(lon, lat) {
    const clamped = Math.max(-90, Math.min(90, lat));
    const absLat = Math.abs(clamped);
    const slot = Math.min(Math.floor(absLat / 5), 17);
    const t = (absLat - slot * 5) / 5;
    const xp = ROBINSON_X[slot] + t * (ROBINSON_X[slot + 1] - ROBINSON_X[slot]);
    const yp = ROBINSON_Y[slot] + t * (ROBINSON_Y[slot + 1] - ROBINSON_Y[slot]);
    return [
      X_SCALE * (lon * Math.PI / 180) * xp,
      clamped === 0 ? 0 : Y_SCALE * yp * Math.sign(clamped),
    ];
  }

  const PAD = 6;
  const MIN_ZOOM = 1;
  const MAX_ZOOM = 16;

  /* Cell size in degrees for a zoom level. Quantised so that panning around at
   * one magnification does not re-bin on every frame. */
  function resolutionForZoom(zoom) {
    if (zoom < 1.6) return 2;
    if (zoom < 3) return 1;
    if (zoom < 6) return 0.5;
    if (zoom < 11) return 0.25;
    return 0.1;
  }

  function style(name, fallback) {
    const value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return value || fallback;
  }

  class Map2D {
    constructor(canvas, tip) {
      this.canvas = canvas;
      this.tip = tip;
      this.context = canvas.getContext('2d');
      this.land = [];
      this.bins = [];
      this.describe = null;
      this.onSelect = null;
      this.scaleByCount = true;
      this.resolution = 2;
      this.zoom = 1;
      this.panX = 0;
      this.panY = 0;
      this.selected = null;
      this.drag = null;
      this.frame = null;

      canvas.addEventListener('mousemove', (e) => this._onMove(e));
      canvas.addEventListener('mouseleave', () => this._onLeave());
      canvas.addEventListener('mousedown', (e) => this._onDown(e));
      canvas.addEventListener('wheel', (e) => this._onWheel(e), { passive: false });
      canvas.addEventListener('dblclick', (e) => this._onDoubleClick(e));
      window.addEventListener('mouseup', (e) => this._onUp(e));
      window.addEventListener('mousemove', (e) => this._onDragMove(e));

      if (typeof ResizeObserver === 'function') {
        this.observer = new ResizeObserver(() => this.redraw());
        this.observer.observe(canvas.parentElement || canvas);
      } else if (typeof window.addEventListener === 'function') {
        window.addEventListener('resize', () => this.redraw());
      }
    }

    setLand(rings) { this.land = rings || []; }

    /* points: { lat, lon, rows } parallel arrays plus the source row index. */
    setPoints(points, options = {}) {
      this.points = points;
      if (options.scaleByCount !== undefined) this.scaleByCount = options.scaleByCount;
      if (options.describe !== undefined) this.describe = options.describe;
      if (options.onSelect !== undefined) this.onSelect = options.onSelect;
      this.selected = null;
      this._bin();
      this.redraw();
    }

    setScaleByCount(on) {
      this.scaleByCount = on;
      this.redraw();
    }

    clearSelection() {
      this.selected = null;
      this.redraw();
    }

    resetView() {
      this.zoom = 1;
      this.panX = 0;
      this.panY = 0;
      this._reBinIfNeeded();
      this.redraw();
    }

    zoomBy(factor) {
      const rect = this.canvas.getBoundingClientRect();
      this._zoomAt(rect.width / 2, rect.height / 2, factor);
    }

    _zoomAt(mx, my, factor) {
      const before = this.zoom;
      const after = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, before * factor));
      if (after === before) return;
      const { scale, cx, cy } = this.transform || {};
      if (scale) {
        // Keep the geography under the pointer where it is.
        const worldX = (mx - cx - this.panX) / (scale * before);
        const worldY = (cy + this.panY - my) / (scale * before);
        this.panX = mx - cx - worldX * scale * after;
        this.panY = my - cy + worldY * scale * after;
      }
      this.zoom = after;
      if (after === MIN_ZOOM) { this.panX = 0; this.panY = 0; }
      this._reBinIfNeeded();
      this.redraw();
    }

    _reBinIfNeeded() {
      const wanted = resolutionForZoom(this.zoom);
      if (wanted !== this.resolution) {
        this.resolution = wanted;
        const keep = this.selected;
        this._bin();
        // The old cell no longer exists at the new grid size; drop the
        // highlight rather than leave it pointing at nothing.
        if (keep) this.selected = null;
      }
    }

    _clampPan() {
      const { scale } = this.transform || {};
      if (!scale) return;
      const rect = this.canvas.getBoundingClientRect();
      const halfWidth = X_MAX * scale * this.zoom;
      const halfHeight = Y_MAX * scale * this.zoom;
      const slackX = Math.max(0, halfWidth - rect.width / 2);
      const slackY = Math.max(0, halfHeight - rect.height / 2);
      this.panX = Math.max(-slackX, Math.min(slackX, this.panX));
      this.panY = Math.max(-slackY, Math.min(slackY, this.panY));
    }

    _bin() {
      const points = this.points;
      this.bins = [];
      this.plotted = 0;
      this.offMap = 0;
      if (!points) return;

      const step = this.resolution;
      const columns = Math.ceil(360 / step) + 2;
      const buckets = new Map();
      const { lat, lon, rows } = points;
      for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const y = lat[row];
        const x = lon[row];
        if (Number.isNaN(x) || Number.isNaN(y)) { this.offMap++; continue; }
        // Snap to bin centre so a cell's marker never drifts with its members.
        const by = Math.floor((y + 90) / step);
        const bx = Math.floor((x + 180) / step);
        const key = by * columns + bx;
        let bin = buckets.get(key);
        if (!bin) {
          bin = { lat: (by + 0.5) * step - 90, lon: (bx + 0.5) * step - 180, count: 0, rows: [] };
          buckets.set(key, bin);
        }
        bin.count++;
        bin.rows.push(row);
        this.plotted++;
      }
      this.bins = [...buckets.values()].sort((a, b) => a.count - b.count);
      this.maxCount = this.bins.reduce((m, b) => Math.max(m, b.count), 1);
    }

    redraw() {
      const canvas = this.canvas;
      const parent = canvas.parentElement;
      const width = Math.max(320, (parent ? parent.clientWidth : canvas.clientWidth) || 640);
      const height = canvas.clientHeight || 420;
      const ratio = window.devicePixelRatio || 1;

      canvas.width = Math.round(width * ratio);
      canvas.height = Math.round(height * ratio);
      const context = this.context;
      context.setTransform(ratio, 0, 0, ratio, 0, 0);
      context.clearRect(0, 0, width, height);

      const scale = Math.min((width - 2 * PAD) / (2 * X_MAX), (height - 2 * PAD) / (2 * Y_MAX));
      this.transform = { scale, cx: width / 2, cy: height / 2 };
      this._clampPan();

      const zoom = this.zoom;
      const cx = width / 2 + this.panX;
      const cy = height / 2 + this.panY;
      const toScreen = (lon, lat) => {
        const [x, y] = project(lon, lat);
        return [cx + x * scale * zoom, cy - y * scale * zoom];
      };

      context.fillStyle = style('--ocean', '#ddecef');
      context.beginPath();
      for (let lat = -90; lat <= 90; lat += 2) { const [x, y] = toScreen(180, lat); lat === -90 ? context.moveTo(x, y) : context.lineTo(x, y); }
      for (let lat = 90; lat >= -90; lat -= 2) { const [x, y] = toScreen(-180, lat); context.lineTo(x, y); }
      context.closePath();
      context.fill();
      context.save();
      context.clip();

      context.fillStyle = style('--land', '#f3efe6');
      context.strokeStyle = style('--coast', '#a89e8c');
      context.lineWidth = 0.6;
      for (const ring of this.land) {
        context.beginPath();
        for (let i = 0; i < ring.length; i += 2) {
          const [x, y] = toScreen(ring[i], ring[i + 1]);
          i === 0 ? context.moveTo(x, y) : context.lineTo(x, y);
        }
        context.closePath();
        context.fill();
        context.stroke();
      }

      // Graticule spacing follows the zoom, so it stays a guide and never a grid.
      const gridStep = zoom < 2 ? 30 : zoom < 5 ? 15 : zoom < 10 ? 5 : 2;
      context.strokeStyle = style('--coast', '#a89e8c');
      context.globalAlpha = 0.28;
      context.lineWidth = 0.5;
      for (let lon = -180 + gridStep; lon < 180; lon += gridStep) {
        context.beginPath();
        for (let lat = -90; lat <= 90; lat += 2) {
          const [x, y] = toScreen(lon, lat);
          lat === -90 ? context.moveTo(x, y) : context.lineTo(x, y);
        }
        context.stroke();
      }
      for (let lat = -90 + gridStep; lat < 90; lat += gridStep) {
        context.beginPath();
        for (let lon = -180; lon <= 180; lon += 5) {
          const [x, y] = toScreen(lon, lat);
          lon === -180 ? context.moveTo(x, y) : context.lineTo(x, y);
        }
        context.stroke();
      }
      context.globalAlpha = 1;

      const fill = style('--blue-4', '#1c5cab');
      const edge = style('--blue-5', '#0d366b');
      const highlight = style('--accent', '#c2622d');
      const maxRadius = Math.max(5, Math.min(14, scale * 0.055)) * Math.min(1 + (zoom - 1) * 0.35, 2.6);
      for (const bin of this.bins) {
        const [x, y] = toScreen(bin.lon, bin.lat);
        bin.x = x;
        bin.y = y;
        bin.r = this.scaleByCount
          ? 1.7 + (maxRadius - 1.7) * Math.sqrt(bin.count / this.maxCount)
          : 2.6 * Math.min(1 + (zoom - 1) * 0.3, 2.2);
        const picked = bin === this.selected;
        context.beginPath();
        context.arc(x, y, bin.r, 0, Math.PI * 2);
        context.globalAlpha = picked ? 0.95 : 0.62;
        context.fillStyle = picked ? highlight : fill;
        context.fill();
        context.globalAlpha = 0.9;
        context.lineWidth = picked ? 1.8 : 0.55;
        context.strokeStyle = picked ? highlight : edge;
        context.stroke();
        if (picked) {
          context.beginPath();
          context.arc(x, y, bin.r + 4.5, 0, Math.PI * 2);
          context.globalAlpha = 0.85;
          context.lineWidth = 1.2;
          context.stroke();
        }
      }
      context.globalAlpha = 1;
      context.restore();
    }

    _scheduleRedraw() {
      if (this.frame) return;
      this.frame = requestAnimationFrame(() => { this.frame = null; this.redraw(); });
    }

    _at(event) {
      const rect = this.canvas.getBoundingClientRect();
      return [event.clientX - rect.left, event.clientY - rect.top, rect];
    }

    _hit(mx, my) {
      let best = null;
      let bestDistance = Infinity;
      for (const bin of this.bins) {
        if (bin.x === undefined) continue;
        const d = Math.hypot(bin.x - mx, bin.y - my);
        if (d < Math.max(bin.r + 4, 7) && d < bestDistance) { best = bin; bestDistance = d; }
      }
      return best;
    }

    _onWheel(event) {
      event.preventDefault();
      const [mx, my] = this._at(event);
      // Trackpads send many small deltas; the exponent keeps them proportional.
      this._zoomAt(mx, my, Math.pow(1.0015, -event.deltaY));
    }

    _onDoubleClick(event) {
      event.preventDefault();
      const [mx, my] = this._at(event);
      this._zoomAt(mx, my, 1.9);
    }

    _onDown(event) {
      if (event.button !== 0) return;
      const [mx, my] = this._at(event);
      this.drag = { startX: event.clientX, startY: event.clientY, panX: this.panX, panY: this.panY, moved: false, mx, my };
      this.canvas.style.cursor = 'grabbing';
    }

    _onDragMove(event) {
      if (!this.drag) return;
      const dx = event.clientX - this.drag.startX;
      const dy = event.clientY - this.drag.startY;
      if (!this.drag.moved && Math.hypot(dx, dy) < 4) return;   // still a click
      this.drag.moved = true;
      this.panX = this.drag.panX + dx;
      this.panY = this.drag.panY + dy;
      this._onLeave();
      this._scheduleRedraw();
    }

    _onUp(event) {
      if (!this.drag) return;
      const drag = this.drag;
      this.drag = null;
      this.canvas.style.cursor = '';
      if (drag.moved) return;
      const [mx, my] = this._at(event);
      const bin = this._hit(mx, my);
      this.selected = bin || null;
      this.redraw();
      if (this.onSelect) this.onSelect(bin);
    }

    _onMove(event) {
      if (this.drag) return;
      if (!this.bins.length) return;
      const [mx, my, rect] = this._at(event);
      const best = this._hit(mx, my);
      this.canvas.style.cursor = best ? 'pointer' : 'grab';
      if (!best) return this._onLeave();

      const lines = [
        `<b>${best.count.toLocaleString()}</b> ${best.count === 1 ? 'record' : 'records'}`,
        `${Math.abs(best.lat).toFixed(1)}°${best.lat < 0 ? 'S' : 'N'}, ` +
        `${Math.abs(best.lon).toFixed(1)}°${best.lon < 0 ? 'W' : 'E'}`,
      ];
      if (this.describe) {
        const extra = this.describe(best.rows.slice(0, 6), best.count);
        if (extra) lines.push(extra);
      }
      lines.push('<i>Click for details</i>');
      this.tip.innerHTML = lines.join('<br>');
      this.tip.hidden = false;
      const tipRect = this.tip.getBoundingClientRect();
      const left = Math.min(Math.max(4, mx + 12), rect.width - tipRect.width - 4);
      const top = Math.max(4, my - tipRect.height - 10);
      this.tip.style.left = `${left}px`;
      this.tip.style.top = `${top}px`;
    }

    _onLeave() { this.tip.hidden = true; }
  }

  return { Map2D, project, resolutionForZoom };
})();
