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
      Y_SCALE * yp * Math.sign(clamped || 1) * (clamped === 0 ? 0 : 1),
    ];
  }

  const PAD = 6;

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
      this.scaleByCount = true;
      this.resolution = 2;

      this._onMove = this._onMove.bind(this);
      this._onLeave = this._onLeave.bind(this);
      canvas.addEventListener('mousemove', this._onMove);
      canvas.addEventListener('mouseleave', this._onLeave);

      if (typeof ResizeObserver === 'function') {
        this.observer = new ResizeObserver(() => this.redraw());
        this.observer.observe(canvas.parentElement || canvas);
      } else {
        window.addEventListener('resize', () => this.redraw());
      }
    }

    setLand(rings) { this.land = rings || []; }

    /* points: { lat, lon, rows } parallel arrays plus the source row index. */
    setPoints(points, options = {}) {
      this.points = points;
      if (options.scaleByCount !== undefined) this.scaleByCount = options.scaleByCount;
      if (options.describe !== undefined) this.describe = options.describe;
      this._bin();
      this.redraw();
    }

    setScaleByCount(on) {
      this.scaleByCount = on;
      this.redraw();
    }

    _bin() {
      const points = this.points;
      this.bins = [];
      this.plotted = 0;
      this.offMap = 0;
      if (!points) return;

      const step = this.resolution;
      const buckets = new Map();
      const { lat, lon, rows } = points;
      for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const y = lat[row];
        const x = lon[row];
        if (Number.isNaN(x) || Number.isNaN(y)) { this.offMap++; continue; }
        // Snap to bin centre so a bin's marker never drifts with its members.
        const by = Math.floor((y + 90) / step);
        const bx = Math.floor((x + 180) / step);
        const key = by * 100000 + bx;
        let bin = buckets.get(key);
        if (!bin) {
          bin = { lat: (by + 0.5) * step - 90, lon: (bx + 0.5) * step - 180, count: 0, rows: [] };
          buckets.set(key, bin);
        }
        bin.count++;
        if (bin.rows.length < 6) bin.rows.push(row);
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
      const cx = width / 2;
      const cy = height / 2;
      this.transform = { scale, cx, cy };

      const toScreen = (lon, lat) => {
        const [x, y] = project(lon, lat);
        return [cx + x * scale, cy - y * scale];
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

      context.strokeStyle = style('--coast', '#a89e8c');
      context.globalAlpha = 0.28;
      context.lineWidth = 0.5;
      for (let lon = -150; lon <= 150; lon += 30) {
        context.beginPath();
        for (let lat = -90; lat <= 90; lat += 2) {
          const [x, y] = toScreen(lon, lat);
          lat === -90 ? context.moveTo(x, y) : context.lineTo(x, y);
        }
        context.stroke();
      }
      for (let lat = -60; lat <= 60; lat += 30) {
        context.beginPath();
        for (let lon = -180; lon <= 180; lon += 5) {
          const [x, y] = toScreen(lon, lat);
          lon === -180 ? context.moveTo(x, y) : context.lineTo(x, y);
        }
        context.stroke();
      }
      context.globalAlpha = 1;
      context.restore();

      const fill = style('--blue-4', '#1c5cab');
      const edge = style('--blue-5', '#0d366b');
      const maxRadius = Math.max(5, Math.min(13, scale * 0.055));
      for (const bin of this.bins) {
        const [x, y] = toScreen(bin.lon, bin.lat);
        bin.x = x;
        bin.y = y;
        bin.r = this.scaleByCount
          ? 1.7 + (maxRadius - 1.7) * Math.sqrt(bin.count / this.maxCount)
          : 2.6;
        context.beginPath();
        context.arc(x, y, bin.r, 0, Math.PI * 2);
        context.globalAlpha = 0.62;
        context.fillStyle = fill;
        context.fill();
        context.globalAlpha = 0.9;
        context.lineWidth = 0.55;
        context.strokeStyle = edge;
        context.stroke();
      }
      context.globalAlpha = 1;
    }

    _onMove(event) {
      if (!this.bins.length) return;
      const rect = this.canvas.getBoundingClientRect();
      const mx = event.clientX - rect.left;
      const my = event.clientY - rect.top;
      let best = null;
      let bestDistance = Infinity;
      for (const bin of this.bins) {
        if (bin.x === undefined) continue;
        const d = Math.hypot(bin.x - mx, bin.y - my);
        if (d < Math.max(bin.r + 4, 7) && d < bestDistance) { best = bin; bestDistance = d; }
      }
      if (!best) return this._onLeave();

      const lines = [
        `<b>${best.count.toLocaleString()}</b> ${best.count === 1 ? 'record' : 'records'}`,
        `${Math.abs(best.lat).toFixed(1)}°${best.lat < 0 ? 'S' : 'N'}, ` +
        `${Math.abs(best.lon).toFixed(1)}°${best.lon < 0 ? 'W' : 'E'}`,
      ];
      if (this.describe) {
        const extra = this.describe(best.rows, best.count);
        if (extra) lines.push(extra);
      }
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

  return { Map2D, project };
})();
