"""Convert the Natural Earth 110 m land shapefile into a compact outline for the
web map.

pyshp is not installed in the environment used for these builds, and the
shapefile polygon format is simple enough to read directly: a 100-byte header,
then records of (number, length, shape type, bounding box, part offsets,
points). Reading it here keeps the site self-contained -- no tile server, no
runtime GIS dependency, and no third-party map service to depend on or overload.

Rings are emitted as flat [lon, lat, lon, lat, ...] arrays rounded to two
decimals, with consecutive duplicate points dropped. At 110 m source resolution
two decimals is below the data's own precision, and flat arrays avoid one JSON
sub-array per coordinate.

Natural Earth is public domain.

Usage:
    python site_build/build_land_outline.py
"""
import json
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SHP = ROOT / 'figures/assets/naturalearth/ne_110m_land.shp'
OUT = ROOT / 'docs/data/land.json'

POLYGON = 5
MIN_RING_POINTS = 4

raw = SHP.read_bytes()
magic, = struct.unpack_from('>i', raw, 0)
if magic != 9994:
    raise SystemExit(f'{SHP} is not a shapefile (magic {magic})')
file_length = struct.unpack_from('>i', raw, 24)[0] * 2
shape_type = struct.unpack_from('<i', raw, 32)[0]
print(f'shape type {shape_type}, {file_length:,} bytes')

rings, offset, skipped = [], 100, 0
while offset < file_length:
    _number, content_words = struct.unpack_from('>ii', raw, offset)
    content = offset + 8
    kind, = struct.unpack_from('<i', raw, content)
    if kind != POLYGON:
        offset = content + content_words * 2
        continue
    n_parts, n_points = struct.unpack_from('<ii', raw, content + 36)
    parts = struct.unpack_from(f'<{n_parts}i', raw, content + 44)
    coords_at = content + 44 + n_parts * 4
    points = struct.unpack_from(f'<{n_points * 2}d', raw, coords_at)
    bounds = list(parts) + [n_points]
    for start, end in zip(bounds[:-1], bounds[1:]):
        flat, last = [], None
        for i in range(start, end):
            lon = round(points[i * 2], 2)
            lat = round(points[i * 2 + 1], 2)
            if (lon, lat) == last:
                continue
            flat += [lon, lat]
            last = (lon, lat)
        if len(flat) >= MIN_RING_POINTS * 2:
            rings.append(flat)
        else:
            skipped += 1
    offset = content + content_words * 2

rings.sort(key=len, reverse=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(rings, separators=(',', ':')))
total = sum(len(r) // 2 for r in rings)
print(f'{len(rings)} rings, {total:,} points, {skipped} too small to draw')
print(f'wrote {OUT} ({OUT.stat().st_size / 1e3:.0f} kB)')
