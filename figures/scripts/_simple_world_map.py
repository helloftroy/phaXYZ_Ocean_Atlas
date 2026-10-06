"""Lightweight Robinson world-map drawing without Cartopy.

The local laptop environment used for figure polishing does not always have
compiled GIS dependencies available. This helper keeps the paper figures
regenerable from a small Natural Earth shapefile plus pyshp.
"""
import math
from functools import lru_cache
from pathlib import Path

import matplotlib.patches as mpatches
import numpy as np
import shapefile

OUT = Path(__file__).resolve().parent.parent
LAND_SHP = OUT / "assets" / "naturalearth" / "ne_110m_land.shp"

OCEAN = "#DDECEF"
LAND = "#F3EFE6"
COAST = "#A89E8C"
GRID = "#78909A"

# Robinson projection lookup table, every 5 degrees from equator to pole.
_ROBIN_X = np.array([
    1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
    0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
    0.7186, 0.6732, 0.6213, 0.5722, 0.5322,
])
_ROBIN_Y = np.array([
    0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
    0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
    0.8435, 0.8936, 0.9394, 0.9761, 1.0000,
])
_X_SCALE = 0.8487
_Y_SCALE = 1.3523
XMAX = _X_SCALE * math.pi
YMAX = _Y_SCALE


def project(lon, lat):
    """Project lon/lat degrees to Robinson x/y."""
    lon_arr = np.asarray(lon, dtype=float)
    lat_arr = np.asarray(lat, dtype=float)
    clipped = np.clip(lat_arr, -90, 90)
    abs_lat = np.abs(clipped)
    xp = np.interp(abs_lat, np.arange(0, 91, 5), _ROBIN_X)
    yp = np.interp(abs_lat, np.arange(0, 91, 5), _ROBIN_Y)
    x = _X_SCALE * np.radians(lon_arr) * xp
    y = _Y_SCALE * yp * np.sign(clipped)
    return x, y


@lru_cache(maxsize=1)
def land_parts():
    reader = shapefile.Reader(str(LAND_SHP))
    parts = []
    for shape in reader.shapes():
        idxs = list(shape.parts) + [len(shape.points)]
        for start, end in zip(idxs[:-1], idxs[1:]):
            pts = shape.points[start:end]
            if len(pts) < 3:
                continue
            current = [pts[0]]
            for prev, point in zip(pts, pts[1:]):
                if abs(point[0] - prev[0]) > 180:
                    if len(current) >= 3:
                        parts.append(current)
                    current = [point]
                else:
                    current.append(point)
            if len(current) >= 3:
                parts.append(current)
    return parts


@lru_cache(maxsize=1)
def boundary_xy():
    """The true Robinson outline: the projected 180 deg meridians.

    This used to be approximated by an ellipse, which is wrong in a way that
    silently loses data rather than merely looking off. Robinson's edge is
    flatter than an ellipse at mid latitudes, so the ellipse cuts inside it
    there -- at 35 deg S the ellipse half-width is 2.402 against the
    projection's own 2.502. Any point beyond about 172 deg E/W at those
    latitudes was clipped away with no warning: it sat inside the axis limits,
    inside the drawn ocean, and simply did not render. Found when a cluster
    with 11 of 12 genomes at one vent at 179.1 deg E, 34.9 deg S vanished from
    the map (figures/scripts/plot_phac_cluster_biogeography.py).
    """
    lats = np.linspace(-90, 90, 361)
    x_east, y_east = project(np.full_like(lats, 180.0), lats)
    x_west, y_west = project(np.full_like(lats, -180.0), lats[::-1])
    return np.column_stack([np.concatenate([x_east, x_west]),
                            np.concatenate([y_east, y_west])])


def setup_ax(ax):
    ax.set_aspect("equal")
    ax.set_xlim(-XMAX * 1.03, XMAX * 1.03)
    ax.set_ylim(-YMAX * 1.04, YMAX * 1.04)
    ax.axis("off")
    ocean = mpatches.Polygon(boundary_xy(), closed=True,
                             facecolor=OCEAN, edgecolor="#8FA1A5",
                             linewidth=0.7, zorder=0)
    ax.add_patch(ocean)
    for lon in range(-120, 181, 60):
        lats = np.linspace(-88, 88, 260)
        lons = np.full_like(lats, lon)
        x, y = project(lons, lats)
        ax.plot(x, y, color=GRID, alpha=0.24, linewidth=0.45, zorder=1, clip_path=ocean)
    for lat in range(-60, 61, 30):
        lons = np.linspace(-180, 180, 361)
        lats = np.full_like(lons, lat)
        x, y = project(lons, lats)
        ax.plot(x, y, color=GRID, alpha=0.24, linewidth=0.45, zorder=1, clip_path=ocean)
    return ocean


def draw_land(ax, clip_path):
    for pts in land_parts():
        lons = [p[0] for p in pts]
        lats = [p[1] for p in pts]
        x, y = project(lons, lats)
        patch = mpatches.Polygon(np.column_stack([x, y]), closed=True,
                                 facecolor=LAND, edgecolor=COAST,
                                 linewidth=0.28, zorder=2)
        patch.set_clip_path(clip_path)
        ax.add_patch(patch)


def scatter(ax, lon, lat, clip_path=None, **kwargs):
    x, y = project(lon, lat)
    artist = ax.scatter(x, y, **kwargs)
    if clip_path is not None:
        artist.set_clip_path(clip_path)
    return artist
