"""Render small structure cartoons with PyMOL, if PyMOL is available.

This script finds ColabFold PDB/CIF outputs for the 12 regional-specialist
PhaC proteins, writes one PyMOL script per structure, renders transparent PNGs,
and makes a simple contact sheet for figure drafting.
"""
import csv
import math
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent.parent
BASE = ROOT / "figures" / "structures" / "regional_specialists"
META = BASE / "phaC_regional_specialists_structure_metadata.tsv"
PRED = BASE / "colabfold_single_sequence"
RENDER = BASE / "cartoons"

COLORS = [
    "teal", "orange", "purple", "forest", "yelloworange", "raspberry",
    "marine", "brightorange", "brown", "palegreen", "tv_red", "gray60",
]


def find_model(structure_id):
    patterns = [
        f"{structure_id}*rank_001*.pdb",
        f"{structure_id}*rank_001*.cif",
        f"{structure_id}*.pdb",
        f"{structure_id}*.cif",
    ]
    for pattern in patterns:
        hits = sorted(PRED.glob(pattern))
        if hits:
            return hits[0]
    return None


def render_one(pymol, row, color):
    sid = row["structure_id"]
    model = find_model(sid)
    if model is None:
        print(f"missing model for {sid}")
        return None

    png = RENDER / f"{sid}.png"
    pml = RENDER / f"{sid}.pml"
    pml.write_text(f"""
load {model.resolve()}, model
hide everything
show cartoon, model
color {color}, model
bg_color white
set ray_opaque_background, off
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
orient model
zoom model, 3
ray 1200, 900
png {png.resolve()}, dpi=300
quit
""".strip() + "\n")
    subprocess.run([pymol, "-cq", str(pml)], check=True)
    return png


def contact_sheet(rows, pngs):
    valid = [(r, p) for r, p in zip(rows, pngs) if p and p.exists()]
    if not valid:
        return
    ncols = 4
    nrows = math.ceil(len(valid) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 2.7 * nrows), dpi=300)
    axes = list(axes.flat if hasattr(axes, "flat") else [axes])
    for ax in axes:
        ax.axis("off")
    for ax, (row, png) in zip(axes, valid):
        ax.imshow(mpimg.imread(png))
        ax.set_title(
            f"{row['dominant_genus']}\\n{row['structure_id'][:12]}",
            fontsize=8,
            fontstyle="italic",
        )
        ax.axis("off")
    out = RENDER / "phaC_regional_specialists_structure_contact_sheet.png"
    fig.tight_layout(pad=0.4)
    fig.savefig(out, transparent=False, facecolor="white")
    print(f"wrote {out}")


def main():
    pymol = shutil.which("pymol") or shutil.which("pymol-open-source")
    if pymol is None:
        raise SystemExit("PyMOL not found. Load a pymol module, install pymol-open-source, or render locally.")

    RENDER.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(META.open(), delimiter="\t"))
    pngs = [render_one(pymol, row, COLORS[i % len(COLORS)]) for i, row in enumerate(rows)]
    contact_sheet(rows, pngs)


if __name__ == "__main__":
    main()
