"""Lay out the PyMOL renders of the novel-phaC superposition as one figure.

PyMOL ray-traces the structures (render_novel_phac_superposition.pml, run that
first); this places the three views, sets the type in the same face as the rest
of the manuscript figures, and writes the measured numbers underneath. Keeping
the labelling here rather than in PyMOL is deliberate -- PyMOL's label fonts do
not match the other panels and cannot be positioned as precisely.

Each render is cropped to its own alpha bounding box before placement, because
PyMOL pads a ray trace to the requested canvas and the three views fill
different fractions of it; placing them uncropped would have made the same
molecule appear at three different scales.

Usage:
    python structure_prediction/compose_superposition_figure.py

Input:  figures/_superposition_raw/view{1,2,3}.png
Output: figures/novel_phac_superposition.png / .pdf
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / 'figures/_superposition_raw'

CANDIDATE_COLOR, REFERENCE_COLOR, TRIAD_COLOR = '#1C5CAB', '#C2622D', '#0D9488'
TEXT_DARK, TEXT_MUTED = '#1E2630', '#5B6670'

# Measured by render_novel_phac_superposition.pml and
# superpose_novel_phac_example.py; see PHA_CLEAN_RESULTS.md section 13.4.
CEALIGN_RMSD, CEALIGN_N = 4.14, 328
KABSCH_TM, KABSCH_N, KABSCH_RMSD, KABSCH_CLOSE = 0.860, 339, 2.32, 309
ALIGNMENT_IDENTITY = 15.9
TRIAD_DISTANCES = [('Cys', 0.49), ('Asp', 1.67), ('His', 1.23)]
CYS_SG_DISTANCE = 0.52

PANELS = [
    ('view1.png', 'A', 'Superposed, full length'),
    ('view2.png', 'B', 'Rotated 90°'),
    ('view3.png', 'C', 'Catalytic triad'),
]


def cropped(path):
    """Trim a ray trace to its own content, so the three views share a scale."""
    image = Image.open(path).convert('RGBA')
    alpha = np.asarray(image)[:, :, 3]
    rows, cols = np.where(alpha > 8)
    if not len(rows):
        return image
    image = image.crop((cols.min(), rows.min(), cols.max() + 1, rows.max() + 1))
    flat = Image.new('RGBA', image.size, 'white')
    flat.alpha_composite(image)
    return flat.convert('RGB')


images = [(cropped(RAW / name), letter, caption) for name, letter, caption in PANELS]
for (image, letter, _) in images:
    print(f'panel {letter}: {image.size[0]} x {image.size[1]} px after cropping')

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig = plt.figure(figsize=(11.0, 4.9), dpi=300)
gs = fig.add_gridspec(1, 3, wspace=0.03, left=0.012, right=0.988, top=0.955, bottom=0.215)

for column, (image, letter, caption) in enumerate(images):
    ax = fig.add_subplot(gs[0, column])
    ax.imshow(image)
    ax.set_axis_off()
    ax.text(0.0, 1.0, letter, transform=ax.transAxes, ha='left', va='top',
            fontsize=14, fontweight='bold', color=TEXT_DARK)
    ax.text(0.5, -0.035, caption, transform=ax.transAxes, ha='center', va='top',
            fontsize=10, color=TEXT_DARK)

handles = [
    plt.Line2D([0], [0], color=CANDIDATE_COLOR, linewidth=4,
               label='Candidate, 466 aa — $\\it{Haliangium}$ (Myxococcota)'),
    plt.Line2D([0], [0], color=REFERENCE_COLOR, linewidth=4,
               label='A0A1H8CR08, 340 aa — $\\it{Mesobacillus\\ persicus}$ (Bacillota)'),
    plt.Line2D([0], [0], color=TRIAD_COLOR, linewidth=4,
               label='Catalytic triad, both structures'),
]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8.9,
           bbox_to_anchor=(0.5, 0.085), handlelength=1.5, columnspacing=1.5, handletextpad=0.6)

triad_text = ', '.join(f'{name} {d:.2f} Å' for name, d in TRIAD_DISTANCES)
fig.text(0.5, 0.062,
         f'TM-score {KABSCH_TM:.3f} over {KABSCH_N} aligned residues  ·  '
         f'{ALIGNMENT_IDENTITY:.1f}% identity within the alignment  ·  '
         f'backbone RMSD {KABSCH_RMSD:.2f} Å over the {KABSCH_CLOSE} pairs within 5 Å '
         f'(PyMOL cealign: {CEALIGN_RMSD:.2f} Å over {CEALIGN_N})',
         ha='center', va='top', fontsize=9.3, color=TEXT_MUTED)
fig.text(0.5, 0.030,
         f'Triad Cα separation after superposition — {triad_text}; catalytic Cys Sγ–Sγ {CYS_SG_DISTANCE:.2f} Å',
         ha='center', va='top', fontsize=9.3, color=TEXT_MUTED)

fig.savefig(ROOT / 'figures/novel_phac_superposition.png', dpi=300, facecolor='white')
fig.savefig(ROOT / 'figures/novel_phac_superposition.pdf', facecolor='white')
print('saved', ROOT / 'figures/novel_phac_superposition.png')
