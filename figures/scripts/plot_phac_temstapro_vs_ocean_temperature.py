"""Plot TemStaPro predicted phaC stability class against WOA ocean temperature.

QC: removes target_ids whose best phaC reference query is one of the
confirmed mislabeled/non-phaC entries listed in _phac_qc.py.

Length QC: also removes proteins shorter than MIN_LENGTH aa. Checked
directly (2026-09): across all 128,000 TemStaPro-scored phaC proteins,
"hyperthermophilic"-labeled ones have a median length of 117aa vs 386aa
for the dataset overall, and 86.6% of all hyperthermophilic calls are
<150aa -- short fragments (common at MAG contig breaks) produce noisy,
non-monotonic per-threshold probabilities (confidence dips below 0.5 at
mid-range thresholds then jumps back up at the highest ones -- a pattern
a genuinely stable full-length protein shouldn't show) that TemStaPro's
own "clash" flag only partly catches. Applied to the WHOLE plotted
population here (not just hyperthermophilic calls) to avoid selection
bias -- a short-fragment mesophilic call is exactly as unreliable as a
short-fragment hyperthermophilic one, just less visually dramatic.
300aa is a safe floor even for the shorter class III/IV PhaC (~320-370aa
typical) vs class I/II (~560-600aa). This filter is deliberately scoped
to TemStaPro-derived stability analysis only -- presence/absence and
architecture calls elsewhere in this project rely on homology, which
remains valid on partial sequences, and are NOT length-filtered.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np


ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "figures" / "scripts"))
import _phac_qc

IN = ROOT / "temstapro" / "phaC_temstapro_predictions_with_metadata.tsv"
PHAC_META_CANDIDATES = [
    ROOT / "PHA_bioprospecting" / "omdb_search" / "results" / "phaC_unique_targets_with_metadata_depth.tsv",
    ROOT.parent / "fair_ocean_agent" / "phaC_unique_targets_with_metadata_depth.tsv",
]
OUT = ROOT / "figures"

LABEL_TO_TEMP = {
    "<40": 37.5,
    "[40-45)": 42.5,
    "[45-50)": 47.5,
    "[50-55)": 52.5,
    "[55-60)": 57.5,
    "[60-65)": 62.5,
    "[65-70)": 67.5,
    "[70-75)": 72.5,
    "[75-80)": 77.5,
    "80<=": 82.5,
}
TEMP_ORDER = [37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5]
MIN_LENGTH = 300
THERM_COLORS = {
    "mesophilic": "#638A8C",
    "thermophilic": "#C5653B",
    "hyperthermophilic": "#8B3F67",
    "undetermined": "#A8A8A0",
}


def parse_float(raw: str):
    if raw is None or raw == "":
        return None
    try:
        val = float(raw)
    except ValueError:
        return None
    if math.isnan(val):
        return None
    return val


def median(vals):
    vals = sorted(vals)
    if not vals:
        return math.nan
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def load_rows():
    phac_meta = next((p for p in PHAC_META_CANDIDATES if p.exists()), None)
    if phac_meta is None:
        raise SystemExit(
            "Could not find phaC_unique_targets_with_metadata_depth.tsv for QC filtering. "
            "Checked: " + ", ".join(str(p) for p in PHAC_META_CANDIDATES)
        )
    bad_targets = _phac_qc.load_bad_targets(phac_meta)
    rows = []
    counts = Counter()
    with IN.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            counts["all_temstapro_rows"] += 1
            if row.get("target_id") in bad_targets:
                counts["qc_removed_bad_target"] += 1
                continue
            length = parse_float(row.get("clean_length", ""))
            if length is None or length < MIN_LENGTH:
                counts["qc_removed_short_fragment"] += 1
                continue
            ocean_temp = parse_float(row.get("woa23_temp_annual_degC", ""))
            if ocean_temp is None:
                counts["missing_woa_temperature"] += 1
                continue
            label = row.get("right_hand_label", "")
            pred_temp = LABEL_TO_TEMP.get(label)
            if pred_temp is None:
                counts["missing_temstapro_label"] += 1
                continue
            rows.append({
                "target_id": row.get("target_id", ""),
                "ocean_temp": ocean_temp,
                "pred_temp": pred_temp,
                "label": label,
                "thermophilicity": row.get("thermophilicity", "undetermined") or "undetermined",
                "genus": row.get("gtdb_genus", "") or "Unknown",
                "species": row.get("gtdb_species", "") or "Unknown",
                "family": row.get("gtdb_family", "") or "Unknown",
                "study_id": row.get("study_id", "") or "Unknown",
                "genome": row.get("genome", ""),
                "phylum": row.get("gtdb_phylum", "") or "Unknown",
                "depth": parse_float(row.get("depth_m", "")),
                "depth_zone": row.get("depth_zone", "") or "Unknown",
                "length": length,
            })
    counts["plotted_rows"] = len(rows)
    return rows, counts


def binned_medians(rows, bin_width=5):
    bins = defaultdict(list)
    for r in rows:
        b0 = math.floor(r["ocean_temp"] / bin_width) * bin_width
        bins[b0].append(r["pred_temp"])
    out = []
    for b0 in sorted(bins):
        vals = bins[b0]
        if len(vals) >= 25:
            out.append((b0 + bin_width / 2, median(vals), len(vals)))
    return out


def write_summary(rows, counts):
    out = OUT / "phaC_temstapro_vs_ocean_temperature_summary.tsv"
    therm_counts = Counter(r["thermophilicity"] for r in rows)
    label_counts = Counter(r["label"] for r in rows)
    with out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for key in sorted(counts):
            w.writerow([key, counts[key]])
        w.writerow(["median_ocean_temp_C", f"{median([r['ocean_temp'] for r in rows]):.3f}"])
        w.writerow(["median_temstapro_class_midpoint_C", f"{median([r['pred_temp'] for r in rows]):.3f}"])
        for key, val in therm_counts.most_common():
            w.writerow([f"thermophilicity:{key}", val])
        for key, val in label_counts.most_common():
            w.writerow([f"right_hand_label:{key}", val])
    print(f"wrote {out}")


def write_cold_high_stability_tables(rows):
    cold_hot = [r for r in rows if r["ocean_temp"] <= 5 and r["pred_temp"] >= 82.5]
    detail_out = OUT / "phaC_temstapro_cold_ocean_predicted_80C.tsv"
    with detail_out.open("w", newline="") as f:
        fields = [
            "target_id", "genome", "woa23_temp_annual_degC", "predicted_bin_midpoint_C",
            "temstapro_label", "thermophilicity", "length_aa", "gtdb_phylum", "gtdb_family",
            "gtdb_genus", "gtdb_species", "study_id", "depth_m", "depth_zone",
        ]
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in sorted(cold_hot, key=lambda x: (x["genus"], x["species"], x["ocean_temp"])):
            w.writerow({
                "target_id": r["target_id"],
                "genome": r["genome"],
                "woa23_temp_annual_degC": f"{r['ocean_temp']:.3f}",
                "predicted_bin_midpoint_C": f"{r['pred_temp']:.1f}",
                "temstapro_label": r["label"],
                "thermophilicity": r["thermophilicity"],
                "length_aa": f"{r['length']:.0f}",
                "gtdb_phylum": r["phylum"],
                "gtdb_family": r["family"],
                "gtdb_genus": r["genus"],
                "gtdb_species": r["species"],
                "study_id": r["study_id"],
                "depth_m": "" if r["depth"] is None else f"{r['depth']:.1f}",
                "depth_zone": r["depth_zone"],
            })

    summary_out = OUT / "phaC_temstapro_cold_ocean_predicted_80C_taxa.tsv"
    grouped = defaultdict(list)
    for r in cold_hot:
        grouped[(r["genus"], r["species"], r["family"], r["phylum"])].append(r)
    with summary_out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "gtdb_genus", "gtdb_species", "gtdb_family", "gtdb_phylum", "n_proteins",
            "median_ocean_temp_C", "median_depth_m", "studies",
        ])
        for key, vals in sorted(grouped.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            depths = [r["depth"] for r in vals if r["depth"] is not None]
            w.writerow([
                key[0], key[1], key[2], key[3], len(vals),
                f"{median([r['ocean_temp'] for r in vals]):.3f}",
                "" if not depths else f"{median(depths):.1f}",
                ",".join(sorted({r["study_id"] for r in vals})),
            ])
    print(f"wrote {detail_out}")
    print(f"wrote {summary_out}")
    return cold_hot


def main():
    rows, counts = load_rows()
    if not rows:
        raise SystemExit("No QC-passed rows with both TemStaPro label and WOA temperature.")
    write_summary(rows, counts)
    cold_hot = write_cold_high_stability_tables(rows)

    rng = np.random.default_rng(7)
    x = np.array([r["ocean_temp"] for r in rows])
    y = np.array([r["pred_temp"] for r in rows])
    # Light vertical jitter reveals density within the discrete TemStaPro bins.
    y_jitter = y + rng.uniform(-0.85, 0.85, size=len(y))

    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig = plt.figure(figsize=(15.0, 7.4), dpi=300)
    ax = fig.add_axes([0.07, 0.17, 0.58, 0.70])
    ax_hist = fig.add_axes([0.73, 0.33, 0.18, 0.34], sharey=ax)

    # Plot the dominant mesophilic cloud first, then warmer classes over it.
    for group in ["mesophilic", "undetermined", "thermophilic", "hyperthermophilic"]:
        idx = [i for i, r in enumerate(rows) if r["thermophilicity"] == group]
        if not idx:
            continue
        ax.scatter(
            x[idx], y_jitter[idx],
            s=9 if group == "mesophilic" else 13,
            c=THERM_COLORS.get(group, "#888888"),
            alpha=0.55 if group == "mesophilic" else 0.80,
            edgecolors="none",
            rasterized=True,
            label=f"{group} (n={len(idx):,})",
        )

    # The binned-median overlay that used to sit here is gone. It drew the median
    # predicted class within each 5 deg C WOA bin, which is 37.5 in every bin
    # because 3,827 of 4,224 proteins fall in the single <40 class -- a flat black
    # line across the bottom of the panel carrying no information, and read as an
    # unexplained annotation. binned_medians() is kept; the stats table still
    # reports those numbers, where a reader can see the bin sizes alongside them.

    ax.set_xlabel("WOA23 annual ocean temperature at genome sample depth (deg C)")
    ax.set_ylabel("TemStaPro predicted stability class (deg C bin midpoint)")
    ax.set_xlim(-2.5, 31)
    ax.set_ylim(34.5, 85.5)
    ax.set_yticks(TEMP_ORDER)
    ax.set_yticklabels(["<40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", ">=80"])
    ax.grid(True, color="#D8E0DD", linewidth=0.6, alpha=0.85)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles, labels,
        frameon=False, fontsize=8.5,
        loc="upper left", bbox_to_anchor=(0.72, 0.88),
        handletextpad=0.25, borderpad=0.2,
    )

    # Side histogram by predicted stability bin.
    label_counts = Counter(r["pred_temp"] for r in rows)
    hist_y = TEMP_ORDER
    hist_x = [label_counts.get(t, 0) for t in hist_y]
    ax_hist.barh(hist_y, hist_x, height=3.4, color="#9CB2AF", edgecolor="white", linewidth=0.5)
    ax_hist.set_xlabel("n proteins")
    ax_hist.tick_params(axis="y", labelleft=False)
    ax_hist.grid(True, axis="x", color="#D8E0DD", linewidth=0.6, alpha=0.85)
    ax_hist.set_title("Predicted class count", fontsize=10.5, fontweight="bold", pad=8)

    # Call out the cold-ocean subset.
    cold = [r for r in rows if r["ocean_temp"] < 0]
    if cold:
        # Report how many of the sub-zero proteins are predicted thermostable
        # rather than their median predicted bin. That median is 37.5 for every
        # subset of this data, for the same reason the binned-median line was
        # removed above, so it looked like a measurement while saying nothing.
        cold_hot_n = sum(1 for r in cold if r["pred_temp"] >= 55)
        fig.text(
            0.72, 0.125,
            f"WOA < 0 deg C: n={len(cold):,}\n"
            f"{cold_hot_n} predicted >=55 deg C ({100 * cold_hot_n / len(cold):.1f}%)",
            fontsize=8.8, color="#324346",
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor="#B9C7C4", linewidth=0.6, alpha=0.88),
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
        )
    if cold_hot:
        top_genera = Counter(r["genus"] for r in cold_hot).most_common(5)
        fig.text(
            0.72, 0.155,
            "WOA <= 5 deg C and predicted >=80 deg C\n"
            + "\n".join(
                [
                    "; ".join(f"{g} ({n})" for g, n in top_genera[:3]),
                    "; ".join(f"{g} ({n})" for g, n in top_genera[3:]),
                ]
            ),
            fontsize=8.4, color="#324346", ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor="#B9C7C4", linewidth=0.6, alpha=0.88),
        )

    fig.text(
        0.07, 0.072,
        "TemStaPro class is shown as the midpoint of right_hand_label bins (<40 plotted at 37.5 deg C; >=80 at 82.5 deg C).\n"
        "Points are QC-filtered with figures/scripts/_phac_qc.py to remove known mislabeled phaC reference hits, and length-filtered "
        f"to >={MIN_LENGTH}aa (short fragments produce unreliable, non-monotonic stability calls).",
        fontsize=8.2, color="#5B6B69",
    )
    fig.text(
        0.07, 0.025,
        f"Rows: {counts['all_temstapro_rows']:,} TemStaPro predictions; "
        f"{counts['qc_removed_bad_target']:,} removed by phaC QC; "
        f"{counts['qc_removed_short_fragment']:,} removed as <{MIN_LENGTH}aa fragments; "
        f"{counts['plotted_rows']:,} passed predictions with WOA temperature plotted.",
        fontsize=8.2, color="#5B6B69",
    )

    png = OUT / "phaC_temstapro_vs_ocean_temperature.png"
    pdf = OUT / "phaC_temstapro_vs_ocean_temperature.pdf"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    print(f"plotted_rows={len(rows)} qc_removed={counts['qc_removed_bad_target']}")


if __name__ == "__main__":
    main()
