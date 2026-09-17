"""Export the coldest WOA-mapped phaC genomes for quick inspection."""
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent
IN = ROOT / "data" / "temstapro_inputs" / "phaC_genomes_woa23_annual_temperature.tsv"
OUT = ROOT / "data" / "temstapro_inputs" / "phaC_genomes_woa23_coldest.tsv"


def main():
    with IN.open(newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r.get("woa23_temp_annual_degC")]
    rows.sort(key=lambda r: float(r["woa23_temp_annual_degC"]))
    cols = [
        "genome", "biosample", "study_id", "sample_id", "gtdb_genus",
        "gtdb_species", "latitude_degN", "longitude_degE", "depth_m",
        "depth_zone", "woa23_temp_annual_degC", "woa23_match_status",
        "woa23_depth_method",
    ]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=cols)
        w.writeheader()
        w.writerows({c: r.get(c, "") for c in cols} for r in rows)
    print(f"wrote {OUT}")
    for threshold in [-1, 0, 1, 2, 4]:
        print(f"<= {threshold} C", sum(float(r["woa23_temp_annual_degC"]) <= threshold for r in rows))


if __name__ == "__main__":
    main()
