#!/usr/bin/env python3
import csv
import math
import os
import re
import sys


PAT = re.compile(r"particles_(\d+)_all\.csv$")


def collect_case(case_dir, zmin=-2e-6, zmax=0.0):
    files = []
    for name in os.listdir(case_dir):
        m = PAT.match(name)
        if m:
            files.append((int(m.group(1)), os.path.join(case_dir, name)))
    files.sort()

    raw = []
    means = []
    for frame, path in files:
        pts = []
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    typ = int(float(row["rank_0/type"]))
                    z = float(row["Points:2"])
                    ed = float(row["rank_0/plastic_strain_rate"])
                    ys = float(row["rank_0/yield_stress"])
                except Exception:
                    continue
                if typ != 0:
                    continue
                if not (zmin <= z <= zmax):
                    continue
                if not (math.isfinite(ed) and math.isfinite(ys)):
                    continue
                if ed <= 0.0 or ys <= 0.0:
                    continue
                pts.append((ed, ys))
                raw.append((frame, ed, ys))
        if pts:
            means.append(
                (
                    frame,
                    sum(p[0] for p in pts) / len(pts),
                    sum(p[1] for p in pts) / len(pts),
                    len(pts),
                )
            )
    return raw, means


def main() -> int:
    if len(sys.argv) != 4:
        print(
            "usage: extract_case_interface_history.py RUNS_ROOT SUMMARY_FILE OUT_DIR",
            file=sys.stderr,
        )
        return 2

    runs_root, summary_file, out_dir = sys.argv[1:]
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    with open(summary_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            case_name = line.split()[0]
            case_dir = os.path.join(runs_root, f"run_{case_name}")
            if not os.path.isdir(case_dir):
                continue
            raw, means = collect_case(case_dir)
            if not raw or not means:
                continue

            raw_tsv = os.path.join(out_dir, f"{case_name}_raw.tsv")
            mean_tsv = os.path.join(out_dir, f"{case_name}_frame_mean.tsv")
            png = os.path.join(out_dir, f"{case_name}_stress_rate.png")

            with open(raw_tsv, "w", newline="", encoding="utf-8") as rf:
                writer = csv.writer(rf, delimiter="\t")
                writer.writerow(["frame", "plastic_strain_rate_1ps", "yield_stress_Pa"])
                writer.writerows(raw)

            with open(mean_tsv, "w", newline="", encoding="utf-8") as mf:
                writer = csv.writer(mf, delimiter="\t")
                writer.writerow(
                    ["frame", "mean_plastic_strain_rate_1ps", "mean_yield_stress_Pa", "count"]
                )
                writer.writerows(means)

            rows.append((case_name, raw_tsv, mean_tsv, png, len(raw), len(means)))

    manifest = os.path.join(out_dir, "case_interface_history_manifest.tsv")
    with open(manifest, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            ["case_name", "raw_tsv", "mean_tsv", "out_png", "raw_points", "frames_with_points"]
        )
        writer.writerows(rows)

    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
