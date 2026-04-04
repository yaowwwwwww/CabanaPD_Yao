#!/usr/bin/env python3
import csv
import math
import os
import re
import sys
from collections import defaultdict


PAT = re.compile(r"particles_(\d+)_all\.csv$")


def collect_run(run_dir, zmin=-2e-6, zmax=0.0):
    files = []
    for name in os.listdir(run_dir):
        m = PAT.match(name)
        if m:
            files.append((int(m.group(1)), os.path.join(run_dir, name)))
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
    if len(sys.argv) < 4:
        print(
            "usage: plot_interface_history_overlay.py OUT_DIR OUT_TAG label::RUN_DIR [label::RUN_DIR ...]",
            file=sys.stderr,
        )
        return 2

    out_dir = sys.argv[1]
    out_tag = sys.argv[2]
    specs = sys.argv[3:]
    os.makedirs(out_dir, exist_ok=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"ERROR: matplotlib unavailable: {exc}", file=sys.stderr)
        return 1

    fig, ax = plt.subplots(figsize=(9, 6), dpi=180)
    wrote_any = False
    summary_rows = []

    for spec in specs:
        if "::" not in spec:
            print(f"skip malformed spec: {spec}", file=sys.stderr)
            continue
        label, run_dir = spec.split("::", 1)
        if not os.path.isdir(run_dir):
            print(f"skip missing run dir: {run_dir}", file=sys.stderr)
            continue
        raw, means = collect_run(run_dir)
        raw_tsv = os.path.join(out_dir, f"{out_tag}_{label}_raw.tsv")
        mean_tsv = os.path.join(out_dir, f"{out_tag}_{label}_frame_mean.tsv")

        with open(raw_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["frame", "plastic_strain_rate_1ps", "yield_stress_Pa"])
            writer.writerows(raw)

        with open(mean_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(
                ["frame", "mean_plastic_strain_rate_1ps", "mean_yield_stress_Pa", "count"]
            )
            writer.writerows(means)

        if raw:
            # thin scatter to keep output readable
            stride = max(1, len(raw) // 3000)
            xs = [r[1] for r in raw[::stride]]
            ys = [r[2] / 1e9 for r in raw[::stride]]
            ax.scatter(xs, ys, s=4, alpha=0.08)
        if means:
            mx = [m[1] for m in means]
            my = [m[2] / 1e9 for m in means]
            ax.plot(mx, my, linewidth=2.2, label=label)
            wrote_any = True
            summary_rows.append((label, len(raw), len(means), raw_tsv, mean_tsv))

    ax.set_xscale("log")
    ax.set_xlabel("Plastic strain rate (1/s)")
    ax.set_ylabel("Yield stress (GPa)")
    ax.set_title("Interface point stress-rate history through impact")
    ax.grid(True, which="both", alpha=0.25)
    if wrote_any:
        ax.legend(fontsize=8)
    fig.tight_layout()
    out_png = os.path.join(out_dir, f"{out_tag}.png")
    fig.savefig(out_png)

    summary_tsv = os.path.join(out_dir, f"{out_tag}_sources.tsv")
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["label", "raw_points", "frames_with_points", "raw_tsv", "mean_tsv"])
        writer.writerows(summary_rows)

    print(out_png)
    print(summary_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
