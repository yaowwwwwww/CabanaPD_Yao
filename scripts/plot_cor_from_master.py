#!/usr/bin/env python3
import csv
import math
import os
import sys
from collections import defaultdict


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: plot_cor_from_master.py MASTER_SUMMARY_TSV OUT_PNG", file=sys.stderr)
        return 2

    summary_path, out_png = sys.argv[1:]
    out_tsv = os.path.splitext(out_png)[0] + ".tsv"

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"ERROR: matplotlib unavailable: {exc}", file=sys.stderr)
        return 1

    series = defaultdict(list)
    with open(summary_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                k0 = float(row["drag_K0"])
                vin = abs(float(row["vin_mps"]))
                cor = float(row["CoR"])
            except Exception:
                continue
            if not (math.isfinite(vin) and math.isfinite(cor)):
                continue
            series[k0].append((vin, cor))

    rows = []
    fig, ax = plt.subplots(figsize=(9, 6), dpi=180)
    for k0, pts in sorted(series.items()):
        pts.sort()
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, "-o", linewidth=2, markersize=5, label=f"K0={k0:g} Pa")
        for vin, cor in pts:
            rows.append((k0, vin, cor))

    ax.set_xlabel("Velocity (m/s)")
    ax.set_ylabel("CoR")
    ax.set_title("CoR vs velocity for log-drag K0 scan")
    ax.grid(True, alpha=0.3)
    if series:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png)

    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["K0_Pa", "velocity_mps", "CoR"])
        writer.writerows(sorted(rows))

    print(out_png)
    print(out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
