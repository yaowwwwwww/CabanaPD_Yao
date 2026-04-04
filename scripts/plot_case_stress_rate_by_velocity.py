#!/usr/bin/env python3
import csv
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_RE_V = re.compile(r"^(v\d+)_")
CASE_RE_K0 = re.compile(r"_K0_([^_]+)_")


def read_manifest(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def read_raw_points(path):
    pts = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            try:
                x = float(row[1])
                y = float(row[2])
            except ValueError:
                continue
            if x > 0.0 and y > 0.0:
                pts.append((x, y))
    return pts


def binned_mean_curve(points):
    if len(points) < 6:
        return [], []
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    logx = [math.log10(x) for x in xs]
    nbins = min(30, max(12, len(xs) // 20))
    lo = min(logx)
    hi = max(logx)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return [], []
    edges = [lo + (hi - lo) * i / nbins for i in range(nbins + 1)]
    bx = []
    by = []
    for b in range(nbins):
        if b < nbins - 1:
            idx = [i for i, lx in enumerate(logx) if edges[b] <= lx < edges[b + 1]]
        else:
            idx = [i for i, lx in enumerate(logx) if edges[b] <= lx <= edges[b + 1]]
        if len(idx) < 3:
            continue
        mean_logx = sum(logx[i] for i in idx) / len(idx)
        mean_y = sum(ys[i] for i in idx) / len(idx)
        bx.append(10 ** mean_logx)
        by.append(mean_y)
    return bx, by


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: plot_case_stress_rate_by_velocity.py <manifest.tsv>")

    manifest = sys.argv[1]
    out_dir = os.path.dirname(os.path.abspath(manifest))
    rows = read_manifest(manifest)

    grouped = {}
    for row in rows:
        case_name = row["case_name"]
        m_v = CASE_RE_V.search(case_name)
        m_k0 = CASE_RE_K0.search(case_name)
        if not m_v or not m_k0:
            continue
        vel = m_v.group(1)
        k0 = m_k0.group(1)
        grouped.setdefault(vel, []).append((k0, row["raw_tsv"]))

    for vel, entries in sorted(grouped.items()):
        fig, ax = plt.subplots(figsize=(8.6, 6.2), dpi=180)
        plotted = 0
        for k0, raw_tsv in sorted(entries, key=lambda kv: float(kv[0].replace("e", "E"))):
            points = read_raw_points(raw_tsv)
            bx, by = binned_mean_curve(points)
            if len(bx) < 2:
                continue
            ax.semilogx(
                bx,
                [y / 1e9 for y in by],
                "-o",
                markersize=4.5,
                markerfacecolor="none",
                linewidth=1.8,
                label=f"K0={k0}",
            )
            plotted += 1
        if plotted == 0:
            plt.close(fig)
            continue
        ax.set_xlabel("Plastic strain rate (1/s)")
        ax.set_ylabel("Flow stress (GPa)")
        ax.set_title(f"{vel}: mean flow stress by strain-rate bin")
        ax.grid(True, which="both", alpha=0.35)
        ax.legend(loc="best")
        fig.tight_layout()
        out_png = os.path.join(out_dir, f"{vel}_stress_rate_by_k0.png")
        fig.savefig(out_png)
        plt.close(fig)


if __name__ == "__main__":
    main()
