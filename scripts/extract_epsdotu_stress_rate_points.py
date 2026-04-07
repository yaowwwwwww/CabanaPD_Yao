#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


RUN_DIR_RE = re.compile(r"run_.*_U_(?P<epsu>[^_]+)_tq_")
FRAME_CSV_RE = re.compile(r"particles_(\d+)_all\.csv$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract all active near-impact substrate points for stress-rate plotting."
    )
    p.add_argument("run_root", type=Path, help="Run root containing run_* case directories")
    p.add_argument(
        "--epsdot-u",
        nargs="*",
        default=["6.8e6", "4e7", "6.8e5", "5.8e4"],
        help="epsdot_u values to include, in output order",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Defaults to run_root/epsdotu_pointcloud_r5e-6",
    )
    p.add_argument("--zmin", type=float, default=-2.0e-6)
    p.add_argument("--zmax", type=float, default=0.0)
    p.add_argument("--rmax", type=float, default=5.0e-6)
    return p.parse_args()


def iter_case_frames(case_dir: Path):
    frames = []
    for path in case_dir.iterdir():
        match = FRAME_CSV_RE.match(path.name)
        if match:
            frames.append((int(match.group(1)), path))
    return sorted(frames)


def extract_case_points(case_dir: Path, epsu: str, zmin: float, zmax: float, rmax: float | None):
    band_thickness = abs(zmax - zmin)
    rmax2 = None if rmax is None else rmax * rmax
    out_rows = []

    for frame, csv_path in iter_case_frames(case_dir):
        pts = []
        top_z = None
        with csv_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                try:
                    typ = int(float(row["rank_0/type"]))
                    x = float(row["Points:0"])
                    y = float(row["Points:1"])
                    z = float(row["Points:2"])
                    eps_p = float(row["rank_0/plastic_strain"])
                    eps_dot = float(row["rank_0/plastic_strain_rate"])
                    ys = float(row["rank_0/yield_stress"])
                except Exception:
                    continue
                if typ != 0:
                    continue
                if not (math.isfinite(eps_p) and math.isfinite(eps_dot) and math.isfinite(ys)):
                    continue
                pts.append((x, y, z, eps_p, eps_dot, ys))
                top_z = z if top_z is None else max(top_z, z)

        if top_z is None:
            continue

        z_lo = top_z - band_thickness
        z_hi = top_z

        for x, y, z, eps_p, eps_dot, ys in pts:
            if not (z_lo <= z <= z_hi):
                continue
            if rmax2 is not None and (x * x + y * y) > rmax2:
                continue
            if eps_p <= 0.0 and eps_dot <= 0.0:
                continue
            out_rows.append((epsu, frame, eps_p, eps_dot, ys))

    return out_rows


def worker(job):
    case_dir, epsu, zmin, zmax, rmax = job
    return extract_case_points(case_dir, epsu, zmin, zmax, rmax)


def main() -> int:
    args = parse_args()
    run_root = args.run_root.resolve()
    out_dir = (args.output_dir or (run_root / "epsdotu_pointcloud_r5e-6")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = list(args.epsdot_u)
    jobs = []
    found = set()
    for case_dir in sorted(run_root.iterdir()):
        if not case_dir.is_dir() or not case_dir.name.startswith("run_"):
            continue
        m = RUN_DIR_RE.search(case_dir.name)
        if not m:
            continue
        epsu = m.group("epsu")
        if epsu in wanted and epsu not in found:
            jobs.append((case_dir, epsu, args.zmin, args.zmax, args.rmax))
            found.add(epsu)

    if not jobs:
        raise SystemExit(f"No matching run_* case directories found in {run_root}")

    all_rows = []
    with ProcessPoolExecutor() as ex:
        for rows in ex.map(worker, jobs):
            all_rows.extend(rows)

    wanted_order = {epsu: i for i, epsu in enumerate(wanted)}
    all_rows.sort(key=lambda r: (wanted_order.get(r[0], 999), r[1], r[3], r[4]))

    out_tsv = out_dir / "stress_rate_points_by_epsdot_u.tsv"
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["epsdot_u_label", "frame", "plastic_strain", "plastic_strain_rate_1ps", "yield_stress_Pa"]
        )
        w.writerows(all_rows)

    print(out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
