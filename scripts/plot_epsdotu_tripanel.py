#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


RUN_DIR_RE = re.compile(r"run_.*_U_(?P<epsu>[^_]+)_tq_")
FRAME_CSV_RE = re.compile(r"particles_(\d+)_all\.csv$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot 3-panel mean substrate response vs time for different jc_epsdot_u values."
    )
    p.add_argument("run_root", type=Path, help="Run root containing run_* case directories")
    p.add_argument(
        "--epsdot-u",
        nargs="*",
        default=["6.8e6", "4e7", "6.8e5", "5.8e4"],
        help="epsdot_u values to include, in plot order",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for output TSV/PNG files. Defaults to run_root/epsdotu_tripanel",
    )
    p.add_argument("--zmin", type=float, default=-2.0e-6)
    p.add_argument("--zmax", type=float, default=0.0)
    p.add_argument("--rmax", type=float, default=None)
    return p.parse_args()


def read_case_time_scale(case_dir: Path) -> float:
    input_json = case_dir / "input.json"
    if not input_json.exists():
        return 2.0e-9
    data = json.loads(input_json.read_text(encoding="utf-8"))
    return float(data["timestep"]["value"]) * float(data["output_frequency"]["value"])


def iter_case_frames(case_dir: Path):
    frames = []
    for path in case_dir.iterdir():
        match = FRAME_CSV_RE.match(path.name)
        if match:
            frames.append((int(match.group(1)), path))
    return sorted(frames)


def extract_case_history(case_dir: Path, zmin: float, zmax: float, rmax: float | None):
    band_thickness = abs(zmax - zmin)
    rmax2 = None if rmax is None else rmax * rmax
    rows = {}

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
        region_count = 0
        active_count = 0
        eps_sum = 0.0
        epsdot_sum = 0.0
        ys_sum = 0.0

        for x, y, z, eps_p, eps_dot, ys in pts:
            if not (z_lo <= z <= z_hi):
                continue
            if rmax2 is not None and (x * x + y * y) > rmax2:
                continue
            region_count += 1
            if eps_p > 0.0 or eps_dot > 0.0:
                active_count += 1
                eps_sum += eps_p
                epsdot_sum += eps_dot
                ys_sum += ys

        if region_count == 0:
            continue

        if active_count:
            rows[frame] = (
                eps_sum / active_count,
                epsdot_sum / active_count,
                ys_sum / active_count,
                region_count,
                active_count,
            )
        else:
            rows[frame] = (0.0, 0.0, 0.0, region_count, 0)

    return rows


def collect_case_job(job):
    case_dir, epsu, zmin, zmax, rmax = job
    hist = extract_case_history(case_dir, zmin, zmax, rmax)
    dt_frame = read_case_time_scale(case_dir)
    return epsu, hist, dt_frame


def write_tsv(out_tsv: Path, data, dt_map, ordered_labels):
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "epsdot_u_label",
                "frame",
                "time_s",
                "mean_plastic_strain",
                "mean_plastic_strain_rate_1ps",
                "mean_yield_stress_Pa",
                "region_count",
                "active_count",
            ]
        )
        for epsu in ordered_labels:
            case = data[epsu]
            dt_frame = dt_map[epsu]
            max_frame = max(case) if case else 0
            for frame in range(max_frame + 1):
                eps, epsdot, ys, region_count, active_count = case.get(frame, (0.0, 0.0, 0.0, 0, 0))
                w.writerow([epsu, frame, frame * dt_frame, eps, epsdot, ys, region_count, active_count])


def main() -> int:
    args = parse_args()
    run_root = args.run_root.resolve()
    out_dir = (args.output_dir or (run_root / "epsdotu_tripanel")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = list(args.epsdot_u)
    jobs = []
    seen = OrderedDict()
    for case_dir in sorted(run_root.iterdir()):
        if not case_dir.is_dir() or not case_dir.name.startswith("run_"):
            continue
        m = RUN_DIR_RE.search(case_dir.name)
        if not m:
            continue
        epsu = m.group("epsu")
        if epsu in wanted and epsu not in seen:
            seen[epsu] = case_dir

    missing = [e for e in wanted if e not in seen]
    if missing:
        print("Missing epsdot_u cases:", " ".join(missing))

    if not seen:
        raise SystemExit(f"No matching run_* case directories found in {run_root}")

    with ProcessPoolExecutor() as ex:
        results = list(
            ex.map(
                collect_case_job,
                [(case_dir, epsu, args.zmin, args.zmax, args.rmax) for epsu, case_dir in seen.items()],
            )
        )

    data = OrderedDict()
    dt_map = {}
    for epsu, hist, dt_frame in results:
        if hist:
            data[epsu] = hist
            dt_map[epsu] = dt_frame

    ordered_labels = [e for e in wanted if e in data]
    out_tsv = out_dir / "mean_response_vs_time_by_epsdot_u.tsv"
    write_tsv(out_tsv, data, dt_map, ordered_labels)
    print(out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
