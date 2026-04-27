#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


def read_value(node):
    if isinstance(node, dict) and "value" in node:
        return node["value"]
    return node


def scalar(node) -> float:
    value = read_value(node)
    if isinstance(value, list):
        return float(value[0])
    return float(value)


@dataclass
class RunCase:
    run_dir: Path
    label: str
    input_path: Path
    dx: float
    dt: float
    output_frequency: int
    alpha1: float
    beta: float
    r0: float
    radius: float
    volume: float
    bulk_modulus: float
    cap: float
    fmax: float
    lj_linearize_r: float
    lj_linearize_slope: float

    def raw_lj_force_density(self, r: float) -> float:
        alpha_eff = (
            self.bulk_modulus
            * 18.0
            / (math.pi * (self.delta**4) * self.alpha1)
            * self.r0
            * self.r0
            * self.volume
            * self.volume
            / 72.0
            / (self.beta ** (7.0 / 3.0))
        )
        term13 = (self.r0 / r) ** 13.0
        term7 = (self.r0 / r) ** 7.0
        return ((12.0 * alpha_eff) / self.r0) * (self.beta * term7 - term13) / self.volume

    @property
    def delta(self) -> float:
        return self._delta

    @delta.setter
    def delta(self, value: float) -> None:
        self._delta = value

    def lj_force_density(self, r: float) -> float:
        if r <= 1.0e-14 or r > self.radius:
            return 0.0
        raw = self.raw_lj_force_density(r)
        if (
            self.cap <= 0.0
            or self.lj_linearize_r <= 0.0
            or r >= self.lj_linearize_r
            or raw <= self.cap
        ):
            return raw
        return self.fmax + self.lj_linearize_slope * r


def build_case(run_dir: Path, label: str) -> RunCase:
    input_path = run_dir / "input.json"
    if not input_path.exists():
        input_path = run_dir / "_current_input.json"
    data = json.loads(input_path.read_text(encoding="utf-8"))

    low = read_value(data["low_corner"])
    high = read_value(data["high_corner"])
    ncells = read_value(data["num_cells"])
    dx = (float(high[0]) - float(low[0])) / float(ncells[0])

    elastic = read_value(data["elastic_modulus"])
    poisson = read_value(data["Poisson's_ratio"])
    elastic0 = float(elastic[0] if isinstance(elastic, list) else elastic)
    poisson0 = float(poisson[0] if isinstance(poisson, list) else poisson)
    bulk_modulus = elastic0 / (3.0 * (1.0 - 2.0 * poisson0))

    delta = scalar(data["horizon"])
    alpha1 = scalar(data["LJalpha"])
    beta = scalar(data["LJbeta"])
    r0 = scalar(data["LJr0"]) * dx
    radius = scalar(data["contact_horizon_factor"]) * r0
    dt = scalar(data["timestep"])
    output_frequency = int(read_value(data["output_frequency"]))
    cap = scalar(data.get("LJ_linearize_force_density_cap", {"value": 0.0}))
    fmax = scalar(data.get("LJ_linearize_force_density_max", {"value": 0.0}))
    volume = dx**3

    case = RunCase(
        run_dir=run_dir,
        label=label,
        input_path=input_path,
        dx=dx,
        dt=dt,
        output_frequency=output_frequency,
        alpha1=alpha1,
        beta=beta,
        r0=r0,
        radius=radius,
        volume=volume,
        bulk_modulus=bulk_modulus,
        cap=cap,
        fmax=fmax,
        lj_linearize_r=0.0,
        lj_linearize_slope=0.0,
    )
    case.delta = delta
    prepare_linearization(case)
    return case


def prepare_linearization(case: RunCase) -> None:
    if case.cap <= 0.0 or case.fmax <= 0.0:
        case.lj_linearize_r = 0.0
        case.lj_linearize_slope = 0.0
        return

    lo = max(case.r0 * 1.0e-6, 1.0e-14)
    hi = case.r0
    raw_hi = case.raw_lj_force_density(hi)
    raw_lo = case.raw_lj_force_density(lo)

    if raw_hi > case.cap:
        case.lj_linearize_r = hi
    elif raw_lo <= case.cap:
        case.lj_linearize_r = 0.0
        case.lj_linearize_slope = 0.0
        return
    else:
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if case.raw_lj_force_density(mid) > case.cap:
                lo = mid
            else:
                hi = mid
        case.lj_linearize_r = hi

    case.lj_linearize_slope = (case.cap - case.fmax) / case.lj_linearize_r


def particle_csvs(run_dir: Path) -> list[Path]:
    def frame_id(path: Path) -> int:
        match = re.search(r"particles_(\d+)_all\.csv$", path.name)
        return int(match.group(1)) if match else -1

    return sorted(run_dir.glob("particles_*_all.csv"), key=frame_id)


def frame_id(path: Path) -> int:
    match = re.search(r"particles_(\d+)_all\.csv$", path.name)
    if not match:
        raise ValueError(f"Unexpected csv filename: {path}")
    return int(match.group(1))


def minimum_pair_distance(
    csv_path: Path, projectile_type: int = 0, substrate_type: int = 1
) -> tuple[float, tuple[float, float, float], tuple[float, float, float]]:
    projectile_points = []
    substrate_points = []

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ptype = int(float(row["rank_0/type"]))
            point = (
                float(row["Points:0"]),
                float(row["Points:1"]),
                float(row["Points:2"]),
            )
            if ptype == projectile_type:
                projectile_points.append(point)
            elif ptype == substrate_type:
                substrate_points.append(point)

    if not projectile_points or not substrate_points:
        raise RuntimeError(f"Could not find projectile/substrate particles in {csv_path}")

    best_projectile = None
    best_substrate = None
    best_d2 = None
    for px, py, pz in projectile_points:
        for sx, sy, sz in substrate_points:
            d2 = (sx - px) ** 2 + (sy - py) ** 2 + (sz - pz) ** 2
            if best_d2 is None or d2 < best_d2:
                best_d2 = d2
                best_projectile = (px, py, pz)
                best_substrate = (sx, sy, sz)

    return math.sqrt(best_d2), best_projectile, best_substrate


def write_single_case_outputs(case: RunCase) -> tuple[Path, list[tuple[float, float, float]]]:
    rows = []
    for csv_path in particle_csvs(case.run_dir):
        frame = frame_id(csv_path)
        time_ns = frame * case.output_frequency * case.dt * 1.0e9
        distance_m, projectile_point, substrate_point = minimum_pair_distance(csv_path)
        force_density = case.lj_force_density(distance_m)
        force_n = force_density * case.volume
        rows.append(
            (
                frame,
                time_ns,
                distance_m,
                force_n,
                projectile_point[0],
                projectile_point[1],
                projectile_point[2],
                substrate_point[0],
                substrate_point[1],
                substrate_point[2],
            )
        )

    out_tsv = case.run_dir / "minimum_distance_and_equiv_force_vs_time.tsv"
    with out_tsv.open("w", encoding="utf-8") as handle:
        handle.write(
            "frame\ttime_ns\tdistance_m\tequiv_force_N\t"
            "proj_x\tproj_y\tproj_z\tnearest_sub_x\tnearest_sub_y\tnearest_sub_z\n"
        )
        for row in rows:
            handle.write(
                f"{row[0]}\t{row[1]:.9g}\t{row[2]:.9g}\t{row[3]:.9g}\t"
                f"{row[4]:.9g}\t{row[5]:.9g}\t{row[6]:.9g}\t"
                f"{row[7]:.9g}\t{row[8]:.9g}\t{row[9]:.9g}\n"
            )
    return out_tsv, rows


def run_gnuplot(script: str) -> None:
    subprocess.run(["gnuplot"], input=script, text=True, check=True)


def plot_overlay(cases: list[RunCase], output_png: Path, title: str) -> Path:
    generated = []
    for case in cases:
        tsv, rows = write_single_case_outputs(case)
        generated.append((case, tsv, rows))

    with tempfile.TemporaryDirectory(prefix="minimum_distance_force_") as tmpdir:
        legend_lines = []
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
        idx = 1
        for (case, tsv, _), color in zip(generated, colors):
            q = str(tsv).replace("\\", "\\\\").replace("'", "\\'")
            label = case.label.replace("'", "\\'")
            legend_lines.append(
                f"'{q}' using 2:($3*1e6) with lines lw 2 lc rgb '{color}' title '{label} distance', \\"
            )
            legend_lines.append(
                f"'{q}' using 2:4 axes x1y2 with lines dt 2 lw 2 lc rgb '{color}' title '{label} force', \\"
            )
            idx += 1

        if legend_lines:
            legend_lines[-1] = legend_lines[-1].rstrip(" ,\\")

        script = f"""
set terminal pngcairo size 1400,900 enhanced
set output '{str(output_png).replace("\\\\", "\\\\\\\\").replace("'", "\\\\'")}'
set title '{title.replace("'", "\\\\'")}'
set xlabel 'Time (ns)'
set ylabel 'Minimum projectile-substrate distance (um)'
set y2label 'Equivalent contact force from minimum distance (N)'
set ytics nomirror
set y2tics
set grid
set key outside right center
plot \\
{chr(10).join(legend_lines)}
"""
        run_gnuplot(script)

    combined_tsv = output_png.with_suffix(".tsv")
    with combined_tsv.open("w", encoding="utf-8") as handle:
        handle.write("series\tframe\ttime_ns\tdistance_m\tequiv_force_N\n")
        for case, _, rows in generated:
            for row in rows:
                handle.write(
                    f"{case.label}\t{row[0]}\t{row[1]:.9g}\t{row[2]:.9g}\t{row[3]:.9g}\n"
                )
    return combined_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For each frame, compute the minimum projectile-substrate distance and "
            "map that distance through the current "
            "LJ contact law."
        )
    )
    parser.add_argument("run_dirs", nargs="+", help="One or more run directories.")
    parser.add_argument("--labels", nargs="*", help="Optional labels matching run_dirs.")
    parser.add_argument("--output", help="Output PNG for overlay plot.")
    parser.add_argument("--title", default="Minimum Distance and LJ Force")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.labels and len(args.labels) != len(args.run_dirs):
        raise SystemExit("--labels must have the same length as run_dirs")

    cases = []
    for i, run_dir_text in enumerate(args.run_dirs):
        run_dir = Path(run_dir_text)
        label = args.labels[i] if args.labels else run_dir.name
        cases.append(build_case(run_dir, label))

    if args.output:
        output_png = Path(args.output)
    elif len(cases) == 1:
        output_png = cases[0].run_dir / "minimum_distance_and_equiv_force_vs_time.png"
    else:
        output_png = Path.cwd() / "minimum_distance_and_equiv_force_overlay.png"

    combined_tsv = plot_overlay(cases, output_png, args.title)
    print(output_png)
    print(combined_tsv)


if __name__ == "__main__":
    main()
