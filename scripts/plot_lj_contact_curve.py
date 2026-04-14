#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ContactCase:
    input_path: Path
    alpha_input: float
    beta: float
    r0: float
    radius: float
    delta: float
    dx: float
    volume: float
    bulk_modulus: float
    czm_scale: float
    czm_yield_stretch: float
    czm_decay_rate: float
    alpha_effective: float

    def lj_pair_force(self, r: float) -> float:
        term13 = (self.r0 / r) ** 13.0
        term7 = (self.r0 / r) ** 7.0
        return ((12.0 * self.alpha_effective) / self.r0) * (
            term13 - self.beta * term7
        )

    def czm_force(self, r: float) -> float:
        s = (r - 2.0e-6) / 2.0e-6
        if s < 0.0 and s >= -self.czm_yield_stretch:
            return self.czm_scale * (-s)
        if s < -self.czm_yield_stretch:
            return (
                self.czm_scale
                * self.czm_yield_stretch
                * math.exp(-self.czm_decay_rate * (-s - self.czm_yield_stretch))
            )
        return 0.0

    def total_pair_force(self, r: float) -> tuple[float, float, float]:
        lj = self.lj_pair_force(r)
        czm = self.czm_force(r)
        return lj - czm, lj, czm


def shell_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace("'", "\\'")


def scientific_label(value: float) -> str:
    return f"{value:.1e}"


def linspace(start: float, stop: float, count: int) -> list[float]:
    if count <= 1:
        return [start]
    step = (stop - start) / (count - 1)
    return [start + i * step for i in range(count)]


def run_gnuplot(script_text: str) -> None:
    subprocess.run(
        ["gnuplot"],
        input=script_text,
        text=True,
        check=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot the implemented LJ contact force curve from CabanaPD input.json "
            "files or from a scan directory containing multiple runs."
        )
    )
    parser.add_argument(
        "target",
        help="Path to an input.json file or to a directory containing run subdirectories.",
    )
    parser.add_argument(
        "--output",
        help="Output PNG path. Defaults to a file next to the target.",
    )
    parser.add_argument(
        "--tsv-output",
        help="Optional output TSV path for sampled curve data.",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=500,
        help="Number of r-samples per curve.",
    )
    parser.add_argument(
        "--r-min-factor",
        type=float,
        default=0.9,
        help="Minimum plotted separation as a multiple of r0.",
    )
    return parser.parse_args()


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_value(node):
    if isinstance(node, dict) and "value" in node:
        return node["value"]
    return node


def ensure_list(value) -> list[float]:
    if isinstance(value, list):
        return value
    return [value]


def input_paths_from_target(target: Path) -> list[Path]:
    if target.is_file():
        return [target]
    paths = sorted(p for p in target.rglob("input.json") if p.is_file())
    if not paths:
        raise FileNotFoundError(f"No input.json files found under {target}")
    return paths


def build_case(input_path: Path) -> ContactCase:
    data = read_json(input_path)
    low_corner = read_value(data["low_corner"])
    high_corner = read_value(data["high_corner"])
    num_cells = read_value(data["num_cells"])
    dx = (high_corner[0] - low_corner[0]) / num_cells[0]

    delta = float(read_value(data["horizon"])) + 1.0e-10
    elastic_modulus = ensure_list(read_value(data["elastic_modulus"]))
    poisson_ratio = ensure_list(read_value(data["Poisson's_ratio"]))
    bulk_moduli = [
        e / (3.0 * (1.0 - 2.0 * nu)) for e, nu in zip(elastic_modulus, poisson_ratio)
    ]
    bulk_modulus = sum(bulk_moduli) / len(bulk_moduli)

    r0 = float(read_value(data["LJr0"])) * dx
    alpha_input = float(read_value(data["LJalpha"]))
    beta = float(read_value(data["LJbeta"]))
    radius = float(read_value(data["contact_horizon_factor"])) * r0

    volume = dx**3
    alpha_effective = (
        (18.0 * bulk_modulus / (math.pi * delta**4 * alpha_input))
        * r0
        * r0
        * volume
        * volume
        / 72.0
        / (beta ** (7.0 / 3.0))
    )

    return ContactCase(
        input_path=input_path,
        alpha_input=alpha_input,
        beta=beta,
        r0=r0,
        radius=radius,
        delta=delta,
        dx=dx,
        volume=volume,
        bulk_modulus=bulk_modulus,
        czm_scale=float(read_value(data.get("CZM_cohesive_scaling", {"value": 0.0}))),
        czm_yield_stretch=float(read_value(data.get("CZM_yield_stretch", {"value": 0.0}))),
        czm_decay_rate=float(read_value(data.get("CZM_degradation_rate", {"value": 0.0}))),
        alpha_effective=alpha_effective,
    )


def default_output_path(target: Path) -> Path:
    if target.is_file():
        return target.with_name("lj_contact_model_curve.png")
    return target / "lj_contact_model_curves.png"


def default_tsv_path(output_path: Path) -> Path:
    return output_path.with_suffix(".tsv")


def plot_cases(
    cases: list[ContactCase],
    output_path: Path,
    tsv_path: Path,
    samples: int,
    r_min_factor: float,
) -> None:
    cases = sorted(cases, key=lambda case: case.alpha_input)
    rows = [
        "alpha_input\tr_m\tr_um\tr_over_r0\tpair_force_total_N\tpair_force_lj_N\tpair_force_czm_N"
    ]

    tsv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    reference = cases[0]
    zero_crossing = reference.r0 * reference.beta ** (-1.0 / 6.0)

    with tempfile.TemporaryDirectory(prefix="lj_contact_plot_") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        curve_paths: list[Path] = []
        for case in cases:
            r_min = max(r_min_factor * case.r0, 1.0e-14)
            curve_path = tmp_dir / f"curve_alpha_{case.alpha_input:.3e}.dat"
            curve_paths.append(curve_path)
            with curve_path.open("w", encoding="utf-8") as handle:
                handle.write("# r_um\tpair_force_total_N\tpair_force_lj_N\tpair_force_czm_N\tr_over_r0\n")
                for r in linspace(r_min, case.radius, samples):
                    total_force, lj_force, czm_force = case.total_pair_force(r)
                    handle.write(
                        f"{r * 1.0e6:.12e}\t{total_force:.12e}\t{lj_force:.12e}\t"
                        f"{czm_force:.12e}\t{r / case.r0:.12e}\n"
                    )
                    rows.append(
                        f"{case.alpha_input:.12e}\t{r:.12e}\t{r * 1.0e6:.12e}\t{r / case.r0:.12e}\t"
                        f"{total_force:.12e}\t{lj_force:.12e}\t{czm_force:.12e}"
                    )

        tsv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        label_text = (
            f"dx={reference.dx * 1.0e6:.3f} um; "
            f"r0={reference.r0 * 1.0e6:.3f} um; "
            f"cutoff={reference.radius * 1.0e6:.3f} um; "
            f"beta={reference.beta:g}; "
            f"CZM scale={reference.czm_scale:g}"
        )
        plot_items = ", \\\n".join(
            [
                f"'{shell_quote(path)}' using 1:2 with lines lw 2 title 'LJalpha={scientific_label(case.alpha_input)}'"
                for case, path in zip(cases, curve_paths)
            ]
        )
        gnuplot_script = f"""
set terminal pngcairo size 1400,900 enhanced
set output '{shell_quote(output_path)}'
set title 'CabanaPD LJ Contact Force Curve'
set xlabel 'Separation r [um]'
set ylabel 'Pair contact force F_c(r) [N]'
set grid
set key top right
set zeroaxis
set arrow 1 from {reference.r0 * 1.0e6:.12e}, graph 0 to {reference.r0 * 1.0e6:.12e}, graph 1 nohead dt 2 lc rgb '#666666'
set arrow 2 from {zero_crossing * 1.0e6:.12e}, graph 0 to {zero_crossing * 1.0e6:.12e}, graph 1 nohead dt 3 lc rgb '#999999'
set label 1 '{label_text}' at graph 0.02, graph 0.95 front
plot {plot_items}
"""
        run_gnuplot(gnuplot_script)


def main() -> None:
    args = parse_args()
    target = Path(args.target).resolve()
    if not target.exists():
        raise FileNotFoundError(f"Target does not exist: {target}")

    input_paths = input_paths_from_target(target)
    cases = [build_case(path) for path in input_paths]

    output_path = Path(args.output).resolve() if args.output else default_output_path(target)
    tsv_path = Path(args.tsv_output).resolve() if args.tsv_output else default_tsv_path(output_path)

    plot_cases(
        cases=cases,
        output_path=output_path,
        tsv_path=tsv_path,
        samples=args.samples,
        r_min_factor=args.r_min_factor,
    )

    print(f"Wrote plot: {output_path}")
    print(f"Wrote data: {tsv_path}")
    print("Input files:")
    for case in sorted(cases, key=lambda item: item.alpha_input):
        print(
            f"  {case.input_path} | LJalpha={case.alpha_input:.3e}, "
            f"alpha_effective={case.alpha_effective:.3e}, r0={case.r0:.3e} m"
        )


if __name__ == "__main__":
    main()
