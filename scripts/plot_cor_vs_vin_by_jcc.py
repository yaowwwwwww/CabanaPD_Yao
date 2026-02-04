#!/usr/bin/env python3
"""Plot CoR vs impact speed with one curve per jc_c."""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


CASE_RE = re.compile(
    r"^vin_(?P<vin>[^_]+)_alpha_(?P<alpha>[^_]+)_beta_(?P<beta>[^_]+)_A_(?P<A>[^_]+)_B_(?P<B>[^_]+)_C_(?P<C>[^_]+)$"
)


def parse_sort_key(value: str):
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def sanitize_filename(value: str) -> str:
    return (
        value.replace("+", "p")
        .replace("-", "m")
        .replace(".", "d")
        .replace("/", "_")
        .replace(" ", "_")
    )


def parse_file(path: Path) -> Dict[Tuple[str, str], Dict[str, List[Tuple[float, float]]]]:
    grouped: Dict[Tuple[str, str], Dict[str, List[Tuple[float, float]]]] = defaultdict(
        lambda: defaultdict(list)
    )

    with path.open("r", encoding="utf-8") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 4:
                print(f"[WARN] line {line_no}: expected >=4 columns, got {len(parts)}")
                continue

            case_name = parts[0]
            m = CASE_RE.match(case_name)
            if not m:
                print(f"[WARN] line {line_no}: case_name not matched: {case_name}")
                continue

            try:
                vin_col = float(parts[1])
                cor = float(parts[3])
            except ValueError:
                print(f"[WARN] line {line_no}: failed to parse vin/cor")
                continue

            impact_speed = abs(vin_col)
            if math.isnan(impact_speed):
                try:
                    impact_speed = abs(float(m.group("vin")))
                except ValueError:
                    print(f"[WARN] line {line_no}: failed to parse impact speed")
                    continue

            jc_a = m.group("A")
            jc_b = m.group("B")
            jc_c = m.group("C")
            grouped[(jc_a, jc_b)][jc_c].append((impact_speed, cor))

    return grouped


def _write_data_file(data_file: Path, points: List[Tuple[float, float]]) -> None:
    data_file.parent.mkdir(parents=True, exist_ok=True)
    with data_file.open("w", encoding="utf-8") as f:
        for speed, cor in points:
            if math.isnan(speed) or math.isnan(cor):
                continue
            f.write(f"{speed:.12g} {cor:.12g}\n")


def _run_gnuplot(script_text: str) -> None:
    subprocess.run(["gnuplot"], input=script_text, text=True, check=True)


def plot_grouped_data(
    grouped: Dict[Tuple[str, str], Dict[str, List[Tuple[float, float]]]], out_dir: Path
) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: List[Path] = []
    data_root = out_dir / "_curve_data"

    sorted_keys = sorted(
        grouped.keys(), key=lambda key: (parse_sort_key(key[0]), parse_sort_key(key[1]))
    )

    for jc_a, jc_b in sorted_keys:
        curve_specs: List[Tuple[str, Path]] = []
        for jc_c in sorted(grouped[(jc_a, jc_b)].keys(), key=parse_sort_key):
            points = sorted(grouped[(jc_a, jc_b)][jc_c], key=lambda item: item[0])
            data_file = (
                data_root
                / f"jcA_{sanitize_filename(jc_a)}"
                / f"jcB_{sanitize_filename(jc_b)}"
                / f"jcC_{sanitize_filename(jc_c)}.dat"
            )
            _write_data_file(data_file, points)
            if data_file.stat().st_size > 0:
                curve_specs.append((jc_c, data_file))

        if not curve_specs:
            continue

        out_file = (
            out_dir
            / f"cor_vs_vin_jcA_{sanitize_filename(jc_a)}_jcB_{sanitize_filename(jc_b)}_by_jcC.png"
        )
        gnuplot_curves = ", \\\n  ".join(
            f"'{path}' using 1:2 with linespoints lw 2 pt 7 ps 1 title 'jc_c={jc_c}'"
            for jc_c, path in curve_specs
        )
        gnuplot_script = f"""
set terminal pngcairo size 1400,900 enhanced font 'Arial,15'
set output '{out_file}'
set xlabel 'Impact speed (m/s)'
set ylabel 'CoR'
set title 'CoR vs impact speed (jc_a={jc_a}, jc_b={jc_b})'
set grid
set key left top
set yrange [0:*]
plot {gnuplot_curves}
"""
        _run_gnuplot(gnuplot_script)

        outputs.append(out_file)

    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot CoR-vs-impact-speed with curves for different jc_c."
    )
    parser.add_argument("summary_file", type=Path, help="Summary text file path")
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <summary_dir>/plots_cor_vs_vin_by_jcC)",
    )
    args = parser.parse_args()

    summary_file = args.summary_file.resolve()
    if not summary_file.is_file():
        raise FileNotFoundError(f"Input file not found: {summary_file}")

    out_dir = (
        args.out_dir.resolve()
        if args.out_dir is not None
        else summary_file.parent / "plots_cor_vs_vin_by_jcC"
    )

    grouped = parse_file(summary_file)
    if not grouped:
        raise RuntimeError("No valid data parsed from summary file.")

    if shutil.which("gnuplot") is None:
        raise RuntimeError("gnuplot not found. Please install gnuplot first.")

    outputs = plot_grouped_data(grouped, out_dir)
    print(f"Generated {len(outputs)} figure(s) in: {out_dir}")
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
