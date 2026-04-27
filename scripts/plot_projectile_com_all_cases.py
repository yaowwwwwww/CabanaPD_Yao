#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import subprocess
from pathlib import Path


COLORS = [
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#9467bd",
    "#ff7f0e",
    "#17becf",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
]


def frame_id(path: Path) -> int:
    match = re.search(r"particles_(\d+)_all\.csv$", path.name)
    if not match:
        raise ValueError(f"Unexpected particle CSV filename: {path}")
    return int(match.group(1))


def read_value(node):
    if isinstance(node, dict) and "value" in node:
        return node["value"]
    return node


def scalar(node) -> float:
    value = read_value(node)
    if isinstance(value, list):
        return float(value[0])
    return float(value)


def parse_case_token(case_name: str, key: str) -> str | None:
    match = re.search(rf"{re.escape(key)}_(.+?)(?:_[A-Za-z]+_|$)", case_name)
    if match:
        return match.group(1)
    return None


def decode_token(token: str | None) -> str | None:
    if token is None:
        return None
    return token.replace("p", ".").replace("_", "-")


def parse_float_token(token: str | None) -> float | None:
    decoded = decode_token(token)
    if decoded is None:
        return None
    try:
        return float(decoded)
    except ValueError:
        return None


def projectile_com(csv_path: Path, projectile_type: int) -> tuple[float, float, float, int]:
    sx = sy = sz = 0.0
    count = 0
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if int(float(row["rank_0/type"])) == projectile_type:
                sx += float(row["Points:0"])
                sy += float(row["Points:1"])
                sz += float(row["Points:2"])
                count += 1
    if count == 0:
        raise RuntimeError(f"No projectile particles of type {projectile_type} in {csv_path}")
    return sx / count, sy / count, sz / count, count


def load_case(case_dir: Path, projectile_type: int) -> list[dict[str, float | int | str]]:
    input_path = case_dir / "input.json"
    if not input_path.exists():
        input_path = case_dir / "_current_input.json"
    if not input_path.exists():
        raise FileNotFoundError(f"No input.json or _current_input.json in {case_dir}")

    data = json.loads(input_path.read_text(encoding="utf-8"))
    dt_out = scalar(data["timestep"]) * int(read_value(data["output_frequency"]))

    files = sorted(case_dir.glob("particles_*_all.csv"), key=frame_id)
    if not files:
        raise FileNotFoundError(f"No particles_*_all.csv files in {case_dir}")

    coms = [projectile_com(path, projectile_type) for path in files]
    x0, y0, z0, _ = coms[0]
    previous = None
    path_nm = 0.0
    rows: list[dict[str, float | int | str]] = []

    for path, (x, y, z, count) in zip(files, coms):
        frame = frame_id(path)
        time_ns = frame * dt_out * 1.0e9
        dx_nm = (x - x0) * 1.0e9
        dy_nm = (y - y0) * 1.0e9
        dz_nm = (z - z0) * 1.0e9
        displacement_nm = math.sqrt(dx_nm * dx_nm + dy_nm * dy_nm + dz_nm * dz_nm)
        if previous is not None:
            sx_nm = (x - previous[0]) * 1.0e9
            sy_nm = (y - previous[1]) * 1.0e9
            sz_nm = (z - previous[2]) * 1.0e9
            path_nm += math.sqrt(sx_nm * sx_nm + sy_nm * sy_nm + sz_nm * sz_nm)

        rows.append(
            {
                "frame": frame,
                "time_ns": time_ns,
                "com_x_m": x,
                "com_y_m": y,
                "com_z_m": z,
                "displacement_nm": displacement_nm,
                "path_distance_nm": path_nm,
                "n_projectile": count,
                "case": case_dir.name,
            }
        )
        previous = (x, y, z)
    return rows


def choose_labels(case_dirs: list[Path]) -> tuple[list[str], list[float]]:
    case_names = [path.name.removeprefix("run_") for path in case_dirs]

    vin_tokens = [parse_case_token(name, "vin") for name in case_names]
    alpha_tokens = [parse_case_token(name, "alpha") for name in case_names]
    beta_tokens = [parse_case_token(name, "beta") for name in case_names]

    if len({token for token in vin_tokens if token is not None}) > 1:
        labels = [f"v{decode_token(token)}" for token in vin_tokens]
        sort_values = [parse_float_token(token) or float("inf") for token in vin_tokens]
        return labels, sort_values

    if len({token for token in alpha_tokens if token is not None}) > 1:
        labels = [f"alpha={decode_token(token)}" for token in alpha_tokens]
        sort_values = [parse_float_token(token) or float("inf") for token in alpha_tokens]
        return labels, sort_values

    if len({token for token in beta_tokens if token is not None}) > 1:
        labels = [f"beta={decode_token(token)}" for token in beta_tokens]
        sort_values = [parse_float_token(token) or float("inf") for token in beta_tokens]
        return labels, sort_values

    return case_names, [float(idx) for idx in range(len(case_names))]


def load_cor_values(summary_path: Path | None) -> dict[str, float]:
    if summary_path is None or not summary_path.exists():
        return {}

    cor_values: dict[str, float] = {}
    with summary_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 4:
                continue
            try:
                cor = float(parts[3])
            except ValueError:
                continue

            case_name = parts[0]
            cor_values[case_name] = cor
            cor_values[f"run_{case_name}"] = cor

    return cor_values


def find_cor_summary(run_root: Path) -> Path | None:
    matches = sorted(run_root.glob("summary_cor_hmax_recomputed*.txt"))
    if matches:
        return matches[-1]
    return None


def append_cor_to_labels(
    case_dirs: list[Path], labels: list[str], cor_values: dict[str, float]
) -> list[str]:
    if not cor_values:
        return labels
    out: list[str] = []
    for case_dir, label in zip(case_dirs, labels):
        cor = cor_values.get(case_dir.name)
        if cor is None:
            cor = cor_values.get(case_dir.name.removeprefix("run_"))
        if cor is None:
            out.append(label)
        else:
            out.append(f"{label} CoR={cor:.4g}")
    return out


def format_tick(value: float, decimals: int = 1) -> str:
    if abs(value) >= 1000:
        return f"{value:.0f}"
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.{decimals}f}"


def tick_values(lo: float, hi: float, count: int) -> list[float]:
    if count <= 1:
        return [lo]
    if abs(hi - lo) < 1.0e-12:
        return [lo for _ in range(count)]
    return [lo + idx * (hi - lo) / (count - 1) for idx in range(count)]


def filter_rows(
    rows: list[dict[str, float | int | str]],
    time_min: float,
    inclusive: bool,
) -> list[dict[str, float | int | str]]:
    if time_min <= 0.0:
        return rows
    if inclusive:
        return [row for row in rows if float(row["time_ns"]) >= time_min]
    return [row for row in rows if float(row["time_ns"]) > time_min]


def write_tsv(
    cases: list[tuple[Path, str, list[dict[str, float | int | str]]]],
    path: Path,
    time_min: float,
    inclusive: bool,
) -> int:
    fieldnames = [
        "case",
        "label",
        "frame",
        "time_ns",
        "com_x_m",
        "com_y_m",
        "com_z_m",
        "displacement_nm",
        "path_distance_nm",
        "n_projectile",
    ]
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for case_dir, label, rows in cases:
            for row in filter_rows(rows, time_min, inclusive):
                out = dict(row)
                out["case"] = case_dir.name
                out["label"] = label
                writer.writerow(out)
                count += 1
    return count


def svg_polyline_points(
    times: list[float],
    values: list[float],
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    left: float,
    plot_width: float,
    top: float,
    plot_height: float,
) -> str:
    def scale_x(value: float) -> float:
        if x_max == x_min:
            return left + plot_width / 2.0
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def scale_y(value: float) -> float:
        if y_max == y_min:
            return top + plot_height / 2.0
        return top + plot_height - (value - y_min) / (y_max - y_min) * plot_height

    return " ".join(
        f"{scale_x(t):.2f},{scale_y(v):.2f}" for t, v in zip(times, values)
    )


def render_png(svg_path: Path, png_path: Path) -> None:
    commands = [
        ["inkscape", str(svg_path), "--export-filename", str(png_path)],
        ["convert", str(svg_path), str(png_path)],
    ]
    for cmd in commands:
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return
        except Exception:
            continue
    raise RuntimeError(f"Failed to render PNG for {svg_path}")


def write_svg(
    cases: list[tuple[Path, str, list[dict[str, float | int | str]]]],
    run_root: Path,
    path: Path,
    projectile_type: int,
    time_min: float,
    inclusive: bool,
) -> int:
    width = 1220
    height = 697
    left = 105
    right = 45
    top = 168
    bottom = 74
    plot_width = width - left - right
    plot_height = height - top - bottom

    filtered_cases: list[tuple[str, list[dict[str, float | int | str]], str]] = []
    for idx, (_, label, rows) in enumerate(cases):
        selected = filter_rows(rows, time_min, inclusive)
        if selected:
            filtered_cases.append((label, selected, COLORS[idx % len(COLORS)]))

    if not filtered_cases:
        raise RuntimeError("No case data available for requested time window")

    all_times = [float(row["time_ns"]) for _, rows, _ in filtered_cases for row in rows]
    all_moves = [float(row["displacement_nm"]) for _, rows, _ in filtered_cases for row in rows]

    x_min = min(all_times)
    x_max = max(all_times)
    y_min = 0.0
    y_max = max(all_moves)
    if abs(x_max - x_min) < 1.0e-12:
        x_max = x_min + 1.0
    if y_max <= y_min:
        y_max = y_min + 1.0
    y_max *= 1.05

    def scale_x(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def scale_y(value: float) -> float:
        return top + plot_height - (value - y_min) / (y_max - y_min) * plot_height

    if time_min > 0.0:
        title = f"Projectile COM move vs time after {time_min:g} ns"
        cmp_text = ">=" if inclusive else ">"
        subtitle = (
            f"{run_root.name} | {len(filtered_cases)} plotted cases | "
            f"time_ns {cmp_text} {time_min:g} | y auto-scaled to post-cutoff data"
        )
    else:
        title = "Projectile COM move vs time"
        subtitle = (
            f"{run_root.name} | {len(filtered_cases)} cases | projectile type {projectile_type} | "
            f"move = net COM displacement from first frame"
        )

    x_ticks = tick_values(x_min, x_max, 6 if time_min > 0.0 else 7)
    y_ticks = tick_values(y_min, y_max, 6)

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        """<style>
  .bg { fill: #fff; }
  .title { font-family: Arial, sans-serif; font-size: 21px; font-weight: 700; fill: #111; }
  .subtitle { font-family: Arial, sans-serif; font-size: 13px; font-weight: 400; fill: #444; }
  .tick, .legend, .note { font-family: Arial, sans-serif; font-size: 12px; font-weight: 400; fill: #333; }
  .axis-label { font-family: Arial, sans-serif; font-size: 14px; font-weight: 400; fill: #222; }
  .grid { stroke: #dddddd; stroke-width: 1; }
  .frame { fill: none; stroke: #222; stroke-width: 1.2; }
</style>""",
        '<rect width="100%" height="100%" class="bg"/>',
        f'<text x="{width / 2}" y="32" text-anchor="middle" class="title">{html.escape(title)}</text>',
        f'<text x="{width / 2}" y="56" text-anchor="middle" class="subtitle">{html.escape(subtitle)}</text>',
    ]

    legend_cols = 4
    legend_x0 = 105.0
    legend_dx = 265.0
    legend_y0 = 96.0 if time_min > 0.0 else 84.0
    for idx, (label, _, color) in enumerate(filtered_cases):
        col = idx % legend_cols
        row = idx // legend_cols
        x0 = legend_x0 + col * legend_dx
        y0 = legend_y0 + row * 24.0
        svg.append(
            f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x0 + 28:.1f}" y2="{y0:.1f}" '
            f'stroke="{color}" stroke-width="3"/>'
        )
        svg.append(
            f'<text x="{x0 + 35:.1f}" y="{y0 + 4:.1f}" class="legend">{html.escape(label)}</text>'
        )

    for tick in x_ticks:
        x = scale_x(tick)
        svg.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" class="grid"/>')
        svg.append(
            f'<text x="{x:.2f}" y="{top + plot_height + 25}" text-anchor="middle" class="tick">'
            f"{format_tick(tick, 1)}</text>"
        )

    for tick in y_ticks:
        y = scale_y(tick)
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" class="grid"/>')
        svg.append(
            f'<text x="{left - 14}" y="{y + 4:.2f}" text-anchor="end" class="tick">'
            f"{format_tick(tick, 0)}</text>"
        )

    svg.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" class="frame"/>'
    )

    for label, rows, color in filtered_cases:
        times = [float(row["time_ns"]) for row in rows]
        values = [float(row["displacement_nm"]) for row in rows]
        svg.append(
            f'<polyline points="{svg_polyline_points(times, values, x_min, x_max, y_min, y_max, left, plot_width, top, plot_height)}" '
            f'fill="none" stroke="{color}" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round"/>'
        )
        for time, value in [(times[0], values[0]), (times[-1], values[-1])]:
            svg.append(
                f'<circle cx="{scale_x(time):.2f}" cy="{scale_y(value):.2f}" r="3.4" '
                f'fill="white" stroke="{color}" stroke-width="1.5"/>'
            )

    svg.append(
        f'<text x="{left + plot_width / 2}" y="{height - 30}" text-anchor="middle" class="axis-label">time (ns)</text>'
    )
    svg.append(
        f'<text x="34" y="{top + plot_height / 2}" text-anchor="middle" class="axis-label" '
        f'transform="rotate(-90 34 {top + plot_height / 2})">COM move / net displacement (nm)</text>'
    )
    svg.append("</svg>")

    path.write_text("\n".join(svg) + "\n", encoding="utf-8")
    return len(filtered_cases)


def default_prefix(run_root: Path, time_min: float) -> Path:
    if time_min <= 0.0:
        return run_root / "projectile_com_all_cases_move_vs_time"
    if abs(time_min - round(time_min)) < 1.0e-9:
        return run_root / f"projectile_com_all_cases_move_vs_time_after{int(round(time_min))}ns"
    safe = str(time_min).replace(".", "p")
    return run_root / f"projectile_com_all_cases_move_vs_time_after{safe}ns"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot projectile COM move vs time for all run_* cases under a run root. "
            "Writes SVG, PNG, and TSV."
        )
    )
    parser.add_argument("run_root", type=Path, help="Directory containing run_* case directories")
    parser.add_argument("--projectile-type", type=int, default=0, help="Projectile rank_0/type value")
    parser.add_argument("--time-min", type=float, default=0.0, help="Only include rows after this time (ns)")
    parser.add_argument(
        "--inclusive",
        action="store_true",
        help="Use time_ns >= time_min instead of time_ns > time_min",
    )
    parser.add_argument("--output-prefix", default=None, help="Optional explicit output path prefix")
    parser.add_argument(
        "--label-cor",
        action="store_true",
        help="Append CoR values from summary_cor_hmax_recomputed*.txt to legend labels",
    )
    parser.add_argument(
        "--cor-summary",
        type=Path,
        default=None,
        help="Explicit summary file to use for --label-cor",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_root = args.run_root.resolve()

    case_dirs = sorted(
        [path for path in run_root.iterdir() if path.is_dir() and path.name.startswith("run_")]
    )
    if not case_dirs:
        raise FileNotFoundError(f"No run_* directories found in {run_root}")

    labels, sort_values = choose_labels(case_dirs)
    if args.label_cor:
        summary_path = args.cor_summary.resolve() if args.cor_summary else find_cor_summary(run_root)
        labels = append_cor_to_labels(case_dirs, labels, load_cor_values(summary_path))
    combined = sorted(zip(sort_values, case_dirs, labels), key=lambda item: item[0])

    cases: list[tuple[Path, str, list[dict[str, float | int | str]]]] = []
    for _, case_dir, label in combined:
        rows = load_case(case_dir, args.projectile_type)
        cases.append((case_dir, label, rows))

    prefix = Path(args.output_prefix).resolve() if args.output_prefix else default_prefix(run_root, args.time_min)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    svg_path = prefix.with_suffix(".svg")
    png_path = prefix.with_suffix(".png")
    tsv_path = prefix.with_suffix(".tsv")

    write_tsv(cases, tsv_path, args.time_min, args.inclusive)
    plotted_cases = write_svg(cases, run_root, svg_path, args.projectile_type, args.time_min, args.inclusive)
    render_png(svg_path, png_path)

    print(f"cases={plotted_cases}")
    print(svg_path)
    print(png_path)
    print(tsv_path)


if __name__ == "__main__":
    main()
