#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from pathlib import Path


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


def load_rows(
    case_dir: Path, projectile_type: int, last: int | None
) -> tuple[list[dict[str, float | int]], float]:
    input_path = case_dir / "input.json"
    if not input_path.exists():
        input_path = case_dir / "_current_input.json"
    if not input_path.exists():
        raise FileNotFoundError(f"No input.json or _current_input.json in {case_dir}")

    data = json.loads(input_path.read_text(encoding="utf-8"))
    dt_out = scalar(data["timestep"]) * int(read_value(data["output_frequency"]))

    files = sorted(case_dir.glob("particles_*_all.csv"), key=frame_id)
    if last is not None:
        if last <= 1:
            raise ValueError("--last must be greater than 1")
        files = files[-last:]
    if not files:
        raise FileNotFoundError(f"No particles_*_all.csv files in {case_dir}")

    frames = [frame_id(path) for path in files]
    coms = [projectile_com(path, projectile_type) for path in files]

    x0, y0, z0, _ = coms[0]
    path_nm = 0.0
    previous = None
    rows: list[dict[str, float | int]] = []

    for frame, (x, y, z, count) in zip(frames, coms):
        time_ns = frame * dt_out * 1.0e9
        dx_nm = (x - x0) * 1.0e9
        dy_nm = (y - y0) * 1.0e9
        dz_nm = (z - z0) * 1.0e9
        net_nm = math.sqrt(dx_nm * dx_nm + dy_nm * dy_nm + dz_nm * dz_nm)
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
                "displacement_nm": net_nm,
                "path_distance_nm": path_nm,
                "n_projectile": count,
            }
        )
        previous = (x, y, z)

    return rows, dt_out


def write_tsv(rows: list[dict[str, float | int]], path: Path) -> None:
    fieldnames = [
        "frame",
        "time_ns",
        "com_x_m",
        "com_y_m",
        "com_z_m",
        "displacement_nm",
        "path_distance_nm",
        "n_projectile",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def nice_bounds(values: list[float], include_zero: bool = True) -> tuple[float, float]:
    lo = min(values)
    hi = max(values)
    if include_zero:
        lo = min(lo, 0.0)
        hi = max(hi, 0.0)
    if abs(hi - lo) < 1.0e-12:
        lo -= 1.0
        hi += 1.0
    pad = 0.10 * (hi - lo)
    return lo - pad, hi + pad


def tick_values(lo: float, hi: float, count: int) -> list[float]:
    return [lo + idx * (hi - lo) / (count - 1) for idx in range(count)]


def format_tick(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def write_svg(
    rows: list[dict[str, float | int]],
    dt_out: float,
    case_label: str,
    path: Path,
) -> None:
    width, height = 1100, 500
    left, right = 105, 55
    top_margin, panel_height = 90, 300
    plot_width = width - left - right
    panel_top = top_margin

    frames = [int(row["frame"]) for row in rows]
    times = [float(row["time_ns"]) for row in rows]
    series = {
        "displacement": [float(row["displacement_nm"]) for row in rows],
        "path": [float(row["path_distance_nm"]) for row in rows],
    }
    colors = {
        "displacement": "#111111",
        "path": "#9467bd",
    }

    x_min, x_max = min(times), max(times)
    y_min, y_max = nice_bounds([v for values in series.values() for v in values])

    def scale_x(value: float) -> float:
        if x_max == x_min:
            return left + plot_width / 2.0
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def scale_y(value: float, y_min: float, y_max: float, panel_top: float) -> float:
        return panel_top + panel_height - (value - y_min) / (y_max - y_min) * panel_height

    def polyline(values: list[float], y_min: float, y_max: float, panel_top: float) -> str:
        return " ".join(
            f"{scale_x(x):.2f},{scale_y(y, y_min, y_max, panel_top):.2f}"
            for x, y in zip(times, values)
        )

    def add_panel(
        svg: list[str],
        panel_top: float,
        y_min: float,
        y_max: float,
        title: str,
        ylabel: str,
        series: dict[str, list[float]],
        labels: list[str],
    ) -> None:
        x0, y0 = left, panel_top
        x1, y1 = left + plot_width, panel_top + panel_height
        svg.append(
            f'<text x="{left}" y="{panel_top - 34}" class="panel-title">'
            f"{html.escape(title)}</text>"
        )
        for tick in tick_values(x_min, x_max, 6):
            x = scale_x(tick)
            svg.append(
                f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y1}" class="grid"/>'
            )
            svg.append(
                f'<text x="{x:.2f}" y="{y1 + 24}" class="tick" text-anchor="middle">'
                f"{format_tick(tick)}</text>"
            )
        for tick in tick_values(y_max, y_min, 5):
            y = scale_y(tick, y_min, y_max, panel_top)
            svg.append(
                f'<line x1="{x0}" y1="{y:.2f}" x2="{x1}" y2="{y:.2f}" class="grid"/>'
            )
            svg.append(
                f'<text x="{left - 14}" y="{y + 4:.2f}" class="tick" text-anchor="end">'
                f"{format_tick(tick)}</text>"
            )
        if y_min <= 0.0 <= y_max:
            zero_y = scale_y(0.0, y_min, y_max, panel_top)
            svg.append(
                f'<line x1="{x0}" y1="{zero_y:.2f}" x2="{x1}" y2="{zero_y:.2f}" '
                f'class="zero"/>'
            )
        svg.append(
            f'<rect x="{x0}" y="{y0}" width="{plot_width}" height="{panel_height}" '
            f'class="frame"/>'
        )
        svg.append(
            f'<text x="30" y="{panel_top + panel_height / 2}" class="ylabel" '
            f'transform="rotate(-90 30 {panel_top + panel_height / 2})" '
            f'text-anchor="middle">{html.escape(ylabel)}</text>'
        )

        legend_spacing = 145 if any(len(label) > 6 for label in labels) else 82
        legend_x = left + plot_width - legend_spacing * len(labels)
        legend_y = panel_top - 38
        for idx, key in enumerate(labels):
            x = legend_x + idx * legend_spacing
            svg.append(
                f'<line x1="{x}" y1="{legend_y}" x2="{x + 25}" y2="{legend_y}" '
                f'stroke="{colors[key]}" stroke-width="3"/>'
            )
            svg.append(
                f'<text x="{x + 32}" y="{legend_y + 4}" class="legend">'
                f"{html.escape(key)}</text>"
            )

        for key in labels:
            values = series[key]
            svg.append(
                f'<polyline points="{polyline(values, y_min, y_max, panel_top)}" '
                f'fill="none" stroke="{colors[key]}" stroke-width="2.6" '
                f'stroke-linejoin="round" stroke-linecap="round"/>'
            )
            for idx, (frame, time, value) in enumerate(zip(frames, times, values)):
                if idx == 0 or idx == len(frames) - 1 or frame % 5 == 0:
                    svg.append(
                        f'<circle cx="{scale_x(time):.2f}" '
                        f'cy="{scale_y(value, y_min, y_max, panel_top):.2f}" '
                        f'r="3.3" fill="white" stroke="{colors[key]}" stroke-width="1.8"/>'
                    )

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        """<style>
  .bg { fill: #fff; }
  .title { font: 700 21px Arial, sans-serif; fill: #111; }
  .subtitle { font: 13px Arial, sans-serif; fill: #444; }
  .panel-title { font: 700 15px Arial, sans-serif; fill: #111; }
  .tick, .legend { font: 12px Arial, sans-serif; fill: #333; }
  .ylabel, .xlabel { font: 13px Arial, sans-serif; fill: #222; }
  .grid { stroke: #d9d9d9; stroke-width: 1; }
  .zero { stroke: #555; stroke-width: 1.3; stroke-dasharray: 5 4; }
  .frame { fill: none; stroke: #222; stroke-width: 1.2; }
</style>""",
        f'<rect width="{width}" height="{height}" class="bg"/>',
        f'<text x="{width / 2}" y="34" class="title" text-anchor="middle">'
        f"Projectile COM displacement vs time</text>",
        f'<text x="{width / 2}" y="58" class="subtitle" text-anchor="middle">'
        f"case: {html.escape(case_label)}, frames {frames[0]}-{frames[-1]} "
        f"({len(frames)} frames), dt = {dt_out * 1.0e9:.3g} ns, "
        f"origin = COM at frame {frames[0]}</text>",
    ]
    add_panel(
        svg,
        panel_top,
        y_min,
        y_max,
        f"Combined COM displacement from frame {frames[0]}",
        "displacement (nm)",
        series,
        ["displacement", "path"],
    )
    svg.append(
        f'<text x="{left + plot_width / 2}" y="{height - 26}" class="xlabel" '
        f'text-anchor="middle">time (ns)</text>'
    )
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def default_prefix(last: int | None) -> str:
    if last is None:
        return "projectile_com_all_frames_move_vs_time"
    return f"projectile_com_last{last}_move_vs_time"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute projectile center-of-mass motion from particles_*_all.csv "
            "and write an SVG plot plus a TSV table."
        )
    )
    parser.add_argument("case_dir", type=Path, help="Case directory containing input.json and particle CSVs")
    parser.add_argument("--last", type=int, default=None, help="Plot only the last N frames")
    parser.add_argument("--projectile-type", type=int, default=0, help="Projectile rank_0/type value")
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output path prefix. Defaults to a file inside case_dir.",
    )
    parser.add_argument("--label", default=None, help="Case label for the plot title")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    rows, dt_out = load_rows(case_dir, args.projectile_type, args.last)

    if args.output_prefix is None:
        prefix = case_dir / default_prefix(args.last)
    else:
        prefix = Path(args.output_prefix)
        if not prefix.is_absolute():
            prefix = Path.cwd() / prefix

    prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = prefix.with_suffix(".tsv")
    svg_path = prefix.with_suffix(".svg")

    write_tsv(rows, tsv_path)
    write_svg(rows, dt_out, args.label or case_dir.name, svg_path)

    frames = [int(row["frame"]) for row in rows]
    print(f"frames={frames[0]}-{frames[-1]} n={len(frames)}")
    print(svg_path)
    print(tsv_path)


if __name__ == "__main__":
    main()
