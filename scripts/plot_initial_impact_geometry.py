#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import json
import math
import shutil
import subprocess
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


def vector(node) -> list[float]:
    value = read_value(node)
    if not isinstance(value, list):
        raise TypeError(f"Expected list value, got {type(value)!r}")
    return [float(item) for item in value]


def fmt_um(value_m: float) -> str:
    return f"{value_m * 1.0e6:.2f} um"


def fmt_ms(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value:.0f} m/s"


def ellipsize(text: str, max_len: int = 92) -> str:
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def tick_values(lo: float, hi: float, count: int) -> list[float]:
    if count <= 1 or abs(hi - lo) < 1.0e-15:
        return [lo]
    return [lo + idx * (hi - lo) / (count - 1) for idx in range(count)]


def fmt_tick_um(value_m: float) -> str:
    value_um = value_m * 1.0e6
    if abs(value_um) >= 10:
        return f"{value_um:.0f}"
    return f"{value_um:.1f}"


def draw_dimension(svg: list[str], x0: float, y0: float, x1: float, y1: float, label: str) -> None:
    dx = x1 - x0
    dy = y1 - y0
    length = math.hypot(dx, dy)
    if length < 1.0e-12:
        return
    ux = dx / length
    uy = dy / length
    px = -uy
    py = ux
    arrow = 8.0
    offset = 6.0
    sx = x0 + px * offset
    sy = y0 + py * offset
    ex = x1 + px * offset
    ey = y1 + py * offset
    svg.append(
        f'<line x1="{sx:.2f}" y1="{sy:.2f}" x2="{ex:.2f}" y2="{ey:.2f}" class="dim"/>'
    )
    for ax, ay, sign in ((sx, sy, 1.0), (ex, ey, -1.0)):
        bx = ax + sign * ux * arrow + px * 4.0
        by = ay + sign * uy * arrow + py * 4.0
        cx = ax + sign * ux * arrow - px * 4.0
        cy = ay + sign * uy * arrow - py * 4.0
        svg.append(
            f'<polygon points="{ax:.2f},{ay:.2f} {bx:.2f},{by:.2f} {cx:.2f},{cy:.2f}" class="dim-fill"/>'
        )
    lx = (sx + ex) / 2.0 + px * 14.0
    ly = (sy + ey) / 2.0 + py * 14.0
    svg.append(
        f'<text x="{lx:.2f}" y="{ly:.2f}" text-anchor="middle" class="note">'
        f"{html.escape(label)}</text>"
    )


def panel_frame(svg: list[str], x: float, y: float, w: float, h: float, title: str) -> None:
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" class="panel"/>')
    svg.append(
        f'<text x="{x}" y="{y - 14}" class="panel-title">{html.escape(title)}</text>'
    )


def draw_axes_and_grid(
    svg: list[str],
    x: float,
    y: float,
    w: float,
    h: float,
    x_min: float,
    x_max: float,
    z_min: float,
    z_max: float,
    sx,
    sz,
) -> None:
    for xv in tick_values(x_min, x_max, 7):
        px = sx(xv)
        svg.append(f'<line x1="{px:.2f}" y1="{y}" x2="{px:.2f}" y2="{y + h}" class="grid"/>')
        svg.append(
            f'<text x="{px:.2f}" y="{y + h + 22}" text-anchor="middle" class="tick">'
            f"{fmt_tick_um(xv)}</text>"
        )
    for zv in tick_values(z_min, z_max, 7):
        pz = sz(zv)
        svg.append(f'<line x1="{x}" y1="{pz:.2f}" x2="{x + w}" y2="{pz:.2f}" class="grid"/>')
        svg.append(
            f'<text x="{x - 12}" y="{pz + 4:.2f}" text-anchor="end" class="tick">'
            f"{fmt_tick_um(zv)}</text>"
        )
    svg.append(
        f'<text x="{x + w / 2:.2f}" y="{y + h + 44}" text-anchor="middle" class="axis-label">x (um)</text>'
    )
    svg.append(
        f'<text x="{x - 58:.2f}" y="{y + h / 2:.2f}" text-anchor="middle" class="axis-label" '
        f'transform="rotate(-90 {x - 58:.2f} {y + h / 2:.2f})">z (um)</text>'
    )


def draw_geometry(svg: list[str], sx, sz, x_min, x_max, plate_z_min, plate_z_max, cx, cz, radius) -> None:
    left = sx(x_min)
    right = sx(x_max)
    plate_top = sz(plate_z_max)
    plate_bottom = sz(plate_z_min)
    svg.append(
        f'<rect x="{left:.2f}" y="{plate_top:.2f}" width="{right - left:.2f}" '
        f'height="{plate_bottom - plate_top:.2f}" class="plate"/>'
    )
    svg.append(
        f'<circle cx="{sx(cx):.2f}" cy="{sz(cz):.2f}" r="{abs(sx(cx + radius) - sx(cx)):.2f}" class="projectile"/>'
    )
    svg.append(
        f'<line x1="{sx(cx):.2f}" y1="{sz(plate_z_min):.2f}" x2="{sx(cx):.2f}" y2="{sz(cz + radius):.2f}" class="centerline"/>'
    )


def write_svg(
    input_path: Path,
    out_path: Path,
    title: str,
    subtitle: str,
    low_corner: list[float],
    high_corner: list[float],
    num_cells: list[float],
    ball_center: list[float],
    ball_radius: float,
    plate_z_min: float,
    plate_z_max: float,
    horizon: float,
    vin: float | None,
) -> None:
    width = 1380
    height = 860
    left_margin = 105
    top_margin = 120
    panel_size = 520
    gap_between = 90
    right_x = left_margin + panel_size + gap_between
    bottom_y = top_margin + panel_size

    x_min, x_max = low_corner[0], high_corner[0]
    z_min, z_max = low_corner[2], high_corner[2]
    cx, _, cz = ball_center
    radius = ball_radius
    gap = cz - radius - plate_z_max
    diameter = 2.0 * radius
    plate_thickness = plate_z_max - plate_z_min
    domain_dx = high_corner[0] - low_corner[0]
    domain_dz = high_corner[2] - low_corner[2]
    cell_dx = domain_dx / num_cells[0]

    zoom_x_min = cx - 18.0e-6
    zoom_x_max = cx + 18.0e-6
    zoom_z_min = plate_z_min - 2.0e-6
    zoom_z_max = cz + radius + 6.0e-6

    def make_scale(px0, py0, pw, ph, ax_min, ax_max, az_min, az_max):
        def sx(value: float) -> float:
            return px0 + (value - ax_min) / (ax_max - ax_min) * pw

        def sz(value: float) -> float:
            return py0 + ph - (value - az_min) / (az_max - az_min) * ph

        return sx, sz

    full_sx, full_sz = make_scale(left_margin, top_margin, panel_size, panel_size, x_min, x_max, z_min, z_max)
    zoom_sx, zoom_sz = make_scale(right_x, top_margin, panel_size, panel_size, zoom_x_min, zoom_x_max, zoom_z_min, zoom_z_max)

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        """<style>
  .bg { fill: #ffffff; }
  .title { font-family: Arial, sans-serif; font-size: 24px; font-weight: 700; fill: #111; }
  .subtitle { font-family: Arial, sans-serif; font-size: 13px; font-weight: 400; fill: #444; }
  .panel-title { font-family: Arial, sans-serif; font-size: 16px; font-weight: 700; fill: #111; }
  .tick, .note, .legend, .meta { font-family: Arial, sans-serif; font-size: 12px; font-weight: 400; fill: #333; }
  .axis-label { font-family: Arial, sans-serif; font-size: 14px; font-weight: 400; fill: #222; }
  .panel { fill: none; stroke: #222; stroke-width: 1.15; }
  .grid { stroke: #dddddd; stroke-width: 1; }
  .plate { fill: #a6cee3; stroke: #3b6d8e; stroke-width: 1.4; }
  .projectile { fill: #f9c784; stroke: #b5621d; stroke-width: 1.6; }
  .centerline { stroke: #666; stroke-width: 1.2; stroke-dasharray: 5 4; }
  .arrow { stroke: #cc2b2b; stroke-width: 2.8; }
  .arrow-fill { fill: #cc2b2b; }
  .dim { stroke: #444; stroke-width: 1.4; }
  .dim-fill { fill: #444; }
  .meta-box { fill: #f7f7f7; stroke: #c8c8c8; stroke-width: 1; }
</style>""",
        f'<rect width="{width}" height="{height}" class="bg"/>',
        f'<text x="{width / 2:.1f}" y="36" text-anchor="middle" class="title">{html.escape(title)}</text>',
        f'<text x="{width / 2:.1f}" y="60" text-anchor="middle" class="subtitle">'
        f"{html.escape(subtitle)}</text>",
    ]

    panel_frame(svg, left_margin, top_margin, panel_size, panel_size, "Full Domain (x-z cross-section)")
    panel_frame(svg, right_x, top_margin, panel_size, panel_size, "Zoom Near Impact Region")

    draw_axes_and_grid(
        svg,
        left_margin,
        top_margin,
        panel_size,
        panel_size,
        x_min,
        x_max,
        z_min,
        z_max,
        full_sx,
        full_sz,
    )
    draw_axes_and_grid(
        svg,
        right_x,
        top_margin,
        panel_size,
        panel_size,
        zoom_x_min,
        zoom_x_max,
        zoom_z_min,
        zoom_z_max,
        zoom_sx,
        zoom_sz,
    )

    draw_geometry(svg, full_sx, full_sz, x_min, x_max, plate_z_min, plate_z_max, cx, cz, radius)
    draw_geometry(svg, zoom_sx, zoom_sz, zoom_x_min, zoom_x_max, plate_z_min, plate_z_max, cx, cz, radius)

    # Initial velocity arrow in zoom panel.
    arrow_x = zoom_sx(cx)
    arrow_y0 = zoom_sz(cz + radius + 4.5e-6)
    arrow_y1 = zoom_sz(cz + radius + 0.8e-6)
    svg.append(f'<line x1="{arrow_x:.2f}" y1="{arrow_y0:.2f}" x2="{arrow_x:.2f}" y2="{arrow_y1:.2f}" class="arrow"/>')
    svg.append(
        f'<polygon points="{arrow_x:.2f},{arrow_y1:.2f} {arrow_x - 7:.2f},{arrow_y1 - 12:.2f} {arrow_x + 7:.2f},{arrow_y1 - 12:.2f}" class="arrow-fill"/>'
    )
    svg.append(
        f'<text x="{arrow_x + 14:.2f}" y="{arrow_y0 - 4:.2f}" class="note">initial velocity = {html.escape(fmt_ms(vin))}</text>'
    )

    # Labels in zoom panel.
    svg.append(
        f'<text x="{zoom_sx(cx + radius) + 16:.2f}" y="{zoom_sz(cz) - 6:.2f}" class="note">projectile</text>'
    )
    svg.append(
        f'<text x="{zoom_sx(cx + radius) + 16:.2f}" y="{zoom_sz(-8.0e-6):.2f}" class="note">substrate plate</text>'
    )

    # Dimensions in zoom panel.
    draw_dimension(
        svg,
        zoom_sx(cx),
        zoom_sz(cz - radius),
        zoom_sx(cx),
        zoom_sz(cz + radius),
        f"diameter = {fmt_um(diameter)}",
    )
    draw_dimension(
        svg,
        zoom_sx(cx - 12.0e-6),
        zoom_sz(plate_z_max),
        zoom_sx(cx - 12.0e-6),
        zoom_sz(cz - radius),
        f"initial gap = {fmt_um(gap)}",
    )
    draw_dimension(
        svg,
        zoom_sx(cx + 14.0e-6),
        zoom_sz(plate_z_min),
        zoom_sx(cx + 14.0e-6),
        zoom_sz(plate_z_max),
        f"plate thickness = {fmt_um(plate_thickness)}",
    )

    # Info box.
    info_x = left_margin
    info_y = bottom_y + 90
    info_w = width - 2 * left_margin
    info_h = 150
    svg.append(f'<rect x="{info_x}" y="{info_y}" width="{info_w}" height="{info_h}" rx="6" ry="6" class="meta-box"/>')

    rows = [
        ("domain x", f"{fmt_um(x_min)} to {fmt_um(x_max)}"),
        ("domain z", f"{fmt_um(z_min)} to {fmt_um(z_max)}"),
        ("grid", f"{int(num_cells[0])} x {int(num_cells[1])} x {int(num_cells[2])} cells"),
        ("cell size", fmt_um(cell_dx)),
        ("ball center", f"(x,z) = ({fmt_um(cx)}, {fmt_um(cz)})"),
        ("ball radius", fmt_um(radius)),
        ("plate top", fmt_um(plate_z_max)),
        ("plate bottom", fmt_um(plate_z_min)),
        ("horizon", fmt_um(horizon)),
        ("applies to scan", "geometry same across all cases in this run set"),
    ]

    col_x = [info_x + 24, info_x + 430, info_x + 820]
    row_y = [info_y + 34, info_y + 60, info_y + 86, info_y + 112]
    for idx, (key, value) in enumerate(rows):
        cx_pos = col_x[idx % 3]
        cy_pos = row_y[idx // 3]
        svg.append(f'<text x="{cx_pos}" y="{cy_pos}" class="meta">{html.escape(key)}:</text>')
        svg.append(f'<text x="{cx_pos + 110}" y="{cy_pos}" class="meta">{html.escape(value)}</text>')

    svg.append("</svg>")
    out_path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw an initial geometry schematic from an impact-case input.json."
    )
    parser.add_argument("input_json", type=Path, help="Path to input.json")
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=None,
        help="Output file prefix. Defaults next to input_json.",
    )
    parser.add_argument(
        "--title",
        default="Initial Impact Geometry Setup",
        help="Figure title",
    )
    parser.add_argument(
        "--subtitle",
        default=None,
        help="Optional subtitle. Defaults to a short representative-case description.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input_json.resolve()
    data = json.loads(input_path.read_text(encoding="utf-8"))

    low_corner = vector(data["low_corner"])
    high_corner = vector(data["high_corner"])
    num_cells = vector(data["num_cells"])
    ball_center = vector(data["ball_center"])
    ball_radius = scalar(data["ball_radius"])
    plate_z_min = scalar(data["plate_z_min"])
    plate_z_max = scalar(data["plate_z_max"])
    horizon = scalar(data["horizon"])
    vin = scalar(data["ball_initial_velocity"]) if "ball_initial_velocity" in data else None

    if args.output_prefix is None:
        prefix = input_path.parent / "initial_impact_geometry_schematic"
    else:
        prefix = args.output_prefix.resolve()

    if args.subtitle is None:
        scan_name = input_path.parent.parent.name if input_path.parent.parent != input_path.parent else input_path.parent.name
        subtitle = (
            f"Representative case: v0 = {fmt_ms(vin)} | scan: {scan_name} | geometry identical across cases"
        )
    else:
        subtitle = args.subtitle

    prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = prefix.with_suffix(".svg")
    png_path = prefix.with_suffix(".png")

    write_svg(
        input_path,
        svg_path,
        args.title,
        ellipsize(subtitle),
        low_corner,
        high_corner,
        num_cells,
        ball_center,
        ball_radius,
        plate_z_min,
        plate_z_max,
        horizon,
        vin,
    )

    if shutil.which("convert"):
        subprocess.run(["convert", str(svg_path), str(png_path)], check=True)
        print(png_path)
    print(svg_path)


if __name__ == "__main__":
    main()
