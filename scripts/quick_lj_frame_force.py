#!/usr/bin/env python3
import argparse
import csv
import heapq
import json
import math
from pathlib import Path


def load_input(input_path: Path):
    with input_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def frame_csv_path(run_dir: Path, frame: int) -> Path:
    return run_dir / f"particles_{frame}_all.csv"


def compute_constants(inp):
    low = inp["low_corner"]["value"]
    high = inp["high_corner"]["value"]
    nc = inp["num_cells"]["value"]

    dx = (high[0] - low[0]) / nc[0]
    vol = dx ** 3

    E = float(inp["elastic_modulus"]["value"][0])
    nu = float(inp["Poisson's_ratio"]["value"][0])
    K = E / (3.0 * (1.0 - 2.0 * nu))
    delta = float(inp["horizon"]["value"])

    beta = float(inp["LJbeta"]["value"])
    alpha1 = float(inp["LJalpha"]["value"])
    r0 = float(inp["LJr0"]["value"]) * dx
    radius = float(inp["contact_horizon_factor"]["value"]) * r0

    c_czm = float(inp["CZM_cohesive_scaling"]["value"])
    sy = float(inp["CZM_yield_stretch"]["value"])
    m_czm = float(inp["CZM_degradation_rate"]["value"])

    c = 18.0 * K / (math.pi * (delta ** 4) * alpha1)
    alpha_eff = c * r0 * r0 * vol * vol / 72.0 / (beta ** (7.0 / 3.0))

    return {
        "low": low,
        "high": high,
        "vol": vol,
        "beta": beta,
        "r0": r0,
        "radius": radius,
        "alpha_eff": alpha_eff,
        "c_czm": c_czm,
        "sy": sy,
        "m_czm": m_czm,
    }


def lj_coeff_and_fc_total(r, const):
    if r <= 1e-14 or r > const["radius"]:
        return 0.0, 0.0

    r0 = const["r0"]
    beta = const["beta"]

    term13 = (r0 / r) ** 13
    term7 = (r0 / r) ** 7
    fc = ((12.0 * const["alpha_eff"]) / r0) * (term13 - beta * term7)

    s = (r - 2.0e-6) / 2.0e-6
    f_czm = 0.0
    if s < 0.0 and s >= -const["sy"]:
        f_czm = const["c_czm"] * (-s)
    elif s < -const["sy"]:
        f_czm = const["c_czm"] * const["sy"] * math.exp(-const["m_czm"] * (-s - const["sy"]))

    fc_total = fc - f_czm
    coeff = fc_total / const["vol"]
    return coeff, fc_total


def load_particles(csv_path: Path):
    # Keep only what we need: row_index, type, x, y, z
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rd = csv.DictReader(f)
        rows = []
        for idx, r in enumerate(rd):
            rows.append(
                (
                    idx,
                    int(float(r["rank_0/type"])),
                    float(r["Points:0"]),
                    float(r["Points:1"]),
                    float(r["Points:2"]),
                )
            )
    return rows


def bin_index(x, y, z, low, h):
    ix = int(math.floor((x - low[0]) / h))
    iy = int(math.floor((y - low[1]) / h))
    iz = int(math.floor((z - low[2]) / h))
    return ix, iy, iz


def analyze_frame(run_dir: Path, frame: int, topk: int):
    inp = load_input(run_dir / "input.json")
    const = compute_constants(inp)
    csv_path = frame_csv_path(run_dir, frame)
    pts = load_particles(csv_path)

    low = const["low"]
    h = const["radius"]
    r2cut = h * h

    substrate = []
    projectile = []
    for row in pts:
        if row[1] == 0:
            substrate.append(row)
        elif row[1] == 1:
            projectile.append(row)

    bins = {}
    for row in substrate:
        _, _, x, y, z = row
        key = bin_index(x, y, z, low, h)
        bins.setdefault(key, []).append(row)

    n_active = 0
    n_pos = 0
    n_neg = 0
    coeff_sum = 0.0
    max_coeff = -float("inf")
    min_coeff = float("inf")
    max_abs_coeff = 0.0
    r_min = float("inf")
    r_max = 0.0

    # Heap of (abs_coeff, row_i, row_j, r, coeff, fc_total)
    heap = []

    for irow in projectile:
        i_idx, _, xi, yi, zi = irow
        bix, biy, biz = bin_index(xi, yi, zi, low, h)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = (bix + dx, biy + dy, biz + dz)
                    cand = bins.get(key)
                    if not cand:
                        continue
                    for jrow in cand:
                        j_idx, _, xj, yj, zj = jrow
                        rx = xj - xi
                        ry = yj - yi
                        rz = zj - zi
                        rr2 = rx * rx + ry * ry + rz * rz
                        if rr2 <= 1e-28 or rr2 >= r2cut:
                            continue

                        r = math.sqrt(rr2)
                        coeff, fc_total = lj_coeff_and_fc_total(r, const)
                        n_active += 1
                        coeff_sum += coeff
                        if coeff > 0.0:
                            n_pos += 1
                        elif coeff < 0.0:
                            n_neg += 1

                        if coeff > max_coeff:
                            max_coeff = coeff
                        if coeff < min_coeff:
                            min_coeff = coeff
                        ac = abs(coeff)
                        if ac > max_abs_coeff:
                            max_abs_coeff = ac
                        if r < r_min:
                            r_min = r
                        if r > r_max:
                            r_max = r

                        item = (ac, i_idx, j_idx, r, coeff, fc_total)
                        if len(heap) < topk:
                            heapq.heappush(heap, item)
                        elif ac > heap[0][0]:
                            heapq.heapreplace(heap, item)

    if n_active == 0:
        avg_coeff = 0.0
        max_coeff = 0.0
        min_coeff = 0.0
        r_min = float("nan")
        r_max = 0.0
    else:
        avg_coeff = coeff_sum / n_active

    top_pairs = sorted(heap, key=lambda x: x[0], reverse=True)
    return {
        "frame": frame,
        "csv_path": str(csv_path),
        "n_active": n_active,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "avg_coeff": avg_coeff,
        "max_coeff": max_coeff,
        "min_coeff": min_coeff,
        "max_abs_coeff": max_abs_coeff,
        "r_min": r_min,
        "r_max": r_max,
        "radius": const["radius"],
        "r0": const["r0"],
        "r_star": const["r0"] / (const["beta"] ** (1.0 / 6.0)),
        "top_pairs": top_pairs,
    }


def write_top_pairs_tsv(out_path: Path, result):
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["frame", "row_i(type1)", "row_j(type0)", "r_m", "coeff_signed", "fc_total_signed", "abs_coeff"])
        for ac, i_idx, j_idx, r, coeff, fc_total in result["top_pairs"]:
            w.writerow([result["frame"], i_idx, j_idx, f"{r:.8e}", f"{coeff:.8e}", f"{fc_total:.8e}", f"{ac:.8e}"])


def main():
    ap = argparse.ArgumentParser(description="Quick LJ contact-force check for one frame.")
    ap.add_argument("--run-dir", required=True, help="Case directory containing input.json and particles_<frame>_all.csv")
    ap.add_argument("--frame", required=True, type=int, help="Frame index N for particles_N_all.csv")
    ap.add_argument("--topk", type=int, default=20, help="Number of strongest pairs to keep")
    ap.add_argument("--out-tsv", default="", help="Optional output TSV path for strongest pairs")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    result = analyze_frame(run_dir, args.frame, args.topk)

    print(f"frame={result['frame']} file={result['csv_path']}")
    print(f"n_active={result['n_active']} n_pos={result['n_pos']} n_neg={result['n_neg']}")
    print(
        "avg_coeff={:.8e} max_coeff={:.8e} min_coeff={:.8e} max_abs_coeff={:.8e}".format(
            result["avg_coeff"], result["max_coeff"], result["min_coeff"], result["max_abs_coeff"]
        )
    )
    print("r_min={:.8e} r_max={:.8e} radius={:.8e} r0={:.8e} r_star={:.8e}".format(
        result["r_min"], result["r_max"], result["radius"], result["r0"], result["r_star"]
    ))

    if args.out_tsv:
        out_path = Path(args.out_tsv)
    else:
        out_path = run_dir / f"lj_top_pairs_frame_{args.frame}.tsv"
    write_top_pairs_tsv(out_path, result)
    print(f"wrote_top_pairs={out_path}")


if __name__ == "__main__":
    main()
