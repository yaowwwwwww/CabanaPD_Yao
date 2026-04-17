#!/usr/bin/env python3
import csv
import json
import math
import re
import sys
from pathlib import Path


def frame_id(path: Path) -> int:
    m = re.match(r"particles_(\d+)_all\.csv$", path.name)
    if not m:
        raise ValueError(f"unexpected file name: {path.name}")
    return int(m.group(1))


def mean_or_nan(values):
    if not values:
        return float("nan")
    return sum(values) / len(values)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: extract_interface_time_history.py <case_dir> <out_tsv>", file=sys.stderr)
        return 2

    case_dir = Path(sys.argv[1])
    out_tsv = Path(sys.argv[2])

    with (case_dir / "input.json").open("r", encoding="utf-8") as f:
        inp = json.load(f)

    dt = float(inp["timestep"]["value"])
    out_freq = int(inp["output_frequency"]["value"])

    csv_files = sorted(case_dir.glob("particles_*_all.csv"), key=frame_id)

    with out_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "frame",
                "time_s",
                "count",
                "mean_yield_stress_Pa",
                "mean_plastic_strain",
                "mean_plastic_strain_rate_1_per_s",
                "mean_temperature_K",
            ]
        )

        for csv_path in csv_files:
            frame = frame_id(csv_path)
            ys_vals = []
            eps_vals = []
            edot_vals = []
            temp_vals = []

            with csv_path.open("r", encoding="utf-8", newline="") as fp:
                rd = csv.DictReader(fp)
                for row in rd:
                    try:
                        typ = int(float(row["rank_0/type"]))
                        z = float(row["Points:2"])
                    except Exception:
                        continue

                    if typ != 0:
                        continue
                    if not (-2.0e-6 <= z <= 0.0):
                        continue

                    try:
                        ys = float(row["rank_0/yield_stress"])
                        eps = float(row["rank_0/plastic_strain"])
                        edot = float(row["rank_0/plastic_strain_rate"])
                        temp = float(row["rank_0/temperature"])
                    except Exception:
                        continue

                    if math.isfinite(ys):
                        ys_vals.append(ys)
                    if math.isfinite(eps):
                        eps_vals.append(eps)
                    if math.isfinite(edot):
                        edot_vals.append(edot)
                    if math.isfinite(temp):
                        temp_vals.append(temp)

            w.writerow(
                [
                    frame,
                    f"{frame * dt * out_freq:.8e}",
                    len(ys_vals),
                    f"{mean_or_nan(ys_vals):.8e}",
                    f"{mean_or_nan(eps_vals):.8e}",
                    f"{mean_or_nan(edot_vals):.8e}",
                    f"{mean_or_nan(temp_vals):.8e}",
                ]
            )

    print(out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
