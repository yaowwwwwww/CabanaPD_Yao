import math
import re
import subprocess
from pathlib import Path

RUN_SCRIPT = "/home/wuwen/program/CabanaPD_Yao/scripts/automaticallyscanljparameter.sh"
LOG_FILE = Path("coarse_lj_scan_log.csv")

ALPHAS = [
    1.0e-6,
]

BETAS = [
    0.5,
]


def run_case(alpha, beta, vin):
    cmd = [
        RUN_SCRIPT,
        f"{alpha:.2e}",
        f"{beta:.5f}",
        str(vin),
    ]

    result = subprocess.run(
        cmd,
        text=True,
        capture_output=True,
        check=False,
    )

    match = re.search(r"CoR\s*=\s*([0-9.eE+-]+|nan|NaN)", result.stdout)

    if not match:
        return float("nan")

    return float(match.group(1))


def main():
    LOG_FILE.write_text("alpha,beta,cor100,cor600,status\n")

    for alpha in ALPHAS:
        for beta in BETAS:
            print(f"\nRunning alpha={alpha:.2e}, beta={beta:.5f}", flush=True)

            cor100 = run_case(alpha, beta, 100)
            cor600 = run_case(alpha, beta, 600)

            if math.isnan(cor100) and math.isnan(cor600):
                status = "both_nan"
            elif math.isnan(cor100):
                status = "cor100_nan"
            elif math.isnan(cor600):
                status = "cor600_nan"
            else:
                status = "ok"

            line = f"{alpha:.3e},{beta:.5f},{cor100},{cor600},{status}\n"

            with LOG_FILE.open("a") as f:
                f.write(line)

            print(line.strip(), flush=True)


if __name__ == "__main__":
    main()