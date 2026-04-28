import math
import re
import subprocess
from pathlib import Path

try:
    from skopt import gp_minimize
    from skopt.space import Real
except ImportError:
    raise SystemExit(
        "Missing package. Install with:\n"
        "pip install scikit-optimize numpy"
    )

RUN_SCRIPT = "/home/wuwen/program/CabanaPD_Yao/scripts/automaticallyscanljparameter.sh"
LOG_FILE = Path("optimization_lj_log.csv")

BAD_LOSS = 1.0
FIXED_ALPHA = 1.5e-6

FAILED_BETA = []

def run_case(alpha: float, beta: float, vin: int) -> float:
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
        print("STDOUT:")
        print(result.stdout)
        print("STDERR:")
        print(result.stderr)
        return float("nan")

    return float(match.group(1))


def objective(x):
    #log10_alpha, beta = x

    #log10_alpha = round(log10_alpha, 1) 

    beta = x[0]
 
    beta = max(0.0001, round(beta, 5))
    alpha= FIXED_ALPHA
    #alpha = 10.0 ** log10_alpha
    for b_fail in FAILED_BETA:
        if abs(beta - b_fail) < 0.0002:
            print(f"Skip unstable beta={beta:.5f}, near failed beta={b_fail:.5f}")
            return BAD_LOSS
    cor100 = run_case(alpha, beta, 100)
    if math.isnan(cor100):
        FAILED_BETA.append(beta)
        return BAD_LOSS

    cor600 = run_case(alpha, beta, 600)

    if math.isnan(cor600):
        FAILED_BETA.append(beta)
        loss = 0.8
    else:
        loss_100 = max(0.0, 0.14 + cor100) ** 2
        loss_600 = min(0.0, -cor600 - 0.03) ** 2
        loss = loss_100 + loss_600

    line = (
        f"{alpha:.3e},{beta:.5f},"
        f"{cor100},{cor600},{loss:.5e}\n"
    )

    with LOG_FILE.open("a") as f:
        f.write(line)

    print(
        f"alpha={alpha:.3e}, beta={beta:.5f}, "
        f"CoR100={cor100}, CoR600={cor600}, "
        f"loss={loss:.5e}",
        flush=True,
    )

    return loss


def main():
    LOG_FILE.write_text("alpha,beta,cor100,cor600,loss\n")

    space = [
        #Real(-8.0, -5.0, name="log10_alpha"),  # alpha = 1e-8 ~ 1e-5
        Real(0.0001, 1.0, name="beta"),
    ]

    result = gp_minimize(
        objective,
        space,
        n_calls=30,
        n_initial_points=8,
        random_state=1,
    )

 #   best_log_alpha, best_beta = result.x
 #   best_log_alpha = round(best_log_alpha, 1)
    best_beta = round(result.x[0], 5)
    best_alpha = FIXED_ALPHA

    print("\nBest result:")
    print(f"alpha = {best_alpha:.3e}")
    print(f"beta  = {best_beta:.5f}")
    print(f"loss  = {result.fun:.5e}")
    print(f"log   = {LOG_FILE}")


if __name__ == "__main__":
    main()