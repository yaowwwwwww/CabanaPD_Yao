#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"
BASE_DIR="${ROOT_DIR}/build"
DRIVER="${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"

VIN="${VIN:-200}"
ALPHA="${ALPHA:-1.5e-6}"
JC_A="${JC_A:-9e7}"
JC_B="${JC_B:-2.92e8}"
JC_C="${JC_C:-0.025}"
JC_M="${JC_M:-1.09}"
FINAL_TIME="${FINAL_TIME:-8e-8}"
TQ="${TQ:-0}"
OUT_FREQ="${OUT_FREQ:-200}"
DT="${DT:-1e-11}"

LOWER_NAN_BETA="${LOWER_NAN_BETA:-0.11}"
UPPER_ZERO_BETA="${UPPER_ZERO_BETA:-0.115}"
INITIAL_BETAS="${INITIAL_BETAS:-0.111 0.112 0.113 0.114}"
MAX_ITER="${MAX_ITER:-20}"
MIN_INTERVAL="${MIN_INTERVAL:-1e-6}"

run_scan() {
  local tag="$1"
  local betas="$2"
  echo "=== run_scan tag=${tag} betas=${betas} ==="
  (
    cd "${ROOT_DIR}"
    RUN_TAG_BASE="${tag}" \
    CASE_MATRIX_STR="betaScan ${VIN} ${ALPHA} 0.5 ${JC_A} ${JC_B} ${JC_C} ${FINAL_TIME}" \
    BETA_LIST_STR="${betas}" \
    JC_M_LIST_STR="${JC_M}" \
    TQ_LIST_STR="${TQ}" \
    OUTPUT_FREQUENCY_FIXED="${OUT_FREQ}" \
    TIMESTEP_FIXED="${DT}" \
    bash "${DRIVER}"
  )
}

latest_summary_for_tag() {
  local tag="$1"
  local run_dir
  run_dir="$(ls -dt "${BASE_DIR}/runs_${tag}_"* 2>/dev/null | head -n1 || true)"
  if [[ -z "${run_dir}" ]]; then
    return 1
  fi
  ls "${run_dir}"/summary_cor_hmax_recomputed_*.txt 2>/dev/null | head -n1
}

analyze_summary() {
  python3 - "$@" <<'PY'
import math, sys

lower=float(sys.argv[1])
upper=float(sys.argv[2])
paths=sys.argv[3:]

states={lower:"nan", upper:"zero"}
positive=[]

for path in paths:
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            parts=line.split()
            name=parts[0]
            beta=None
            fields=name.split("_")
            for i,tok in enumerate(fields):
                if tok=="beta" and i+1 < len(fields):
                    beta=float(fields[i+1])
                    break
            if beta is None:
                continue
            cor=parts[3]
            if cor == "NaN":
                states[beta]="nan"
            else:
                val=float(cor)
                if math.isclose(val, 0.0, abs_tol=1e-15):
                    states[beta]="zero"
                elif val > 0.0:
                    states[beta]="positive"
                    positive.append((beta,val))

if positive:
    positive.sort()
    beta,val=positive[0]
    print(f"FOUND {beta:.12g} {val:.12g}")
    sys.exit(0)

nans=sorted(b for b,s in states.items() if s=="nan")
zeros=sorted(b for b,s in states.items() if s=="zero")

new_lower=None
new_upper=None
for b in reversed(nans):
    larger=[z for z in zeros if z>b]
    if larger:
        new_lower=b
        new_upper=min(larger)
        break

if new_lower is None or new_upper is None:
    print("NO_BRACKET")
else:
    print(f"BRACKET {new_lower:.12g} {new_upper:.12g}")
PY
}

make_refined_betas() {
  python3 - "$1" "$2" <<'PY'
import sys
lo=float(sys.argv[1]); hi=float(sys.argv[2])
vals=[lo + (hi-lo)*i/5.0 for i in range(1,5)]
print(" ".join(f"{v:.12g}" for v in vals))
PY
}

echo "Adaptive beta-until-positive scan starts"
echo "Initial bracket: (${LOWER_NAN_BETA}, ${UPPER_ZERO_BETA}]"
echo "Initial betas: ${INITIAL_BETAS}"

all_summaries=()
current_betas="${INITIAL_BETAS}"
iter=1
while (( iter <= MAX_ITER )); do
  tag=$(printf 'eq16_scan_beta_until_positive_iter%02d_v200_alpha1p5e-6' "${iter}")
  run_scan "${tag}" "${current_betas}"
  summary="$(latest_summary_for_tag "${tag}")"
  if [[ -z "${summary:-}" ]]; then
    echo "ERROR: summary not found for ${tag}"
    exit 1
  fi
  all_summaries+=("${summary}")
  result="$(analyze_summary "${LOWER_NAN_BETA}" "${UPPER_ZERO_BETA}" "${all_summaries[@]}")"
  echo "analysis[iter=${iter}]: ${result}"

  kind="$(awk '{print $1}' <<<"${result}")"
  if [[ "${kind}" == "FOUND" ]]; then
    echo "Found positive CoR: ${result}"
    exit 0
  fi
  if [[ "${kind}" != "BRACKET" ]]; then
    echo "No usable bracket left"
    exit 0
  fi

  new_lower="$(awk '{print $2}' <<<"${result}")"
  new_upper="$(awk '{print $3}' <<<"${result}")"

  if python3 - "${new_lower}" "${new_upper}" "${MIN_INTERVAL}" <<'PY'
import sys
lo=float(sys.argv[1]); hi=float(sys.argv[2]); tol=float(sys.argv[3])
sys.exit(0 if (hi-lo) <= tol else 1)
PY
  then
    echo "Bracket width <= MIN_INTERVAL; stop at (${new_lower}, ${new_upper}]"
    exit 0
  fi

  LOWER_NAN_BETA="${new_lower}"
  UPPER_ZERO_BETA="${new_upper}"
  current_betas="$(make_refined_betas "${LOWER_NAN_BETA}" "${UPPER_ZERO_BETA}")"
  echo "next betas: ${current_betas}"
  iter=$((iter+1))
done

echo "Reached MAX_ITER=${MAX_ITER} without finding positive CoR"
