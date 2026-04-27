#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"
BASE_DIR="${ROOT_DIR}/build"
DRIVER="${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"

ALPHA="${ALPHA:-1.5e-6}"
VIN="${VIN:-200}"
JC_A="${JC_A:-9e7}"
JC_B="${JC_B:-2.92e8}"
JC_C="${JC_C:-0.025}"
JC_M="${JC_M:-1.09}"
FINAL_TIME="${FINAL_TIME:-8e-8}"
TQ="${TQ:-0}"
OUT_FREQ="${OUT_FREQ:-200}"
DT="${DT:-1e-11}"

ANCHOR_SUMMARY="${ANCHOR_SUMMARY:-${BASE_DIR}/2026-04-15_20-02_runs_eq16_scan_beta_0p1_0p5_0p9_v200_alpha1p5e-6/summary_cor_hmax_recomputed_2026Apr15-20-02.txt}"

PRIMARY_BETAS="${PRIMARY_BETAS:-0.2 0.3 0.4}"

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

primary_all_zero_or_nan() {
  python3 - "$1" <<'PY'
import math, sys
path = sys.argv[1]
states = []
with open(path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        cor = parts[3]
        if cor == "NaN":
            states.append("nan")
        else:
            val = float(cor)
            if math.isclose(val, 0.0, abs_tol=1e-15):
                states.append("zero")
            else:
                states.append("positive")
print("yes" if states and all(s in {"nan", "zero"} for s in states) else "no")
PY
}

build_refined_betas() {
  python3 - "$@" <<'PY'
import math, sys
anchors = sys.argv[1]
primary = sys.argv[2]

def parse(path):
    out = {}
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            parts=line.split()
            name=parts[0]
            beta=None
            for token in name.split("_"):
                pass
            fields=name.split("_")
            for i, tok in enumerate(fields):
                if tok == "beta" and i + 1 < len(fields):
                    beta=float(fields[i+1])
                    break
            if beta is None:
                continue
            cor=parts[3]
            if cor == "NaN":
                state="nan"
            else:
                val=float(cor)
                if math.isclose(val, 0.0, abs_tol=1e-15):
                    state="zero"
                elif val > 0.0:
                    state="positive"
                else:
                    state="other"
            out[beta]=state
    return out

data={}
for path in [anchors, primary]:
    data.update(parse(path))

nan_betas=sorted([b for b,s in data.items() if s=="nan"])
zero_betas=sorted([b for b,s in data.items() if s=="zero"])

if not nan_betas or not zero_betas:
    sys.exit(2)

left=max([b for b in nan_betas if any(z>b for z in zero_betas)], default=None)
right=min([b for b in zero_betas if any(n<b for n in nan_betas)], default=None)

if left is None or right is None or not (left < right):
    sys.exit(3)

vals=[left + (right-left)*0.25, left + (right-left)*0.5, left + (right-left)*0.75]
print(" ".join(f"{v:.6g}" for v in vals))
PY
}

echo "Adaptive beta scan starts"
echo "Anchor summary: ${ANCHOR_SUMMARY}"
echo "Primary betas: ${PRIMARY_BETAS}"

PRIMARY_TAG="eq16_scan_beta_0p2_0p3_0p4_v200_alpha1p5e-6"
run_scan "${PRIMARY_TAG}" "${PRIMARY_BETAS}"
PRIMARY_SUMMARY="$(latest_summary_for_tag "${PRIMARY_TAG}")"

if [[ -z "${PRIMARY_SUMMARY:-}" ]]; then
  echo "ERROR: primary summary not found"
  exit 1
fi

echo "Primary summary: ${PRIMARY_SUMMARY}"

if [[ "$(primary_all_zero_or_nan "${PRIMARY_SUMMARY}")" != "yes" ]]; then
  echo "Primary scan already contains positive CoR; no refined scan launched."
  exit 0
fi

REFINED_BETAS="$(build_refined_betas "${ANCHOR_SUMMARY}" "${PRIMARY_SUMMARY}" || true)"

if [[ -z "${REFINED_BETAS}" ]]; then
  echo "No valid NaN/zero bracket found after primary scan; no refined scan launched."
  exit 0
fi

echo "Refined betas: ${REFINED_BETAS}"

REFINED_TAG="eq16_scan_beta_refined_v200_alpha1p5e-6"
run_scan "${REFINED_TAG}" "${REFINED_BETAS}"

REFINED_SUMMARY="$(latest_summary_for_tag "${REFINED_TAG}" || true)"
echo "Refined summary: ${REFINED_SUMMARY}"
echo "Adaptive beta scan done"
