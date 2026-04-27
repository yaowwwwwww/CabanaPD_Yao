#!/usr/bin/env bash
# ============================================================
# Generic and resumable scan driver for CabanaPD.
# - Exact resume/skip based on summary case_name.
# - Clear staging outputs only when a case fails.
# - Does not overwrite the template input JSON.
# ============================================================

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"
INPUT_JSON_TEMPLATE="${ROOT_DIR}/examples/mechanics/inputs/simple_impact_thermal.json"
CABANAPD_EXE="${ROOT_DIR}/build/examples/mechanics/ColdSprayImpactThermal"
PY_SCRIPT="${ROOT_DIR}/scripts/silo2csv.py"
OCT_SCRIPT="${ROOT_DIR}/scripts/avg-velocity.m"
BASE_DIR="${ROOT_DIR}/build"

# ========= Run naming =========
# Output directory:
#   ${BASE_DIR}/${DATE_TAG}_runs_${RUN_TAG_BASE}
# You can override full output directory with RUNS_ROOT_OVERRIDE.
# Keep tag short by default (date + scan name).
RUN_TAG_BASE="${RUN_TAG_BASE:-scan_jcc}"
# Optional: provide a summary file path to replay exact case_name entries.
CASE_LIST_FILE="${CASE_LIST_FILE:-}"

# ========= Fixed params =========
CZM_SCALE_FIXED=0
CZM_DECAY_FIXED=1.0
CZM_YIELD_FIXED=0.05
JC_N_FIXED=0.31

# ========= Scan params (edit these lists) =========
BALL_VIN_LIST=(400 500 600)
LJ_ALPHA_LIST=(1e-6)
LJ_BETA_LIST=(0.5)
JC_A_LIST=(2e9)
JC_B_LIST=(2e9)
JC_C_LIST=(0.025 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 25000)

require_cmd() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: command not found: ${cmd}" >&2
    exit 1
  fi
}

append_nan() {
  local case_name="$1"
  printf "%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" "${case_name}" >> "${SUMMARY_FILE}"
}

clear_staging_outputs() {
  rm -f "${BASE_DIR}"/*.silo "${BASE_DIR}"/*.csv 2>/dev/null || true
}

run_recompute() {
  local case_dir="$1"
  {
    echo "===== $(date +"%F %T") recompute: ${case_dir} ====="
    SUMMARY_FILE="${SUMMARY_FILE}" SINGLE_CASE_DIR="${case_dir}" \
      "${OCTAVE_CMD}" --no-gui --quiet --eval "source('${OCT_SCRIPT}'); fflush(stdout);"
  } >> "${RUNS_ROOT}/octave_analysis_output.txt" 2>&1
}

case_exists() {
  local case_name="$1"
  awk -v c="${case_name}" '$1==c{found=1; exit} END{exit (found?0:1)}' "${SUMMARY_FILE}" 2>/dev/null
}

run_one_case() {
  local BALL_VIN="$1"
  local LJ_ALPHA="$2"
  local LJ_BETA="$3"
  local JC_A="$4"
  local JC_B="$5"
  local JC_C="$6"
  local CASE_NAME="vin_${BALL_VIN}_alpha_${LJ_ALPHA}_beta_${LJ_BETA}_A_${JC_A}_B_${JC_B}_C_${JC_C}"
  local CASE_TAG="run_${CASE_NAME}"
  local CASE_DIR="${RUNS_ROOT}/${CASE_TAG}"

  if case_exists "${CASE_NAME}"; then
    echo "SKIP (already in summary): ${CASE_NAME}"
    return 0
  fi

  echo "-----------------------------------------------------"
  echo "RUN: ${CASE_NAME}"

  if ! jq --indent 2 \
    --argjson vin       "${BALL_VIN}" \
    --argjson LJalpha   "${LJ_ALPHA}" \
    --argjson LJbeta    "${LJ_BETA}" \
    --argjson jc_a      "${JC_A}" \
    --argjson jc_b      "${JC_B}" \
    --argjson jc_n      "${JC_N_FIXED}" \
    --argjson jc_c      "${JC_C}" \
    --argjson czm_scale "${CZM_SCALE_FIXED}" \
    --argjson czm_yield "${CZM_YIELD_FIXED}" \
    --argjson czm_decay "${CZM_DECAY_FIXED}" \
    '
    .ball_initial_velocity.value   = $vin       |
    .ball_center.value[2]          = 1.25e-5   |
    .contact_horizon_extend_factor.value = 0.05 |
    .contact_horizon_factor.value  = ((1.5e-6 - (0.05 * ((.high_corner.value[0] - .low_corner.value[0]) / .num_cells.value[0]))) / (2 * .LJr0.value * ((.high_corner.value[0] - .low_corner.value[0]) / .num_cells.value[0]))) |
    .LJalpha.value                 = $LJalpha   |
    .LJbeta.value                  = $LJbeta    |
    .yield_stress.value            = [ $jc_a, $jc_a ] |
    .jc_B.value                    = [ $jc_b, $jc_b ] |
    .jc_n.value                    = [ $jc_n, $jc_n ] |
    .jc_C.value                    = [ $jc_c, $jc_c ] |
    .CZM_cohesive_scaling.value    = $czm_scale |
    .CZM_yield_stretch.value       = $czm_yield |
    .CZM_degradation_rate.value    = $czm_decay
    ' "${INPUT_JSON_TEMPLATE}" > "${WORK_INPUT_JSON}"; then
    echo "WARN: jq failed for ${CASE_NAME}"
    append_nan "${CASE_NAME}"
    return 0
  fi

  cd "${BASE_DIR}" || exit 1

  "${CABANAPD_EXE}" "${WORK_INPUT_JSON}" 2>&1 | tee "${RUNS_ROOT}/cabana_last_output.txt"
  CABANA_STATUS=${PIPESTATUS[0]}
  if [ "${CABANA_STATUS}" -ne 0 ]; then
    echo "WARN: CabanaPD failed (${CABANA_STATUS}) for ${CASE_NAME}"
    clear_staging_outputs
    append_nan "${CASE_NAME}"
    return 0
  fi

  if ! pvpython "${PY_SCRIPT}" >/dev/null 2>&1; then
    echo "WARN: CSV conversion failed for ${CASE_NAME}"
    clear_staging_outputs
    append_nan "${CASE_NAME}"
    return 0
  fi

  mkdir -p "${CASE_DIR}"
  cp "${WORK_INPUT_JSON}" "${CASE_DIR}/input.json"

  shopt -s nullglob
  mv -f "${BASE_DIR}"/*.silo "${CASE_DIR}/" 2>/dev/null || true
  mv -f "${BASE_DIR}"/*.csv  "${CASE_DIR}/" 2>/dev/null || true
  shopt -u nullglob

  if ! run_recompute "${CASE_DIR}"; then
    echo "WARN: Octave recompute failed for ${CASE_NAME}"
    append_nan "${CASE_NAME}"
    return 0
  fi

  return 0
}

require_cmd jq
require_cmd pvpython
require_cmd awk
require_cmd tee

if command -v octave-cli >/dev/null 2>&1; then
  OCTAVE_CMD="octave-cli"
elif command -v octave >/dev/null 2>&1; then
  OCTAVE_CMD="octave"
else
  echo "ERROR: octave or octave-cli is required." >&2
  exit 1
fi

if [ ! -f "${INPUT_JSON_TEMPLATE}" ]; then
  echo "ERROR: input JSON not found: ${INPUT_JSON_TEMPLATE}" >&2
  exit 1
fi
if [ ! -x "${CABANAPD_EXE}" ]; then
  echo "ERROR: executable not found or not executable: ${CABANAPD_EXE}" >&2
  exit 1
fi
if [ ! -f "${PY_SCRIPT}" ]; then
  echo "ERROR: script not found: ${PY_SCRIPT}" >&2
  exit 1
fi
if [ ! -f "${OCT_SCRIPT}" ]; then
  echo "ERROR: script not found: ${OCT_SCRIPT}" >&2
  exit 1
fi

DATE_TAG=$(date +"%Y-%m-%d_%H-%M")
DATE_TAG_SHORT=$(date +"%Y%b%d-%H-%M")

if [ -n "${RUNS_ROOT_OVERRIDE:-}" ]; then
  RUNS_ROOT="${RUNS_ROOT_OVERRIDE}"
  mkdir -p "${RUNS_ROOT}"
  SUMMARY_FILE=$(ls -1t "${RUNS_ROOT}"/summary_cor_hmax_recomputed_*.txt 2>/dev/null | head -n 1)
  if [ -z "${SUMMARY_FILE}" ]; then
    SUMMARY_FILE="${RUNS_ROOT}/summary_cor_hmax_recomputed_${DATE_TAG_SHORT}.txt"
  fi
else
  RUNS_ROOT="${BASE_DIR}/${DATE_TAG}_runs_${RUN_TAG_BASE}"
  mkdir -p "${RUNS_ROOT}"
  SUMMARY_FILE="${RUNS_ROOT}/summary_cor_hmax_recomputed_${DATE_TAG_SHORT}.txt"
fi

if [ ! -f "${SUMMARY_FILE}" ]; then
  {
    printf "# params: BALL_VIN_LIST=%s; LJ_ALPHA_LIST=%s; LJ_BETA_LIST=%s; JC_A_LIST=%s; JC_B_LIST=%s; JC_C_LIST=%s; JC_N_FIXED=%s; CZM=(%s,%s,%s)\n" \
      "${BALL_VIN_LIST[*]}" \
      "${LJ_ALPHA_LIST[*]}" \
      "${LJ_BETA_LIST[*]}" \
      "${JC_A_LIST[*]}" \
      "${JC_B_LIST[*]}" \
      "${JC_C_LIST[*]}" \
      "${JC_N_FIXED}" \
      "${CZM_SCALE_FIXED}" "${CZM_YIELD_FIXED}" "${CZM_DECAY_FIXED}"
    printf "# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_max(m)    best_frame    h_residual(m)    A_residual(m2)    h_mean_residual(m)\n"
  } > "${SUMMARY_FILE}"
fi

WORK_INPUT_JSON="${RUNS_ROOT}/_current_input.json"
cleanup() {
  rm -f "${WORK_INPUT_JSON}" "${BASE_DIR}/cabana_last_output.txt" 2>/dev/null || true
}
trap cleanup EXIT

if [ -n "${CASE_LIST_FILE}" ]; then
  if [ ! -f "${CASE_LIST_FILE}" ]; then
    echo "ERROR: CASE_LIST_FILE not found: ${CASE_LIST_FILE}" >&2
    exit 1
  fi
  TOTAL_CASES=$(awk '!/^#/ && NF{print $1}' "${CASE_LIST_FILE}" | awk '!seen[$1]++' | wc -l)
else
  TOTAL_CASES=$(( ${#BALL_VIN_LIST[@]} * ${#LJ_ALPHA_LIST[@]} * ${#LJ_BETA_LIST[@]} * ${#JC_A_LIST[@]} * ${#JC_B_LIST[@]} * ${#JC_C_LIST[@]} ))
fi

echo "==== Generic scan start ($(date)) ===="
echo "Run root: ${RUNS_ROOT}"
echo "Summary: ${SUMMARY_FILE}"
echo "Total planned cases: ${TOTAL_CASES}"
echo "Fixed: jc_n=${JC_N_FIXED}, czm_scale=${CZM_SCALE_FIXED}, czm_yield=${CZM_YIELD_FIXED}, czm_decay=${CZM_DECAY_FIXED}"
if [ -n "${CASE_LIST_FILE}" ]; then
  echo "Replay mode from CASE_LIST_FILE: ${CASE_LIST_FILE}"
else
  echo "Scan vin    = ${BALL_VIN_LIST[*]}"
  echo "Scan alpha  = ${LJ_ALPHA_LIST[*]}"
  echo "Scan beta   = ${LJ_BETA_LIST[*]}"
  echo "Scan JC_A   = ${JC_A_LIST[*]}"
  echo "Scan JC_B   = ${JC_B_LIST[*]}"
  echo "Scan JC_C   = ${JC_C_LIST[*]}"
fi
echo

if [ -n "${CASE_LIST_FILE}" ]; then
  while IFS= read -r CASE_NAME_FROM_FILE; do
    if [[ "${CASE_NAME_FROM_FILE}" =~ ^vin_([^_]+)_alpha_([^_]+)_beta_([^_]+)_A_([^_]+)_B_([^_]+)_C_([^_]+)$ ]]; then
      run_one_case \
        "${BASH_REMATCH[1]}" \
        "${BASH_REMATCH[2]}" \
        "${BASH_REMATCH[3]}" \
        "${BASH_REMATCH[4]}" \
        "${BASH_REMATCH[5]}" \
        "${BASH_REMATCH[6]}"
    else
      echo "WARN: cannot parse case_name, skip: ${CASE_NAME_FROM_FILE}"
    fi
  done < <(awk '!/^#/ && NF{print $1}' "${CASE_LIST_FILE}" | awk '!seen[$1]++')
else
  for BALL_VIN in "${BALL_VIN_LIST[@]}"; do
    for LJ_ALPHA in "${LJ_ALPHA_LIST[@]}"; do
      for LJ_BETA in "${LJ_BETA_LIST[@]}"; do
        for JC_A in "${JC_A_LIST[@]}"; do
          for JC_B in "${JC_B_LIST[@]}"; do
            for JC_C in "${JC_C_LIST[@]}"; do
              run_one_case "${BALL_VIN}" "${LJ_ALPHA}" "${LJ_BETA}" "${JC_A}" "${JC_B}" "${JC_C}"
            done
          done
        done
      done
    done
  done
fi

echo
echo "Scan complete."
echo "Summary: ${SUMMARY_FILE}"
echo "Cases archived under: ${RUNS_ROOT}"
