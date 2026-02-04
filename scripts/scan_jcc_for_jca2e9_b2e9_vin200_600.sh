#!/bin/bash
# ============================================================
# Scan: LJ_alpha with fixed JC(A/B/n/C) at vin=300 m/s
# Fixed: LJ beta, JC params, CZM params
# ============================================================

set -u

INPUT_JSON="/home/wuwen/program/CabanaPD_Yao/examples/mechanics/inputs/simple_impact_thermal.json"
CABANAPD_EXE="/home/wuwen/program/CabanaPD_Yao/build/examples/mechanics/ColdSprayImpactThermal"
PY_SCRIPT="/home/wuwen/program/CabanaPD_Yao/scripts/silo2csv.py"
OCT_SCRIPT="/home/wuwen/program/CabanaPD_Yao/scripts/avg-velocity.m"
BASE_DIR="/home/wuwen/program/CabanaPD_Yao/build"

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
  RUNS_ROOT="${BASE_DIR}/runs_ljalpha_scan_A9e7_B2p92e8_C0p025_v300_${DATE_TAG}"
  mkdir -p "${RUNS_ROOT}"
  SUMMARY_FILE="${RUNS_ROOT}/summary_cor_hmax_recomputed_${DATE_TAG_SHORT}.txt"
fi

# ========= Fixed params =========
CZM_SCALE_FIXED=0
CZM_DECAY_FIXED=1.0
CZM_YIELD_FIXED=0.05
LJ_BETA_FIXED=0.5
JC_A_FIXED=9e7
JC_B_FIXED=2.92e8
JC_N_FIXED=0.31
JC_C_FIXED=0.025
# ========= Scan params =========
BALL_VIN_FIXED=300
LJ_ALPHA_LIST=(
  1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10
  1e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-18 1e-19 1e-20
)

{
  echo "==== Scan: LJ_alpha with fixed JC and vin=300 (timestamp = ${DATE_TAG}) ===="
  date
  echo "Fixed LJ: beta=${LJ_BETA_FIXED}"
  echo "Fixed vin: ${BALL_VIN_FIXED}"
  echo "Fixed JC: A=${JC_A_FIXED}, B=${JC_B_FIXED}, n=${JC_N_FIXED}, C=${JC_C_FIXED}"
  echo "Fixed CZM: scale=${CZM_SCALE_FIXED}, decay=${CZM_DECAY_FIXED}, yield=${CZM_YIELD_FIXED}"
  echo "Scanning LJ_alpha = ${LJ_ALPHA_LIST[*]}"
  echo ""
}

if [ ! -f "${SUMMARY_FILE}" ]; then
  printf "# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_max(m)    best_frame    h_residual(m)    A_residual(m2)    h_mean_residual(m)\n" \
    > "${SUMMARY_FILE}"
fi

run_recompute() {
  local case_dir="$1"
  local octave_cmd="octave"
  if command -v octave-cli >/dev/null 2>&1; then
    octave_cmd="octave-cli"
  fi

  {
    echo "===== $(date +"%F %T") recompute: ${case_dir} ====="
    SUMMARY_FILE="${SUMMARY_FILE}" SINGLE_CASE_DIR="${case_dir}" \
      "${octave_cmd}" --no-gui --quiet --eval "source('${OCT_SCRIPT}'); fflush(stdout);"
  } >> "${RUNS_ROOT}/octave_analysis_output.txt" 2>&1
}

for LJ_ALPHA in "${LJ_ALPHA_LIST[@]}"; do
    CASE_NAME="vin_${BALL_VIN_FIXED}_alpha_${LJ_ALPHA}_beta_${LJ_BETA_FIXED}_A_${JC_A_FIXED}_B_${JC_B_FIXED}_C_${JC_C_FIXED}"
    CASE_TAG="run_vin_${BALL_VIN_FIXED}_alpha_${LJ_ALPHA}_beta_${LJ_BETA_FIXED}_A_${JC_A_FIXED}_B_${JC_B_FIXED}_C_${JC_C_FIXED}"
    CASE_DIR="${RUNS_ROOT}/${CASE_TAG}"

    if awk -v c="${CASE_NAME}" '$1==c{found=1; exit} END{exit (found?0:1)}' "${SUMMARY_FILE}" 2>/dev/null; then
      echo "⏭️  Skip existing case in summary: ${CASE_NAME}"
      continue
    fi

    echo "-----------------------------------------------------"
    echo "🚀 Running case: V_in=${BALL_VIN_FIXED}, LJ_alpha=${LJ_ALPHA}, A=${JC_A_FIXED}, B=${JC_B_FIXED}, C=${JC_C_FIXED}"

    jq --indent 2 \
      --argjson vin       "$BALL_VIN_FIXED" \
      --argjson LJalpha   "$LJ_ALPHA" \
      --argjson LJbeta    "$LJ_BETA_FIXED" \
      --argjson jc_a      "$JC_A_FIXED" \
      --argjson jc_b      "$JC_B_FIXED" \
      --argjson jc_n      "$JC_N_FIXED" \
      --argjson jc_c      "$JC_C_FIXED" \
      --argjson czm_scale "$CZM_SCALE_FIXED" \
      --argjson czm_yield "$CZM_YIELD_FIXED" \
      --argjson czm_decay "$CZM_DECAY_FIXED" \
      '
      .ball_initial_velocity.value   = $vin       |
      .LJalpha.value                 = $LJalpha   |
      .LJbeta.value                  = $LJbeta    |
      .yield_stress.value            = [ $jc_a, $jc_a ] |
      .jc_B.value                    = [ $jc_b, $jc_b ] |
      .jc_n.value                    = [ $jc_n, $jc_n ] |
      .jc_C.value                    = [ $jc_c, $jc_c ] |
      .CZM_cohesive_scaling.value    = $czm_scale |
      .CZM_yield_stretch.value       = $czm_yield |
      .CZM_degradation_rate.value    = $czm_decay
      ' "${INPUT_JSON}" > tmp.json && mv tmp.json "${INPUT_JSON}"

    cd "${BASE_DIR}"
    rm -f "${BASE_DIR}"/*.silo "${BASE_DIR}"/*.csv 2>/dev/null || true
    ${CABANAPD_EXE} "${INPUT_JSON}" 2>&1 | tee cabana_last_output.txt
    CABANA_STATUS=${PIPESTATUS[0]}
    if [ $CABANA_STATUS -ne 0 ]; then
      echo "❌ CabanaPD failed (code $CABANA_STATUS)"
      printf "vin_%s_alpha_%s_beta_%s_A_%s_B_%s_C_%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" \
        "$BALL_VIN_FIXED" "$LJ_ALPHA" "$LJ_BETA_FIXED" "$JC_A_FIXED" "$JC_B_FIXED" "$JC_C_FIXED" >> "${SUMMARY_FILE}"
      continue
    fi

    pvpython "${PY_SCRIPT}" >/dev/null 2>&1
    CSV_STATUS=$?
    if [ $CSV_STATUS -ne 0 ]; then
      echo "⚠️ CSV conversion failed"
      printf "vin_%s_alpha_%s_beta_%s_A_%s_B_%s_C_%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" \
        "$BALL_VIN_FIXED" "$LJ_ALPHA" "$LJ_BETA_FIXED" "$JC_A_FIXED" "$JC_B_FIXED" "$JC_C_FIXED" >> "${SUMMARY_FILE}"
      continue
    fi

    mkdir -p "${CASE_DIR}"
    cp "${INPUT_JSON}" "${CASE_DIR}/input.json"

    shopt -s nullglob
    mv -f "${BASE_DIR}"/*.silo "${CASE_DIR}/" 2>/dev/null || true
    mv -f "${BASE_DIR}"/*.csv  "${CASE_DIR}/" 2>/dev/null || true
    shopt -u nullglob
    rm -f "${BASE_DIR}/cabana_last_output.txt"

    if ! run_recompute "${CASE_DIR}"; then
      echo "⚠️ Octave recompute failed for ${CASE_NAME}"
      printf "vin_%s_alpha_%s_beta_%s_A_%s_B_%s_C_%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" \
        "$BALL_VIN_FIXED" "$LJ_ALPHA" "$LJ_BETA_FIXED" "$JC_A_FIXED" "$JC_B_FIXED" "$JC_C_FIXED" >> "${SUMMARY_FILE}"
      continue
    fi
done

echo ""
echo "🎯 LJ_alpha scan complete."
echo "Summary: ${SUMMARY_FILE}"
echo "Cases archived under: ${RUNS_ROOT}"
