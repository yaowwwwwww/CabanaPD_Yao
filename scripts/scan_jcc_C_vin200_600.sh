#!/bin/bash
# ============================================================
# Scan: JC_A × JC_B for JC hardening influence
# Fixed: LJ params, JC n/C, CZM params
# Scan vin in [200, 600] m/s
# ============================================================

set -u

INPUT_JSON="/home/wuwen/program/CabanaPD_Yao/examples/mechanics/inputs/simple_impact_thermal.json"
CABANAPD_EXE="/home/wuwen/program/CabanaPD_Yao/build/examples/mechanics/ColdSprayImpactThermal"
PY_SCRIPT="/home/wuwen/program/CabanaPD_Yao/scripts/silo2csv.py"
OCT_SCRIPT="/home/wuwen/program/CabanaPD_Yao/scripts/avg-velocity.m"
BASE_DIR="/home/wuwen/program/CabanaPD_Yao/build"

DATE_TAG=$(date +"%Y-%m-%d_%H-%M")
DATE_TAG_SHORT=$(date +"%Y%b%d-%H-%M")
RUNS_ROOT="${BASE_DIR}/runs_jcc_AB_scan_v200_600_${DATE_TAG}"
mkdir -p "${RUNS_ROOT}"
SUMMARY_FILE="${RUNS_ROOT}/summary_cor_hmax_recomputed_${DATE_TAG_SHORT}.txt"

# ========= Fixed params =========
CZM_SCALE_FIXED=0
CZM_DECAY_FIXED=1.0
CZM_YIELD_FIXED=0.05
LJ_ALPHA_FIXED=1e-12
LJ_BETA_FIXED=0.5
JC_N_FIXED=0.31
JC_C_FIXED=0.0

# ========= Scan params =========
BALL_VIN_LIST=(600 500 400 300 200)
JC_A_LIST=(8e7 1e8 5e8 1e9 2e9)
JC_B_LIST=(0 1e6 8e6 3.65e9)

{
  echo "==== Scan: JC_A × JC_B with vin={200..600} (timestamp = ${DATE_TAG}) ===="
  date
  echo "Fixed LJ: alpha=${LJ_ALPHA_FIXED}, beta=${LJ_BETA_FIXED}"
  echo "Fixed JC: n=${JC_N_FIXED}, C=${JC_C_FIXED}"
  echo "Fixed CZM: scale=${CZM_SCALE_FIXED}, decay=${CZM_DECAY_FIXED}, yield=${CZM_YIELD_FIXED}"
  echo "Scanning vin = ${BALL_VIN_LIST[*]}"
  echo "Scanning JC_A = ${JC_A_LIST[*]}"
  echo "Scanning JC_B = ${JC_B_LIST[*]}"
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

for BALL_VIN in "${BALL_VIN_LIST[@]}"; do
  for JC_A in "${JC_A_LIST[@]}"; do
    for JC_B in "${JC_B_LIST[@]}"; do
      echo "-----------------------------------------------------"
      echo "🚀 Running case: V_in=${BALL_VIN}, A=${JC_A}, B=${JC_B}, C=${JC_C_FIXED}"

      jq --indent 2 \
        --argjson vin       "$BALL_VIN" \
        --argjson LJalpha   "$LJ_ALPHA_FIXED" \
        --argjson LJbeta    "$LJ_BETA_FIXED" \
        --argjson jc_a      "$JC_A" \
        --argjson jc_b      "$JC_B" \
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
      ${CABANAPD_EXE} "${INPUT_JSON}" 2>&1 | tee cabana_last_output.txt
      CABANA_STATUS=${PIPESTATUS[0]}
      if [ $CABANA_STATUS -ne 0 ]; then
        echo "❌ CabanaPD failed (code $CABANA_STATUS)"
        printf "vin_%s_alpha_%s_beta_%s_A_%s_B_%s_C_%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" \
          "$BALL_VIN" "$LJ_ALPHA_FIXED" "$LJ_BETA_FIXED" "$JC_A" "$JC_B" "$JC_C_FIXED" >> "${SUMMARY_FILE}"
        continue
      fi

      pvpython "${PY_SCRIPT}" >/dev/null 2>&1
      CSV_STATUS=$?
      if [ $CSV_STATUS -ne 0 ]; then
        echo "⚠️ CSV conversion failed"
        printf "vin_%s_alpha_%s_beta_%s_A_%s_B_%s_C_%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" \
          "$BALL_VIN" "$LJ_ALPHA_FIXED" "$LJ_BETA_FIXED" "$JC_A" "$JC_B" "$JC_C_FIXED" >> "${SUMMARY_FILE}"
        continue
      fi

      CASE_TAG="run_vin_${BALL_VIN}_alpha_${LJ_ALPHA_FIXED}_beta_${LJ_BETA_FIXED}_A_${JC_A}_B_${JC_B}_C_${JC_C_FIXED}"
      CASE_DIR="${RUNS_ROOT}/${CASE_TAG}"
      mkdir -p "${CASE_DIR}"
      cp "${INPUT_JSON}" "${CASE_DIR}/input.json"

      shopt -s nullglob
      mv -f "${BASE_DIR}"/*.silo "${CASE_DIR}/" 2>/dev/null || true
      mv -f "${BASE_DIR}"/*.csv  "${CASE_DIR}/" 2>/dev/null || true
      shopt -u nullglob
      rm -f "${BASE_DIR}/cabana_last_output.txt"

      if ! run_recompute "${CASE_DIR}"; then
        echo "⚠️ Octave recompute failed for vin_${BALL_VIN}_alpha_${LJ_ALPHA_FIXED}_beta_${LJ_BETA_FIXED}_A_${JC_A}_B_${JC_B}_C_${JC_C_FIXED}"
        printf "vin_%s_alpha_%s_beta_%s_A_%s_B_%s_C_%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n" \
          "$BALL_VIN" "$LJ_ALPHA_FIXED" "$LJ_BETA_FIXED" "$JC_A" "$JC_B" "$JC_C_FIXED" >> "${SUMMARY_FILE}"
        continue
      fi
    done
  done
done

echo ""
echo "🎯 JC_A × JC_B scan complete."
echo "Summary: ${SUMMARY_FILE}"
echo "Cases archived under: ${RUNS_ROOT}"
