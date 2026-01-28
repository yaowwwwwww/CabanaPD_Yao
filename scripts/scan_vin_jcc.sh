#!/bin/bash 
# ============================================================
# Scan: VIN × JC_C
# Fixed: LJ params, CZM params, yield_stress
# Record: vout, CoR, Lateralmax, h_residual, A_residual, h_mean_residual
# ============================================================

set -u

INPUT_JSON="/home/wuwen/program/CabanaPD_Yao/examples/mechanics/inputs/simple_impact_thermal.json"
CABANAPD_EXE="/home/wuwen/program/CabanaPD_Yao/build/examples/mechanics/ColdSprayImpactThermal"
PY_SCRIPT="/home/wuwen/program/CabanaPD_Yao/build/silo2csv.py"
OCT_SCRIPT="/home/wuwen/program/CabanaPD_Yao/build/avg-velocity.m"
BASE_DIR="/home/wuwen/program/CabanaPD_Yao/build"

DATE_TAG=$(date +"%Y-%m-%d_%H-%M")
RUNS_ROOT="${BASE_DIR}/runs_jcc_scan_${DATE_TAG}"
mkdir -p "${RUNS_ROOT}"
LOG_FILE="${BASE_DIR}/summary_VIN_JCC_scan_${DATE_TAG}.txt"

# ========= Fixed params =========
CZM_SCALE_FIXED=0
CZM_DECAY_FIXED=1.0
CZM_YIELD_FIXED=0.05

YIELD_STRESS_FIXED=80E+6

LJbeta_FIXED=0.5
LJalpha_FIXED=2.5e-12   # <-- FIXED (was 2.5e9e-12 typo)

# ========= Scan params =========
BALL_VIN_LIST=(600 500 400 300 200)
JC_C_LIST=(0.0 0.01 0.03 0.05)

# ========= Log header =========
{
  echo "==== Scan: VIN × JC_C (timestamp = ${DATE_TAG}) ===="
  date
  echo "Fixed CZM: scale=${CZM_SCALE_FIXED}, decay=${CZM_DECAY_FIXED}, yield=${CZM_YIELD_FIXED}"
  echo "Fixed yield_stress (Cu[1]) = ${YIELD_STRESS_FIXED} Pa"
  echo "Fixed LJ: beta=${LJbeta_FIXED}, alpha=${LJalpha_FIXED}"
  echo "Scanning JC_C = ${JC_C_LIST[*]}"
  echo "Scanning vin  = ${BALL_VIN_LIST[*]}"
  echo ""
} > "$LOG_FILE"

# ========= Main loop =========
for JC_C in "${JC_C_LIST[@]}"; do
  for BALL_VIN in "${BALL_VIN_LIST[@]}"; do

    echo "-----------------------------------------------------"
    echo "🚀 Running case: V_in=${BALL_VIN} m/s, JC_C=${JC_C} (beta=${LJbeta_FIXED})"

    # ====== Update JSON ======
    jq --indent 2 \
      --argjson LJalpha   "$LJalpha_FIXED" \
      --argjson LJbeta    "$LJbeta_FIXED" \
      --argjson vin       "$BALL_VIN" \
      --argjson jc_c      "$JC_C" \
      --argjson czm_scale "$CZM_SCALE_FIXED" \
      --argjson czm_yield "$CZM_YIELD_FIXED" \
      --argjson czm_decay "$CZM_DECAY_FIXED" \
      --argjson ys        "$YIELD_STRESS_FIXED" \
      '
      .LJalpha.value                 = $LJalpha   |
      .LJbeta.value                  = $LJbeta    |
      .ball_initial_velocity.value   = $vin       |
      .jc_C.value                    = [ $jc_c, $jc_c ] |
      .CZM_cohesive_scaling.value    = $czm_scale |
      .CZM_yield_stretch.value       = $czm_yield |
      .CZM_degradation_rate.value    = $czm_decay |
      .yield_stress.value[1]         = $ys
      ' "${INPUT_JSON}" > tmp.json && mv tmp.json "${INPUT_JSON}"

    # ====== Run CabanaPD ======
    cd "${BASE_DIR}"
    echo ">>> Running CabanaPD simulation..."
    ${CABANAPD_EXE} "${INPUT_JSON}" 2>&1 | tee cabana_last_output.txt
    CABANA_STATUS=${PIPESTATUS[0]}
    if [ $CABANA_STATUS -ne 0 ]; then
      echo "❌ CabanaPD failed (code $CABANA_STATUS)"
      printf "vin=%-8s  JC_C=%-6s  ❌ CabanaPD failed (code %d)\n" \
        "$BALL_VIN" "$JC_C" "$CABANA_STATUS" >> "$LOG_FILE"
      continue
    fi

    # ====== Silo → CSV ======
    pvpython "${PY_SCRIPT}" >/dev/null 2>&1
    CSV_STATUS=$?
    if [ $CSV_STATUS -ne 0 ]; then
      echo "⚠️ CSV conversion failed"
      printf "vin=%-8s  JC_C=%-6s  ⚠️ CSV conversion failed\n" \
        "$BALL_VIN" "$JC_C" >> "$LOG_FILE"
      continue
    fi

    # ====== Octave analysis ======
    octave --no-gui --quiet --eval "source('${OCT_SCRIPT}'); fflush(stdout);" \
      > octave_analysis_output.txt 2>&1
    OCT_STATUS=$?
    if [ $OCT_STATUS -ne 0 ]; then
      echo "⚠️ Octave analysis failed"
      printf "vin=%-8s  JC_C=%-6s  ⚠️ Octave analysis failed\n" \
        "$BALL_VIN" "$JC_C" >> "$LOG_FILE"
      continue
    fi

    # ====== Extract results ======
    VIN=$(grep -a "vin ="  octave_analysis_output.txt | sed -E 's/.*vin *= *([0-9.+-eE]+).*/\1/')
    VOUT=$(grep -a "vout =" octave_analysis_output.txt | sed -E 's/.*vout *= *([0-9.+-eE]+).*/\1/')
    COR=$(grep -a "CoR"    octave_analysis_output.txt | sed -E 's/.*= *([0-9.+-eE]+).*/\1/')
    DEFORM=$(grep -a "Lateralmax" octave_analysis_output.txt | sed -E 's/.*Lateralmax *= *([0-9.+-eE]+).*/\1/')
    HRES=$(grep -a "h_residual" octave_analysis_output.txt | sed -E 's/.*h_residual *= *([0-9.+-eE]+).*/\1/')
    ARES=$(grep -a "A_residual" octave_analysis_output.txt | sed -E 's/.*A_residual *= *([0-9.+-eE]+).*/\1/')
    HMEAN_RES=$(grep -a "h_mean_residual" octave_analysis_output.txt | sed -E 's/.*h_mean_residual *= *([0-9.+-eE]+).*/\1/')

    # ====== Write summary ======
    if [ -n "${COR:-}" ] && [ "${COR:-NaN}" != "NaN" ]; then
      printf "vin=%-8s  JC_C=%-6s  vout=%8s  CoR=%8s  Lmax=%8s  h_res=%8s  A_res=%8s  h_mean_res=%8s\n" \
        "$BALL_VIN" "$JC_C" \
        "${VOUT:-N/A}" "${COR:-N/A}" "${DEFORM:-N/A}" \
        "${HRES:-N/A}" "${ARES:-N/A}" "${HMEAN_RES:-N/A}" >> "$LOG_FILE"
    else
      printf "vin=%-8s  JC_C=%-6s  ❌ No valid CoR result (h_res=%s A_res=%s h_mean_res=%s)\n" \
        "$BALL_VIN" "$JC_C" \
        "${HRES:-N/A}" "${ARES:-N/A}" "${HMEAN_RES:-N/A}" >> "$LOG_FILE"
    fi

    # ====== Archive outputs ======
    CASE_TAG="run_vin_${BALL_VIN}_jcc_${JC_C}"
    CASE_DIR="${RUNS_ROOT}/${CASE_TAG}"
    mkdir -p "${CASE_DIR}"

    cp "${INPUT_JSON}" "${CASE_DIR}/input.json"

    shopt -s nullglob
    mv -f "${BASE_DIR}"/*.silo "${CASE_DIR}/" 2>/dev/null || true
    mv -f "${BASE_DIR}"/*.csv  "${CASE_DIR}/" 2>/dev/null || true
    shopt -u nullglob

    rm -f "${BASE_DIR}/cabana_last_output.txt" "${BASE_DIR}/octave_analysis_output.txt"

  done
done

echo ""
echo "🎯 JC_C scan complete. Summary: ${LOG_FILE}"
echo "Cases archived under: ${RUNS_ROOT}"
