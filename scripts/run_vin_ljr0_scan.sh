#!/bin/bash 
# ============================================================
# Scan: VIN × LJr0
# Fixed: LJbeta, LJalpha, CZM params, yield_stress
# Record: vout, CoR, Lateralmax, h_max
# ============================================================

set -u

INPUT_JSON="/home/wuwen/program/CabanaPD_Yao/examples/mechanics/inputs/simple_impact.json"
CABANAPD_EXE="/home/wuwen/program/CabanaPD_Yao/build/examples/mechanics/ColdSprayImpact"
PY_SCRIPT="/home/wuwen/program/CabanaPD_Yao/scripts/silo2csv.py"
OCT_SCRIPT="/home/wuwen/program/CabanaPD_Yao/scripts/avg-velocity.m"
BASE_DIR="/home/wuwen/program/CabanaPD_Yao/build"

DATE_TAG=$(date +"%Y-%m-%d_%H-%M")
RUNS_ROOT="${BASE_DIR}/runs_ljr0_scan_${DATE_TAG}"
mkdir -p "${RUNS_ROOT}"
LOG_FILE="${RUNS_ROOT}/scan_status_${DATE_TAG}.txt"
SUMMARY_FILE="${RUNS_ROOT}/summary_cor_hmax_recomputed.txt"

# ========= Fixed params =========
CZM_SCALE_FIXED=0
CZM_DECAY_FIXED=1.0
CZM_YIELD_FIXED=0.05
YIELD_STRESS_FIXED=1e9

LJbeta_FIXED=0.5
LJalpha_FIXED=2.5e-12   # <-- FIXED (was 2.5e9e-12 typo)

# ========= Scan params =========
BALL_VIN_LIST=(600 500 400 300 200 800 700 650)
LJ_R0_LIST=(0.5 0.6 0.7 0.8 0.9 1.01)

# ========= Status log header =========
{
  echo "==== Scan status: VIN × LJr0 (timestamp = ${DATE_TAG}) ===="
  date
  echo "Fixed CZM: scale=${CZM_SCALE_FIXED}, decay=${CZM_DECAY_FIXED}, yield=${CZM_YIELD_FIXED}"
  echo "Fixed yield_stress (Cu[1]) = ${YIELD_STRESS_FIXED} Pa"
  echo "Fixed LJ: beta=${LJbeta_FIXED}, alpha=${LJalpha_FIXED}"
  echo "Scanning LJr0 = ${LJ_R0_LIST[*]}"
  echo "Scanning vin  = ${BALL_VIN_LIST[*]}"
  echo ""
} > "$LOG_FILE"

# ========= Main loop =========
for LJ_R0 in "${LJ_R0_LIST[@]}"; do
  for BALL_VIN in "${BALL_VIN_LIST[@]}"; do

    echo "-----------------------------------------------------"
    echo "🚀 Running case: V_in=${BALL_VIN} m/s, LJ_r0=${LJ_R0} (beta=${LJbeta_FIXED})"

    # ====== Update JSON ======
    jq --indent 2 \
      --argjson LJr0      "$LJ_R0" \
      --argjson LJalpha   "$LJalpha_FIXED" \
      --argjson LJbeta    "$LJbeta_FIXED" \
      --argjson vin       "$BALL_VIN" \
      --argjson czm_scale "$CZM_SCALE_FIXED" \
      --argjson czm_yield "$CZM_YIELD_FIXED" \
      --argjson czm_decay "$CZM_DECAY_FIXED" \
      --argjson ys        "$YIELD_STRESS_FIXED" \
      '
      .LJr0.value                    = $LJr0      |
      .LJalpha.value                 = $LJalpha   |
      .LJbeta.value                  = $LJbeta    |
      .ball_initial_velocity.value   = $vin       |
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
      printf "vin=%-8s  LJr0=%-6s  ❌ CabanaPD failed (code %d)\n" \
        "$BALL_VIN" "$LJ_R0" "$CABANA_STATUS" >> "$LOG_FILE"
      continue
    fi

    # ====== Silo → CSV ======
    pvpython "${PY_SCRIPT}" >/dev/null 2>&1
    CSV_STATUS=$?
    if [ $CSV_STATUS -ne 0 ]; then
      echo "⚠️ CSV conversion failed"
      printf "vin=%-8s  LJr0=%-6s  ⚠️ CSV conversion failed\n" \
        "$BALL_VIN" "$LJ_R0" >> "$LOG_FILE"
      continue
    fi

    printf "vin=%-8s  LJr0=%-6s  ✅ simulation + CSV ok\n" \
      "$BALL_VIN" "$LJ_R0" >> "$LOG_FILE"

    # ====== Archive outputs ======
    CASE_TAG="run_vin_${BALL_VIN}_r0_${LJ_R0}"
    CASE_DIR="${RUNS_ROOT}/${CASE_TAG}"
    mkdir -p "${CASE_DIR}"

    cp "${INPUT_JSON}" "${CASE_DIR}/input.json"

    shopt -s nullglob
    mv -f "${BASE_DIR}"/*.silo "${CASE_DIR}/" 2>/dev/null || true
    mv -f "${BASE_DIR}"/*.csv  "${CASE_DIR}/" 2>/dev/null || true
    shopt -u nullglob

    rm -f "${BASE_DIR}/cabana_last_output.txt"

  done
done

echo "Running recompute analysis in ${RUNS_ROOT} ..."
(
  cd "${RUNS_ROOT}" || exit 1
  octave --no-gui --quiet --eval "source('${OCT_SCRIPT}'); fflush(stdout);" \
    > octave_analysis_output.txt 2>&1
)
OCT_STATUS=$?
if [ $OCT_STATUS -ne 0 ]; then
  echo "⚠️ Octave recompute failed. See ${RUNS_ROOT}/octave_analysis_output.txt"
else
  echo "✅ Octave recompute done: ${SUMMARY_FILE}"
fi

echo ""
echo "🎯 LJr0 scan complete."
echo "Status log: ${LOG_FILE}"
echo "Cases archived under: ${RUNS_ROOT}"
