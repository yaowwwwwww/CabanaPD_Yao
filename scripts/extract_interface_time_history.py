#!/usr/bin/env bash

set -e
set -u
set -o pipefail

ROOT="/home/wuwen/program/CabanaPD_Yao"
BUILD="${ROOT}/build"

INPUT_JSON="${ROOT}/examples/mechanics/inputs/simple_impact_thermal.json"
EXE="${BUILD}/examples/mechanics/ColdSprayImpactThermal"
SILO2CSV="${ROOT}/scripts/silo2csv.py"
AVG_M="${ROOT}/scripts/avg-velocity.m"

DATE_TAG="$(date +"%Y-%m-%d_%H-%M")"
RUN_DIR="${BUILD}/${DATE_TAG}_runs_json_target_v100_600"
SUMMARY_FILE="${RUN_DIR}/summary_velocity_scan.txt"

mkdir -p "${RUN_DIR}"

cat > "${SUMMARY_FILE}" <<EOF
# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_max(m)    best_frame    h_residual(m)    A_residual(m2)    h_mean_residual(m)
EOF

if command -v octave-cli >/dev/null 2>&1; then
    OCTAVE_CMD="octave-cli"
elif command -v octave >/dev/null 2>&1; then
    OCTAVE_CMD="octave"
else
    echo "ERROR: octave or octave-cli not found"
    exit 1
fi

for v in $(seq 100 50 600); do
    case_name="vin_${v}"
    case_dir="${RUN_DIR}/${case_name}"
    work_json="${case_dir}/input.json"

    echo
    echo "========== Running ${case_name} =========="

    mkdir -p "${case_dir}"

    rm -f "${BUILD}"/*.silo "${BUILD}"/*.csv

    jq --indent 2 \
        --argjson vin "${v}" \
        '.ball_initial_velocity.value = $vin' \
        "${INPUT_JSON}" > "${work_json}"

    cd "${BUILD}"

    "${EXE}" "${work_json}" 2>&1 | tee "${case_dir}/cabana_output.log"

    pvpython "${SILO2CSV}" 2>&1 | tee "${case_dir}/csv_conversion.log"

    mv -f "${BUILD}"/*.silo "${case_dir}/" 2>/dev/null || true
    mv -f "${BUILD}"/*.csv  "${case_dir}/" 2>/dev/null || true

    SUMMARY_FILE="${SUMMARY_FILE}" SINGLE_CASE_DIR="${case_dir}" \
        "${OCTAVE_CMD}" --no-gui --quiet --eval "source('${AVG_M}'); fflush(stdout);" \
        2>&1 | tee "${case_dir}/avg_velocity.log"

    echo "Done: ${case_dir}"
done

echo
echo "All done."
echo "Run dir: ${RUN_DIR}"
echo "Summary: ${SUMMARY_FILE}"