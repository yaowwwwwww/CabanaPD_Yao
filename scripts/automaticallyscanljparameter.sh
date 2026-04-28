#!/usr/bin/env bash
set -uo pipefail

ROOT="/home/wuwen/program/CabanaPD_Yao"
BUILD="${ROOT}/build"

INPUT_JSON="${ROOT}/examples/mechanics/inputs/simple_impact_thermal.json"
EXE="${BUILD}/examples/mechanics/ColdSprayImpactThermal"
SILO2CSV="${ROOT}/scripts/silo2csv.py"
AVG_M="${ROOT}/scripts/avg-velocity.m"

alpha="$1"
beta="$2"
v="$3"

RUN_ROOT="${BUILD}/opt_lj_runs"
case_name="v_${v}_a_${alpha}_b_${beta}"
case_dir="${RUN_ROOT}/${case_name}"
work_json="${case_dir}/input.json"
summary_tmp="${case_dir}/summary.txt"

mkdir -p "${case_dir}"

rm -f "${BUILD}"/*.silo "${BUILD}"/*.csv

jq --indent 2 \
    --argjson vin "${v}" \
    --argjson alpha "${alpha}" \
    --argjson beta "${beta}" \
    '
    .ball_initial_velocity.value = $vin
    | .LJalpha.value = $alpha
    | .LJbeta.value  = $beta
    ' \
    "${INPUT_JSON}" > "${work_json}"

cd "${BUILD}" || {
    echo "CoR=nan" | tee "${case_dir}/result.txt"
    exit 0
}

"${EXE}" "${work_json}" > "${case_dir}/cabana_output.log" 2>&1
if [ $? -ne 0 ]; then
    echo "CoR=nan" | tee "${case_dir}/result.txt"
    exit 0
fi

pvpython "${SILO2CSV}" > "${case_dir}/csv_conversion.log" 2>&1
if [ $? -ne 0 ]; then
    echo "CoR=nan" | tee "${case_dir}/result.txt"
    exit 0
fi

mv -f "${BUILD}"/*.silo "${case_dir}/" 2>/dev/null || true
mv -f "${BUILD}"/*.csv  "${case_dir}/" 2>/dev/null || true

SUMMARY_FILE="${summary_tmp}" SINGLE_CASE_DIR="${case_dir}" \
    octave-cli --no-gui --quiet --eval "source('${AVG_M}'); fflush(stdout);" \
    > "${case_dir}/avg_velocity.log" 2>&1

if [ $? -ne 0 ] || [ ! -f "${summary_tmp}" ]; then
    echo "CoR=nan" | tee "${case_dir}/result.txt"
    exit 0
fi

COR=$(awk 'NF && $1 !~ /^#/ {print $4; exit}' "${summary_tmp}")

if [ -z "${COR}" ]; then
    COR="nan"
fi

echo "CoR=${COR}" | tee "${case_dir}/result.txt"