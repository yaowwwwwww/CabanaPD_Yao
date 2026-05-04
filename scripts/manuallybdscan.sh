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

# Keep the same manual-scan style as manuallyvelscan.sh.
# If the caller does not set threads, use a conservative default.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"

if [ ! -f "${AVG_M}" ]; then
    AVG_M="${ROOT}/src/avg-velocity.m"
fi

DATE_TAG="$(date +"%Y-%m-%d_%H-%M")"
RUN_DIR="${BUILD}/${DATE_TAG}_runs_manual_bd_scan"
SUMMARY_FILE="${RUN_DIR}/summary_bd_scan.txt"

mkdir -p "${RUN_DIR}"

echo "Using OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "Using OMP_PROC_BIND=${OMP_PROC_BIND}"
echo "Using OMP_PLACES=${OMP_PLACES}"

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

# Fixed baseline parameters for Bd scan. Override from the command line if needed:
#   V_LIST_STR="300 600" BD_LIST_STR="0 1e-5 2e-5" bash scripts/manuallybdscan.sh
ALPHA_LIST_STR="${ALPHA_LIST_STR:-6.67e-7}"
BETA_LIST_STR="${BETA_LIST_STR:-0.5}"
V_LIST_STR="${V_LIST_STR:-300 600}"
BD_LIST_STR="${BD_LIST_STR:-0 2e-6 5e-6 1e-5 2e-5 5e-5}"

JC_EPSDOT_U="${JC_EPSDOT_U:-6.8e6}"
TAYLOR_QUINNEY="${TAYLOR_QUINNEY:-0.9}"
DRAG_RHO_MOBILE_AL="${DRAG_RHO_MOBILE_AL:-1.0e13}"
DRAG_RHO_MOBILE_CU="${DRAG_RHO_MOBILE_CU:-1.0e13}"
DRAG_BURGERS_AL="${DRAG_BURGERS_AL:-2.86e-10}"
DRAG_BURGERS_CU="${DRAG_BURGERS_CU:-2.56e-10}"
OUTPUT_FREQUENCY="${OUTPUT_FREQUENCY:-200}"

read -r -a ALPHA_LIST <<< "${ALPHA_LIST_STR}"
read -r -a BETA_LIST <<< "${BETA_LIST_STR}"
read -r -a V_LIST <<< "${V_LIST_STR}"
read -r -a BD_LIST <<< "${BD_LIST_STR}"

echo "alpha list: ${ALPHA_LIST[*]}"
echo "beta list: ${BETA_LIST[*]}"
echo "velocity list: ${V_LIST[*]}"
echo "drag_Bd list: ${BD_LIST[*]}"
echo "jc_epsdot_u=${JC_EPSDOT_U}, taylor_quinney=${TAYLOR_QUINNEY}"
echo "drag rho mobile: Al=${DRAG_RHO_MOBILE_AL}, Cu=${DRAG_RHO_MOBILE_CU}"
echo "drag burgers: Al=${DRAG_BURGERS_AL}, Cu=${DRAG_BURGERS_CU}"
echo

for alpha in "${ALPHA_LIST[@]}"; do
for beta in "${BETA_LIST[@]}"; do
for bd in "${BD_LIST[@]}"; do
for v in "${V_LIST[@]}"; do

    case_name="v_${v}_a_${alpha}_b_${beta}_Bd_${bd}_U_${JC_EPSDOT_U}_tq_${TAYLOR_QUINNEY}"
    case_dir="${RUN_DIR}/${case_name}"
    work_json="${case_dir}/input.json"

    echo
    echo "========== Running ${case_name} =========="

    mkdir -p "${case_dir}"

    # Clean only this case directory.
    rm -f "${case_dir}"/*.silo
    rm -f "${case_dir}"/*.csv
    rm -f "${case_dir}"/*.log

    jq --indent 2 \
        --argjson vin "${v}" \
        --argjson alpha "${alpha}" \
        --argjson beta "${beta}" \
        --argjson bd "${bd}" \
        --argjson epsdot_u "${JC_EPSDOT_U}" \
        --argjson tq "${TAYLOR_QUINNEY}" \
        --argjson drag_rho_al "${DRAG_RHO_MOBILE_AL}" \
        --argjson drag_rho_cu "${DRAG_RHO_MOBILE_CU}" \
        --argjson drag_b_al "${DRAG_BURGERS_AL}" \
        --argjson drag_b_cu "${DRAG_BURGERS_CU}" \
        --argjson out_freq "${OUTPUT_FREQUENCY}" \
        '
        .ball_initial_velocity.value = $vin
        | .LJalpha.value = $alpha
        | .LJbeta.value  = $beta
        | .drag_Bd.value = [ $bd, $bd ]
        | .drag_mobile_dislocation_density.value = [ $drag_rho_al, $drag_rho_cu ]
        | .drag_burgers_vector.value = [ $drag_b_al, $drag_b_cu ]
        | .jc_epsdot_u.value = [ $epsdot_u, $epsdot_u ]
        | .taylor_quinney.value = $tq
        | .output_frequency.value = $out_freq
        ' \
        "${INPUT_JSON}" > "${work_json}"

    cd "${BUILD}"

    # Clean staging files before each case.
    rm -f particles_*.silo
    rm -f particles_*_all.csv
    rm -f avg_vmag.csv
    rm -f csv_conversion.log
    rm -f cabana_output.log

    set +e
    "${EXE}" "${work_json}" 2>&1 | tee "${case_dir}/cabana_output.log"
    run_status=${PIPESTATUS[0]}
    set -e

    if [ "${run_status}" -ne 0 ]; then
        echo "FAILED: ${case_name}, exit=${run_status}"
        echo "${case_name}  ${v}  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN" >> "${SUMMARY_FILE}"

        # Save any partial outputs for debugging.
        mv -f particles_*.silo "${case_dir}/" 2>/dev/null || true
        mv -f particles_*_all.csv "${case_dir}/" 2>/dev/null || true
        continue
    fi

    pvpython "${SILO2CSV}" 2>&1 | tee "${case_dir}/csv_conversion.log"

    mv -f particles_*.silo "${case_dir}/" 2>/dev/null || true
    mv -f particles_*_all.csv "${case_dir}/" 2>/dev/null || true

    SUMMARY_FILE="${SUMMARY_FILE}" SINGLE_CASE_DIR="${case_dir}" \
        "${OCTAVE_CMD}" --no-gui --quiet --eval "source('${AVG_M}'); fflush(stdout);" \
        2>&1 | tee "${case_dir}/avg_velocity.log"

    cd "${BUILD}"
    rm -f particles_*.silo
    rm -f particles_*_all.csv
    rm -f avg_vmag.csv

    echo "Done: ${case_dir}"

done
done
done
done

echo
echo "All done."
echo "Run dir: ${RUN_DIR}"
echo "Summary: ${SUMMARY_FILE}"
