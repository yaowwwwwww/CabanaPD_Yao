#!/usr/bin/env bash

set -e
set -u
set -o pipefail

ROOT="/home/wuwen/program/CabanaPD_Yao"
BUILD="${ROOT}/build"

INPUT_JSON="${ROOT}/examples/mechanics/inputs/Lucas-2021-IJMS.json"
EXE="${BUILD}/examples/mechanics/ColdSprayImpactThermal"
SILO2CSV="${ROOT}/scripts/silo2csv.py"
AVG_M="${ROOT}/scripts/avg-velocity.m"

# 不在脚本里强制写死线程；如果外面没设置，默认用 8
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-60}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"

if [ ! -f "${AVG_M}" ]; then
    AVG_M="${ROOT}/src/avg-velocity.m"
fi

DATE_TAG="$(date +"%Y-%m-%d_%H-%M")"
RUN_DIR="${BUILD}/${DATE_TAG}_runs_json_target_v100_600"
SUMMARY_FILE="${RUN_DIR}/summary_velocity_scan.txt"

mkdir -p "${RUN_DIR}"

echo "Using OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "Using OMP_PROC_BIND=${OMP_PROC_BIND}"
echo "Using OMP_PLACES=${OMP_PLACES}"

cat > "${SUMMARY_FILE}" <<EOF
# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_surface(m)    h_penetration(m)    A_residual(m2)
EOF

if command -v octave-cli >/dev/null 2>&1; then
    OCTAVE_CMD="octave-cli"
elif command -v octave >/dev/null 2>&1; then
    OCTAVE_CMD="octave"
else
    echo "ERROR: octave or octave-cli not found"
    exit 1
fi

for alpha in 3e-5; do
for beta in 0.5; do

for v in 200; do

    case_name="v_${v}_a_${alpha}_b_${beta}"
    case_dir="${RUN_DIR}/${case_name}"
    work_json="${case_dir}/input.json"

    echo
    echo "========== Running ${case_name} =========="

    mkdir -p "${case_dir}"

    # 清理 case_dir
    rm -f "${case_dir}"/*.silo
    rm -f "${case_dir}"/*.csv
    rm -f "${case_dir}"/*.log

    # 生成 input
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

    # 回到 BUILD，保持和以前稳定运行方式一致
    cd "${BUILD}"

    # 关键：每个 case 前彻底清理 BUILD 中旧输出
    rm -f particles_*.silo
    rm -f particles_*_all.csv
    rm -f avg_vmag.csv
    rm -f csv_conversion.log
    rm -f cabana_output.log

    # 运行 EXE；失败也不要中断整个扫描
    set +e
    "${EXE}" "${work_json}" 2>&1 | tee "${case_dir}/cabana_output.log"
    run_status=${PIPESTATUS[0]}
    set -e

    if [ "${run_status}" -ne 0 ]; then
        echo "FAILED: ${case_name}, exit=${run_status}"
        echo "${case_name}  ${v}  NaN  NaN  NaN  NaN  NaN  NaN" >> "${SUMMARY_FILE}"

        # 保存可能残留的部分输出，方便检查
        mv -f particles_*.silo "${case_dir}/" 2>/dev/null || true
        mv -f particles_*_all.csv "${case_dir}/" 2>/dev/null || true
        continue
    fi

    # 在 BUILD 里转换当前 case 的 silo
    pvpython "${SILO2CSV}" 2>&1 | tee "${case_dir}/csv_conversion.log"

    # 转换完成后，马上把文件移到 case_dir
    mv -f particles_*.silo "${case_dir}/" 2>/dev/null || true
    mv -f particles_*_all.csv "${case_dir}/" 2>/dev/null || true

    # Octave 只分析这个 case_dir
    SUMMARY_FILE="${SUMMARY_FILE}" SINGLE_CASE_DIR="${case_dir}" \
        "${OCTAVE_CMD}" --no-gui --quiet --eval "source('${AVG_M}'); fflush(stdout);" \
        2>&1 | tee "${case_dir}/avg_velocity.log"

    # 再清理 BUILD，避免影响下一个 case
    cd "${BUILD}"
    rm -f particles_*.silo
    rm -f particles_*_all.csv
    rm -f avg_vmag.csv

    echo "Done: ${case_dir}"

done
done
done

echo
echo "All done."
echo "Run dir: ${RUN_DIR}"
echo "Summary: ${SUMMARY_FILE}"