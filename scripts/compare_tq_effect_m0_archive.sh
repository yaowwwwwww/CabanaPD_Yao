#!/usr/bin/env bash
# Compare tq effect / thermal-softening combinations.
# For each selected parameter group, run M x tq combinations, then compare CoR/vout.

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"
INPUT_JSON_TEMPLATE="${ROOT_DIR}/examples/mechanics/inputs/simple_impact_thermal.json"
CABANAPD_EXE="${ROOT_DIR}/build/examples/mechanics/ColdSprayImpactThermal"
PY_SCRIPT="${ROOT_DIR}/scripts/silo2csv.py"
OCT_SCRIPT="${ROOT_DIR}/scripts/avg-velocity.m"
BASE_DIR="${ROOT_DIR}/build"

RUN_TAG_BASE="${RUN_TAG_BASE:-tq_effect_m0_archive_compare}"
DATE_TAG="$(date +"%Y-%m-%d_%H-%M")"
DATE_TAG_SHORT="$(date +"%Y%b%d-%H-%M")"
RUNS_ROOT="${BASE_DIR}/runs_${RUN_TAG_BASE}_${DATE_TAG}"
SUMMARY_FILE="${RUNS_ROOT}/summary_cor_hmax_recomputed_${DATE_TAG_SHORT}.txt"

# Fixed mechanical settings for this validation.
# Use list to scan M values (default: 0 and 1.09).
JC_M_LIST_STR="${JC_M_LIST_STR:-0 1.09}"
read -r -a JC_M_LIST <<< "${JC_M_LIST_STR}"
JC_N_FIXED=0.31
CZM_SCALE_FIXED=0
CZM_DECAY_FIXED=1.0
CZM_YIELD_FIXED=0.05
TQ_LIST_STR="${TQ_LIST_STR:-0 0.9}"
read -r -a TQ_LIST <<< "${TQ_LIST_STR}"
OUTPUT_FREQUENCY_FIXED="${OUTPUT_FREQUENCY_FIXED:-200}"

# Each case: short_name vin alpha beta A B C final_time
CASE_MATRIX=(
  "refB 300 6.67e-7 0.5 9e7 2.92e8 0.025 8e-8"
)

# Optional override for CASE_MATRIX.
# Format:
#   CASE_MATRIX_STR="name vin alpha beta A B C final_time;name2 vin alpha beta A B C final_time"
CASE_MATRIX_STR="${CASE_MATRIX_STR:-}"
if [ -n "${CASE_MATRIX_STR}" ]; then
  IFS=';' read -r -a CASE_MATRIX <<< "${CASE_MATRIX_STR}"
fi

# Optional alpha sweep override (space-separated), e.g.
# ALPHA_LIST_STR="1e-6 1e-5 1e-4 ... 1e+6"
ALPHA_LIST_STR="${ALPHA_LIST_STR:-}"
if [ -n "${ALPHA_LIST_STR}" ]; then
  read -r -a ALPHA_LIST <<< "${ALPHA_LIST_STR}"
fi

# Optional beta sweep override (space-separated), e.g.
# BETA_LIST_STR="0.0005 0.005 0.05 0.5 5"
BETA_LIST_STR="${BETA_LIST_STR:-}"
if [ -n "${BETA_LIST_STR}" ]; then
  read -r -a BETA_LIST <<< "${BETA_LIST_STR}"
fi

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
  local short_name="$1"
  local vin="$2"
  local lj_alpha="$3"
  local lj_beta="$4"
  local jc_a="$5"
  local jc_b="$6"
  local jc_c="$7"
  local final_time="$8"
  local tq="$9"
  local jc_m="${10}"
  local case_name="${short_name}_vin_${vin}_alpha_${lj_alpha}_beta_${lj_beta}_A_${jc_a}_B_${jc_b}_C_${jc_c}_M_${jc_m}_tq_${tq}_ft_${final_time}"
  local case_tag="run_${case_name}"
  local case_dir="${RUNS_ROOT}/${case_tag}"

  if case_exists "${case_name}"; then
    echo "SKIP (already in summary): ${case_name}"
    return 0
  fi

  echo "RUN: ${case_name}"

  if ! jq --indent 2 \
    --argjson vin       "${vin}" \
    --argjson lj_alpha  "${lj_alpha}" \
    --argjson lj_beta   "${lj_beta}" \
    --argjson jc_a      "${jc_a}" \
    --argjson jc_b      "${jc_b}" \
    --argjson jc_n      "${JC_N_FIXED}" \
    --argjson jc_c      "${jc_c}" \
    --argjson jc_m      "${jc_m}" \
    --argjson tq        "${tq}" \
    --argjson final_t   "${final_time}" \
    --argjson out_freq  "${OUTPUT_FREQUENCY_FIXED}" \
    --argjson czm_scale "${CZM_SCALE_FIXED}" \
    --argjson czm_yield "${CZM_YIELD_FIXED}" \
    --argjson czm_decay "${CZM_DECAY_FIXED}" \
    '
    .ball_initial_velocity.value   = $vin       |
    .LJalpha.value                 = $lj_alpha  |
    .LJbeta.value                  = $lj_beta   |
    .yield_stress.value            = [ $jc_a, $jc_a ] |
    .jc_B.value                    = [ $jc_b, $jc_b ] |
    .jc_n.value                    = [ $jc_n, $jc_n ] |
    .jc_C.value                    = [ $jc_c, $jc_c ] |
    .jc_m.value                    = [ $jc_m, $jc_m ] |
    .taylor_quinney.value          = $tq        |
    .final_time.value              = $final_t   |
    .output_frequency.value        = $out_freq  |
    .CZM_cohesive_scaling.value    = $czm_scale |
    .CZM_yield_stretch.value       = $czm_yield |
    .CZM_degradation_rate.value    = $czm_decay
    ' "${INPUT_JSON_TEMPLATE}" > "${WORK_INPUT_JSON}"; then
    echo "WARN: jq failed for ${case_name}"
    append_nan "${case_name}"
    return 0
  fi

  cd "${BASE_DIR}" || exit 1
  "${CABANAPD_EXE}" "${WORK_INPUT_JSON}" >> "${RUNS_ROOT}/cabana_output.log" 2>&1
  local cabana_status=$?
  if [ "${cabana_status}" -ne 0 ]; then
    echo "WARN: CabanaPD failed (${cabana_status}) for ${case_name}"
    clear_staging_outputs
    append_nan "${case_name}"
    return 0
  fi

  if ! pvpython "${PY_SCRIPT}" >/dev/null 2>&1; then
    echo "WARN: CSV conversion failed for ${case_name}"
    clear_staging_outputs
    append_nan "${case_name}"
    return 0
  fi

  mkdir -p "${case_dir}"
  cp "${WORK_INPUT_JSON}" "${case_dir}/input.json"

  shopt -s nullglob
  mv -f "${BASE_DIR}"/*.silo "${case_dir}/" 2>/dev/null || true
  mv -f "${BASE_DIR}"/*.csv  "${case_dir}/" 2>/dev/null || true
  shopt -u nullglob

  if ! run_recompute "${case_dir}"; then
    echo "WARN: Octave recompute failed for ${case_name}"
    append_nan "${case_name}"
    return 0
  fi

  return 0
}

generate_reports() {
  local compare_file="${RUNS_ROOT}/compare_tq_effect_m0.csv"
  local archive_ref_file="${RUNS_ROOT}/archive_reference_subset.csv"
  local archive_cmp_file="${RUNS_ROOT}/compare_to_archive.csv"
  export SUMMARY_FILE RUNS_ROOT compare_file archive_ref_file archive_cmp_file

  python3 - <<'PY'
import csv
import os
import re
from pathlib import Path

summary = Path(os.environ["SUMMARY_FILE"])
compare_file = Path(os.environ["compare_file"])
archive_ref_file = Path(os.environ["archive_ref_file"])
archive_cmp_file = Path(os.environ["archive_cmp_file"])

rows = []
with summary.open() as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        case_name, vin, vout, cor = p[0], float(p[1]), float(p[2]), float(p[3])
        tq_m = re.search(r"_tq_([^_]+)_ft_", case_name)
        grp_m = re.search(r"^(.*)_alpha_", case_name)
        rows.append({
            "case_name": case_name,
            "group": grp_m.group(1) if grp_m else case_name,
            "vin_abs_mps": int(abs(vin)),
            "vout": vout,
            "cor": cor,
            "tq": float(tq_m.group(1)) if tq_m else None,
        })

by_key = {}
for r in rows:
    key = (r["group"], r["vin_abs_mps"])
    by_key.setdefault(key, {})[r["tq"]] = r

pair_rows = []
for (group, vin_abs), d in sorted(by_key.items()):
    r0 = d.get(0.0)
    r9 = d.get(0.9)
    if not r0 or not r9:
        continue
    pair_rows.append({
        "group": group,
        "vin_abs_mps": vin_abs,
        "cor_tq0": r0["cor"],
        "cor_tq09": r9["cor"],
        "delta_cor_tq09_minus_tq0": r9["cor"] - r0["cor"],
        "vout_tq0": r0["vout"],
        "vout_tq09": r9["vout"],
        "delta_vout_tq09_minus_tq0": r9["vout"] - r0["vout"],
        "case_tq0": r0["case_name"],
        "case_tq09": r9["case_name"],
    })

with compare_file.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "group", "vin_abs_mps",
        "cor_tq0", "cor_tq09", "delta_cor_tq09_minus_tq0",
        "vout_tq0", "vout_tq09", "delta_vout_tq09_minus_tq0",
        "case_tq0", "case_tq09",
    ])
    w.writeheader()
    w.writerows(pair_rows)

archive_rows = []
archive_rows.append({
    "group": "refA_v400",
    "source_desc": "archive runs_adiabatic_vin400_tq0_2026-02-08 (csv-recomputed)",
    "vin_abs_mps": 400,
    "tq": 0.0,
    "cor": 0.105192,
    "vout": 42.0767,
})
archive_rows.append({
    "group": "refA_v400",
    "source_desc": "archive runs_adiabatic_vin400_tq09_2026-02-08 (csv-recomputed)",
    "vin_abs_mps": 400,
    "tq": 0.9,
    "cor": 0.106696,
    "vout": 42.6785,
})
archive_rows.append({
    "group": "refA_v400",
    "source_desc": "archive runs_recheck_vin400_tq09_2026-02-08 (csv-recomputed)",
    "vin_abs_mps": 400,
    "tq": 0.9,
    "cor": 0.0812451,
    "vout": 32.4980,
})

cor_compare = Path("/home/wuwen/program/CabanaPD_Yao/build/runs_cor_compare_tq09_final2e7_2026-02-08/cor_vs_vin_compare_tq09_final2e7.dat")
if cor_compare.exists():
    with cor_compare.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            vin, cor_no, cor_with, delta = line.split()
            vin_abs = int(abs(float(vin)))
            if vin_abs in (300, 400, 600):
                archive_rows.append({
                    "group": f"refB_v{vin_abs}",
                    "source_desc": "archive runs_cor_compare_tq09_final2e7_2026-02-08:no_soft",
                    "vin_abs_mps": vin_abs,
                    "tq": 0.9,
                    "cor": float(cor_no),
                    "vout": "",
                })

with archive_ref_file.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["group", "source_desc", "vin_abs_mps", "tq", "cor", "vout"])
    w.writeheader()
    w.writerows(archive_rows)

archive_by_key = {}
for r in archive_rows:
    key = (r["group"], r["vin_abs_mps"], r["tq"])
    archive_by_key.setdefault(key, []).append(r)

cmp_rows = []
for r in rows:
    key = (r["group"], r["vin_abs_mps"], r["tq"])
    refs = archive_by_key.get(key, [])
    if not refs:
        continue
    for ref in refs:
        ref_cor = ref["cor"]
        ref_vout = ref["vout"]
        cmp_rows.append({
            "group": r["group"],
            "vin_abs_mps": r["vin_abs_mps"],
            "tq": r["tq"],
            "new_cor": r["cor"],
            "archive_cor": ref_cor,
            "delta_cor_new_minus_archive": r["cor"] - ref_cor if ref_cor != "" else "",
            "new_vout": r["vout"],
            "archive_vout": ref_vout,
            "delta_vout_new_minus_archive": (r["vout"] - ref_vout) if ref_vout != "" else "",
            "archive_source": ref["source_desc"],
            "new_case_name": r["case_name"],
        })

with archive_cmp_file.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "group", "vin_abs_mps", "tq",
        "new_cor", "archive_cor", "delta_cor_new_minus_archive",
        "new_vout", "archive_vout", "delta_vout_new_minus_archive",
        "archive_source", "new_case_name",
    ])
    w.writeheader()
    w.writerows(cmp_rows)

print(compare_file)
print(archive_ref_file)
print(archive_cmp_file)
PY
}

require_cmd jq
require_cmd pvpython
require_cmd awk
require_cmd python3

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

mkdir -p "${RUNS_ROOT}"
{
  printf "# params: CASE_MATRIX=%s; JC_M_LIST=%s; TQ_LIST=%s; JC_N_FIXED=%s; CZM=(%s,%s,%s); OUTPUT_FREQUENCY_FIXED=%s; ALPHA_LIST_STR=%s; BETA_LIST_STR=%s\n" \
    "${CASE_MATRIX[*]}" \
    "${JC_M_LIST[*]}" \
    "${TQ_LIST[*]}" \
    "${JC_N_FIXED}" \
    "${CZM_SCALE_FIXED}" "${CZM_YIELD_FIXED}" "${CZM_DECAY_FIXED}" \
    "${OUTPUT_FREQUENCY_FIXED}" \
    "${ALPHA_LIST_STR:-<default>}" \
    "${BETA_LIST_STR:-<default>}"
  printf "# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_max(m)    best_frame    h_residual(m)    A_residual(m2)    h_mean_residual(m)\n"
} > "${SUMMARY_FILE}"

WORK_INPUT_JSON="${RUNS_ROOT}/_current_input.json"
cleanup() {
  rm -f "${WORK_INPUT_JSON}" "${BASE_DIR}/cabana_last_output.txt" 2>/dev/null || true
}
trap cleanup EXIT

echo "==== tq effect validation start ($(date)) ===="
echo "Run root: ${RUNS_ROOT}"
echo "Summary: ${SUMMARY_FILE}"
echo "M list: ${JC_M_LIST[*]}, jc_n=${JC_N_FIXED}, CZM=(${CZM_SCALE_FIXED},${CZM_YIELD_FIXED},${CZM_DECAY_FIXED})"
echo "tq list: ${TQ_LIST[*]}"
echo

for case_spec in "${CASE_MATRIX[@]}"; do
  # shellcheck disable=SC2086
  set -- ${case_spec}
  short_name="$1"
  vin="$2"
  lj_alpha="$3"
  lj_beta="$4"
  jc_a="$5"
  jc_b="$6"
  jc_c="$7"
  final_time="$8"
  if [ -n "${ALPHA_LIST_STR}" ]; then
    CURRENT_ALPHAS=("${ALPHA_LIST[@]}")
  else
    CURRENT_ALPHAS=("${lj_alpha}")
  fi
  if [ -n "${BETA_LIST_STR}" ]; then
    CURRENT_BETAS=("${BETA_LIST[@]}")
  else
    CURRENT_BETAS=("${lj_beta}")
  fi
  for alpha_val in "${CURRENT_ALPHAS[@]}"; do
    for beta_val in "${CURRENT_BETAS[@]}"; do
      for jc_m in "${JC_M_LIST[@]}"; do
        for tq in "${TQ_LIST[@]}"; do
          run_one_case "${short_name}" "${vin}" "${alpha_val}" "${beta_val}" "${jc_a}" "${jc_b}" "${jc_c}" "${final_time}" "${tq}" "${jc_m}"
        done
      done
    done
  done
done

generate_reports

echo
echo "Validation complete."
echo "Summary: ${SUMMARY_FILE}"
echo "Pair comparison: ${RUNS_ROOT}/compare_tq_effect_m0.csv"
echo "Archive refs: ${RUNS_ROOT}/archive_reference_subset.csv"
echo "New vs archive: ${RUNS_ROOT}/compare_to_archive.csv"
echo "Cases archived under: ${RUNS_ROOT}"
