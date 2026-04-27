#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/home/wuwen/program/CabanaPD_Yao}"
DRIVER="${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"

export ROOT_DIR
export BASE_DIR="${BASE_DIR:-${ROOT_DIR}/build}"
export CABANAPD_EXE="${CABANAPD_EXE:-${ROOT_DIR}/install/bin/ColdSprayImpactThermal}"
export RUN_TAG_BASE="${RUN_TAG_BASE:-scan_alpha_beta_v300}"

# Single-velocity Cu/Cu thermal impact baseline taken from simple_impact_thermal.json.
export CASE_MATRIX_STR="${CASE_MATRIX_STR:-ab300 300 1.5e-6 0.5 9e7 2.92e8 0.025 8e-8}"

# Moderate 2D scan around the current baseline values.
export ALPHA_LIST_STR="${ALPHA_LIST_STR:-7e-7 1e-6 1.5e-6}"
export BETA_LIST_STR="${BETA_LIST_STR:-0.1 0.3 0.5 0.8}"

export JC_M_LIST_STR="${JC_M_LIST_STR:-1.09}"
export TQ_LIST_STR="${TQ_LIST_STR:-0}"
export JC_N_FIXED="${JC_N_FIXED:-0.31}"
export JC_C2_FIXED="${JC_C2_FIXED:-0.908}"
export JC_EPSDOT_U_FIXED="${JC_EPSDOT_U_FIXED:-680000.0}"
export DRAG_BD_FIXED="${DRAG_BD_FIXED:-1.0e-5}"
export DRAG_RHO_MOBILE_FIXED="${DRAG_RHO_MOBILE_FIXED:-1.0e13}"
export DRAG_BURGERS_FIXED="${DRAG_BURGERS_FIXED:-2.56e-10}"
export OUTPUT_FREQUENCY_FIXED="${OUTPUT_FREQUENCY_FIXED:-200}"
export TIMESTEP_FIXED="${TIMESTEP_FIXED:-1e-11}"

exec bash "${DRIVER}"
