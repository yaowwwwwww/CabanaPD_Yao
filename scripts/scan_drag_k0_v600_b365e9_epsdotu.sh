#!/usr/bin/env bash
# Two-case drag-Bd scan at vin=600 m/s for selected epsdot_u values.

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"

export RUN_TAG_BASE="${RUN_TAG_BASE:-scan_drag_bd_v600_b365e9_epsdotu}"

export CASE_MATRIX_STR="${CASE_MATRIX_STR:-\
v600 600 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8}"

export DRAG_BD_LIST_STR="${DRAG_BD_LIST_STR:-1e-5}"
export JC_EPSDOT_U_LIST_STR="${JC_EPSDOT_U_LIST_STR:-6.8e6 4e7}"

export JC_M_LIST_STR="${JC_M_LIST_STR:-1.09}"
export TQ_LIST_STR="${TQ_LIST_STR:-0.9}"
export JC_N_FIXED="${JC_N_FIXED:-0.31}"
export JC_C2_FIXED="${JC_C2_FIXED:-0.908}"
export DRAG_RHO_MOBILE_FIXED="${DRAG_RHO_MOBILE_FIXED:-1.0e13}"
export DRAG_BURGERS_FIXED="${DRAG_BURGERS_FIXED:-2.56e-10}"
export OUTPUT_FREQUENCY_FIXED="${OUTPUT_FREQUENCY_FIXED:-200}"

exec bash "${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"
