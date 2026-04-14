#!/usr/bin/env bash
# Run the current Eq.16 drag model at 200 and 600 m/s using the active
# constitutive inputs from simple_impact_thermal.json.

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"

export RUN_TAG_BASE="${RUN_TAG_BASE:-eq16_two_case_v200_v600}"
export CASE_MATRIX_STR="${CASE_MATRIX_STR:-\
v200 200 1e-6 0.5 9e7 2.92e8 0.025 8e-8;\
v600 600 1e-6 0.5 9e7 2.92e8 0.025 8e-8}"

export JC_M_LIST_STR="${JC_M_LIST_STR:-1.09}"
export TQ_LIST_STR="${TQ_LIST_STR:-0}"
export JC_C2_FIXED="${JC_C2_FIXED:-0.908}"
export JC_EPSDOT_U_FIXED="${JC_EPSDOT_U_FIXED:-680000.0}"
export DRAG_BD_LIST_STR="${DRAG_BD_LIST_STR:-1.0e-5}"
export DRAG_BD_AL_FIXED="${DRAG_BD_AL_FIXED:-1.0e-5}"
export DRAG_BD_CU_FIXED="${DRAG_BD_CU_FIXED:-1.0e-5}"
export DRAG_RHO_MOBILE_AL_FIXED="${DRAG_RHO_MOBILE_AL_FIXED:-1.0e13}"
export DRAG_RHO_MOBILE_CU_FIXED="${DRAG_RHO_MOBILE_CU_FIXED:-1.0e13}"
export DRAG_BURGERS_AL_FIXED="${DRAG_BURGERS_AL_FIXED:-2.56e-10}"
export DRAG_BURGERS_CU_FIXED="${DRAG_BURGERS_CU_FIXED:-2.56e-10}"
export OUTPUT_FREQUENCY_FIXED="${OUTPUT_FREQUENCY_FIXED:-200}"

exec bash "${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"
