#!/usr/bin/env bash
# Thin old-style wrapper for scanning Taylor-Quinney coefficient values.
# It fixes the requested velocity list and tq list, then hands control to the
# existing single-driver scan script.

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"

export RUN_TAG_BASE="${RUN_TAG_BASE:-scan_tq_coeff_v100_700}"

# Requested velocity list.
export CASE_MATRIX_STR="${CASE_MATRIX_STR:-\
v100 100 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v300 300 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v450 450 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v500 500 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v550 550 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v600 600 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v650 650 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v700 700 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8}"

# Requested tq scan values.
export TQ_LIST_STR="${TQ_LIST_STR:-0.1 0.3 0.5 0.8}"

# Keep other model settings fixed by default.
export JC_M_LIST_STR="${JC_M_LIST_STR:-1.09}"
export DRAG_K0_LIST_STR="${DRAG_K0_LIST_STR:-1e8}"
export JC_N_FIXED="${JC_N_FIXED:-0.31}"
export JC_C2_FIXED="${JC_C2_FIXED:-0.908}"
export JC_EPSDOT_U_FIXED="${JC_EPSDOT_U_FIXED:-680000.0}"
export DRAG_M_FIXED="${DRAG_M_FIXED:-0.008}"
export DRAG_A_FIXED="${DRAG_A_FIXED:-1.0}"
export DRAG_BETA_G_FIXED="${DRAG_BETA_G_FIXED:-0.9}"
export OUTPUT_FREQUENCY_FIXED="${OUTPUT_FREQUENCY_FIXED:-200}"

exec bash "${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"
