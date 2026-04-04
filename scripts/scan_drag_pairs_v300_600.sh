#!/usr/bin/env bash
# Thin old-style wrapper for scanning drag K0 values.
# It only exports the desired parameter lists and then hands control to the
# existing single-driver scan script.

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"

export RUN_TAG_BASE="${RUN_TAG_BASE:-scan_drag_k0_v100_700}"

# Current baseline for the log-drag model.
export CASE_MATRIX_STR="${CASE_MATRIX_STR:-\
v100 100 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v150 150 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v200 200 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v250 250 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v300 300 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v350 350 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v400 400 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v450 450 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v500 500 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v550 550 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v600 600 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v650 650 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8;\
v700 700 6.67e-7 0.5 9e7 3.65e9 0.025 8e-8}"

# Current model: scan K0 only.
export DRAG_K0_LIST_STR="${DRAG_K0_LIST_STR:-1e8 3e8 5e8 8e8}"

# Fixed settings.
export JC_M_LIST_STR="${JC_M_LIST_STR:-1.09}"
export TQ_LIST_STR="${TQ_LIST_STR:-0.9}"
export JC_N_FIXED="${JC_N_FIXED:-0.31}"
export JC_C2_FIXED="${JC_C2_FIXED:-0.908}"
export JC_EPSDOT_U_FIXED="${JC_EPSDOT_U_FIXED:-680000.0}"
export DRAG_M_FIXED="${DRAG_M_FIXED:-0.008}"
export DRAG_A_FIXED="${DRAG_A_FIXED:-1.0}"
export DRAG_BETA_G_FIXED="${DRAG_BETA_G_FIXED:-0.9}"
export OUTPUT_FREQUENCY_FIXED="${OUTPUT_FREQUENCY_FIXED:-200}"

exec bash "${ROOT_DIR}/scripts/compare_tq_effect_m0_archive.sh"
