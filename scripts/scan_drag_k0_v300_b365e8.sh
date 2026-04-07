#!/usr/bin/env bash
# Single-velocity drag-K0 scan on the drag-enabled thermal impact path.

set -u
set -o pipefail

ROOT_DIR="/home/wuwen/program/CabanaPD_Yao"

export RUN_TAG_BASE="${RUN_TAG_BASE:-scan_drag_k0_v300_b365e8}"

export CASE_MATRIX_STR="${CASE_MATRIX_STR:-\
v300 300 6.67e-7 0.5 9e7 3.65e8 0.025 8e-8}"

export DRAG_K0_LIST_STR="${DRAG_K0_LIST_STR:-1e7 1e8 1e9}"

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
