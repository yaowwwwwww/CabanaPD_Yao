#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 RUN_ROOT SUMMARY_FILE [POLL_SECONDS]" >&2
  exit 2
fi

RUN_ROOT="$1"
SUMMARY_FILE="$2"
POLL_SECONDS="${3:-60}"
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SCAN_WRAPPER="$ROOT_DIR/scripts/scan_drag_pairs_v300_600.sh"
WATCH_LOG="$RUN_ROOT/watch_resume.log"
LAUNCH_LOG="$RUN_ROOT/watch_resume_launch.log"
LOCKDIR="$RUN_ROOT/.watch_resume_lock"

if [[ -n "${TOTAL_CASES:-}" ]]; then
  TOTAL_CASES="${TOTAL_CASES}"
else
  case_count=13
  k0_count=4
  TOTAL_CASES=$(( case_count * k0_count ))
fi

mkdir -p "$RUN_ROOT"

timestamp() { date '+%F %T'; }
count_rows() {
  awk '!/^#/ && NF {c++} END{print c+0}' "$SUMMARY_FILE" 2>/dev/null || echo 0
}
scan_active() {
  pgrep -af "ColdSprayImpactThermal .*${RUN_ROOT}/_current_input\.json" >/dev/null || \
  pgrep -af "compare_tq_effect_m0_archive\.sh" | grep -F "$RUN_ROOT" >/dev/null
}

{
  echo "[$(timestamp)] watchdog start run_root=$RUN_ROOT summary=$SUMMARY_FILE poll=${POLL_SECONDS}s total_cases=$TOTAL_CASES"
  while true; do
    rows="$(count_rows)"
    if [[ "$rows" -ge "$TOTAL_CASES" ]]; then
      echo "[$(timestamp)] complete rows=$rows/$TOTAL_CASES; watchdog exit"
      exit 0
    fi

    if scan_active; then
      echo "[$(timestamp)] active rows=$rows/$TOTAL_CASES"
      sleep "$POLL_SECONDS"
      continue
    fi

    echo "[$(timestamp)] inactive rows=$rows/$TOTAL_CASES; attempting resume"
    if mkdir "$LOCKDIR" 2>/dev/null; then
      {
        echo "[$(timestamp)] resume launch start"
        RUNS_ROOT_OVERRIDE="$RUN_ROOT" \
        SUMMARY_FILE_OVERRIDE="$SUMMARY_FILE" \
        bash -x "$SCAN_WRAPPER"
        rc=$?
        echo "[$(timestamp)] resume launch exit rc=$rc"
      } >>"$LAUNCH_LOG" 2>&1 || true
      rmdir "$LOCKDIR" 2>/dev/null || true
    else
      echo "[$(timestamp)] resume skipped; lock exists"
    fi
    sleep "$POLL_SECONDS"
  done
} >>"$WATCH_LOG" 2>&1
