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
WATCH_SCRIPT="$ROOT_DIR/scripts/watch_resume_scan.sh"

mkdir -p "$RUN_ROOT"

SCAN_PID_FILE="$RUN_ROOT/daemon_scan.pid"
WATCH_PID_FILE="$RUN_ROOT/daemon_watch.pid"
DAEMON_LOG="$RUN_ROOT/daemon_resume_scan.log"

# Stop stale pidfile references.
for f in "$SCAN_PID_FILE" "$WATCH_PID_FILE"; do
  if [[ -f "$f" ]]; then
    pid="$(cat "$f" 2>/dev/null || true)"
    if [[ -n "${pid:-}" ]] && ! kill -0 "$pid" 2>/dev/null; then
      rm -f "$f"
    fi
  fi
done

# Launch detached watchdog if not already alive.
if [[ ! -f "$WATCH_PID_FILE" ]] || ! kill -0 "$(cat "$WATCH_PID_FILE")" 2>/dev/null; then
  setsid bash -lc "exec '$WATCH_SCRIPT' '$RUN_ROOT' '$SUMMARY_FILE' '$POLL_SECONDS'" \
    </dev/null >>"$DAEMON_LOG" 2>&1 &
  echo $! > "$WATCH_PID_FILE"
fi

# If no active scan for this run root, trigger one detached resume immediately.
if ! pgrep -af "ColdSprayImpactThermal .*${RUN_ROOT}/_current_input\\.json" >/dev/null && \
   ! pgrep -af "compare_tq_effect_m0_archive\\.sh" | grep -F "$RUN_ROOT" >/dev/null; then
  setsid bash -lc "cd '$ROOT_DIR' && export RUNS_ROOT_OVERRIDE='$RUN_ROOT' SUMMARY_FILE_OVERRIDE='$SUMMARY_FILE'; exec bash scripts/scan_drag_pairs_v300_600.sh" \
    </dev/null >>"$DAEMON_LOG" 2>&1 &
  echo $! > "$SCAN_PID_FILE"
fi

printf 'run_root=%s\n' "$RUN_ROOT"
printf 'watch_pid=%s\n' "$(cat "$WATCH_PID_FILE")"
if [[ -f "$SCAN_PID_FILE" ]]; then
  printf 'scan_pid=%s\n' "$(cat "$SCAN_PID_FILE")"
fi
printf 'log=%s\n' "$DAEMON_LOG"
