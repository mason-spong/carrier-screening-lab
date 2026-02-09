#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PID_FILE="logs/mason_local_runner.pid"
LOG_FILE="logs/mason_local_runner.log"

mkdir -p logs

if [[ -f "$PID_FILE" ]]; then
  old_pid="$(cat "$PID_FILE" 2>/dev/null || true)"
  if [[ -n "${old_pid:-}" ]] && ps -p "$old_pid" >/dev/null 2>&1; then
    echo "Mason orchestrator is already running (pid: $old_pid)."
    echo "Log: $LOG_FILE"
    exit 0
  fi
fi

nohup caffeinate -dimsu bash scripts/run_mason_orchestrator.sh > "$LOG_FILE" 2>&1 &
echo $! > "$PID_FILE"

sleep 1
new_pid="$(cat "$PID_FILE")"
if ps -p "$new_pid" >/dev/null 2>&1; then
  echo "Started Mason orchestrator (pid: $new_pid)."
  echo "Log: $LOG_FILE"
  echo "Check: ps -p $new_pid -o pid,etime,command"
  echo "Tail:  tail -n 40 logs/mason_orchestrator.log"
else
  echo "Failed to start orchestrator. Check $LOG_FILE"
  exit 1
fi
