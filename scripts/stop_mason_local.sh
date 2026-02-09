#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PID_FILE="logs/mason_local_runner.pid"

stop_pid() {
  local pid="$1"
  if [[ -n "${pid:-}" ]] && ps -p "$pid" >/dev/null 2>&1; then
    kill "$pid" >/dev/null 2>&1 || true
  fi
}

if [[ -f "$PID_FILE" ]]; then
  pid="$(cat "$PID_FILE" 2>/dev/null || true)"
  stop_pid "$pid"
fi

pkill -f "run_mason_orchestrator.sh" >/dev/null 2>&1 || true
pkill -f "run_mason_full_analysis.sh" >/dev/null 2>&1 || true
pkill -f "bwa mem" >/dev/null 2>&1 || true
pkill -f "samtools view" >/dev/null 2>&1 || true
pkill -f "samtools sort" >/dev/null 2>&1 || true
pkill -f "samtools fixmate" >/dev/null 2>&1 || true
pkill -f "samtools markdup" >/dev/null 2>&1 || true
pkill -f "bcftools" >/dev/null 2>&1 || true
pkill -f "snpEff" >/dev/null 2>&1 || true

sleep 1

if ps -Ao command | rg -q "run_mason_orchestrator.sh|run_mason_full_analysis.sh|bwa mem|samtools view|samtools sort|samtools fixmate|samtools markdup|bcftools|snpEff"; then
  pkill -9 -f "run_mason_orchestrator.sh|run_mason_full_analysis.sh|bwa mem|samtools view|samtools sort|samtools fixmate|samtools markdup|bcftools|snpEff" >/dev/null 2>&1 || true
fi

rm -f "$PID_FILE"

echo "Stopped Mason local run processes (if any were active)."
