#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

LOG="logs/mason_orchestrator.log"
REF="data/refs/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SA="${REF}.sa"
STATE_DIR="data/interim/mason/state"
LOCK_DIR="$STATE_DIR/orchestrator.lock"
mkdir -p logs data/refs/grch38 "$STATE_DIR"

release_orchestrator_lock() {
  rm -rf "$LOCK_DIR"
}

acquire_orchestrator_lock() {
  if mkdir "$LOCK_DIR" 2>/dev/null; then
    echo "$$" > "$LOCK_DIR/pid"
    return
  fi

  local existing_pid=""
  if [[ -f "$LOCK_DIR/pid" ]]; then
    existing_pid="$(cat "$LOCK_DIR/pid" 2>/dev/null || true)"
  fi

  if [[ -n "$existing_pid" ]] && ps -p "$existing_pid" >/dev/null 2>&1; then
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Orchestrator already running (pid: $existing_pid); skipping duplicate launch"
    exit 0
  fi

  rm -rf "$LOCK_DIR"
  if ! mkdir "$LOCK_DIR" 2>/dev/null; then
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] ERROR: unable to acquire orchestrator lock at $LOCK_DIR"
    exit 1
  fi
  echo "$$" > "$LOCK_DIR/pid"
}

exec >>"$LOG" 2>&1
trap release_orchestrator_lock EXIT INT TERM
acquire_orchestrator_lock

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Orchestrator started"

if [[ ! -s "$REF" ]]; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] ERROR: reference FASTA missing: $REF"
  exit 1
fi

while [[ ! -s "$SA" ]]; do
  if pgrep -f "bwa index $REF" >/dev/null 2>&1; then
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Waiting for active bwa index to finish"
    sleep 60
    continue
  fi

  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Launching/repairing bwa index"
  bwa index "$REF"
done

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Reference index complete"
echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting full Mason pipeline under caffeinate"

caffeinate -dimsu bash scripts/run_mason_full_analysis.sh

status=$?
echo "[$(date +"%Y-%m-%d %H:%M:%S")] Full pipeline exit status: $status"
exit $status
