#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

LOG="logs/mason_orchestrator.log"
REF="data/refs/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SA="${REF}.sa"
mkdir -p logs data/refs/grch38

exec >>"$LOG" 2>&1

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
