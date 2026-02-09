#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PORT="${1:-8788}"

echo "Starting Mason status monitor on http://127.0.0.1:${PORT}/"
echo "Press Ctrl+C to stop."
python3 scripts/serve_status_ui_maximalist.py "$PORT"
