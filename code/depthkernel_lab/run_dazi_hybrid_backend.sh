#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 CASE_DIR [OMP_THREADS]" >&2
  exit 2
fi

case_dir="$1"
threads="${2:-${DAZI_THREADS:-4}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec python3 "$script_dir/run_dazi_configured.py" "$case_dir" --threads "$threads"
