#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 CASE_DIR [threads]" >&2
  echo "Example: $0 examples/test4_Yunnan/voro 4" >&2
  exit 2
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$1"
THREADS="${2:-${DAZI_THREADS:-4}}"
EXE="${DAZI_EXE:-$ROOT/code/bin/DAzimSurfTomo}"

if [ ! -x "$EXE" ]; then
  echo "ERROR: executable not found: $EXE" >&2
  echo "Run ./build.sh first, or set DAZI_EXE=/path/to/DAzimSurfTomo" >&2
  exit 2
fi

if [ ! -f "$CASE_DIR/para.in" ]; then
  echo "ERROR: missing para.in in $CASE_DIR" >&2
  exit 2
fi

export OMP_NUM_THREADS="$THREADS"
export DAZI_THREADS="$THREADS"

cd "$CASE_DIR"
"$EXE" para.in | tee para.in_inv.log

