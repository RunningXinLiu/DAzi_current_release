#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$ROOT/code/bin"

echo "Building forward solver..."
make -C "$ROOT/code/src_forward" -B

echo "Building depth-kernel worker..."
make -C "$ROOT/code/depthkernel_lab" -B

echo "Building inversion solver..."
make -C "$ROOT/code/src_inv_iso_joint" -B

echo "Build finished."
echo "Inversion executable: $ROOT/code/bin/DAzimSurfTomo"
echo "Forward executable:   $ROOT/code/bin/SurfAAForward"

