#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

cd -- "$TMP_DIR"

gfortran -cpp -DNDIM=3 -ffree-line-length-none -fcheck=all \
  -J "$TMP_DIR" -I "$TMP_DIR" \
  "$SCRIPT_DIR/test_sidm_cell_classification.f90" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_cell_geometry.f90" \
  -o ./test_sidm_cell_classification
./test_sidm_cell_classification
