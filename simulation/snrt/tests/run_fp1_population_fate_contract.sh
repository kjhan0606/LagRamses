#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/kjhan/LRD_JWST
SNRT_ROOT="$ROOT/simulation/snrt"
PYTHON="$SNRT_ROOT/.venv/bin/python"
BUILD_DIR="$ROOT/build/fp1_population_fate"
DATA_DIR="$SNRT_ROOT/data"

mkdir -p "$BUILD_DIR" "$DATA_DIR"
"$PYTHON" "$SNRT_ROOT/tests/fp1_population_fate.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_population_fate.py" \
  --json-out "$DATA_DIR/fp1_population_fate_audit.json"
"$PYTHON" "$SNRT_ROOT/tests/fp1_fate_admission.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_fate_admission.py" \
  --json-out "$DATA_DIR/fp1_fate_admission_audit.json"
