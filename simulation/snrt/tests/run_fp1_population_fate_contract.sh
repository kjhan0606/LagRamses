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
"$PYTHON" "$SNRT_ROOT/tests/fp1_source_node_contract.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_source_node_projection.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_source_node_contract.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_terminal_deposition_contract.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_terminal_deposition_contract.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_physical_package_admission.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_physical_package_admission.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_low_mass_seam.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_low_mass_seam.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_high_mass_seam.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_high_mass_seam.py"
echo "FP1_POPULATION_FATE_CONTRACT_OK"
