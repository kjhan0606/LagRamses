#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/kjhan/LRD_JWST
SNRT_ROOT="$ROOT/simulation/snrt"
PYTHON="$SNRT_ROOT/.venv/bin/python"
BUILD_DIR="$ROOT/build/fp1_population_fate"
DATA_DIR="$SNRT_ROOT/data"
HIGH_MASS_REVIEW="$DATA_DIR/fp1_high_mass_seam_review.json"

mkdir -p "$BUILD_DIR" "$DATA_DIR"

# Regenerate the high-mass evidence before either admission consumer can read
# it.  The byte-identical result is intentional; a changed result must fail
# the code-owned admission lock rather than make a stale report look current.
"$PYTHON" "$SNRT_ROOT/tests/fp1_high_mass_seam.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_high_mass_seam.py" \
  --json-out "$HIGH_MASS_REVIEW"
high_mass_review_sha256="$(sha256sum "$HIGH_MASS_REVIEW" | awk '{print $1}')"

"$PYTHON" "$SNRT_ROOT/tests/fp1_population_fate.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_population_fate.py" \
  --json-out "$DATA_DIR/fp1_population_fate_audit.json"
"$PYTHON" "$SNRT_ROOT/tests/fp1_source_identity_rights.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_lc18_failed_wind_crosscheck.py" \
  --json-out "$DATA_DIR/fp1_lc18_failed_wind_crosscheck.json"
"$PYTHON" "$SNRT_ROOT/tests/fp1_source_node_contract.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_source_node_projection.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_source_node_contract.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_terminal_deposition_contract.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_terminal_deposition_contract.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_physical_package_admission.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_physical_package_admission.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_fate_admission.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_fate_admission.py" \
  --json-out "$DATA_DIR/fp1_fate_admission_audit.json"
"$PYTHON" "$SNRT_ROOT/tests/fp1_publication_rights.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_low_mass_seam.py"
"$PYTHON" "$SNRT_ROOT/tools/audit_fp1_low_mass_seam.py"
"$PYTHON" "$SNRT_ROOT/tests/fp1_high_mass_freshness.py" \
  --high-mass-review "$HIGH_MASS_REVIEW" \
  --physical-package "$DATA_DIR/fp1_physical_package_admission_audit.json" \
  --fate-admission "$DATA_DIR/fp1_fate_admission_audit.json" \
  --expected-sha256 "$high_mass_review_sha256"
echo "FP1_POPULATION_FATE_CONTRACT_OK"
