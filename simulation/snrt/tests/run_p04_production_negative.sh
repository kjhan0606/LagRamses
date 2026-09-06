#!/usr/bin/env bash
set -euo pipefail

root=/gpfs/kjhan/LRD_JWST
binary="${P04_BINARY:-$root/bin/ramses_final3d}"
namelist="$root/namelist/phase0_validation_pilot.nml"
evidence_dir="$root/build/p04_fail_closed"
baseline_log="$evidence_dir/baseline_gate.log"

mkdir -p "$evidence_dir"
[[ -x "$binary" ]] || {
  echo "P04_BINARY_MISSING binary=$binary" >&2
  exit 1
}

# The current top-level production gate is the unresolved F-P1 fate policy, so
# yield-table-specific failures are intentionally unreachable in this fixture.
# Keep this test aligned with the effective binary rather than asserting stale
# downstream error codes.
set +e
env -u PHASE0_YIELD_TABLE -u PHASE0_SNIA_RUNTIME_CONTRACT \
  timeout 60s "$binary" "$namelist" >"$baseline_log" 2>&1
baseline_status=$?
set -e

if [[ "$baseline_status" -ne 3 ]]; then
  echo "P04_BASELINE_GATE_STATUS_BAD status=$baseline_status" >&2
  exit 1
fi
if ! grep -q 'source model is not implemented for production' "$baseline_log"; then
  echo "P04_BASELINE_GATE_MESSAGE_BAD" >&2
  exit 1
fi

SNRT_PRODUCTION_BINARY="$binary" \
  python3 "$root/simulation/snrt/tests/fp2_snia_production_runtime_negative.py"

echo "P04_PRODUCTION_NEGATIVE_OK baseline=$baseline_status snia_fail_closed=pass"
