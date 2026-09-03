#!/usr/bin/env bash
set -euo pipefail

root=/gpfs/kjhan/LRD_JWST
binary="$root/bin/ramses_final3d"
namelist="$root/namelist/phase0_validation_pilot.nml"
fixture="$root/patch/lagRamses/phase0_validation_yields.dat"
evidence_dir="$root/build/p04_fail_closed"
missing_log="$evidence_dir/missing_table.log"
coverage_log="$evidence_dir/incomplete_coverage.log"

mkdir -p "$evidence_dir"

set +e
env -u PHASE0_YIELD_TABLE timeout 30s "$binary" "$namelist" \
  >"$missing_log" 2>&1
missing_status=$?
env PHASE0_YIELD_TABLE="$fixture" timeout 30s "$binary" "$namelist" \
  >"$coverage_log" 2>&1
coverage_status=$?
set -e

if [[ "$missing_status" -ne 1 ]]; then
  echo "P04_MISSING_TABLE_STATUS_BAD status=$missing_status" >&2
  exit 1
fi
if ! grep -q 'embedded fallback is disabled' "$missing_log"; then
  echo "P04_MISSING_TABLE_MESSAGE_BAD" >&2
  exit 1
fi
if [[ "$coverage_status" -ne 121 ]]; then
  echo "P04_COVERAGE_STATUS_BAD status=$coverage_status" >&2
  exit 1
fi
if ! grep -q 'does not cover enabled channel' "$coverage_log"; then
  echo "P04_COVERAGE_MESSAGE_BAD" >&2
  exit 1
fi

echo "P04_PRODUCTION_NEGATIVE_OK missing=$missing_status coverage=$coverage_status"
