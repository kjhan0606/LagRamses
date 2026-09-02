#!/usr/bin/env bash
set -euo pipefail

# P0.1 production-linked build harness.  The canonical source is selected by
# bin/Makefile (PATCH=../patch/lagRamses).  The native G1 runner remains a
# differential oracle; it is not allowed to stand in for this build.
ROOT=/gpfs/kjhan/LRD_JWST
BUILD_DIR="$ROOT/build/p0_production_linked"
DATA_DIR="$ROOT/simulation/snrt/data"
STATIC_EVIDENCE="$BUILD_DIR/static_source_parity.json"
BUILD_LOG="$BUILD_DIR/make.log"
SMOKE_LOG="$BUILD_DIR/smoke.log"
BINARY="$ROOT/bin/ramses_final3d"
SMOKE_EXPECTED_EXIT_CODE=$(python3 - "$ROOT/simulation/snrt/config/stellar_source_identity_v1.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    print(json.load(handle)["production_linked_harness"]["smoke_expected_exit_code"])
PY
)

mkdir -p "$BUILD_DIR" "$DATA_DIR"

# This preflight is read-only and is intentionally allowed to report BLOCK:
# the final evidence file does not exist until the forced production build
# below has completed successfully.
python3 "$ROOT/simulation/snrt/tools/validate_stellar_source_parity.py" \
  --json > "$STATIC_EVIDENCE"

if [[ "${P0_BUILD:-0}" != 1 ]]; then
  echo "P0_PRODUCTION_LINKED_BUILD_NOT_RUN set P0_BUILD=1 to authorize the forced bin/Makefile build" >&2
  exit 2
fi

# -B is required: an old object in bin/ is not evidence of this source tree.
# Keep the command text stable because the parity validator audits the target.
printf '%s\n' 'P0_BUILD_COMMAND make -C "$ROOT/bin" -B ramses' > "$BUILD_LOG"
make -C "$ROOT/bin" -B ramses 2>&1 | tee -a "$BUILD_LOG"
BINARY_SHA256=$(sha256sum "$BINARY" | awk '{print $1}')
printf 'P0_BINARY_SHA256=%s\n' "$BINARY_SHA256" >> "$BUILD_LOG"

printf '%s\n' 'P0_SMOKE_COMMAND "$ROOT/bin/ramses_final3d"' > "$SMOKE_LOG"
printf 'P0_BINARY_SHA256=%s\n' "$BINARY_SHA256" >> "$SMOKE_LOG"
set +e
"$BINARY" >> "$SMOKE_LOG" 2>&1
SMOKE_STATUS=$?
set -e
printf 'P0_SMOKE_EXIT_CODE=%s\n' "$SMOKE_STATUS" >> "$SMOKE_LOG"
if [[ "$SMOKE_STATUS" -ne "$SMOKE_EXPECTED_EXIT_CODE" ]]; then
  echo "P0_PRODUCTION_LINKED_SMOKE_FAILED status=$SMOKE_STATUS expected=$SMOKE_EXPECTED_EXIT_CODE" >&2
  exit 1
fi

python3 "$ROOT/simulation/snrt/tools/record_p0_production_linked_build.py" \
  --binary "$BINARY" \
  --build-log "$BUILD_LOG" \
  --output "$ROOT/simulation/snrt/data/p0_production_linked_build_evidence.json"

python3 "$ROOT/simulation/snrt/tools/validate_stellar_source_parity.py" --require-pass
echo "P0_PRODUCTION_LINKED_CONTRACT_OK"
