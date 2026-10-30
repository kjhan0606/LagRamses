#!/usr/bin/env bash
set -euo pipefail

# P0.1 production-linked build harness.  The canonical source is selected by
# bin/Makefile (PATCH=../patch/lagRamses).  Build in a detached worktree so
# this gate never overwrites the operator's active bin/ objects or executable.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/../../.." && pwd)"
BUILD_DIR="$ROOT/build/p0_production_linked"
DATA_DIR="$ROOT/simulation/snrt/data"
STATIC_EVIDENCE="$BUILD_DIR/static_source_parity.json"
BUILD_LOG="$BUILD_DIR/make.log"
SMOKE_LOG="$BUILD_DIR/smoke.log"
BINARY="$ROOT/build/p0_production_linked/ramses_final3d"
BUILD_COMMAND='make -C "$SOURCE_WORKTREE/bin" -B ramses'
SMOKE_COMMAND='mpirun -np 1 "$ROOT/build/p0_production_linked/ramses_final3d"'
SMOKE_EXPECTED_EXIT_CODE=$(python3 - "$ROOT/simulation/snrt/config/stellar_source_identity_v1.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    print(json.load(handle)["production_linked_harness"]["smoke_expected_exit_code"])
PY
)

if [[ "${P0_BUILD:-0}" != 1 ]]; then
  echo "P0_PRODUCTION_LINKED_BUILD_NOT_RUN set P0_BUILD=1 to authorize the forced isolated production build" >&2
  exit 2
fi

# Never reuse a prior build directory: its contents are part of the recorded
# evidence and must all come from this one forced build.
if [[ -e "$BUILD_DIR" ]]; then
  echo "P0_BUILD_OUTPUT_EXISTS refusing to overwrite $BUILD_DIR" >&2
  exit 2
fi
if [[ -e "$BINARY" ]]; then
  echo "P0_BUILD_BINARY_EXISTS refusing to overwrite $BINARY" >&2
  exit 2
fi
if [[ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no)" ]]; then
  echo "P0_SOURCE_TREE_HAS_TRACKED_CHANGES refusing to build an uncommitted source tree" >&2
  exit 2
fi

SOURCE_COMMIT=$(git -C "$ROOT" rev-parse HEAD)
WORKTREE_PARENT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
SOURCE_WORKTREE="$WORKTREE_PARENT/lagramses-p0-production-linked-${SLURM_JOB_ID:-$$}"
if [[ ! -d "$WORKTREE_PARENT" || -e "$SOURCE_WORKTREE" ]]; then
  echo "P0_WORKTREE_PATH_UNAVAILABLE $SOURCE_WORKTREE" >&2
  exit 2
fi

mkdir -p "$BUILD_DIR" "$DATA_DIR"

# This preflight is read-only and may report BLOCK until fresh linked-build
# evidence is recorded below.
python3 "$ROOT/simulation/snrt/tools/validate_stellar_source_parity.py" \
  --json > "$STATIC_EVIDENCE"

WORKTREE_ADDED=0
cleanup_worktree() {
  if [[ "$WORKTREE_ADDED" == 1 ]]; then
    git -C "$ROOT" worktree remove --force "$SOURCE_WORKTREE" || \
      echo "P0_WORKTREE_CLEANUP_FAILED $SOURCE_WORKTREE" >&2
  fi
}
trap cleanup_worktree EXIT

git -C "$ROOT" worktree add --detach "$SOURCE_WORKTREE" "$SOURCE_COMMIT"
WORKTREE_ADDED=1
if [[ "$(git -C "$SOURCE_WORKTREE" rev-parse HEAD)" != "$SOURCE_COMMIT" ]] || \
   [[ -n "$(git -C "$SOURCE_WORKTREE" status --porcelain --untracked-files=normal)" ]]; then
  echo "P0_SOURCE_WORKTREE_NOT_CLEAN $SOURCE_WORKTREE" >&2
  exit 2
fi

# -B is required: stale objects in either checkout are not evidence. Keep
# command markers stable because the parity validator audits the exact target.
printf 'P0_SOURCE_COMMIT=%s\n' "$SOURCE_COMMIT" > "$BUILD_LOG"
printf 'P0_BUILD_COMMAND %s\n' "$BUILD_COMMAND" >> "$BUILD_LOG"
make -C "$SOURCE_WORKTREE/bin" -B ramses 2>&1 | tee -a "$BUILD_LOG"
install -m 0755 "$SOURCE_WORKTREE/bin/ramses_final3d" "$BINARY"
BINARY_SHA256=$(sha256sum "$BINARY" | awk '{print $1}')
printf 'P0_BINARY_SHA256=%s\n' "$BINARY_SHA256" >> "$BUILD_LOG"

printf 'P0_SMOKE_COMMAND %s\n' "$SMOKE_COMMAND" > "$SMOKE_LOG"
printf 'P0_BINARY_SHA256=%s\n' "$BINARY_SHA256" >> "$SMOKE_LOG"
set +e
mpirun -np 1 "$BINARY" >> "$SMOKE_LOG" 2>&1
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
