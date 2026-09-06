#!/usr/bin/env bash
set -euo pipefail

SNRT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON="$SNRT_ROOT/.venv/bin/python"
[[ -x "$PYTHON" ]] || { echo "G3_FAIL missing project Python=$PYTHON" >&2; exit 1; }

export JAX_PLATFORMS=cpu
export PYTHONPATH="$SNRT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON" "$SNRT_ROOT/tests/merge_photon_ledgers.py"
"$PYTHON" "$SNRT_ROOT/tests/merge_photon_source_ledgers.py"
"$PYTHON" "$SNRT_ROOT/tests/source_sed_dust_closure.py"

echo "G3_SOURCE_LEDGER_CLOSURE_PASS tests=3 backend=cpu"
