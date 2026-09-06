#!/usr/bin/env bash
set -euo pipefail

SNRT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON="$SNRT_ROOT/.venv/bin/python"
[[ -x "$PYTHON" ]] || { echo "G4_FAIL missing project Python=$PYTHON" >&2; exit 1; }

export JAX_PLATFORMS=cpu
export PYTHONPATH="$SNRT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON" "$SNRT_ROOT/tests/p4_ingestion.py"
"$PYTHON" "$SNRT_ROOT/tests/source_sed_dust_closure.py"
"$PYTHON" "$SNRT_ROOT/tests/p5_dust_runner.py"

echo "G4_DUST_CLOSURE_PASS tests=3 mapping=explicit thermal=one_pass backend=cpu"
