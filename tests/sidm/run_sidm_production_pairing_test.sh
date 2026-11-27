#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT
cd -- "$TMP_DIR"

# Exercise the production scatter body, with deterministic fixture modules
# and random draws. This is not an MPI or full RAMSES integration test.
gfortran -cpp -E -DNDIM=3 -DWITHOUTMPI \
  -I "$SOURCE_ROOT/patch/lagRamses" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_scatter.f90" > scatter_preprocessed.f90
python3 - <<'PY'
from pathlib import Path
import re

source = Path('scatter_preprocessed.f90').read_text()
matches = re.findall(r'^subroutine sub_sidm_scatter\b.*?^end subroutine sub_sidm_scatter\s*$',
                     source, re.IGNORECASE | re.MULTILINE | re.DOTALL)
if len(matches) != 1:
    raise SystemExit('could not isolate exactly one production scatter routine')
Path('scatter_routine.f90').write_text(matches[0] + '\n')
PY
gfortran -cpp -DNDIM=3 -DWITHOUTMPI -ffree-line-length-none \
  -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -J "$TMP_DIR" -I "$TMP_DIR" \
  "$SCRIPT_DIR/test_sidm_production_pairing.f90" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_cell_geometry.f90" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_estimator.f90" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_angular.f90" \
  scatter_routine.f90 -o ./test_sidm_production_pairing
./test_sidm_production_pairing
for case_name in zero one; do
  ./test_sidm_production_pairing "$case_name"
done
for case_name in above_one negative nonfinite; do
  if ./test_sidm_production_pairing "$case_name" > "$case_name.log" 2>&1; then
    echo "Invalid probability case was accepted: $case_name" >&2
    exit 1
  fi
  if ! grep -Fq 'ERROR: invalid SIDM probability' "$case_name.log"; then
    echo "Invalid probability case failed for the wrong reason: $case_name" >&2
    cat "$case_name.log" >&2
    exit 1
  fi
  echo "Invalid SIDM probability rejected: $case_name"
done
