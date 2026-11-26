#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT
cd -- "$TMP_DIR"

# Isolate the unchanged production routine. MPI reductions are a single-rank
# copy stub, so this checks child selection and exchange, not MPI integration.
gfortran -cpp -E -DNDIM=3 -DWITHOUTMPI \
  -I "$SOURCE_ROOT/patch/lagRamses" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_scatter.f90" > scatter_preprocessed.f90
python3 - <<'PY'
from pathlib import Path
import re

source = Path('scatter_preprocessed.f90').read_text()
matches = re.findall(r'^subroutine sidm_baryon_drag\b.*?^end subroutine sidm_baryon_drag\s*$',
                     source, re.IGNORECASE | re.MULTILINE | re.DOTALL)
if len(matches) != 1:
    raise SystemExit('could not isolate exactly one production drag routine')
Path('drag_routine.f90').write_text(matches[0] + '\n')
PY
cat > mpi_stub.c <<'C'
#include <stddef.h>
#include <string.h>
void mpi_allreduce_(const void *src, void *dst, const int *n, const int *type,
                    const int *op, const int *comm, int *error) {
  (void)op;
  (void)comm;
  size_t width = *type == 1 ? sizeof(int) : sizeof(double);
  memcpy(dst, src, (size_t)*n * width);
  *error = 0;
}
C
cc -c mpi_stub.c -o mpi_stub.o
gfortran -cpp -DNDIM=3 -DWITHOUTMPI -ffree-line-length-none \
  -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -J "$TMP_DIR" -I "$TMP_DIR" \
  "$SCRIPT_DIR/test_sidm_baryon_drag_cells.f90" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_cell_geometry.f90" \
  drag_routine.f90 mpi_stub.o -o ./test_sidm_baryon_drag_cells
./test_sidm_baryon_drag_cells
