#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT
cd -- "$TMP_DIR"

# Compile the production wrapper and RNG against small fixture modules. The
# reduction is a one-rank stub; MPI and full RAMSES integration remain separate.
gfortran -cpp -E -DNDIM=3 -DWITHOUTMPI -I "$SOURCE_ROOT/patch/lagRamses" \
  "$SOURCE_ROOT/patch/cuRamses/sidm_scatter.f90" > scatter_preprocessed.f90
python3 - "$SOURCE_ROOT/patch/lagRamses/kjhan.f90" <<'PY'
from pathlib import Path
import re
import sys
source = Path('scatter_preprocessed.f90').read_text()
matches = re.findall(r'^subroutine sidm_scatter\b.*?^end subroutine sidm_scatter\s*$',
                     source, re.IGNORECASE | re.MULTILINE | re.DOTALL)
if len(matches) != 1:
    raise SystemExit('could not isolate exactly one SIDM wrapper')
Path('sidm_wrapper.f90').write_text(matches[0] + '\n')
source = Path(sys.argv[1]).read_text()
matches = re.findall(r'^subroutine pthreadLinkedList\b.*?^end subroutine pthreadLinkedList\s*$',
                     source, re.IGNORECASE | re.MULTILINE | re.DOTALL)
if len(matches) != 1:
    raise SystemExit('could not isolate exactly one production linked-list partitioner')
Path('linked_list.f90').write_text(matches[0] + '\n')
PY
cat > mpi_stub.c <<'C'
void mpi_allreduce_(void *a, void *b, int *n, int *t, int *op, int *c, int *e) {
  (void)a; (void)b; (void)n; (void)t; (void)op; (void)c; *e=0;
}
C
cc -c mpi_stub.c -o mpi_stub.o
gfortran -cpp -DNDIM=3 -DWITHOUTMPI -fopenmp -ffree-line-length-none \
  -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -c "$SCRIPT_DIR/test_sidm_rng_driver_modules.f90" -o modules.o
gfortran -cpp -DNDIM=3 -DWITHOUTMPI -fopenmp -ffree-line-length-none \
  -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -c "$SOURCE_ROOT/patch/cuRamses/random.f90" -o random.o
gfortran -cpp -DNDIM=3 -DWITHOUTMPI -fopenmp -ffree-line-length-none \
  -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow \
  "$SCRIPT_DIR/test_sidm_rng_driver.f90" sidm_wrapper.f90 linked_list.f90 modules.o random.o mpi_stub.o \
  -o test_sidm_rng_driver
OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 ./test_sidm_rng_driver > one.out
OMP_DYNAMIC=FALSE OMP_NUM_THREADS=2 ./test_sidm_rng_driver > two.out
cmp one.out two.out
cat one.out
echo 'One- and two-thread configurations produce identical SIDM wrapper draws: PASS'
