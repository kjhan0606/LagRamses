#!/bin/bash
# Clean bounded build/regression on the user's requested compute host.
set -euo pipefail
test "$(hostname -s)" = LagEunha || { echo 'Run manually on Lageunha, not syntax.' >&2; exit 1; }
code=/home/kjhan/BACKUP/lagRamses-SIDM
evidence=${1:?Specify an existing fresh evidence root}
test -d "$evidence"
test ! -e "$evidence/source" || { echo 'Refusing to overwrite source snapshot.' >&2; exit 1; }
test "$(git -C "$code" rev-parse --show-toplevel)" = "$code"
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17
module list >"$evidence/build_modules.txt" 2>&1
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export I_MPI_FABRICS=shm PYTHONDONTWRITEBYTECODE=1
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:/home/kjhan/miniconda3/lib:${LD_LIBRARY_PATH:-}
snapshot="$evidence/source"
mkdir "$snapshot"
git -C "$code" rev-parse HEAD >"$evidence/base_revision.txt"
git -C "$code" diff --binary -- patch bin >"$evidence/build_source_changes.diff"
git -C "$code" archive HEAD | tar -xf - -C "$snapshot"
# Carry the existing working source too, without touching or committing it.
while IFS= read -r relative; do
    test -f "$code/$relative"
    cp -p "$code/$relative" "$snapshot/$relative"
done < <(git -C "$code" diff --name-only -- patch bin)
head_sha=$(git -C "$code" rev-parse HEAD)
cd "$snapshot/bin"
TIMEFORMAT='build wall_seconds=%R user_seconds=%U system_seconds=%S'
# The current Makefile does not declare every Fortran module dependency;
# a fresh parallel build can race amr_parameters.mod. Use serial ordering.
{ time nice -n 19 timeout 600s make -j1 HDF5=1 USE_FFTW=1 EXEC=ramses_referee_p0 \
      GITBRANCH=local-p0-snapshot GITHASH="${head_sha}-dirty" \
      GITREPO=local-source-snapshot; } >"$evidence/build.log" 2>&1
sha256sum "$snapshot/bin/ramses_referee_p03d" >"$evidence/BINARY_SHA256SUMS"
cd "$evidence"
sha256sum source/patch/lagRamses/dark_cooling_mod.f90 \
    source/patch/cuRamses/dark_cooling_mod.f90 \
    source/patch/cuRamses/sidm_scatter.f90 \
    source/patch/lagRamses/read_params.jaehyun.f90 \
    source/bin/Makefile >SOURCE_SHA256SUMS
# Regression drivers compile small standalone sources; this full build does
# not launch RAMSES, read ICs or checkpoints, or create a simulation dump.
nice -n 19 timeout 90s python3 "$code/tests/sidm/test_sidm_estimator.py" >estimator.log 2>&1
nice -n 19 timeout 90s python3 "$code/tests/sidm/test_angular_transport_integrals.py" >angular_integrals.log 2>&1
nice -n 19 timeout 90s python3 "$code/tests/adm/test_adm_hpm_wiring.py" >hpm_wiring.log 2>&1
echo 'Clean HDF5+FFTW build and estimator/angular/HPM regressions passed.'
