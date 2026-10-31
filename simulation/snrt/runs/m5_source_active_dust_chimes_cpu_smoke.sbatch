#!/usr/bin/env bash
# CPU-only companion to the bounded M5 + stellar/AGN + dust + CHIMES smoke.
# Submit this script on Grammar's normal queue; it builds from the current
# LagRamses source and retains no RAMSES dump files.
#SBATCH --job-name=m5-cpu-src-smoke
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=none
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-source-dust-chimes-cpu-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly ASSETS=/gpfs/kjhan/LRD_JWST
readonly SUN=$ASSETS/.snrt-performance.jDx9Iz/sundials-release
readonly CHIMES_SOURCE=$ROOT/.cvode-gpu/chimes
readonly HDF5_ROOT=/home/kjhan/local/hdf5
readonly PACKAGE=$ROOT/.m5-kroupa-source-package-1245389
readonly BUILD=$ROOT/.m5-source-dust-chimes-cpu-build-$SLURM_JOB_ID
readonly CHIMES_STAGE=$BUILD/chimes-library-stage
readonly CHIMES_INSTALL=$BUILD/chimes-install
readonly EXEC_PREFIX=$BUILD/bin/ramses_m5_cpu
readonly BINARY=${EXEC_PREFIX}3d

[[ "$SLURM_JOB_PARTITION" == normal && -z "${SLURM_JOB_GPUS:-}" ]]
test ! -e "$BUILD"
test -d "$CHIMES_SOURCE/src"
test -f "$SUN/include/cvode/cvode.h"
test -f "$HDF5_ROOT/include/hdf5.h"
test -f "$PACKAGE/effective-ssp-binding.json"
git -C "$ROOT" diff --quiet
git -C "$ROOT" diff --cached --quiet
readonly SOURCE_HEAD=$(git -C "$ROOT" rev-parse HEAD)
mkdir -p "$BUILD/bin"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17

for source_dir in patch hydro aton pm poisson amr utils simulation tests mhd; do
  ln -s "$ROOT/$source_dir" "$BUILD/$source_dir"
done
cp -- "$ROOT/bin/Makefile" "$BUILD/bin/Makefile"
bash "$ROOT/simulation/snrt/tools/build_chimes_cpu_library.sh" \
  "$CHIMES_SOURCE" "$CHIMES_STAGE" "$SUN" "$HDF5_ROOT" reference \
  > "$BUILD/chimes-build.log" 2>&1
mkdir "$CHIMES_INSTALL"
ln -s "$CHIMES_SOURCE/src" "$CHIMES_INSTALL/src"
cp -- "$CHIMES_STAGE/libchimes.so" "$CHIMES_INSTALL/libchimes.so"

cd "$BUILD/bin"
make -B GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$SOURCE_HEAD" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$EXEC_PREFIX" SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 \
  USE_CUDA=0 USE_FFTW=1 CHIMES_CUDA=0 NENER=0 NVECTOR=32 \
  CHIMES_DIR="$CHIMES_INSTALL" SUNDIALS_DIR="$SUN" ramses \
  > "$BUILD/ramses-build.log" 2>&1
test -x "$BINARY"
grep -aFq 'SNRT_MN_AGN_SOURCE accepted_packets=' "$BINARY"
grep -aFq 'SNRT_MN_STELLAR_SOURCE accepted_stars=' "$BINARY"
export LD_LIBRARY_PATH="$CHIMES_INSTALL:$SUN/lib64:$SUN/lib:$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
if ldd "$BINARY" | grep -E 'lib(cudart|cuda)\.so'; then
  echo 'CPU-only binary unexpectedly links CUDA runtime libraries' >&2
  exit 3
fi

binary_sha=$(sha256sum "$BINARY" | awk '{print $1}')
library_sha=$(sha256sum "$CHIMES_INSTALL/libchimes.so" | awk '{print $1}')
cat > "$BUILD/build-result.txt" <<EOF_BUILD
CHIMES_LIBRARY_MODE=reference
CHIMES_CUDA=0
SOURCE_HEAD=$SOURCE_HEAD
BINARY_SHA256=$binary_sha
CHIMES_LIBRARY_SHA256=$library_sha
EOF_BUILD
sha256sum "$BINARY" "$CHIMES_INSTALL/libchimes.so" "$BUILD/build-result.txt" > "$BUILD/build.sha256"
printf 'CPU_COMBINED_SOURCE_BUILD_PASS\nSOURCE_HEAD=%s\nBINARY=%s\n' \
  "$SOURCE_HEAD" "$BINARY" | tee "$BUILD/launch.txt"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm I_MPI_PIN=0
M5_CHIMES_BINARY="$BINARY" \
M5_CHIMES_LIBRARY_MODE=reference \
M5_CHIMES_BACKEND=cpu \
M5_CHIMES_CELL_SCHEDULER=cohort \
M5_TRANSPORT_BACKEND=auto \
M5_MATERIAL_BACKEND=auto \
M5_SOURCE_PROFILE=effective_population \
M5_EFFECTIVE_AGN_IC_MODE=thermal_seed_reference \
M5_STELLAR_PACKAGE="$PACKAGE" \
M5_EFFECTIVE_NSTEP=4 \
M5_CAPTURE_FINAL_DUMP=0 \
M5_VERIFY_RESTART=0 \
M5_AMR_INTERFACE_CHECK=0 \
M5_RUN_ID="grammar-cpu-$SLURM_JOB_ID" \
bash "$ROOT/simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch"

printf 'CPU_COMBINED_SOURCE_SMOKE_PASS job=%s binary=%s\n' "$SLURM_JOB_ID" "$BINARY" \
  | tee -a "$BUILD/launch.txt"
