#!/usr/bin/env bash
# Re-link only the changed hydro-parameter validator against the archived
# M5+CHIMES build. Keep the original object and executable for provenance.
#SBATCH --job-name=m5-dust-validator
#SBATCH --partition=a100,h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-dust-validator-build-%j.sbatch.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly BUILD=$ROOT/.m5-cosmo-dust-chimes-build-419104
readonly BIN=$BUILD/bin/ramses_m5_chimes3d
readonly BACKUP_BIN=$BUILD/bin/ramses_m5_chimes3d.pre_dust_validator
readonly BACKUP_OBJ=$BUILD/bin/read_hydro_params.o.pre_dust_validator

test -x "$BIN"
test -f "$BUILD/bin/read_hydro_params.o"
test ! -e "$BACKUP_BIN"
test ! -e "$BACKUP_OBJ"
cp -p -- "$BIN" "$BACKUP_BIN"
cp -p -- "$BUILD/bin/read_hydro_params.o" "$BACKUP_OBJ"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

cd "$BUILD/bin"
make -j1 \
  GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$(git -C "$ROOT" rev-parse HEAD)" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$BUILD/bin/ramses_m5_chimes_dust_validator" \
  SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_CUDA=1 USE_FFTW=1 \
  NENER=0 NVECTOR=32 CHIMES_DIR="$BUILD/chimes-install" \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release \
  ramses > "$ROOT/.m5-dust-validator-build-${SLURM_JOB_ID}.log" 2>&1

test -x "$BUILD/bin/ramses_m5_chimes_dust_validator3d"
if ldd "$BUILD/bin/ramses_m5_chimes_dust_validator3d" | rg 'not found'; then
  echo 'Unresolved shared-library dependency in diagnostic binary' >&2
  exit 3
fi
sha256sum "$BUILD/bin/ramses_m5_chimes_dust_validator3d" \
  "$BUILD/bin/read_hydro_params.o" > "$ROOT/.m5-dust-validator-build-${SLURM_JOB_ID}.sha256"
