#!/usr/bin/env bash
# Bounded, source-dark cosmological M5 + grain-mass + CHIMES integration.
# Reference controls are explicitly NONPRODUCTION; no raw output is scheduled.
# Only A100/H100 allocations are accepted at runtime (H200 is prohibited).
#SBATCH --job-name=m5-dust-chimes
#SBATCH --partition=a100,h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-cosmo-dust-chimes-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly ASSETS=/gpfs/kjhan/LRD_JWST
readonly SUN=$ASSETS/.snrt-performance.jDx9Iz/sundials-release
readonly CHIMES_SOURCE=$ROOT/.cvode-gpu/chimes
readonly HDF5_ROOT=/home/kjhan/local/hdf5
readonly TABLE=$ASSETS/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
readonly GROUP_DIR=$ASSETS/.dust-extension.AOz7mU/groups-v7
readonly WORK=$ROOT/.m5-cosmo-dust-chimes-$SLURM_JOB_ID
readonly BUILD=$ROOT/.m5-cosmo-dust-chimes-build-$SLURM_JOB_ID
readonly CONFIG=$ROOT/simulation/snrt/config/m5_cosmo_dust_chimes_preflight.nml
readonly BINARY=$BUILD/bin/ramses_m5_chimes3d
readonly CHIMES_STAGE=$BUILD/chimes-library-stage
readonly CHIMES_INSTALL=$BUILD/chimes-install

test ! -e "$WORK"
test ! -e "$BUILD"
mkdir "$WORK" "$BUILD" "$BUILD/bin"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader | sed -n '1p')
case "$gpu_name" in
  *A100*|*H100*) ;;
  *) echo "Refusing non-A100/H100 GPU allocation: $gpu_name" >&2; exit 2 ;;
esac
printf 'GPU=%s\nJOB_ID=%s\nSOURCE_HEAD=%s\n' \
  "$gpu_name" "$SLURM_JOB_ID" "$(git -C "$ROOT" rev-parse HEAD)" | tee "$WORK/launch.txt"

test -d "$ASSETS/.galaxy-preenriched.IJ8aeX/ics/level_007"
test -f "$ASSETS/.effective-population.gQodLz/inputs/history.nml"
test -f "$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/yields.dat"
test -f "$TABLE" -a -f "$GROUP_DIR/manifest.json"
test -f "$SUN/include/cvode/cvode.h" -a -f "$HDF5_ROOT/include/hdf5.h"
for i in {1..9}; do test -f "$GROUP_DIR/group_$(printf '%02d' "$i").hdf5"; done
cp -- "$CONFIG" "$WORK/effective.nml"
cp -- "$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml" "$WORK/"
cp -- "$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml" "$WORK/"
cp -- "$ROOT/simulation/snrt/config/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml" "$WORK/"
sha256sum "$CONFIG" "$WORK/effective.nml" "$WORK"/*.nml \
  "$TABLE" "$GROUP_DIR/manifest.json" "$GROUP_DIR"/group_*.hdf5 \
  "$ASSETS/.effective-population.gQodLz/inputs/history.nml" \
  "$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/yields.dat" > "$WORK/inputs.sha256"

for source_dir in patch hydro aton pm poisson amr utils simulation tests mhd; do
  ln -s "$ROOT/$source_dir" "$BUILD/$source_dir"
done
cp -- "$ROOT/bin/Makefile" "$BUILD/bin/Makefile"

# Build fresh RAMSES objects from this worktree; reuse no stale executable.
bash "$ROOT/simulation/snrt/tools/build_chimes_cpu_library.sh" \
  "$CHIMES_SOURCE" "$CHIMES_STAGE" "$SUN" "$HDF5_ROOT" reference
mkdir "$CHIMES_INSTALL"
ln -s "$CHIMES_SOURCE/src" "$CHIMES_INSTALL/src"
cp -- "$CHIMES_STAGE/libchimes.so" "$CHIMES_INSTALL/libchimes.so"

cd "$BUILD/bin"
make -B GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$(git -C "$ROOT" rev-parse HEAD)" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$BUILD/bin/ramses_m5_chimes" SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_CUDA=1 \
  USE_FFTW=1 NENER=0 NVECTOR=32 CHIMES_DIR="$CHIMES_INSTALL" \
  SUNDIALS_DIR="$SUN" ramses > "$WORK/build.log" 2>&1
test -x "$BINARY"
if ldd "$BINARY" | rg 'not found'; then
  echo 'Unresolved shared-library dependency in M5+CHIMES binary' >&2
  exit 3
fi
sha256sum "$BINARY" "$CHIMES_STAGE/libchimes.so" > "$WORK/build.sha256"

export LD_LIBRARY_PATH="$CHIMES_INSTALL:$SUN/lib64:$SUN/lib:$HDF5_ROOT/lib:/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm I_MPI_PIN=0
export SNRT_BACKEND=auto SNRT_DUST_BACKEND=auto SNRT_REDUCED_C=0.01
export SNRT_ALLOW_REFERENCE_CONTROL=1
export SNRT_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1
export SNRT_GROUP_CONTRACT="$WORK/snrt_group_contract_reference_control_v1.nml"
export SNRT_SECONDARY_TABLE_CONTRACT="$WORK/snrt_secondary_table_contract_v1.nml"
export SNRT_DUST_CONTRACT="$WORK/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml"
export SNRT_CHIMES_BINARY="$BINARY"
export SNRT_CHIMES_MAIN_DATA="$TABLE"
export SNRT_CHIMES_GROUP_DIR="$GROUP_DIR"
export SNRT_CHIMES_BAND_TABLE="$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5"
export SNRT_CHIMES_MOLECULAR_TABLE="$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5"
export SNRT_CHIMES_COMPUTE_BACKEND=cpu
export SNRT_PERF_DIAGNOSTIC=1
export PHASE0_YIELD_TABLE="$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/yields.dat"
export PHASE0_SNIA_RUNTIME_CONTRACT="$ROOT/simulation/snrt/config/fp2_snia_effective_population_runtime_v1.nml"
unset SNRT_RT_ENABLE SNRT_RT_LEVEL SNRT_AGN_MODEL SNRT_STELLAR_SED \
  SNRT_DRIVER_TEST_SEED_SOURCE SNRT_RT_TX_TEST_FAIL_STAGE SNRT_RT_TX_TEST_FAIL_LEAF \
  SNRT_PERF_DIAGNOSTIC_ALL SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND \
  SNRT_CHIMES_SPECTRAL_MODEL SNRT_CHIMES_CELL_SCHEDULER SNRT_CHIMES_GPU_BROKERS \
  SNRT_CHIMES_GPU_WORKER SNRT_CHIMES_GPU_BATCH_CELLS SNRT_CHIMES_GPU_RESULT SNRT_CHIMES_GPU_COOLING

cd "$WORK"
printf 'PROFILE=cosmological_M5_dust_mass_CHIMES_dark_preflight\n' | tee -a launch.txt
printf 'LIMIT=2_coarse_steps; MPI=1; OMP=8; RT=M5_order5; CHIMES_backend=cpu\n' | tee -a launch.txt
printf 'OUTPUT_POLICY=aout=1.1,tout=1e100,foutput=fbackup=1000000; expected_raw_dumps=0\n' | tee -a launch.txt
sha256sum effective.nml "$BINARY" >> launch.txt

srun --mpi=pmi2 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
  --cpu-bind=none "$BINARY" effective.nml > run.log 2>&1

grep -q 'Run completed' run.log
grep -Eq 'Main step= *2([[:space:]]|$)' run.log
grep -Fq 'Radiation transport selector: snrt_mn' run.log
grep -Fq 'Dust mass/thermal evolution (dust_mass_enabled)=T' run.log
grep -Fq 'CHIMES built in double precision.' run.log
grep -Fq 'SNRT M_N live AMR transport' run.log
grep -Eq 'SNRT_MN_COMMIT level/order/substeps=.*[[:space:]]5[[:space:]]' run.log
grep -Eq 'SNRT PERF (DARK|PHOTO|TILE)' run.log
if rg -qi 'ERROR:|Aborting|forrtl:|failed to converge|nonfinite|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'FAIL: unexpected error/nonfinite/nonconvergence marker; preserving all outputs' >&2
  exit 4
fi
shopt -s nullglob
raw_outputs=(output_[0-9][0-9][0-9][0-9][0-9])
if ((${#raw_outputs[@]})); then
  printf 'FAIL: unexpected raw output retained for evaluation: %s\n' "${raw_outputs[*]}" >&2
  exit 5
fi
printf 'M5_COSMO_DUST_CHIMES_PREFLIGHT_PASS source_population=none reference_controls=NONPRODUCTION raw_outputs=0\n' \
  | tee "$WORK/result.txt"
printf 'cleanup_manifest=No raw simulation outputs were produced; logs, effective input, source/build hashes and compact result retained.\n' \
  > "$WORK/cleanup_manifest.txt"
