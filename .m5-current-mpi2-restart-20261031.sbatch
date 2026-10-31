#!/usr/bin/env bash
# Current-source CPU two-rank continuous/restart equivalence run.
# The three full HDF5 dumps are retained only through evaluation.
#SBATCH --job-name=m5-cpu-mpi2-restart
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --gres=none
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-current-mpi2-restart-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly BINARY=$ROOT/.m5-current-source-effective-cpu-build-1248163/bin/ramses_m5_current_cpu3d
readonly PACKAGE=$ROOT/.m5-kroupa-source-package-1245389

printf 'M5_MPI2_RESTART_WRAPPER_SHA256='
sha256sum "$0" | awk '{print $1}'
printf 'BINARY_SHA256='
sha256sum "$BINARY" | awk '{print $1}'
cat "$(dirname "$(dirname "$BINARY")")/build-result.txt"
printf 'MPI=%s OMP_PER_RANK=%s EXPECTED_TOTAL_CPU=8\n' "${SLURM_NTASKS:?}" "${SLURM_CPUS_PER_TASK:?}"

M5_CHIMES_BINARY="$BINARY" \
M5_CHIMES_LIBRARY_MODE=reference \
M5_CHIMES_BACKEND=cpu \
M5_CHIMES_CELL_SCHEDULER=cohort \
M5_TRANSPORT_BACKEND=openmp \
M5_MATERIAL_BACKEND=openmp \
M5_SOURCE_PROFILE=effective_population \
M5_STELLAR_PACKAGE="$PACKAGE" \
M5_EFFECTIVE_NSTEP=4 \
M5_CAPTURE_FINAL_DUMP=0 \
M5_VERIFY_RESTART=1 \
M5_AMR_INTERFACE_CHECK=0 \
bash "$ROOT/simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch"
