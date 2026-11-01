#!/usr/bin/env bash
#SBATCH --job-name=smbh-cap-mpi2
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --partition=a100_pcie

set -euo pipefail

: "${RUN_DIR:?set RUN_DIR to a fresh prepared output directory}"
: "${BINARY:?set BINARY to the source-bound RAMSES executable}"
: "${SOURCE_TREE:?set SOURCE_TREE to the clean source worktree}"
: "${MANIFEST_SHA256:?set MANIFEST_SHA256 to the prepared manifest hash}"
: "${POSTCHECK_VALIDATOR:?set POSTCHECK_VALIDATOR to the moving-fixture validator}"

module purge
module load intel/tbb/latest intel/compiler-rt/2025.3.0 intel/umf/latest \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1

runner="$SOURCE_TREE/tests/sink/smbh-capture-restart/run_smoke.sh"
bash "$runner" --execute --run-dir "$RUN_DIR" --binary "$BINARY" \
  --manifest-sha256 "$MANIFEST_SHA256" --mpi-ranks 2
python3 "$POSTCHECK_VALIDATOR" postcheck "$RUN_DIR" --expected-mpi-ranks 2
