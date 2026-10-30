#!/usr/bin/env bash
# CPU-only production-linked source parity build on grammar's normal queue.
# The harness builds in a temporary detached worktree and preserves bin/.
#SBATCH --job-name=p0-stellar-source
#SBATCH --partition=normal
#SBATCH --gres=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.p0-production-linked-%j.log

set -euo pipefail

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17

export OMP_NUM_THREADS=1
export P0_BUILD=1
exec bash /gpfs/kjhan/chimes-cvode-gpu-20261001/simulation/snrt/tests/run_p0_production_linked_contract.sh
