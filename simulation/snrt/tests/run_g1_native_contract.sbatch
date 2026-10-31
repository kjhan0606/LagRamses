#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=snrt-g1-jax
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/simulation/snrt/data/g1-native-contract-%j.log

set -euo pipefail
module purge
module load python/3.13
module load intel/compiler-rt intel/tbb intel/umf
module load intel/compiler/2025.3.0 intel/mpi/2021.17

ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
cd "$ROOT"
P0_DIAGNOSTIC=1 bash "$ROOT/simulation/snrt/tests/run_g1_native_contract.sh"
