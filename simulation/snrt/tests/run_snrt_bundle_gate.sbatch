#!/usr/bin/env bash
#SBATCH --job-name=snrt-bundle
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=64G
#SBATCH --time=01:30:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.snrt-bundle-gate-%j.log

set -euo pipefail

ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 \
  intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

export CUDA_ARCH=80
export OMP_NUM_THREADS=1
export SNRT_BUNDLE_GATE_SUMMARY="$ROOT/.snrt-bundle-gate-${SLURM_JOB_ID}.summary.txt"
cd "$ROOT"
"$ROOT/simulation/snrt/tests/run_snrt_bundle_gate.sh"
