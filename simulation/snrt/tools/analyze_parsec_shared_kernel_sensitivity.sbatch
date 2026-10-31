#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=parsec-kernel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.parsec-kernel-sensitivity-20261031/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/chimes-cvode-gpu-20261001/.parsec-kernel-sensitivity-20261031/slurm-%j.err

set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

python3 simulation/snrt/tools/analyze_parsec_shared_kernel_sensitivity.py \
  --source-dir /gpfs/kjhan/LRD_JWST/.parsec-sources.9ZIUlJ \
  --package /gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-kroupa-source-package-1245389/input-v5 \
  --group-edges simulation/snrt/config/p0_photon_group_edges_ev.txt \
  --pilot-contract simulation/snrt/config/snrt_group_contract_reference_control_v1.nml \
  --output /gpfs/kjhan/chimes-cvode-gpu-20261001/.parsec-kernel-sensitivity-20261031/results.json
