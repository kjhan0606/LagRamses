#!/bin/bash
#SBATCH --job-name=chimes-mgpu-a100
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:2
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/multigpu-a100-%j.log
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/12.8.1

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
build=${CHIMES_QUAL_BUILD:?submit with CHIMES_QUAL_BUILD pointing to a completed integrated qualification build}
sun=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
table=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
cd "$build"
export LD_LIBRARY_PATH="$build/library:$sun/lib64:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
export SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated
unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_RHS_VERIFY \
  SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_CVODE_DIAGNOSTICS

# Two independent rank-like processes, each pinned to one A100, exercise the
# production rank->GPU mapping while eight OMP workers share three streams.
srun --ntasks=2 --cpus-per-task=8 --gpu-bind=single:1 bash -c '
  printf "TASK=%s CUDA_VISIBLE_DEVICES=%s GPU=" "$SLURM_PROCID" "${CUDA_VISIBLE_DEVICES:-unset}"
  nvidia-smi --query-gpu=uuid,name --format=csv,noheader
  timeout 240 ./test "$1" check-integrated
' _ "$table"

unset SNRT_CHIMES_COMPUTE_BACKEND
srun --ntasks=2 --cpus-per-task=8 --gpu-bind=single:1 bash -c '
  printf "LU_TASK=%s CUDA_VISIBLE_DEVICES=%s GPU=" "$SLURM_PROCID" "${CUDA_VISIBLE_DEVICES:-unset}"
  nvidia-smi --query-gpu=uuid,name --format=csv,noheader
  timeout 240 ./test "$1" check-lu
' _ "$table"

echo MULTIGPU_INTEGRATED_QUALIFICATION_PASS
