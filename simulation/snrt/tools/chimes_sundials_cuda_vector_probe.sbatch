#!/bin/bash
#SBATCH --job-name=sundials-cuda-vector
#SBATCH --partition=a100,h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/sundials-cuda-vector-%j.log

set -euo pipefail
module purge
module load gnu13 cuda/13.0.2

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly SOURCE=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
readonly WORK=$ROOT/.cvode-gpu/sundials-cuda-vector-$SLURM_JOB_ID
readonly BUILD=$WORK/build
mkdir -p "$BUILD"

printf 'PURPOSE=build SUNDIALS 5.8 CUDA N_Vector and run upstream managed-CUDA CVODE smoke test\n'
printf 'JOB_ID=%s NODE=%s\n' "$SLURM_JOB_ID" "$(hostname)"
nvidia-smi --query-gpu=name,uuid --format=csv,noheader

# Keep the existing CPU-only SUNDIALS install untouched; this trial has its
# own build prefix and enables the CUDA vector module explicitly.
cmake -S "$SOURCE" -B "$BUILD" \
  -DCMAKE_INSTALL_PREFIX="$WORK/install" \
  -DCMAKE_CUDA_ARCHITECTURES='80;90' \
  -DSUNDIALS_PRECISION=DOUBLE -DSUNDIALS_INDEX_SIZE=64 \
  -DENABLE_CUDA=ON -DBUILD_NVECTOR_CUDA=ON -DBUILD_CVODE=ON \
  -DBUILD_STATIC_LIBS=OFF -DBUILD_SHARED_LIBS=ON \
  -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF \
  -DEXAMPLES_ENABLE_CUDA=ON -DEXAMPLES_INSTALL=OFF

cmake --build "$BUILD" --target cvAdvDiff_kry_cuda_managed -j4
rg -n 'SUNDIALS_NVECTOR_CUDA|SUNDIALS_VERSION ' "$BUILD/include/sundials/sundials_config.h"
ldd "$BUILD/examples/cvode/cuda/cvAdvDiff_kry_cuda_managed" | rg 'sundials|cuda|not found'
timeout 120 "$BUILD/examples/cvode/cuda/cvAdvDiff_kry_cuda_managed"
printf 'SUNDIALS_CUDA_NVECTOR_PROBE=PASS\n'
