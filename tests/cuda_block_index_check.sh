#!/bin/bash
# Reusable static/compile gate for CUDA block grid-major cell addressing.
# Set RUN_DEVICE=1 on an allocated CUDA node to execute the device unit.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
NVCC=${NVCC:-/opt/ohpc/pub/cuda/13.0.2/bin/nvcc}
MPIIFX=${MPIIFX:-mpiifx}
TMP=$(mktemp -d /tmp/rz-cuda-index.XXXXXX)
trap 'rm -rf "$TMP"' EXIT

cd "$ROOT"
python3 tests/cuda_block_index_census.py

g++ -std=c++11 -O2 tests/cuda_block_index_unit.cpp -o "$TMP/host_unit"
"$TMP/host_unit"

CUDA_FLAGS=(-O2 -gencode arch=compute_80,code=sm_80
  -DNDIM=3 -DNVAR=11 -DNVECTOR=32 -DN_STREAMS=16 -Ipatch/cuRamses)
for src in poisson_cuda_kernels scalar_cuda_kernels particle_cuda_kernels; do
  "$NVCC" "${CUDA_FLAGS[@]}" -c "patch/cuRamses/$src.cu" -o "$TMP/$src.o"
done
"$NVCC" -O2 -gencode arch=compute_80,code=sm_80 \
  tests/cuda_block_index_device.cu -o "$TMP/device_unit"

"$MPIIFX" -fpp -E -qopenmp -DHYDRO_CUDA -DNDIM=3 -DNVECTOR=32 \
  -Ipatch/lagRamses patch/lagRamses/rho_fine.kjhan.f90 \
  > "$TMP/rho.cuda.f90"
"$MPIIFX" -fpp -E -qopenmp -DNDIM=3 -DNVECTOR=32 \
  -Ipatch/lagRamses patch/lagRamses/rho_fine.kjhan.f90 \
  > "$TMP/rho.cpu.f90"
for symbol in cuda_pm_rho_begin_c pm_gpu_append pm_rho_merge; do
  cuda_count=$(awk -v s="$symbol" 'index($0,s){n++} END{print n+0}' \
    "$TMP/rho.cuda.f90")
  cpu_count=$(awk -v s="$symbol" 'index($0,s){n++} END{print n+0}' \
    "$TMP/rho.cpu.f90")
  if [ "$cuda_count" -le 0 ] || [ "$cpu_count" -ne 0 ]; then
    echo "CUDA BLOCK INDEX CHECK: FAIL: preprocessor $symbol" >&2
    exit 1
  fi
done

if [ "${RUN_DEVICE:-0}" = 1 ]; then
  "$TMP/device_unit"
fi

echo "CUDA BLOCK INDEX CHECK: PASS"
