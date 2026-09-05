#!/usr/bin/env bash
set -euo pipefail

# Native CUDA smoke for the production multigroup cap.  It exercises the C
# ABI and the actual GPU kernel; it does not launch RAMSES or import Python.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-cuda-multigroup.XXXXXX)"
trap 'rm -rf "$build_dir"' EXIT

fc="${FC:-mpiifx}"
fc_name="${fc##*/}"
if [[ "$fc_name" == "ifx" || "$fc_name" == "mpiifx" ]]; then
  flags=(-fpp -DNDIM=3 -DNPRE=8 -DNVAR=18 -module "$build_dir" -I"$build_dir")
else
  flags=(-cpp -ffree-line-length-none -DNDIM=3 -DNPRE=8 -DNVAR=18 -J"$build_dir" -I"$build_dir")
fi

cuda_root="${CUDA_ROOT:-/opt/ohpc/pub/cuda/13.0.2}"
nvcc="${NVCC:-$cuda_root/bin/nvcc}"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_cuda_multigroup_interface.f90" \
  -o "$build_dir/snrt_cuda_multigroup_interface.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_cuda_multigroup_smoke.f90" \
  -o "$build_dir/snrt_cuda_multigroup_smoke.o"
"$nvcc" -O2 -gencode arch=compute_86,code=sm_86 -c \
  "$repo_root/patch/lagRamses/snrt_cuda_kernels.cu" \
  -o "$build_dir/snrt_cuda_kernels.o"
"$fc" "$build_dir/snrt_cuda_multigroup_interface.o" \
  "$build_dir/snrt_cuda_multigroup_smoke.o" "$build_dir/snrt_cuda_kernels.o" \
  -L"$cuda_root/lib64" -lcudart -lcublas -lcublasLt -lstdc++ \
  -Wl,-rpath,"$cuda_root/lib64" -o "$build_dir/snrt_cuda_multigroup_smoke"

CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}" \
  "$build_dir/snrt_cuda_multigroup_smoke"
echo 'SNRT_NATIVE_CUDA_MULTIGROUP_ALL_OK'
