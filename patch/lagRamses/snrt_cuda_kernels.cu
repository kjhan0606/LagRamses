// Tensor-Core angular-block contraction for the lagRamses S_N RT backend.
// Inputs are logically row-major.  cuBLAS sees their transposes as column-major
// matrices and evaluates B^T * A^T, whose memory layout is the desired A * B.
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include "snrt_species_dust_cell.h"

namespace {

int cuda_status(cudaError_t status) {
  return status == cudaSuccess ? 0 : static_cast<int>(status);
}

int cublas_status(cublasStatus_t status) {
  return status == CUBLAS_STATUS_SUCCESS ? 0 : 1000 + static_cast<int>(status);
}

}  // namespace

extern "C" int snrt_cuda_available_c() {
  int device_count = 0;
  return cudaGetDeviceCount(&device_count) == cudaSuccess ? device_count : 0;
}

// Called at a synchronization boundary. Match the existing cuRamses pool's
// local-rank mapping; never reset the device or discard its hydro buffers.
extern "C" int snrt_cuda_prepare_c(int local_rank, char *uuid) {
  int count=0;
  if(cudaGetDeviceCount(&count)!=cudaSuccess || count<=0) return 1;
  const int device=local_rank%count;
  if(cudaSetDevice(device)!=cudaSuccess) return 1;
  cudaDeviceProp prop;
  if(cudaGetDeviceProperties(&prop,device)!=cudaSuccess) return 1;
  for(int i=0;i<16;++i) std::sprintf(uuid+2*i,"%02x",static_cast<unsigned char>(prop.uuid.bytes[i]));
  return 0;
}

extern "C" long long snrt_cuda_free_bytes_c() {
  size_t free_bytes=0,total=0;
  if(cudaMemGetInfo(&free_bytes,&total)!=cudaSuccess) return 0;
  return static_cast<long long>(free_bytes);
}

extern "C" int snrt_cuda_angular_reduce_tf32_c(
    const float* directional_host,
    const float* projection_host,
    float* binned_host,
    int nrow,
    int ndirection,
    int nbin) {
  if (directional_host == nullptr || projection_host == nullptr || binned_host == nullptr ||
      nrow < 1 || ndirection < 1 || nbin < 1) {
    return -1;
  }

  const size_t directional_bytes = static_cast<size_t>(nrow) * ndirection * sizeof(float);
  const size_t projection_bytes = static_cast<size_t>(ndirection) * nbin * sizeof(float);
  const size_t binned_bytes = static_cast<size_t>(nrow) * nbin * sizeof(float);
  float* directional_device = nullptr;
  float* projection_device = nullptr;
  float* binned_device = nullptr;
  cublasHandle_t handle = nullptr;
  int ierr = 0;

  if ((ierr = cuda_status(cudaMalloc(&directional_device, directional_bytes))) != 0) goto cleanup;
  if ((ierr = cuda_status(cudaMalloc(&projection_device, projection_bytes))) != 0) goto cleanup;
  if ((ierr = cuda_status(cudaMalloc(&binned_device, binned_bytes))) != 0) goto cleanup;
  if ((ierr = cuda_status(cudaMemcpy(directional_device, directional_host, directional_bytes,
                                     cudaMemcpyHostToDevice))) != 0) goto cleanup;
  if ((ierr = cuda_status(cudaMemcpy(projection_device, projection_host, projection_bytes,
                                     cudaMemcpyHostToDevice))) != 0) goto cleanup;
  if ((ierr = cublas_status(cublasCreate(&handle))) != 0) goto cleanup;
  if ((ierr = cublas_status(cublasSetMathMode(handle, CUBLAS_TF32_TENSOR_OP_MATH))) != 0) goto cleanup;

  {
    const float alpha = 1.0f;
    const float beta = 0.0f;
    const cublasStatus_t gemm_status = cublasGemmEx(
        handle,
        CUBLAS_OP_N,
        CUBLAS_OP_N,
        nbin,
        nrow,
        ndirection,
        &alpha,
        projection_device,
        CUDA_R_32F,
        nbin,
        directional_device,
        CUDA_R_32F,
        ndirection,
        &beta,
        binned_device,
        CUDA_R_32F,
        nbin,
        CUBLAS_COMPUTE_32F_FAST_TF32,
        CUBLAS_GEMM_DEFAULT_TENSOR_OP);
    if ((ierr = cublas_status(gemm_status)) != 0) goto cleanup;
  }
  if ((ierr = cuda_status(cudaMemcpy(binned_host, binned_device, binned_bytes,
                                     cudaMemcpyDeviceToHost))) != 0) goto cleanup;

cleanup:
  if (handle != nullptr) cublasDestroy(handle);
  if (binned_device != nullptr) cudaFree(binned_device);
  if (projection_device != nullptr) cudaFree(projection_device);
  if (directional_device != nullptr) cudaFree(directional_device);
  return ierr;
}

__global__ void snrt_weighted_sum_fp32_kernel(const float *directional,
                                               const float *weights,
                                               float *scalar,
                                               int nrow,
                                               int ndirection) {
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= nrow) return;

  float total = 0.0f;
  const float *row_values = directional + static_cast<size_t>(row) * ndirection;
  for (int idir = 0; idir < ndirection; ++idir) {
    total = fmaf(row_values[idir], weights[idir], total);
  }
  scalar[row] = total;
}

extern "C" int snrt_cuda_weighted_sum_fp32_c(const float *directional_host,
                                               const float *weights_host,
                                               float *scalar_host,
                                               int nrow,
                                               int ndirection) {
  if (directional_host == nullptr || weights_host == nullptr ||
      scalar_host == nullptr || nrow <= 0 || ndirection <= 0) return 1;

  float *directional_device = nullptr;
  float *weights_device = nullptr;
  float *scalar_device = nullptr;
  int status = 0;
  const size_t directional_bytes = static_cast<size_t>(nrow) * ndirection * sizeof(float);
  const size_t weights_bytes = static_cast<size_t>(ndirection) * sizeof(float);
  const size_t scalar_bytes = static_cast<size_t>(nrow) * sizeof(float);

  if (cudaMalloc(&directional_device, directional_bytes) != cudaSuccess) {
    status = 2;
    goto cleanup;
  }
  if (cudaMalloc(&weights_device, weights_bytes) != cudaSuccess) {
    status = 3;
    goto cleanup;
  }
  if (cudaMalloc(&scalar_device, scalar_bytes) != cudaSuccess) {
    status = 4;
    goto cleanup;
  }
  if (cudaMemcpy(directional_device, directional_host, directional_bytes,
                 cudaMemcpyHostToDevice) != cudaSuccess ||
      cudaMemcpy(weights_device, weights_host, weights_bytes,
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    status = 5;
    goto cleanup;
  }

  {
    const int threads = 128;
    const int blocks = (nrow + threads - 1) / threads;
    snrt_weighted_sum_fp32_kernel<<<blocks, threads>>>(directional_device,
                                                         weights_device,
                                                         scalar_device,
                                                         nrow,
                                                         ndirection);
    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess) {
      status = 6;
      goto cleanup;
    }
  }
  if (cudaMemcpy(scalar_host, scalar_device, scalar_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) status = 7;

cleanup:
  if (scalar_device != nullptr) cudaFree(scalar_device);
  if (weights_device != nullptr) cudaFree(weights_device);
  if (directional_device != nullptr) cudaFree(directional_device);
  return status;
}

namespace {

__global__ void snrt_upwind_periodic_kernel(const float *state, float *next,
                                            const float *direction,
                                            int nx, int ny, int nz,
                                            int ndirection, float cdt_over_dx) {
  const long long ncell = static_cast<long long>(nx) * ny * nz;
  const long long total = ncell * ndirection;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int idir = static_cast<int>(linear / ncell);
    const int cell = static_cast<int>(linear - static_cast<long long>(idir) * ncell);
    const int ix = cell % nx;
    const int iy = (cell / nx) % ny;
    const int iz = cell / (nx * ny);
    const float mux = direction[3 * idir + 0];
    const float muy = direction[3 * idir + 1];
    const float muz = direction[3 * idir + 2];
    const int ix_up = mux >= 0.0f ? (ix + nx - 1) % nx : (ix + 1) % nx;
    const int iy_up = muy >= 0.0f ? (iy + ny - 1) % ny : (iy + 1) % ny;
    const int iz_up = muz >= 0.0f ? (iz + nz - 1) % nz : (iz + 1) % nz;
    const int x_up = ix_up + nx * (iy + ny * iz);
    const int y_up = ix + nx * (iy_up + ny * iz);
    const int z_up = ix + nx * (iy + ny * iz_up);
    const long long base = static_cast<long long>(idir) * ncell;
    const float q = state[linear];
    next[linear] = q - cdt_over_dx *
        (fabsf(mux) * (q - state[base + x_up]) +
         fabsf(muy) * (q - state[base + y_up]) +
         fabsf(muz) * (q - state[base + z_up]));
  }
}

}  // namespace

extern "C" int snrt_cuda_upwind_periodic_c(float *state_host,
                                             const float *direction_host,
                                             int nx, int ny, int nz,
                                             int ndirection,
                                             float cdt_over_dx) {
  if (state_host == nullptr || direction_host == nullptr || nx <= 0 || ny <= 0 ||
      nz <= 0 || ndirection <= 0 || cdt_over_dx < 0.0f) {
    return 1;
  }

  const long long ncell = static_cast<long long>(nx) * ny * nz;
  const long long total = ncell * ndirection;
  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t direction_bytes = static_cast<size_t>(3 * ndirection) * sizeof(float);
  float *state_device = nullptr;
  float *next_device = nullptr;
  float *direction_device = nullptr;
  int status = 1;

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&next_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&direction_device, direction_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_upwind_periodic_kernel<<<blocks, threads>>>(state_device, next_device,
        direction_device, nx, ny, nz, ndirection, cdt_over_dx);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  if (cudaDeviceSynchronize() != cudaSuccess) goto done;
  if (cudaMemcpy(state_host, next_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  cudaFree(direction_device);
  cudaFree(next_device);
  cudaFree(state_device);
  return status;
}

namespace {

__global__ void snrt_upwind_sparse_kernel(const float *state, float *next,
                                          const float *direction,
                                          const int *neighbor,
                                          int ncell, int ndirection,
                                          float cdt_over_dx) {
  const long long total = static_cast<long long>(ncell) * ndirection;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int idir = static_cast<int>(linear / ncell);
    const int cell = static_cast<int>(linear - static_cast<long long>(idir) * ncell);
    const float mux = direction[3 * idir + 0];
    const float muy = direction[3 * idir + 1];
    const float muz = direction[3 * idir + 2];
    const int x_up_fortran = neighbor[6 * cell + (mux >= 0.0f ? 0 : 1)];
    const int y_up_fortran = neighbor[6 * cell + (muy >= 0.0f ? 2 : 3)];
    const int z_up_fortran = neighbor[6 * cell + (muz >= 0.0f ? 4 : 5)];
    // The AMR-facing API uses Fortran's 1..ncell convention; zero denotes
    // a physical boundary whose ghost value is supplied as the cell value.
    const int x_up = x_up_fortran > 0 ? x_up_fortran - 1 : -1;
    const int y_up = y_up_fortran > 0 ? y_up_fortran - 1 : -1;
    const int z_up = z_up_fortran > 0 ? z_up_fortran - 1 : -1;
    const long long base = static_cast<long long>(idir) * ncell;
    const float q = state[linear];
    const float qx = x_up >= 0 ? state[base + x_up] : q;
    const float qy = y_up >= 0 ? state[base + y_up] : q;
    const float qz = z_up >= 0 ? state[base + z_up] : q;
    next[linear] = q - cdt_over_dx *
        (fabsf(mux) * (q - qx) + fabsf(muy) * (q - qy) + fabsf(muz) * (q - qz));
  }
}

}  // namespace

extern "C" int snrt_cuda_upwind_sparse_c(float *state_host,
                                           const float *direction_host,
                                           const int *neighbor_host,
                                           int ncell, int ndirection,
                                           float cdt_over_dx) {
  if (state_host == nullptr || direction_host == nullptr || neighbor_host == nullptr ||
      ncell <= 0 || ndirection <= 0 || cdt_over_dx < 0.0f) {
    return 1;
  }

  const long long total = static_cast<long long>(ncell) * ndirection;
  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t direction_bytes = static_cast<size_t>(3 * ndirection) * sizeof(float);
  const size_t neighbor_bytes = static_cast<size_t>(6 * ncell) * sizeof(int);
  float *state_device = nullptr;
  float *next_device = nullptr;
  float *direction_device = nullptr;
  int *neighbor_device = nullptr;
  int status = 1;

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&next_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&direction_device, direction_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&neighbor_device, neighbor_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neighbor_device, neighbor_host, neighbor_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_upwind_sparse_kernel<<<blocks, threads>>>(state_device, next_device,
        direction_device, neighbor_device, ncell, ndirection, cdt_over_dx);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  if (cudaDeviceSynchronize() != cudaSuccess) goto done;
  if (cudaMemcpy(state_host, next_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  cudaFree(neighbor_device);
  cudaFree(direction_device);
  cudaFree(next_device);
  cudaFree(state_device);
  return status;
}

namespace {

__global__ void snrt_absorb_kernel(float *state, float *absorbed_direction,
                                   const float *optical_depth,
                                   int ncell, int ndirection) {
  const long long total = static_cast<long long>(ncell) * ndirection;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int cell = static_cast<int>(linear % ncell);
    const float q_before = state[linear];
    const float tau = fmaxf(0.0f, optical_depth[cell]);
    const float q_after = q_before * expf(-tau);
    state[linear] = q_after;
    absorbed_direction[linear] = q_before - q_after;
  }
}

__global__ void snrt_reduce_direction_kernel(const float *directional,
                                             float *scalar,
                                             int ncell, int ndirection) {
  const int cell = blockIdx.x;
  if (cell >= ncell) return;
  float partial = 0.0f;
  for (int idir = threadIdx.x; idir < ndirection; idir += blockDim.x) {
    partial += directional[static_cast<long long>(idir) * ncell + cell];
  }
  __shared__ float cache[128];
  cache[threadIdx.x] = partial;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (threadIdx.x < stride) cache[threadIdx.x] += cache[threadIdx.x + stride];
    __syncthreads();
  }
  if (threadIdx.x == 0) scalar[cell] = cache[0];
}

}  // namespace

extern "C" int snrt_cuda_absorb_c(float *state_host,
                                    const float *optical_depth_host,
                                    float *absorbed_host,
                                    int ncell, int ndirection) {
  if (state_host == nullptr || optical_depth_host == nullptr || absorbed_host == nullptr ||
      ncell <= 0 || ndirection <= 0) {
    return 1;
  }

  const long long total = static_cast<long long>(ncell) * ndirection;
  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t cell_bytes = static_cast<size_t>(ncell) * sizeof(float);
  float *state_device = nullptr;
  float *tau_device = nullptr;
  float *absorbed_direction_device = nullptr;
  float *absorbed_device = nullptr;
  int status = 1;

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_direction_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_device, optical_depth_host, cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_absorb_kernel<<<blocks, threads>>>(state_device, absorbed_direction_device,
        tau_device, ncell, ndirection);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  snrt_reduce_direction_kernel<<<ncell, 128>>>(absorbed_direction_device,
      absorbed_device, ncell, ndirection);
  if (cudaGetLastError() != cudaSuccess) goto done;
  if (cudaDeviceSynchronize() != cudaSuccess) goto done;
  if (cudaMemcpy(state_host, state_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_host, absorbed_device, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  cudaFree(absorbed_device);
  cudaFree(absorbed_direction_device);
  cudaFree(tau_device);
  cudaFree(state_device);
  return status;
}

extern "C" int snrt_cuda_transport_absorb_c(float *state_host,
                                              const float *direction_host,
                                              const int *neighbor_host,
                                              const float *optical_depth_host,
                                              float *absorbed_host,
                                              int ncell, int ndirection,
                                              float cdt_over_dx) {
  if (state_host == nullptr || direction_host == nullptr || neighbor_host == nullptr ||
      optical_depth_host == nullptr || absorbed_host == nullptr || ncell <= 0 ||
      ndirection <= 0 || cdt_over_dx < 0.0f) {
    return 1;
  }

  const long long total = static_cast<long long>(ncell) * ndirection;
  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t direction_bytes = static_cast<size_t>(3 * ndirection) * sizeof(float);
  const size_t neighbor_bytes = static_cast<size_t>(6 * ncell) * sizeof(int);
  const size_t cell_bytes = static_cast<size_t>(ncell) * sizeof(float);
  float *state_device = nullptr;
  float *transport_device = nullptr;
  float *direction_device = nullptr;
  int *neighbor_device = nullptr;
  float *tau_device = nullptr;
  float *absorbed_device = nullptr;
  int status = 1;

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&transport_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&direction_device, direction_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&neighbor_device, neighbor_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neighbor_device, neighbor_host, neighbor_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_device, optical_depth_host, cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_upwind_sparse_kernel<<<blocks, threads>>>(state_device, transport_device,
        direction_device, neighbor_device, ncell, ndirection, cdt_over_dx);
    if (cudaGetLastError() != cudaSuccess) goto done;
    snrt_absorb_kernel<<<blocks, threads>>>(transport_device, state_device,
        tau_device, ncell, ndirection);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  snrt_reduce_direction_kernel<<<ncell, 128>>>(state_device, absorbed_device,
      ncell, ndirection);
  if (cudaGetLastError() != cudaSuccess) goto done;
  if (cudaDeviceSynchronize() != cudaSuccess) goto done;
  if (cudaMemcpy(state_host, transport_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_host, absorbed_device, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  cudaFree(absorbed_device);
  cudaFree(tau_device);
  cudaFree(neighbor_device);
  cudaFree(direction_device);
  cudaFree(transport_device);
  cudaFree(state_device);
  return status;
}

namespace {

__global__ void snrt_limit_optical_depth_kernel(const float *photon_total,
                                                const float *optical_depth,
                                                const float *neutral_hydrogen,
                                                float *limited_optical_depth,
                                                int ncell) {
  const int cell = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell >= ncell) return;
  float tau = fmaxf(0.0f, optical_depth[cell]);
  const float photons = photon_total[cell];
  const float neutral = fmaxf(0.0f, neutral_hydrogen[cell]);
  if (photons > 0.0f && neutral < photons) {
    const float fraction = fmaxf(0.0f, neutral / photons);
    tau = fminf(tau, -log1pf(-fraction));
  }
  limited_optical_depth[cell] = tau;
}

__global__ void snrt_cap_absorption_kernel(float *state, float *absorbed_direction,
                                            const float *absorbed_total,
                                            const float *neutral_hydrogen,
                                            int ncell, int ndirection) {
  const long long total = static_cast<long long>(ncell) * ndirection;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int cell = static_cast<int>(linear % ncell);
    const float total_absorbed = absorbed_total[cell];
    const float neutral = fmaxf(0.0f, neutral_hydrogen[cell]);
    float cap = 1.0f;
    if (total_absorbed > neutral && total_absorbed > 0.0f) {
      // Leave a small FP32 guard band after the reduction so the returned
      // chemistry source never consumes more atoms than are available.
      cap = fmaxf(0.0f, fminf(1.0f, (neutral / total_absorbed) * 0.99999f));
    }
    const float removed = absorbed_direction[linear];
    const float limited_removed = removed * cap;
    state[linear] += removed - limited_removed;
    absorbed_direction[linear] = limited_removed;
  }
}

}  // namespace

extern "C" int snrt_cuda_transport_absorb_limited_c(float *state_host,
                                                      const float *direction_host,
                                                      const int *neighbor_host,
                                                      const float *optical_depth_host,
                                                      const float *neutral_hydrogen_host,
                                                      float *absorbed_host,
                                                      int ncell, int ndirection,
                                                      float cdt_over_dx) {
  if (state_host == nullptr || direction_host == nullptr || neighbor_host == nullptr ||
      optical_depth_host == nullptr || neutral_hydrogen_host == nullptr ||
      absorbed_host == nullptr || ncell <= 0 || ndirection <= 0 || cdt_over_dx < 0.0f) {
    return 1;
  }

  const long long total = static_cast<long long>(ncell) * ndirection;
  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t direction_bytes = static_cast<size_t>(3 * ndirection) * sizeof(float);
  const size_t neighbor_bytes = static_cast<size_t>(6 * ncell) * sizeof(int);
  const size_t cell_bytes = static_cast<size_t>(ncell) * sizeof(float);
  float *state_device = nullptr;
  float *transport_device = nullptr;
  float *direction_device = nullptr;
  int *neighbor_device = nullptr;
  float *tau_device = nullptr;
  float *neutral_device = nullptr;
  float *photon_total_device = nullptr;
  float *limited_tau_device = nullptr;
  float *absorbed_device = nullptr;
  int status = 1;

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&transport_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&direction_device, direction_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&neighbor_device, neighbor_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&neutral_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&photon_total_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&limited_tau_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neighbor_device, neighbor_host, neighbor_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_device, optical_depth_host, cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neutral_device, neutral_hydrogen_host, cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    const int cell_blocks = (ncell + threads - 1) / threads;
    snrt_upwind_sparse_kernel<<<blocks, threads>>>(state_device, transport_device,
        direction_device, neighbor_device, ncell, ndirection, cdt_over_dx);
    if (cudaGetLastError() != cudaSuccess) goto done;
    snrt_reduce_direction_kernel<<<ncell, 128>>>(transport_device,
        photon_total_device, ncell, ndirection);
    if (cudaGetLastError() != cudaSuccess) goto done;
    snrt_limit_optical_depth_kernel<<<cell_blocks, threads>>>(photon_total_device,
        tau_device, neutral_device, limited_tau_device, ncell);
    if (cudaGetLastError() != cudaSuccess) goto done;
    snrt_absorb_kernel<<<blocks, threads>>>(transport_device, state_device,
        limited_tau_device, ncell, ndirection);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  snrt_reduce_direction_kernel<<<ncell, 128>>>(state_device, absorbed_device,
      ncell, ndirection);
  if (cudaGetLastError() != cudaSuccess) goto done;
  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_cap_absorption_kernel<<<blocks, threads>>>(transport_device, state_device,
        absorbed_device, neutral_device, ncell, ndirection);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  snrt_reduce_direction_kernel<<<ncell, 128>>>(state_device, absorbed_device,
      ncell, ndirection);
  if (cudaGetLastError() != cudaSuccess) goto done;
  if (cudaDeviceSynchronize() != cudaSuccess) goto done;
  if (cudaMemcpy(state_host, transport_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_host, absorbed_device, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  cudaFree(absorbed_device);
  cudaFree(limited_tau_device);
  cudaFree(photon_total_device);
  cudaFree(neutral_device);
  cudaFree(tau_device);
  cudaFree(neighbor_device);
  cudaFree(direction_device);
  cudaFree(transport_device);
  cudaFree(state_device);
  return status;
}

namespace {

__global__ void snrt_multigroup_upwind_kernel(const float *state, float *next,
                                               const float *direction,
                                               const int *neighbor,
                                               int nowned, int nwork, int ndirection,
                                               int ngroup, float cdt_over_dx) {
  const long long per_group = static_cast<long long>(nwork) * ndirection;
  const long long total = per_group * ngroup;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int idir = static_cast<int>((linear / nwork) % ndirection);
    const int cell = static_cast<int>(linear % nwork);
    if (cell >= nowned) continue;
    const float mux = direction[3 * idir + 0];
    const float muy = direction[3 * idir + 1];
    const float muz = direction[3 * idir + 2];
    const int x_up_fortran = neighbor[6 * cell + (mux >= 0.0f ? 0 : 1)];
    const int y_up_fortran = neighbor[6 * cell + (muy >= 0.0f ? 2 : 3)];
    const int z_up_fortran = neighbor[6 * cell + (muz >= 0.0f ? 4 : 5)];
    const int x_up = x_up_fortran > 0 ? x_up_fortran - 1 : -1;
    const int y_up = y_up_fortran > 0 ? y_up_fortran - 1 : -1;
    const int z_up = z_up_fortran > 0 ? z_up_fortran - 1 : -1;
    const long long base = linear - cell;
    const float q = state[linear];
    const float qx = x_up >= 0 ? state[base + x_up] : q;
    const float qy = y_up >= 0 ? state[base + y_up] : q;
    const float qz = z_up >= 0 ? state[base + z_up] : q;
    next[linear] = q - cdt_over_dx *
        (fabsf(mux) * (q - qx) + fabsf(muy) * (q - qy) + fabsf(muz) * (q - qz));
  }
}

__global__ void snrt_multigroup_absorb_kernel(float *state, float *absorbed_direction,
                                               const float *optical_depth,
                                               int nowned, int nwork, int ndirection,
                                               int ngroup) {
  const long long per_group = static_cast<long long>(nwork) * ndirection;
  const long long total = per_group * ngroup;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int cell = static_cast<int>(linear % nwork);
    if (cell >= nowned) continue;
    const int group = static_cast<int>(linear / per_group);
    const float q_before = state[linear];
    const float q_after = q_before * expf(-fmaxf(0.0f, optical_depth[group * nowned + cell]));
    state[linear] = q_after;
    absorbed_direction[linear] = q_before - q_after;
  }
}

__global__ void snrt_reduce_multigroup_kernel(const float *directional, float *scalar,
                                              int nowned, int nwork,
                                              int ndirection, int ngroup) {
  const int cell = blockIdx.x;
  if (cell >= nowned) return;
  const long long per_group = static_cast<long long>(nwork) * ndirection;
  float partial = 0.0f;
  for (int group = 0; group < ngroup; ++group) {
    const long long group_base = static_cast<long long>(group) * per_group;
    for (int idir = threadIdx.x; idir < ndirection; idir += blockDim.x) {
      partial += directional[group_base + static_cast<long long>(idir) * nwork + cell];
    }
  }
  __shared__ float cache[128];
  cache[threadIdx.x] = partial;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (threadIdx.x < stride) cache[threadIdx.x] += cache[threadIdx.x + stride];
    __syncthreads();
  }
  if (threadIdx.x == 0) scalar[cell] = cache[0];
}

__global__ void snrt_reduce_multigroup_group_kernel(const float *directional,
                                                     float *scalar_group,
                                                     int nowned, int nwork,
                                                     int ndirection,
                                                     int ngroup) {
  const long long total = static_cast<long long>(nowned) * ngroup;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int cell = static_cast<int>(linear % nowned);
    const int group = static_cast<int>(linear / nowned);
    const long long group_base = static_cast<long long>(group) * nwork * ndirection;
    float sum = 0.0f;
    for (int idir = 0; idir < ndirection; ++idir) {
      sum += directional[group_base + static_cast<long long>(idir) * nwork + cell];
    }
    scalar_group[linear] = sum;
  }
}

__global__ void snrt_cap_multigroup_absorption_kernel(float *state,
                                                      float *absorbed_direction,
                                                      const float *absorbed_total,
                                                      const float *neutral_hydrogen,
                                                      int nowned, int nwork,
                                                      int ndirection,
                                                      int ngroup) {
  const long long per_group = static_cast<long long>(nwork) * ndirection;
  const long long total = per_group * ngroup;
  for (long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
       linear < total;
       linear += static_cast<long long>(gridDim.x) * blockDim.x) {
    const int cell = static_cast<int>(linear % nwork);
    if (cell >= nowned) continue;
    const float total_absorbed = absorbed_total[cell];
    const float neutral = fmaxf(0.0f, neutral_hydrogen[cell]);
    float cap = 1.0f;
    if (total_absorbed > neutral && total_absorbed > 0.0f) {
      cap = fmaxf(0.0f, fminf(1.0f, (neutral / total_absorbed) * 0.99999f));
    }
    const float removed = absorbed_direction[linear];
    const float limited_removed = removed * cap;
    state[linear] += removed - limited_removed;
    absorbed_direction[linear] = limited_removed;
  }
}

// Apply the production H/He inventory cap.  The transport solve is still
// multigroup, but the inventory is consumed in group order so that a group
// can only use species with non-zero optical depth in that group.  The
// per-cell thread owns all groups and therefore updates the three species
// inventories without atomics or a cross-group race.
__global__ void snrt_cap_multigroup_species_absorption_kernel(
    float *state, float *absorbed_direction,
    const float *optical_depth_species, float *available_species,
    int nowned, int nwork, int ndirection, int ngroup) {
  const int cell = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell >= nowned) return;

  float available[3];
  for (int species = 0; species < 3; ++species) {
    available[species] = fmaxf(0.0f, available_species[species * nowned + cell]);
  }

  for (int group = 0; group < ngroup; ++group) {
    const long long group_base = static_cast<long long>(group) * nwork * ndirection;
    float raw_absorbed = 0.0f;
    for (int idir = 0; idir < ndirection; ++idir) {
      raw_absorbed += fmaxf(0.0f,
          absorbed_direction[group_base + static_cast<long long>(idir) * nwork + cell]);
    }
    if (raw_absorbed <= 0.0f) continue;

    float opacity[3];
    float opacity_sum = 0.0f;
    float eligible_inventory = 0.0f;
    for (int species = 0; species < 3; ++species) {
      const long long opacity_index =
          (static_cast<long long>(species) * ngroup + group) * nowned + cell;
      opacity[species] = fmaxf(0.0f, optical_depth_species[opacity_index]);
      if (opacity[species] > 0.0f) {
        opacity_sum += opacity[species];
        eligible_inventory += available[species];
      }
    }
    // Use the cell-scale inventory magnitude for the residual guard.  A
    // fixed code-density threshold is not scale invariant and can either
    // discard a resolvable remainder or leave an FP32-sized overrun for the
    // double-precision Fortran partition.
    float inventory_scale = fmaxf(fabsf(raw_absorbed), eligible_inventory);
    for (int species = 0; species < 3; ++species) {
      inventory_scale = fmaxf(inventory_scale, fabsf(available[species]));
    }
    const float remainder_tolerance =
        256.0f * FLT_EPSILON * fmaxf(inventory_scale, FLT_MIN);

    float target_absorbed = 0.0f;
    if (opacity_sum > 0.0f && eligible_inventory > 0.0f) {
      // Keep a deterministic FP32 guard band whenever the raw absorption
      // reaches an inventory boundary.  This prevents a later Fortran
      // double-precision partition from seeing a one-ulp overrun.
      const float safe_inventory = eligible_inventory * 0.99995f;
      target_absorbed = fminf(raw_absorbed, safe_inventory);
    }

    // Mirror snrt_partition_absorption: first allocate by opacity, then
    // redistribute only over species that are eligible in this group.  The
    // active-set loop is bounded by the three species: a saturated species
    // is removed from the opacity-weighted pool, and the remaining target is
    // recomputed over the unsaturated species.  This guarantees that the
    // directional cap and the inventory decrement use the same assigned
    // amount; a partially assigned target is returned to the photon field.
    float assigned[3] = {0.0f, 0.0f, 0.0f};
    float remaining = target_absorbed;
    bool active[3] = {opacity[0] > 0.0f && available[0] > 0.0f,
                      opacity[1] > 0.0f && available[1] > 0.0f,
                      opacity[2] > 0.0f && available[2] > 0.0f};
    float active_weight = 0.0f;
    for (int species = 0; species < 3; ++species) {
      if (active[species]) active_weight += opacity[species];
    }
    for (int pass = 0; pass < 3 && remaining > remainder_tolerance &&
         active_weight > 0.0f; ++pass) {
      const float pass_remaining = remaining;
      const float pass_weight = active_weight;
      bool saturated_any = false;
      bool saturated[3] = {false, false, false};
      for (int species = 0; species < 3; ++species) {
        if (!active[species]) continue;
        const float headroom = fmaxf(0.0f, available[species] - assigned[species]);
        if (headroom <= remainder_tolerance) {
          saturated[species] = true;
          saturated_any = true;
          continue;
        }
        const float requested = pass_remaining * opacity[species] / pass_weight;
        const float addition = fminf(requested, headroom);
        assigned[species] += addition;
        if (addition + remainder_tolerance >= headroom) {
          saturated[species] = true;
          saturated_any = true;
        }
      }
      remaining = target_absorbed - (assigned[0] + assigned[1] + assigned[2]);
      if (remaining < 0.0f) remaining = 0.0f;
      if (!saturated_any) {
        remaining = 0.0f;
        break;
      }
      for (int species = 0; species < 3; ++species) {
        if (saturated[species]) {
          active[species] = false;
          active_weight -= opacity[species];
        }
      }
    }
    // The target is below the total eligible inventory, so any residual here
    // is roundoff or a zero-headroom edge case.  Fill it deterministically
    // without ever exceeding a species inventory; the assigned sum remains
    // authoritative for the photon-field cap.
    if (remaining > remainder_tolerance) {
      for (int species = 0; species < 3 && remaining > 0.0f; ++species) {
        if (!active[species]) continue;
        const float headroom = fmaxf(0.0f, available[species] - assigned[species]);
        const float addition = fminf(headroom, remaining);
        assigned[species] += addition;
        remaining -= addition;
      }
    }
    const float limited_absorbed = assigned[0] + assigned[1] + assigned[2];
    const float cap = limited_absorbed > 0.0f ? limited_absorbed / raw_absorbed : 0.0f;
    for (int idir = 0; idir < ndirection; ++idir) {
      const long long index = group_base + static_cast<long long>(idir) * nwork + cell;
      const float removed = absorbed_direction[index];
      const float limited_removed = removed * cap;
      state[index] += removed - limited_removed;
      absorbed_direction[index] = limited_removed;
    }

    for (int species = 0; species < 3; ++species) {
      available[species] = fmaxf(0.0f, available[species] - assigned[species]);
    }
  }

  for (int species = 0; species < 3; ++species) {
    available_species[species * nowned + cell] = available[species];
  }
}

// Validate the new four-component optical-depth contract before the device
// state is touched.  The total tau is accepted from the caller so the DUST-7
// zero-dust path can be bitwise compared with the legacy ABI; the component
// sum is checked with an FP32-relative tolerance.
__global__ void snrt_validate_species_dust_inputs_kernel(
    const float *state, const float *direction, const int *neighbor,
    const float *optical_depth, const float *optical_depth_species,
    const float *optical_depth_dust, const float *available_species,
    int nowned, int nwork, int ndirection, int ngroup, int *invalid) {
  const long long group_count = static_cast<long long>(nowned) * ngroup;
  const long long species_count = 3 * group_count;
  const long long inventory_count = 3 * static_cast<long long>(nowned);
  const long long state_count = static_cast<long long>(nwork) * ndirection * ngroup;
  const long long direction_count = 3 * static_cast<long long>(ndirection);
  const long long neighbor_count = 6 * static_cast<long long>(nowned);
  long long total = state_count;
  if (direction_count > total) total = direction_count;
  if (neighbor_count > total) total = neighbor_count;
  if (group_count > total) total = group_count;
  if (inventory_count > total) total = inventory_count;
  if (species_count > total) total = species_count;
  const long long linear = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= total) return;

  if (linear < state_count) {
    const float value = state[linear];
    if (!isfinite(value) || value < 0.0f) atomicExch(invalid, 1);
  }
  if (linear < direction_count) {
    if (!isfinite(direction[linear])) atomicExch(invalid, 1);
  }
  if (linear < neighbor_count) {
    const int value = neighbor[linear];
    if (value < 0 || value > nwork) atomicExch(invalid, 1);
  }

  if (linear < group_count) {
    const int cell = static_cast<int>(linear % nowned);
    const int group = static_cast<int>(linear / nowned);
    const float total_tau = optical_depth[linear];
    const float dust_tau = optical_depth_dust[linear];
    const long long hhe_base = static_cast<long long>(group) * nowned + cell;
    const float component_sum = optical_depth_species[hhe_base] +
        optical_depth_species[group_count + hhe_base] +
        optical_depth_species[2 * group_count + hhe_base] + dust_tau;
    const float scale = fmaxf(fmaxf(fabsf(total_tau), fabsf(component_sum)), FLT_MIN);
    if (!isfinite(total_tau) || total_tau < 0.0f || !isfinite(dust_tau) ||
        dust_tau < 0.0f || !isfinite(component_sum) ||
        fabsf(total_tau - component_sum) > 8.0f * FLT_EPSILON * scale) {
      atomicExch(invalid, 1);
    }
  }
  if (linear < species_count) {
    const float value = optical_depth_species[linear];
    if (!isfinite(value) || value < 0.0f) atomicExch(invalid, 1);
  }
  if (linear < inventory_count) {
    const float value = available_species[linear];
    if (!isfinite(value) || value < 0.0f) atomicExch(invalid, 1);
  }
}

// DUST-7 cap and ledger kernel.  This is deliberately separate from the
// legacy three-species kernel so its active-set arithmetic cannot perturb the
// old ABI.  One thread owns all groups of a cell, preserving group-order
// reservoir consumption without atomics.
__global__ void snrt_cap_multigroup_species_dust_absorption_kernel(
    float *state, float *absorbed_direction,
    const float *optical_depth_species, const float *optical_depth_dust,
    float *available_species, float *absorbed_hhe_species,
    float *absorbed_dust_group, float *returned_group, float *raw_group,
    float *absorbed_group, float *absorbed_total,
    int nowned, int nwork, int ndirection, int ngroup) {
  const int cell = blockIdx.x * blockDim.x + threadIdx.x;
  snrt_cap_species_dust_cell(state,absorbed_direction,optical_depth_species,optical_depth_dust,
      available_species,absorbed_hhe_species,absorbed_dust_group,returned_group,raw_group,
      absorbed_group,absorbed_total,nowned,nwork,ndirection,ngroup,cell);
}

}  // namespace

static int snrt_cuda_multigroup_rt_step_impl(float *state_host,
                                               const float *direction_host,
                                               const int *neighbor_host,
                                               const float *optical_depth_host,
                                               const float *neutral_hydrogen_host,
                                               float *absorbed_host,
                                               float *absorbed_group_host,
                                               int nowned, int nwork,
                                               int ndirection, int ngroup,
                                               float cdt_over_dx,
                                               const float *optical_depth_species_host,
                                               float *available_species_host) {
  const bool species_aware = optical_depth_species_host != nullptr &&
                             available_species_host != nullptr;
  if ((optical_depth_species_host == nullptr) !=
      (available_species_host == nullptr)) {
    return 1;
  }
  if (state_host == nullptr || direction_host == nullptr || neighbor_host == nullptr ||
      optical_depth_host == nullptr || (!species_aware && neutral_hydrogen_host == nullptr) ||
      absorbed_host == nullptr || absorbed_group_host == nullptr ||
      nowned <= 0 || nwork < nowned || ndirection <= 0 || ngroup <= 0 ||
      cdt_over_dx < 0.0f) {
    return 1;
  }

  const long long per_group = static_cast<long long>(nwork) * ndirection;
  const long long total = per_group * ngroup;
  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t direction_bytes = static_cast<size_t>(3 * ndirection) * sizeof(float);
  const size_t neighbor_bytes = static_cast<size_t>(6 * nowned) * sizeof(int);
  const size_t group_bytes = static_cast<size_t>(nowned) * ngroup * sizeof(float);
  const size_t cell_bytes = static_cast<size_t>(nowned) * sizeof(float);
  const size_t species_bytes = static_cast<size_t>(3) * ngroup * nowned * sizeof(float);
  const size_t inventory_bytes = static_cast<size_t>(3) * nowned * sizeof(float);
  float *state_device = nullptr;
  float *transport_device = nullptr;
  float *direction_device = nullptr;
  int *neighbor_device = nullptr;
  float *tau_device = nullptr;
  float *neutral_device = nullptr;
  float *absorbed_device = nullptr;
  float *absorbed_group_device = nullptr;
  float *optical_depth_species_device = nullptr;
  float *available_species_device = nullptr;
  cudaEvent_t event_start = nullptr;
  cudaEvent_t event_stop = nullptr;
  const char *trace_env = std::getenv("SNRT_CUDA_TRACE");
  const bool trace = trace_env != nullptr && trace_env[0] == '1' && trace_env[1] == '\0';
  int status = 1;

  if (trace) {
    if (cudaEventCreate(&event_start) != cudaSuccess ||
        cudaEventCreate(&event_stop) != cudaSuccess) {
      if (event_stop != nullptr) cudaEventDestroy(event_stop);
      if (event_start != nullptr) cudaEventDestroy(event_start);
      event_start = nullptr;
      event_stop = nullptr;
    }
  }

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&transport_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&direction_device, direction_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&neighbor_device, neighbor_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_device, group_bytes) != cudaSuccess) goto done;
  if (!species_aware && cudaMalloc(&neutral_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_group_device, group_bytes) != cudaSuccess) goto done;
  if (species_aware && cudaMalloc(&optical_depth_species_device, species_bytes) != cudaSuccess) goto done;
  if (species_aware && cudaMalloc(&available_species_device, inventory_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neighbor_device, neighbor_host, neighbor_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_device, optical_depth_host, group_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (!species_aware && cudaMemcpy(neutral_device, neutral_hydrogen_host, cell_bytes,
                                    cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (species_aware && cudaMemcpy(optical_depth_species_device, optical_depth_species_host,
                                  species_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (species_aware && cudaMemcpy(available_species_device, available_species_host,
                                  inventory_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

  if (event_start != nullptr) cudaEventRecord(event_start);
  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_multigroup_upwind_kernel<<<blocks, threads>>>(state_device, transport_device,
        direction_device, neighbor_device, nowned, nwork, ndirection, ngroup, cdt_over_dx);
    if (cudaGetLastError() != cudaSuccess) goto done;
    snrt_multigroup_absorb_kernel<<<blocks, threads>>>(transport_device, state_device,
        tau_device, nowned, nwork, ndirection, ngroup);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  snrt_reduce_multigroup_kernel<<<nowned, 128>>>(state_device, absorbed_device,
      nowned, nwork, ndirection, ngroup);
  if (cudaGetLastError() != cudaSuccess) goto done;
  {
    const long long group_total = static_cast<long long>(nowned) * ngroup;
    const int threads = 256;
    const int blocks = static_cast<int>((group_total + threads - 1) / threads);
    snrt_reduce_multigroup_group_kernel<<<blocks, threads>>>(state_device,
        absorbed_group_device, nowned, nwork, ndirection, ngroup);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  {
    const int threads = 256;
    if (species_aware) {
      const int blocks = (nowned + threads - 1) / threads;
      snrt_cap_multigroup_species_absorption_kernel<<<blocks, threads>>>(transport_device,
          state_device, optical_depth_species_device, available_species_device,
          nowned, nwork, ndirection, ngroup);
    } else {
      const int blocks = static_cast<int>((total + threads - 1) / threads);
      snrt_cap_multigroup_absorption_kernel<<<blocks, threads>>>(transport_device, state_device,
          absorbed_device, neutral_device, nowned, nwork, ndirection, ngroup);
    }
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  // Recompute the returned absorption from the capped removal field.  In this
  // wrapper state_device is the absorbed_direction argument of the absorb and
  // cap kernels; transport_device is the surviving photon field returned to
  // the caller.  Summing transport_device here would invert the ledger.
  snrt_reduce_multigroup_kernel<<<nowned, 128>>>(state_device, absorbed_device,
      nowned, nwork, ndirection, ngroup);
  if (cudaGetLastError() != cudaSuccess) goto done;
  {
    const long long group_total = static_cast<long long>(nowned) * ngroup;
    const int threads = 256;
    const int blocks = static_cast<int>((group_total + threads - 1) / threads);
    snrt_reduce_multigroup_group_kernel<<<blocks, threads>>>(state_device,
        absorbed_group_device, nowned, nwork, ndirection, ngroup);
  }
  if (cudaGetLastError() != cudaSuccess) goto done;
  if (cudaDeviceSynchronize() != cudaSuccess) goto done;
  if (event_stop != nullptr) {
    float elapsed_ms = 0.0f;
    cudaEventRecord(event_stop);
    if (cudaEventSynchronize(event_stop) == cudaSuccess &&
      cudaEventElapsedTime(&elapsed_ms, event_start, event_stop) == cudaSuccess) {
      std::fprintf(stderr,
                   "SNRT CUDA multigroup kernels nowned=%d nwork=%d ndirection=%d ngroup=%d launches=7 ms=%.3f\n",
                   nowned, nwork, ndirection, ngroup, elapsed_ms);
    }
  }
  if (cudaMemcpy(state_host, transport_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_host, absorbed_device, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_group_host, absorbed_group_device, group_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (species_aware && cudaMemcpy(available_species_host, available_species_device,
                                  inventory_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  if (event_stop != nullptr) cudaEventDestroy(event_stop);
  if (event_start != nullptr) cudaEventDestroy(event_start);
  cudaFree(absorbed_group_device);
  cudaFree(absorbed_device);
  cudaFree(available_species_device);
  cudaFree(optical_depth_species_device);
  cudaFree(neutral_device);
  cudaFree(tau_device);
  cudaFree(neighbor_device);
  cudaFree(direction_device);
  cudaFree(transport_device);
  cudaFree(state_device);
  return status;
}

extern "C" int snrt_cuda_multigroup_rt_step_c(float *state_host,
                                                const float *direction_host,
                                                const int *neighbor_host,
                                                const float *optical_depth_host,
                                                const float *neutral_hydrogen_host,
                                                float *absorbed_host,
                                                float *absorbed_group_host,
                                                int ncell, int ndirection, int ngroup,
                                                float cdt_over_dx) {
  return snrt_cuda_multigroup_rt_step_impl(state_host, direction_host, neighbor_host,
      optical_depth_host, neutral_hydrogen_host, absorbed_host, absorbed_group_host,
      ncell, ncell, ndirection, ngroup, cdt_over_dx, nullptr, nullptr);
}

extern "C" int snrt_cuda_multigroup_rt_step_owned_c(float *state_host,
                                                      const float *direction_host,
                                                      const int *neighbor_host,
                                                      const float *optical_depth_host,
                                                      const float *neutral_hydrogen_host,
                                                      float *absorbed_host,
                                                      float *absorbed_group_host,
                                                      int nowned, int nwork,
                                                      int ndirection, int ngroup,
                                                      float cdt_over_dx) {
  return snrt_cuda_multigroup_rt_step_impl(state_host, direction_host, neighbor_host,
      optical_depth_host, neutral_hydrogen_host, absorbed_host, absorbed_group_host,
      nowned, nwork, ndirection, ngroup, cdt_over_dx, nullptr, nullptr);
}

extern "C" int snrt_cuda_multigroup_rt_step_species_c(
    float *state_host, const float *direction_host, const int *neighbor_host,
    const float *optical_depth_host, const float *optical_depth_species_host,
    float *available_species_host, float *absorbed_host,
    float *absorbed_group_host, int nowned, int nwork, int ndirection,
    int ngroup, float cdt_over_dx) {
  return snrt_cuda_multigroup_rt_step_impl(state_host, direction_host, neighbor_host,
      optical_depth_host, nullptr, absorbed_host, absorbed_group_host,
      nowned, nwork, ndirection, ngroup, cdt_over_dx,
      optical_depth_species_host, available_species_host);
}

// DUST-7 keeps the existing three-species ABI above untouched.  This wrapper
// owns the fourth (dust) component and returns separate ledgers so the host
// partition in DUST-8 does not have to reconstruct raw or returned photons.
extern "C" int snrt_cuda_multigroup_rt_step_species_dust_c(
    float *state_host, const float *direction_host, const int *neighbor_host,
    const float *optical_depth_host, const float *optical_depth_species_host,
    const float *optical_depth_dust_host, float *available_species_host,
    float *absorbed_hhe_species_host, float *absorbed_dust_group_host,
    float *returned_group_host, float *raw_group_host,
    float *absorbed_group_host, float *absorbed_host,
    int nowned, int nwork, int ndirection, int ngroup, float cdt_over_dx) {
  if (state_host == nullptr || direction_host == nullptr || neighbor_host == nullptr ||
      optical_depth_host == nullptr || optical_depth_species_host == nullptr ||
      optical_depth_dust_host == nullptr || available_species_host == nullptr ||
      absorbed_hhe_species_host == nullptr || absorbed_dust_group_host == nullptr ||
      returned_group_host == nullptr || raw_group_host == nullptr ||
      absorbed_group_host == nullptr || absorbed_host == nullptr ||
      nowned <= 0 || nwork < nowned || ndirection <= 0 || ngroup <= 0 ||
      !std::isfinite(cdt_over_dx) || cdt_over_dx < 0.0f) {
    return 1;
  }

  const long long max_long = std::numeric_limits<long long>::max();
  if (static_cast<long long>(nwork) > max_long / ndirection) return 1;
  const long long per_group = static_cast<long long>(nwork) * ndirection;
  if (per_group > max_long / ngroup) return 1;
  const long long total = per_group * ngroup;
  const long long group_count = static_cast<long long>(nowned) * ngroup;
  if (group_count > max_long / 3) return 1;
  const long long species_count = 3 * group_count;
  const long long inventory_count = 3 * static_cast<long long>(nowned);
  const long long direction_count = 3 * static_cast<long long>(ndirection);
  const long long neighbor_count = 6 * static_cast<long long>(nowned);
  if (per_group <= 0 || total <= 0 || group_count <= 0 ||
      total > static_cast<long long>(std::numeric_limits<size_t>::max() / sizeof(float)) ||
      group_count > static_cast<long long>(std::numeric_limits<size_t>::max() / sizeof(float)) ||
      species_count > static_cast<long long>(std::numeric_limits<size_t>::max() / sizeof(float)) ||
      inventory_count > static_cast<long long>(std::numeric_limits<size_t>::max() / sizeof(float))) {
    return 1;
  }

  const size_t state_bytes = static_cast<size_t>(total) * sizeof(float);
  const size_t direction_bytes = static_cast<size_t>(direction_count) * sizeof(float);
  const size_t neighbor_bytes = static_cast<size_t>(neighbor_count) * sizeof(int);
  const size_t group_bytes = static_cast<size_t>(group_count) * sizeof(float);
  const size_t species_bytes = static_cast<size_t>(species_count) * sizeof(float);
  const size_t inventory_bytes = static_cast<size_t>(inventory_count) * sizeof(float);
  float *state_device = nullptr;
  float *transport_device = nullptr;
  float *direction_device = nullptr;
  int *neighbor_device = nullptr;
  float *tau_device = nullptr;
  float *tau_species_device = nullptr;
  float *tau_dust_device = nullptr;
  float *available_species_device = nullptr;
  float *absorbed_hhe_species_device = nullptr;
  float *absorbed_dust_group_device = nullptr;
  float *returned_group_device = nullptr;
  float *raw_group_device = nullptr;
  float *absorbed_group_device = nullptr;
  float *absorbed_device = nullptr;
  int *invalid_device = nullptr;
  int status = 1;

  if (cudaMalloc(&state_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&transport_device, state_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&direction_device, direction_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&neighbor_device, neighbor_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_species_device, species_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&tau_dust_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&available_species_device, inventory_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_hhe_species_device, species_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_dust_group_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&returned_group_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&raw_group_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_group_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_device, static_cast<size_t>(nowned) * sizeof(float)) != cudaSuccess) goto done;
  if (cudaMalloc(&invalid_device, sizeof(int)) != cudaSuccess) goto done;

  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neighbor_device, neighbor_host, neighbor_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_device, optical_depth_host, group_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_species_device, optical_depth_species_host, species_bytes,
                 cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_dust_device, optical_depth_dust_host, group_bytes,
                 cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(available_species_device, available_species_host, inventory_bytes,
                 cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemset(invalid_device, 0, sizeof(int)) != cudaSuccess) goto done;

  {
    long long validation_total = total;
    if (direction_count > validation_total) validation_total = direction_count;
    if (neighbor_count > validation_total) validation_total = neighbor_count;
    if (group_count > validation_total) validation_total = group_count;
    if (species_count > validation_total) validation_total = species_count;
    if (inventory_count > validation_total) validation_total = inventory_count;
    const int threads = 256;
    const long long block_count = (validation_total + threads - 1) / threads;
    if (block_count <= 0 || block_count > 2147483647LL) goto done;
    snrt_validate_species_dust_inputs_kernel<<<static_cast<int>(block_count), threads>>>(
        state_device, direction_device, neighbor_device, tau_device,
        tau_species_device, tau_dust_device, available_species_device,
        nowned, nwork, ndirection, ngroup, invalid_device);
    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess) goto done;
  }
  {
    int invalid = 0;
    if (cudaMemcpy(&invalid, invalid_device, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
    if (invalid != 0) {
      status = 2;
      goto done;
    }
  }

  {
    const int threads = 256;
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_multigroup_upwind_kernel<<<blocks, threads>>>(state_device, transport_device,
        direction_device, neighbor_device, nowned, nwork, ndirection, ngroup, cdt_over_dx);
    if (cudaGetLastError() != cudaSuccess) goto done;
    snrt_multigroup_absorb_kernel<<<blocks, threads>>>(transport_device, state_device,
        tau_device, nowned, nwork, ndirection, ngroup);
    if (cudaGetLastError() != cudaSuccess) goto done;
  }
  {
    const int threads = 256;
    const int blocks = (nowned + threads - 1) / threads;
    snrt_cap_multigroup_species_dust_absorption_kernel<<<blocks, threads>>>(
        transport_device, state_device, tau_species_device, tau_dust_device,
        available_species_device, absorbed_hhe_species_device,
        absorbed_dust_group_device, returned_group_device, raw_group_device,
        absorbed_group_device, absorbed_device, nowned, nwork, ndirection, ngroup);
    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess) goto done;
  }

  if (cudaMemcpy(state_host, transport_device, state_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(available_species_host, available_species_device, inventory_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_hhe_species_host, absorbed_hhe_species_device, species_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_dust_group_host, absorbed_dust_group_device, group_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(returned_group_host, returned_group_device, group_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(raw_group_host, raw_group_device, group_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_group_host, absorbed_group_device, group_bytes,
                 cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  if (cudaMemcpy(absorbed_host, absorbed_device,
                 static_cast<size_t>(nowned) * sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) goto done;
  status = 0;

done:
  cudaFree(invalid_device);
  cudaFree(absorbed_device);
  cudaFree(absorbed_group_device);
  cudaFree(raw_group_device);
  cudaFree(returned_group_device);
  cudaFree(absorbed_dust_group_device);
  cudaFree(absorbed_hhe_species_device);
  cudaFree(available_species_device);
  cudaFree(tau_dust_device);
  cudaFree(tau_species_device);
  cudaFree(tau_device);
  cudaFree(neighbor_device);
  cudaFree(direction_device);
  cudaFree(transport_device);
  cudaFree(state_device);
  return status;
}
