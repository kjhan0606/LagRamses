// Tensor-Core angular-block contraction for the lagRamses S_N RT backend.
// Inputs are logically row-major.  cuBLAS sees their transposes as column-major
// matrices and evaluates B^T * A^T, whose memory layout is the desired A * B.
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

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
                                               float cdt_over_dx) {
  if (state_host == nullptr || direction_host == nullptr || neighbor_host == nullptr ||
      optical_depth_host == nullptr || neutral_hydrogen_host == nullptr ||
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
  float *state_device = nullptr;
  float *transport_device = nullptr;
  float *direction_device = nullptr;
  int *neighbor_device = nullptr;
  float *tau_device = nullptr;
  float *neutral_device = nullptr;
  float *absorbed_device = nullptr;
  float *absorbed_group_device = nullptr;
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
  if (cudaMalloc(&neutral_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_device, cell_bytes) != cudaSuccess) goto done;
  if (cudaMalloc(&absorbed_group_device, group_bytes) != cudaSuccess) goto done;
  if (cudaMemcpy(state_device, state_host, state_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(direction_device, direction_host, direction_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neighbor_device, neighbor_host, neighbor_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(tau_device, optical_depth_host, group_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;
  if (cudaMemcpy(neutral_device, neutral_hydrogen_host, cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) goto done;

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
    const int blocks = static_cast<int>((total + threads - 1) / threads);
    snrt_cap_multigroup_absorption_kernel<<<blocks, threads>>>(transport_device, state_device,
        absorbed_device, neutral_device, nowned, nwork, ndirection, ngroup);
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
  status = 0;

done:
  if (event_stop != nullptr) cudaEventDestroy(event_stop);
  if (event_start != nullptr) cudaEventDestroy(event_start);
  cudaFree(absorbed_group_device);
  cudaFree(absorbed_device);
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
      ncell, ncell, ndirection, ngroup, cdt_over_dx);
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
      nowned, nwork, ndirection, ngroup, cdt_over_dx);
}
