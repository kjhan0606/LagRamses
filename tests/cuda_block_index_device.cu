#include "../patch/cuRamses/amr_cuda_index.cuh"

#include <cstdio>
#include <cuda_runtime.h>

__global__ static void check_index(int* errors)
{
    const long long ncoarse = 17;
    const int child_count = 8;
    const int blocks[] = {1, 64, 128};
    for (int ib = 0; ib < 3; ib++) {
        const int b = blocks[ib];
        const int grids[] = {1, b, b + 1, 2 * b};
        for (int ig = 0; ig < 4; ig++) {
            const int grid = grids[ig];
            const long long g0 = (long long)grid - 1;
            for (int child = 1; child <= child_count; child++) {
                const long long want = ncoarse
                    + (g0 / b) * ((long long)child_count * b)
                    + (long long)(child - 1) * b + (g0 % b) + 1;
                if (amr_cuda_cell_1based(
                        grid, child, ncoarse, b, child_count) != want)
                    (*errors)++;
            }
        }
        const int high_waters[] = {0, 1, b, b + 1, 2 * b};
        const long long ncell = ncoarse + (long long)child_count * 2 * b;
        for (int ih = 0; ih < 5; ih++) {
            const int hw = high_waters[ih];
            const long long prefix = amr_cuda_live_cell_prefix(
                ncoarse, hw, b, child_count);
            const long long g0 = (long long)hw - 1;
            const long long want = hw == 0 ? ncoarse : ncoarse
                + (g0 / b) * ((long long)child_count * b)
                + (long long)(child_count - 1) * b + (g0 % b) + 1;
            if (prefix != want || prefix < ncoarse || prefix > ncell)
                (*errors)++;
        }
    }
    if (!amr_cuda_layout_valid(17 + 8LL * 128, 17, 128, 64, 8, 8) ||
        amr_cuda_layout_valid(17 + 8LL * 128, 17, 128, 0, 8, 8) ||
        amr_cuda_layout_valid(17 + 8LL * 128, 17, 128, 64, 4, 8) ||
        amr_cuda_layout_valid(17 + 8LL * 130, 17, 130, 64, 8, 8))
        (*errors)++;
}

int main()
{
    int* errors = nullptr;
    if (cudaMallocManaged(&errors, sizeof(*errors)) != cudaSuccess) return 2;
    *errors = 0;
    check_index<<<1, 1>>>(errors);
    const cudaError_t err = cudaDeviceSynchronize();
    const int count = *errors;
    cudaFree(errors);
    if (err != cudaSuccess || count != 0) {
        std::fprintf(stderr, "CUDA BLOCK INDEX DEVICE: FAIL errors=%d cuda=%s\n",
                     count, cudaGetErrorString(err));
        return 1;
    }
    std::puts("CUDA BLOCK INDEX DEVICE: PASS");
    return 0;
}
