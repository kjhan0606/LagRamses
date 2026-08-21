#ifndef LAGRAMSES_AMR_CUDA_INDEX_CUH
#define LAGRAMSES_AMR_CUDA_INDEX_CUH

// [RESIZABLE] Keep device-side cell addressing identical to
// amr_index::icell_of.  Grid and child indices are 1-based.
#ifdef __CUDACC__
#define AMR_CUDA_HD __host__ __device__
#define AMR_CUDA_INLINE __forceinline__
#else
#define AMR_CUDA_HD
#define AMR_CUDA_INLINE inline
#endif

static AMR_CUDA_HD AMR_CUDA_INLINE long long amr_cuda_cell_1based(
    int igrid, int ichild, long long ncoarse, int block_size, int child_count)
{
    const long long grid0 = (long long)igrid - 1;
    return ncoarse
         + (grid0 / block_size) * ((long long)child_count * block_size)
         + (long long)(ichild - 1) * block_size
         + (grid0 % block_size) + 1;
}

// The last child of the high-water grid is the highest live address.  The
// resulting prefix is contiguous and can be copied component by component.
static AMR_CUDA_HD AMR_CUDA_INLINE long long amr_cuda_live_cell_prefix(
    long long ncoarse, int grid_high_water, int block_size, int child_count)
{
    if (grid_high_water <= 0) return ncoarse;
    return amr_cuda_cell_1based(
        grid_high_water, child_count, ncoarse, block_size, child_count);
}

static AMR_CUDA_HD AMR_CUDA_INLINE bool amr_cuda_layout_valid(
    long long ncell, long long ncoarse, int ngridmax,
    int block_size, int child_count, int expected_child_count)
{
    return block_size > 0 && child_count == expected_child_count &&
           ngridmax > 0 && ngridmax % block_size == 0 &&
           ncell == ncoarse + (long long)child_count * ngridmax;
}

#undef AMR_CUDA_INLINE
#undef AMR_CUDA_HD

#endif
