#include "../patch/cuRamses/amr_cuda_index.cuh"

#include <cstdio>
#include <set>

static long long expected(int grid, int child, long long ncoarse, int b, int c)
{
    const long long g0 = (long long)grid - 1;
    return ncoarse + (g0 / b) * ((long long)c * b)
         + (long long)(child - 1) * b + (g0 % b) + 1;
}

int main()
{
    const long long ncoarse = 17;
    const int child_count = 8;
    const int blocks[] = {1, 64, 128};
    for (int b : blocks) {
        const int grids[] = {1, b, b + 1, 2 * b};
        for (int grid : grids) {
            for (int child = 1; child <= child_count; child++) {
                const long long got = amr_cuda_cell_1based(
                    grid, child, ncoarse, b, child_count);
                if (got != expected(grid, child, ncoarse, b, child_count)) {
                    std::fprintf(stderr, "cell mismatch B=%d g=%d c=%d\n",
                                 b, grid, child);
                    return 1;
                }
            }
        }

        std::set<long long> seen;
        for (int grid = 1; grid <= 2 * b; grid++) {
            for (int child = 1; child <= child_count; child++) {
                if (!seen.insert(amr_cuda_cell_1based(
                        grid, child, ncoarse, b, child_count)).second) {
                    std::fprintf(stderr, "duplicate cell B=%d g=%d c=%d\n",
                                 b, grid, child);
                    return 1;
                }
            }
        }

        const int high_waters[] = {0, 1, b, b + 1, 2 * b};
        const long long ncell = ncoarse + (long long)child_count * 2 * b;
        for (int hw : high_waters) {
            const long long prefix = amr_cuda_live_cell_prefix(
                ncoarse, hw, b, child_count);
            const long long want = hw == 0 ? ncoarse : expected(
                hw, child_count, ncoarse, b, child_count);
            if (prefix != want || prefix < ncoarse || prefix > ncell) {
                std::fprintf(stderr, "prefix mismatch B=%d hw=%d\n", b, hw);
                return 1;
            }
        }
    }
    if (!amr_cuda_layout_valid(17 + 8LL * 128, 17, 128, 64, 8, 8) ||
        amr_cuda_layout_valid(17 + 8LL * 128, 17, 128, 0, 8, 8) ||
        amr_cuda_layout_valid(17 + 8LL * 128, 17, 128, 64, 4, 8) ||
        amr_cuda_layout_valid(17 + 8LL * 130, 17, 130, 64, 8, 8)) {
        std::fputs("layout validator mismatch\n", stderr);
        return 1;
    }
    std::puts("CUDA BLOCK INDEX HOST UNIT: PASS");
    return 0;
}
