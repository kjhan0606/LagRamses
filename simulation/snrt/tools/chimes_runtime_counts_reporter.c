#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>

typedef void (*counts_fn)(unsigned long long *);
typedef unsigned long long (*mask_fn)(void);

__attribute__((destructor)) static void report_chimes_runtime_counts(void)
{
    counts_fn rhs = (counts_fn)dlsym(RTLD_DEFAULT, "snrt_chimes_rhs_counts");
    mask_fn workers = (mask_fn)dlsym(RTLD_DEFAULT, "snrt_chimes_rhs_worker_mask");
    if (rhs) {
        unsigned long long c[7] = {0};
        rhs(c);
        fprintf(stderr,
                "CHIMES_RHS_COUNTS gpu_calls=%llu cpu_fallback_calls=%llu "
                "gpu_batches=%llu device_errors=%llu warm_table_uploads=%llu "
                "peak_inflight_batches=%llu streams_used=%llu\n",
                c[0], c[1], c[2], c[3], c[4], c[5], c[6]);
    }
    if (workers)
        fprintf(stderr, "CHIMES_GPU_WORKER_MASK value=%llu\n", workers());
    fflush(stderr);
}
