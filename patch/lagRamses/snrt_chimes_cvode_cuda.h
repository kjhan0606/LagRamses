#ifndef SNRT_CHIMES_CVODE_CUDA_H
#define SNRT_CHIMES_CVODE_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Call once before any OpenMP cell work. No CUDA initialization on CPU mode.
 * Environment config: SNRT_CHIMES_CVODE_BACKEND=cpu|cuda_batched_lu.
 * GPU mode is experimental: only dense LU setup is offloaded. */
int snrt_chimes_cvode_configure(void);
void snrt_chimes_cvode_counts(unsigned long long out[6]);
void snrt_chimes_cvode_cuda_finalize(void); /* after all cell calls join */
#ifdef __cplusplus
}
#endif
#endif
