#ifndef SNRT_CHIMES_CVODE_CUDA_H
#define SNRT_CHIMES_CVODE_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Call once before any OpenMP cell work. No CUDA initialization on CPU mode.
 * Environment config: SNRT_CHIMES_CVODE_BACKEND=cpu|cuda_batched_lu|cuda_batched_lu_solve.
 * GPU modes are experimental; adaptive control and RHS remain on the host.
 * cuda_batched_lu_solve is a standalone fixture experiment, not selected by
 * SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated. */
int snrt_chimes_cvode_configure(void);
int snrt_chimes_cvode_configure_mode(int use_cuda);
/* CPU/GPU factors, batches, busy leases, errors, largest batch, peak overlap, streams used. */
void snrt_chimes_cvode_counts(unsigned long long out[10]);
void snrt_chimes_cvode_cuda_finalize(void); /* after all cell calls join */
#ifdef __cplusplus
}
#endif
#endif
