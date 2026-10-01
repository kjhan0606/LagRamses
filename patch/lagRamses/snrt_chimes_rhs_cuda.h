#ifndef SNRT_CHIMES_RHS_CUDA_H
#define SNRT_CHIMES_RHS_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Experimental fused dark reaction stage. Cooling and CVODE stay on host. */
int snrt_chimes_rhs_configure(void);
void snrt_chimes_rhs_counts(unsigned long long out[5]);
void snrt_chimes_rhs_finalize(void);
void snrt_chimes_rhs_verification(double out[3]); /* one-off diagnostic maxima */
void snrt_chimes_rhs_cell_begin(void);
int snrt_chimes_rhs_cell_retry(void); /* reset to CPU, return pending retry */
void snrt_chimes_rhs_cell_end(void);
#ifdef SNRT_CHIMES_RHS_TESTING
void snrt_chimes_rhs_test_drop_after(int callbacks);
unsigned long long snrt_chimes_rhs_test_drops(void);
#endif
#ifdef __cplusplus
}
#endif
#endif
