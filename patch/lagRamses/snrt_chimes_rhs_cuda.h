#ifndef SNRT_CHIMES_RHS_CUDA_H
#define SNRT_CHIMES_RHS_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Experimental fused dark reaction stage. CVODE and the temperature floor
 * stay on the host. A dark compact evaluation may also return the net
 * cooling rate so the host can skip calculate_total_cooling_rate. */
int snrt_chimes_rhs_configure(void);
int snrt_chimes_rhs_configure_mode(int use_cuda);
/* Calls, fallback, batches, errors, uploads, peak concurrent batches, streams used. */
void snrt_chimes_rhs_counts(unsigned long long out[7]);
unsigned long long snrt_chimes_rhs_worker_mask(void);
/* Completed compact/full batches and returned payload bytes. */
void snrt_chimes_rhs_result_counts(unsigned long long out[3]);
void snrt_chimes_rhs_finalize(void);
void snrt_chimes_rhs_verification(double out[9]); /* maxima + non-finite counts */
void snrt_chimes_rhs_cell_begin(void);
int snrt_chimes_rhs_cell_retry(void); /* reset to CPU, return pending retry */
void snrt_chimes_rhs_cell_end(void);
/* Cooperative tile phases: callback owns one cell's private input/output.
 * No callback may publish shared simulation state or invoke MPI. CPU CVODE
 * continuations stay on their owning thread; current RHS requests are staged
 * across cells before GPU execution. A busy pool selects CPU for this tile.
 * Returns dispatcher status, not the callbacks' physical acceptance status. */
typedef void (*snrt_chimes_tile_cell_fn)(int cell,void *context);
int snrt_chimes_rhs_tile_run(int cells,snrt_chimes_tile_cell_fn callback,void *context);
int snrt_chimes_rhs_tile_active(void);
#ifdef SNRT_CHIMES_RHS_TESTING
int snrt_chimes_rhs_tile_gpu_selected(void); /* current tile owns a GPU stream lease */
#endif
/* Level queue: stable callback context/outputs live until all workers join.
 * Selected workers broker current-state RHS events; others use native CPU
 * chemistry. CPU-only uses the same level queue with all workers on CPU. */
int snrt_chimes_rhs_level_requested(void);
int snrt_chimes_rhs_level_run(int cells,snrt_chimes_tile_cell_fn callback,void *context);
int snrt_chimes_rhs_cpu_direct_active(void);
#ifdef SNRT_CHIMES_RHS_TESTING
void snrt_chimes_rhs_test_drop_after(int callbacks);
void snrt_chimes_rhs_test_submit_fail_after(int submissions);
unsigned long long snrt_chimes_rhs_test_drops(void);
/* A pending whole-cell retry must reject further RHS calls before reading
 * the discarded trajectory's UserData. Test-only; preserves thread state. */
int snrt_chimes_rhs_test_retry_sticky(void);
#endif
#ifdef __cplusplus
}
#endif
#endif
