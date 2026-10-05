# CHIMES GPU batch/overlap qualification

Date: 2026-10-05
Slurm job: `413581` (`syn101`, A100-SXM4-80GB)
Work directory: `/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/chimes-batch-overlap-413581`

## Scope

Matched 128^3 real-input CPU and GPU-RHS runs, changing only the GPU RHS batch
capacity (64 vs. 128 cells). Configuration: 2 coarse steps, MPI=1, OpenMP=8,
four GPU brokers/streams. Numerical parity was checked over every HDF5 dataset
against the CPU result with a `1e-8` peak-relative limit. Separate 20-second
Nsight windows checked whether kernels from multiple streams actually overlap.

## Results

| Run | Wall time (s) | Speedup vs CPU | GPU RHS batches |
|---|---:|---:|---:|
| CPU | 1943.531 | 1.0000 | 0 |
| GPU, batch 64 | 1951.941 | 0.9957 | 5,934,028 |
| GPU, batch 128 | 1994.955 | 0.9742 | 4,343,411 |

Both GPU runs matched all 1,395 HDF5 datasets: zero nonfinite pairs, zero
discrete mismatches, and maximum field peak-relative difference `8.33e-11`.
The 128-cell batch reduced the number of launches/batches, but was slower in
wall time than batch 64 and CPU-only. Therefore retain 64 as the default; do
not infer that fewer batches alone improve end-to-end performance.

Nsight confirmed real cross-stream overlap in both profiles. In the sampled
windows, four streams and up to four simultaneous kernels were observed:

| Batch | Kernel instances | Peak overlapping kernels | Integrated overlap |
|---:|---:|---:|---:|
| 64 | 66,457 | 4 | 234,881,166 ns |
| 128 | 37,987 | 4 | 111,766,620 ns |

Thus asynchronous multi-stream execution is present, but it does not make this
GPU-RHS path faster end-to-end. The run also reports `cpu_setups` around 23.5M
and zero GPU CVODE factorizations: CVODE integration remains on CPU, which is
the next targeted work item (#2).

Raw simulation snapshots were removed after parity evaluation. Logs, hashes,
field and wall comparisons, cleanup manifest, and Nsight reports remain in the
work directory above. This commit records the qualification; the experimental
batch selector and one-off Slurm driver remain in the mixed CHIMES worktree and
are intentionally not included here.
