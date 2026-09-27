# M1--M5 live moment CPU/CUDA dispatch efficiency, 2026-09-27

Scope: lagRamses `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). This is
execution-placement work; no moment closure, transport equation, quadrature,
conservation rule, or material/source physics was changed. The production
default moment order remains M3. Reference-control source/dust inputs are
not scientifically calibrated.

## Wiring corrections

- Native moment basis, material interface and M1 stress extraction already
  accepted 4-component M1, but the RAMSES namelist and live payload guard
  rejected it. Both guards now admit M1--M5; `mkrun.py` and the namelist
  generator offer the same set. The CUDA closure and face adapter now admit
  M1's four components, preserving explicit CUDA failure for unsupported
  nondefault quadratures.
- OpenMP closure work previously used one 256-cell chunk per worker; a
  256-cell CPU call could occupy only one thread. All CPU arms now claim
  16-cell chunks. Hybrid CUDA retains 256-cell chunks but sends calls/tails
  under 64 cells to CPU. This threshold is a conservative A100 observation,
  not a universal optimal-tuning claim.
- For the current live 32-direction, at-most-256-face transaction, hybrid
  mode keeps face flux on CPU and leaves the shared device streams for
  nonlinear closure work. Material projection likewise uses per-cell OpenMP
  in hybrid/OpenMP mode. Explicit `SNRT_BACKEND=cuda` still exercises both
  CUDA kernels for comparison and admission checks.

## Same-input native evidence

A100-SXM4-80GB, four CPU workers, FP64, 16x24 angular basis. Jobs 406468,
406469, 406471 and 406478 used one CUDA stream, the unchanged
`--fmad=false` kernel build and same native input per order. The test
artifacts are under `.mn-efficiency-20260927.xeDokg/`. Closure was an
interior, weakly anisotropic positive state; timings include transfer,
allocation and synchronization. They are short single measurements, not
throughput confidence intervals or coupled-run speedups.

| Order | 256-cell closure CPU4 / GPU (ms) | 256-face CPU4 / GPU (ms) | 256-cell projection CPU4 / GPU (ms) |
|---|---:|---:|---:|
| M1 | 10.95 / 1.05 | 0.049 / 2.128 | 0.165 / 0.437 |
| M2 | 18.22 / 1.17 | 0.056 / 1.339 | 0.185 / 0.442 |
| M3 | 35.57 / 1.81 | 0.070 / 1.243 | 0.227 / 0.435 |
| M4 | 46.60 / 3.66 | 0.092 / 1.191 | 0.338 / 0.445 |
| M5 | 80.81 / 5.88 | 0.066 / 1.276 | 0.504 / 0.450 |

The complete M1--M5 CUDA closure/project/face probes passed. Maximum
relative CPU/GPU closure-state difference in the probe was 1.77e-15; face
and project differences stayed below 2e-13 absolute. Job 406471 exercised
the actual hybrid dispatcher: at 256 nonzero cells, all five orders sent
closure to the GPU, while all 256 faces stayed on CPU; its M1 first CUDA
closure included startup/JIT overhead (18.3 vs 16.5 ms CPU), so a one-call
M1 speedup is not claimed. Its 16-cell auto call stayed on CPU. The
unavailable-stream policy test passes with four OpenMP workers and all five
orders, including isotropic stress trace.

`snrt_moment_smoke.f90` intentionally remains M2--M5: its crossing-beam
anisotropic-pressure assertion and fifth harmonic do not apply to M1.
M1 has a lower-order closure and cannot retain the same zero-flux
crossing-beam anisotropy as M2--M5.

## Coupled-run boundary

The first M1 launch (406476) was rejected because its inherited environment
fixed `SNRT_RT_LEVEL=2`; the current IR contract requires all levels. The
second launch (406480) reached the M1 payload check, but the exploratory
binary's `snrt_moment_ramses.o` timestamp predates the `nm>=4` source edit
by about two seconds, so it still contained the former `nm>=9` guard. These
are launch/build-consistency failures, not evidence of M1 physics failure.
Both logs are retained. A clean final binary with the post-edit payload
object has SHA256
`2006c79ed79f62e27d6901d0fc371f0ece6e8aa0a893a248c985d83ba860b2f7`;
its source hashes, CUDA symbols, dynamic libraries and M1 parser string
were checked. The final-build Slurm script exited after successful link
during its trailing verification pipeline; the binary was independently
verified, and the check was made pipe-safe for reproducibility.

Job 406485 completed all four M1/M5 OpenMP-vs-hybrid, two-step 4^3
technical control runs on the same A100 node. Effective inputs live in
`.mn-efficiency-20260927.xeDokg/coupled-final-{m1,m5}-{openmp,hybrid}/effective.nml`.
`noutput=1`, unreachable noncosmological `aout=2`, `tout=1e30`,
`foutput=fbackup=1000000` produced zero full dumps; available GPFS free
space was about 271 TiB at launch. Both orders admitted the nonzero source
in step 2 and completed with source/absorption/energy ledgers. Within each
order, CPU/hybrid source terms match; maximum xHII differed by less than
2e-12 absolute. These are reference-control inputs, not a calibrated
physical comparison.

| Order | Active-step `SNRT_MN_COMMIT` OpenMP / hybrid (s) | Total program OpenMP / hybrid (s) | Hybrid closure CPU / GPU cells |
|---|---:|---:|---:|
| M1 | 65.188 / 69.881 | 68.346 / 73.018 | 5,872 / 18,192 |
| M5 | 63.954 / 89.361 | 67.014 / 93.197 | 4,976 / 19,088 |

The first coupled timing pair does **not** establish end-to-end GPU gain.
The `cooling` timer includes M_N material/CHIMES work and cannot isolate
transport. Job 406503 repeated M5 on the same node in reverse order with
the existing `SNRT_PERF_DIAGNOSTIC` enabled, with no new timing code:

| Mode | Active commit (s) | Material preparation (s) | CHIMES photo sum (two-worker CPU-s) | Commit minus material tile (s) |
|---|---:|---:|---:|---:|
| hybrid | 97.818 | 94.399 | 183.856 | 2.969 |
| OpenMP | 112.356 | 108.033 | 210.722 | 3.781 |

The large run-to-run timing spread is in the photo-chemistry preparation,
not the CUDA closure or face kernel. This pair does not identify whether
adaptive chemistry effort or shared-node CPU variability caused the spread.
The remainder column also contains source/ledger overhead and is not a
pure RT-kernel measurement.

## Warm-start correction and verification

The CPU live closure reused its prior dual coefficients but the CUDA live
adapter silently discarded them. The benchmark's weakly anisotropic fresh
states did not reproduce the repeated SSPRK workload. The CUDA adapter now
accepts the same host dual guess, retains the solver's cold retry and
unchanged convergence criterion, and returns the accepted cache as before.
The pre-existing context GPU warm-target path is unchanged. The Makefile
link guard checks the new CUDA symbol. The first build (406507) stopped at
the old guard; after correcting it, 406511 linked the updated binary with
SHA256 `ec5b1082cc13f6b59e4a9a2fc281491684e6b36de8d8adf30baf389928eab864`.

Native A100 job 406506 passed M1--M5 cold/warm CUDA-vs-CPU parity. For the
same 256-cell M5 batch the first CUDA closure took 6.162 ms and its warmed
second call 1.715 ms; M4 was 3.707 / 1.703 ms. The M1 first call included
CUDA startup (17.550 ms) and its warmed call was 1.533 ms. Across orders
the maximum relative moment difference remained 1.77e-15. These are short
single measurements, not full-simulation speedups. Job 406513 additionally
compared one and four CPU workers on the same 256-cell input. The four-worker
speedup was 1.92, 1.97, 1.89, 1.98 and 1.99 for M1--M5, respectively,
with identical moments. Repetition 406520 returned 1.9--2.8x, showing
timing variability. This is useful but not a four-core scaling claim;
the test did not independently establish physical-core placement.

Job 406520 also measured a like-for-like warmed *second* call on 256
identical cells (same A100 node, four CPU workers, shared CUDA stream):

| Order | Warm CPU / warm hybrid-CUDA (ms) |
|---|---:|
| M1 | 2.375 / 0.693 |
| M2 | 2.792 / 1.495 |
| M3 | 3.629 / 1.591 |
| M4 | 5.004 / 1.711 |
| M5 | 3.565 / 1.755 |

All five retained CPU/CUDA relative state agreement better than 2e-10.
These isolated, weakly anisotropic batches support the CUDA placement choice
at 256 cells, not an end-to-end speedup prediction for a stiff source case.

The updated binary passed a source-active two-step M5 hybrid coupled run
(job 406512) using the same inputs, with `SNRT_PERF_DIAGNOSTIC=1`. Its
active-step commit was 64.511 s, including 61.535 s of material preparation
and 119.147 two-worker CPU-s in CHIMES photo integration. Subtracting the
material tile leaves about 2.586 s for transport, source and ledger work,
versus 2.969 s in the cold-hybrid diagnostic pair. The source and absorbed
photon terms matched the cold-hybrid run to roughly 1e-10 relative; maximum
xHII differed by about 1e-12 absolute, and the run completed without a
full raw dump. This is a correctness and bounded-work check, not a proven
whole-run speedup; the CHIMES phase alone varied by many tens of seconds
between otherwise comparable runs.

Conclusion: M1--M5 now have working OpenMP work sharing, CUDA closure/face
parity and a live CUDA warm-start path; hybrid sends the measured cheap face
and projection calls to CPU. Native CUDA closure throughput is favorable on
A100 for 64/256-cell batches, but production-scale end-to-end speedup and
hardware-specific break-even thresholds remain unmeasured. Do not label the
reference-control 4^3 case a production performance qualification.
