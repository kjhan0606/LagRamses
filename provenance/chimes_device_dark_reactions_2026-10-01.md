# Device-resident dark reaction tables and live CVODE RHS slice

Scope: only `/gpfs/kjhan/chimes-cvode-gpu-20261001`, branch
`agent/chimes-cvode-gpu-20261001`. Original LRD_JWST sources and the other
agent's worktree were not modified. No commit/push. Prior LU-only code,
logs, executables, patches and report are retained. This phase keeps LU on
CPU and evaluates a separate experimental dark-reaction GPU path.

## Cost and implemented boundary

H100 job 409533 used the existing `chimes_cost_probe.c` once on the existing
24-cell fixture, serial, without its optional cache experiment. It measured:

| Component | Aggregate CPU seconds |
|---|---:|
| Rate coefficients | 0.011047960 |
| Reaction flux evaluation | 0.007583994 |
| Species creation/destruction accumulation | 0.016164827 |
| Cooling | 0.018115049 |
| Complete RHS (inclusive) | 0.057692140 |
| Complete network | 0.076322114 |

There were 6110 RHS calls, including 2760 numerical-Jacobian RHS calls,
24 independent CVODE integrations, 2556 accepted steps and 144 error-test
failures. Component timings include their initial non-RHS calls too, so they
are not exclusive parts of the RHS total. Do not compare this instrumented
serial measurement directly to uninstrumented OMP wall times.

The new device kernel fuses all eight dark reaction families: temperature
dependent, constant, case-A/B recombination, grain recombination, primary/
secondary cosmic rays, H2 dust formation, H2 collisional dissociation and
CO cosmic rays. It evaluates coefficients, reaction fluxes and species
creation/destruction. Immutable table arrays and a per-species adjacency
list persist in per-stream device storage. Duplicate reactants/products are
retained; each species accumulates reactions in the original family/order,
without floating-point atomics. Each callback uploads its current temperature,
abundances and controls, and downloads rates and species source terms.

This is connected to the **actual `f()`** used by CHIMES CVODE. CPU cooling
consumes the returned CR and H2 rates, H2 dissociation coefficients and
critical densities. GPU success invalidates the owning thread's exact-T CPU
coefficient cache before any later CPU reuse. Vector unpacking, thermal floor/
root handling, cooling, constraints and CVODE control remain host code.

This is **not a complete device RHS/cooling solver**. Per-cell evolving state
still crosses host/device on every GPU RHS call. Only tables/topology and
scratch persist. No SUNDIALS rebuild was needed: original CVODE 5.8.0,
FP64, **64-bit indices**, serial vectors and CPU dense linear solver remain.
There is no cuSOLVERSP batchQR or mixed index ABI.

## Ownership and fallback

`SNRT_CHIMES_RHS_BACKEND=cpu` (or unset) is the default.
`cuda_dark_reactions` is an explicit **experimental** opt-in, built with
`CHIMES=1 USE_CUDA=1 CHIMES_CUDA=1` and the new private CHIMES hook patch.
Unsupported spectra, noncanonical species mapping, isothermal calls or
small OMP teams select CPU. Direct calls outside the bridge cell scope use CPU.

The initial batching approach allowed individual RHS calls to fall back to
CPU. Even tiny CPU/GPU arithmetic differences can contaminate a finite-
difference Jacobian. The implementation now chooses one arithmetic backend
per cell integration. After GPU selection, even 1–3-cell tail batches stay
on GPU. If a lease/device operation becomes unavailable mid-cell, return a
solver failure, discard that private trial, then restart the entire bridge
call from its original input using CPU. Tolerances and conservation thresholds
are unchanged. Different cells retain independent adaptive steps and CVODE
memory; no shared step/order/error norm is introduced.

One coordinator batches pending callbacks within each MPI rank. It uses the
existing stream-pool device selection and lease, and synchronizes DMA before
releasing ownership. No inter-rank gathering is introduced. Existing Fortran
trial commit and MPI collective rejection code is unchanged. Full MPI
transaction validation has not been run.

## Tests and limitations found

Job 409542 caught a real coupling omission: H2 dissociation coefficients and
critical densities had not been returned to CPU cooling. Thermal errors
exceeded 1e-6. The output contract and cache invalidation were corrected.

Job 409545 performed one-off paired evaluations at identical CVODE trial
states, **returning the CPU results to the integrator**. Maximum relative
differences were reaction flux 1.4216940188718989e-13, species creation/
destruction 2.7841308072272412e-14 and net cooling 6.9607533768443757e-16.
That diagnostic pass is not a GPU-integrated trajectory or speedup result.

With device results actually driving integration, the original 1e10 s
fixture matched the 22 mutually successful cells but did **not** match status
for the two original CPU-rejected cells (indices 5 and 17): CPU returned 41,
GPU succeeded. Backend locking did not remove this difference (409550), so
mixed backend arithmetic is **not demonstrated as its cause**. The precise
origin of this rounding-sensitive convergence/acceptance behavior is open.

An explicitly shortened 1e9 s interval (409547/409551) accepted all 24 cells
in both modes, but a trace species exceeded the pre-existing species parity
threshold. Temperature/energy errors were below 9e-14. The threshold was not
relaxed and the tests correctly failed. Shortening the interval is a
diagnostic, not a replacement for the original failed validation.

Before the final diagnostic print/fault-injection additions, job 409550
measured CPU 0.013262514 s and warm GPU 0.117082594 s: **0.11327x**, about
8.83x slower. It used 6346 GPU RHS calls in the warm run, zero CPU fallback,
zero device errors, one initial table upload and zero warm reuploads. Thus
this is real CHIMES callback work, not a detached kernel benchmark. Cold GPU
was 0.124335467 s and pool initialization separately 0.346676667 s.

The small number of simultaneously active cells (OMP8) entails thousands
of queue/launch/transfer/synchronization round trips. Keeping immutable tables
resident alone does not amortize them. CPU cooling remains a substantial
measured host cost. No performance promotion or production readiness claim
is justified.

## Reproduction and next unblock

Final H100 jobs **409553 and 409554 both failed physics qualification**;
the failure exit statuses are intentionally preserved. Job 409553 used the
original 1e10 s interval and measured CPU OMP8 0.0132506927 s, GPU cold
0.124833099 s and GPU warm 0.118023864 s: **0.112271301x**, or 8.91x slower.
Pool initialization was a separate 0.378087437 s. Cold includes first table
upload, warm has no table upload; neither excludes per-callback transfers.
The warm trajectory made 6346 GPU callbacks, with no CPU fallback or CUDA
errors. These are small-fixture measurements, not a production speed claim.

The 22 mutually successful cells passed the existing species threshold;
maximum temperature/thermal-energy relative errors were 1.22e-15/1.55e-15,
nuclei relative error 2.22e-16 and charge absolute residual 1.86e-20.
Cells 5 and 17 still had CPU status 41 versus GPU status 0, so the whole
fixture failed status parity. Busy-lease and serial fallback reproduced CPU
outputs and statuses exactly. Test-only mid-cell GPU failure injection
triggered eight whole-cell CPU retries, also reproducing CPU exactly.
This verifies restart ownership in the native bridge fixture, not a full
multi-rank Fortran transaction integration run.

Job 409554, the diagnostic 1e9 s interval, had all 24 statuses equal to zero
but failed species parity. Worst was cell 5, species 148 (**CO**, confirmed
in private CHIMES `src/chimes_proto.h`): CPU 2.5296233320840222e-14 versus
GPU 1.6530810069166787e-14, absolute difference 8.7654232516734349e-15.
Temperature/energy errors remained below 9e-14. The CO discrepancy is large
relative to that trace abundance; it is not waived on account of small
absolute abundance or thermal agreement. There is no warm-run performance
result for this failed shorter-interval comparison.

Final disposition: retain the experimental slice and evidence, keep CPU
default, and do not promote either this path or the earlier LU-only path.
Neither complete device-resident RHS+cooling nor scientific trajectory
equivalence is achieved. The next correctness unblock is to isolate CO
formation/destruction along the discrepant cell trajectory and identify the
CPU status-41 acceptance failure before expanding the GPU cooling boundary.
The principal performance limitation exposed here is thousands of small
callback round trips with at most eight active host cells; exact queue,
transfer and kernel contributions have not been separately measured.

```
sbatch --export=ALL,SNRT_CHIMES_TEST_DT_FACTOR=1 simulation/snrt/tools/chimes_rhs_gpu_test.sbatch
sbatch --export=ALL,SNRT_CHIMES_TEST_DT_FACTOR=0.1 simulation/snrt/tools/chimes_rhs_gpu_test.sbatch
```

Both use H100, one GPU, eight CPUs, 16 GiB and a 120 s runtime timeout.
Logs and builds are `.cvode-gpu/rhs-JOB.log` and `.cvode-gpu/rhs-build-JOB/`.
The existing fixture is reused; `SNRT_CHIMES_RHS_VERIFY=1` enables only the
diagnostic CPU-result-returning paired mode. Never use its timings as GPU
trajectory performance. The original LU test script still runs with this
additional backend disabled by default.

To reconstruct the private dependency from CHIMES a58e5c0, apply the existing
`chimes_cvode_gpu_baseline.patch`, then `chimes_cvode_gpu_factory.patch`, then
new `chimes_dark_rhs_factory.patch`. The first patch is inherited work.

Next concrete work is to diagnose the short-interval trace-species difference
and the two long-interval status differences without changing tolerances;
then port the full cooling calculation and its callback contract. A useful
larger batch is the finite-difference Jacobian's independent perturbed states
within **each cell's** CVODE, potentially combined across active cells. That
requires explicit matching of CVODE perturbations/state and complete RHS,
not pooling adaptive controllers. It is not implemented here. No new integrator,
global time-step sharing, tolerance relaxation or 32-bit ABI substitution
is proposed as an implicit workaround.

Changed this phase: new `snrt_chimes_rhs_cuda.{cu,h}`, CHIMES RHS/cache hook
patch, bridge cell-scope/config and CPU restart, Makefile linkage, extensions
to the existing native fixture, one H100 runner, and this report. No raw
simulation outputs were produced; compact evidence and executables remain.
