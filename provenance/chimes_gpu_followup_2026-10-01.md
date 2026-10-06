# CHIMES GPU advisory follow-up — 2026-10-01

Workspace: `/gpfs/kjhan/chimes-cvode-gpu-20261001`, branch
`agent/chimes-cvode-gpu-20261001`, parent `9bec08c`. This is a bounded
native bridge diagnostic, not a production simulation or a performance
promotion. The default/production CHIMES route and its acceptance thresholds
were not changed.

## Paired RHS check

The existing one-state CPU/GPU RHS comparison was extended to (a) explicitly
count non-finite CPU/device values and invalid device outputs, (b) compare net
species derivatives and CO's net derivative, and (c) run with either the CPU
or GPU RHS returned to CVODE. The two return modes therefore sample their own
integrator trajectories. The diagnostic reports status counts separately and
does not turn known parity failures into a passing physics gate.

H100 job **409593** completed on `syn09` in 19 s. It ran the fixture at the
original `dt=1e10 s` and `dt=1e9 s`, in both trajectory modes. All reported
non-finite and invalid-device-output counts were zero.

| Timestep | Returned trajectory | CPU rejects | GPU rejects | Status mismatch | Maximum net-species RHS relative difference | CO net RHS absolute / relative difference |
|---|---|---:|---:|---:|---:|---:|
| 1e10 s | CPU | 2 | 2 | 0 | 6.22e-13 | 9.40e-38 / 1.14e-13 |
| 1e10 s | GPU | 2 | 0 | 2 | 1.01e-12 | 1.41e-37 / 1.98e-13 |
| 1e9 s | CPU | 0 | 0 | 0 | 2.67e-13 | 1.76e-38 / 1.42e-13 |
| 1e9 s | GPU | 0 | 0 | 0 | 1.56e-13 | 1.76e-38 / 1.42e-13 |

These maxima aggregate paired callback states; they are not per-cell maxima.
The status divergence remains specifically cells 5 and 17. At `1e9 s`, the
GPU endpoint CO abundance is `1.6530810069e-14` versus CPU
`2.5296233321e-14` (absolute difference `8.7654e-15`); the pre-existing
species parity criterion fails despite temperature/energy agreement below
`9e-14`. No criterion was relaxed. Paired-run wall times include the added
second RHS evaluation and are not performance measurements. The unpaired
409553 result remains the relevant previous speed sample: warm GPU was about
8.9 times slower than CPU OMP8.

## CPU status-41 behavior

H100 job **409599** completed on `syn09` in 18 s. A compile- and
environment-gated CHIMES diagnostic was enabled only for the fixture's
`nH=0.11 cm^-3`, `T=1e4 K`, `H2=0.05`, `dt>=1e10 s` input. CVODE reported:

```
At t = 9.70709e+08 and h = 98.6282, the error test failed repeatedly or with |h| = hmin.
```

Thus the generic CHIMES “recoverable RHS error” warning is misleading for
this case: the observed rejection is CVODE error-test/minimum-step failure,
not a non-finite RHS or a failed elemental/charge acceptance check. At
`dt=1e9 s`, both CPU and GPU status return zero.

The bridge retries up to three times for post-solve conservation rejection,
tightening tolerances only in that case. A nonzero `chimes_network` status
breaks immediately and maps to bridge status 41. CHIMES has an internal
recovery loop for repeated RHS-function errors, but it rebuilds CVODE toward
the same requested endpoint; it does not halve the interval. The RAMSES
callers fail closed: RT transactions roll back/reject chemistry, and the
dark material path aborts the MPI transaction. No state is silently
committed, but there is no production shortened-interval retry.

Do not add automatic subcycling or loosen tolerances based only on this
synthetic fixture. A solver/retry change must preserve the full requested
interval and the coupled photon/thermal/dust accounting, and should first be
validated against the exact failure state.

## Existing adaptive-SPGMR candidate on the exact fixture

The isolated candidate source was rebuilt against the active SUNDIALS 5.8
prefix. Correction: `sunlinsolspgmr.so.3.8.0` is also a library SONAME in
SUNDIALS 5.8; it does not establish a SUNDIALS 3.8 release or incompatibility.
H100 job **409605** exposed a test-link mismatch (the CPU-only candidate's older CHIMES
updater has no GPU cache-invalidation hook); a test-only no-op hook was
compiled only into that CPU-profile executable. Job **409618** then completed
the paired dense/candidate profile with all GPU selectors unset.

At `dt=1e10 s`, dense CVODE rejected cells 5 and 17 with status 41, whereas
the adaptive candidate returned status 0 for all 24 cells. This shows the
candidate can traverse this particular state, but does not by itself establish
an acceptable replacement solution. At `dt=1e9 s`, both solvers accepted all
cells, yet for cells 5/17 dense produced `CO=2.5296233320840222e-14` and
SPGMR produced `1.7614871338403292e-14` (absolute difference
`7.6813619824e-15`, about 30.4% of the dense result). The existing fixture's
per-species parity scale is `1e-17 + 1e-6*abs(reference)`; the CO difference
alone therefore gives a scaled discrepancy of **766.198**, against the
unchanged pass limit of 1. Temperature and H2 agree closely, but this CO
counterexample is sufficient to reject the candidate as a drop-in solver at
present. No acceptance threshold was changed.

Disposition: keep dense CVODE as the production path, do not promote this
adaptive-SPGMR candidate, and do not relax species parity. CPU chemistry
optimization remains worthwhile given its measured runtime share, but the
next candidate must preserve trace-species results under the existing gate as
well as address status 41. The GPU RHS remains experimental: same-state
arithmetic parity was close, but its measured execution was about 8.9x slower
than CPU OMP8.

## Cost relevance and disposition

An already-recorded bounded source-active gas/dust integration measured
CHIMES at `13.598 s` within `16.906 s` of SNRT coupling (`80.43%`). This is a
small 4^3, four-step control, not a 128^3 cosmological fraction. It does show
that CHIMES is worth optimizing; the low-fraction reason to shelve chemistry
does not apply. However, the current RHS GPU implementation is still about
8.9x slower on its native fixture and covers only reaction evaluation, so
this does **not** justify GPU promotion. Keep both GPU paths experimental and
opt-in.

The CPU robustness issue now has a more specific description: dense BDF
error control stalls at a small internal step for two synthetic low-density
molecular states, while same-state CPU/GPU RHS arithmetic remains close and
finite. The rebuilt SPGMR comparison above completes the proposed bounded
follow-up: it traverses the failed fixture but violates endpoint CO parity.
No new production solver selection is approved by this diagnostic.

Evidence: `.cvode-gpu/rhs-409593.log` and `.cvode-gpu/rhs-409599.log`.

## CPU RHS optimization follow-up

Implementation/testing worktree HEAD was `ad26b255010c1864b79d70611a0eeb3d9fbcd31f`.
H100 job 409629 waited on `QOSMaxGRESPerUser` and was cancelled before running.
The CPU-only comparison ran on the previously authorized LagEunha server,
with GCC 13.2, the same SUNDIALS 5.8 prefix, 1/8 OpenMP threads, and no GPU
initialization. No simulation dumps were produced.

The reusable `tools/build_chimes_cpu_library.sh` builds a fresh library from
an already patched CHIMES tree. `reference` preserves the arithmetic;
`lto` enables cross-function compiler optimization; `inactive` applies two
small RHS/cooling patches. All modes use double precision, no fast math,
and disabled floating-point contraction. No solver tolerances are changed.

LTO's five-sample median was 0.518185 -> 0.503226 s at OMP1 and
0.0820036 -> 0.0801186 s at OMP8: only 2.89%/2.30% reduction. Interpolation
and index helpers were already inlined in CHIMES headers, so this is not a
major new optimization. LTO also prevents reliable LD_PRELOAD interception
of internal calls; it remains an optional experiment.

The `inactive` implementation omits temperature-dependent reaction coefficient
interpolation when a reactant has both an exactly zero elemental budget and
exactly zero abundance. It similarly skips the 1-D cooling coefficient of a
coolant with that property. This applies only with zero radiation spectra.
The corresponding physical contribution is identically zero; positive
budgets, including arbitrarily small ones, retain the original calculation.
Hydrogen and active ions at temporarily zero abundance are not removed.
Only these intermediate coefficients become zero; the network, reaction
ordering, density-dependent cooling, and dense CVODE solver remain intact.

Each timing is 20 repeats of the existing 24-cell, `dt=1e9 s` native bridge
fixture. Five alternating reference/candidate samples gave:

| Threads | Reference median (s) | Inactive median (s) | Reduction |
|---|---:|---:|---:|
| 1 | 0.517193329 | 0.467300834 | 9.65% |
| 8 | 0.0821329728 | 0.0744738542 | 9.33% |

A second run after adding mixed-composition validation gave 9.14%/11.94%
reductions on the same timing fixture; OMP8 absolute times varied between
runs. These are native chemistry timings for a C/O-only mixture, not a
cosmological runtime forecast or a measured gain for fully enriched gas.

All 157 output abundances, temperature and status codes were hashed without
structure padding and matched bit-for-bit between builds and thread counts.
The accepted fixture hash was
`9c8b307e640e4987a24807ccf4caf5a631e0674ad415a47f64917d08520697da`.
The original `dt=1e10 s` fixture still rejected the same two cells with
status 41 and unchanged returned states; this optimization does not fix that
solver failure. A 24-cell mixture of C/O-only, all-metal, metal-free and
tiny-positive (`1e-40`) metal budgets also matched exactly, with zero rejected
cells at OMP1/OMP8 (hash
`e4bd958624151e4016098bdc0b17ca501a77baa140aa1374b80907acb4d4bcc2`).

Evidence: `.cvode-gpu/cpu-lto-lageunha-r2.log`,
`.cvode-gpu/cpu-inactive-lageunha.log`, and
`.cvode-gpu/cpu-inactive-lageunha-mixed.log`. The last run's libraries and
per-sample logs are under `.cvode-gpu/cpu-build-lageunha-20261001-inactive-mixed/`.
No running simulation library or global default was replaced. The real-cell
comparison below supersedes the synthetic timing as evidence for adoption.

### Real-cell replay and disposition

The existing replay driver (SHA256
`0feef185d4155d28865e6320e8ba85ec71115f54aa5ea66fb2ce9ed13bbeae1c`)
replayed captured live cell 18447 from
`.snrt-performance.jDx9Iz/live-cell-20260925/live-403698.log`. Its temperature
is 158.854 K, nH is 0.165361 cm^-3, and dt is 8.307e11 s. Five paired runs
matched the complete printed 158-value state exactly. Median wall times were
11.1830 ms reference and 11.1370 ms `inactive`, a negligible 0.41% change.
Every elemental budget in this captured cell is positive, including the
tiny metal budgets, so the absent-element optimization does not apply.

A second `zero` experiment skipped coefficients for exactly zero current
reactants instead. It disables temperature-coefficient cache reuse whenever
any coefficient was skipped, forcing recomputation if a later RHS/Jacobian
perturbation activates that species. Isothermal initialization and the
coefficient directly consumed by H2 heating are preserved. Cooling coefficients
are recomputed on each call. The same native tests and mixed compositions
matched exactly. However, five real-cell pairs gave 11.1380 ms reference and
11.4861 ms candidate (3.13% slower). Neither skip variant is recommended for
production adoption on this evidence. Both remain isolated build experiments.

### Reducing finite-difference Jacobian work

The existing serial cost probe can now vary CVODE's maximum Jacobian age
using `SNRT_COST_JACOBIAN_AGE`. It does not change tolerances or convergence
failure handling. This is an LD_PRELOAD experiment, not a runtime default.
On the same captured cell:

| Maximum age | Jacobian builds | Linear RHS calls | Total RHS calls | BDF steps | Error / convergence failures |
|---|---:|---:|---:|---:|---:|
| default | 4 | 632 | 838 | 187 | 1 / 0 |
| 100 | 2 | 316 | 522 | 187 | 1 / 0 |
| 200 | 1 | 158 | 364 | 187 | 1 / 0 |

Instrumented network times were 11.317, 9.009 and 7.458 ms, respectively;
these include profiling overhead and are not production speed estimates.
Against the reference, maximum species-scaled differences were 9.78e-11
and 3.62e-10 (limit 1), and temperature differences were 1.12e-15 relative.
Thus this cell supports reducing Jacobian construction cost.

However, the existing 24-cell `dt=1e10 s` comparison rejected global age 200:
the worst CO scaled discrepancy was 1413.27 at cell 1 (1000 K). An optional
cold/low-density restriction (`T<=200 K`, `nH<=0.2`, zero spectra) was also
tested, rather than adopting the global setting. It still failed at cell 3:
CO changed from zero to 4.1280578e-17, giving a scaled discrepancy of 4.128
for both age 100 and 200. Status-41 failures remained the same two cells.
These counterexamples reject Jacobian-age tuning for production under the
existing criteria; no threshold was relaxed or favorable state selected as
the acceptance result.

The fixture now supports loading/saving complete states to apply its existing
species/temperature/energy/conservation comparisons between library builds;
hash equality remains useful for arithmetic-preserving changes. Evidence is
under `.cvode-gpu/cpu-build-lageunha-20261001-zero/`: `replay-*.log`,
`jac-age-*.log`, `age-200-1.log`, `age-cold-1.log`, and `age-cold100-1.log`.

This bounded optimization round is complete without a production promotion.
The useful next algorithmic target is reducing Jacobian construction cost
while retaining its required refreshes (e.g. conservative analytic derivatives
or dependency-aware differentiation), rather than increasing its age or
repeating cache microbenchmarks. The active thermal equation couples chemical
abundances through temperature, so sparsity must include that coupling; a
reaction-only graph must not be assumed to describe the full Jacobian.

### Parallel dense-Jacobian implementation and bounded evaluation

Implemented a CPU column-parallel finite-difference Jacobian experiment in
`simulation/snrt/tests/fixtures/phase0/chimes_parallel_jacobian.c`. It preserves
the SUNDIALS 5.8 increment formula, constraints, dense linear solver, refresh
policy and tolerances, and does not approximate thermal dependencies. This is
an explicitly version-pinned LD_PRELOAD experiment using private SUNDIALS
headers, NOT a production API or a new runtime default.

Each worker owns its abundance vector, species rates, rate buffers, gas state
and a copy of the bridge's cell context. The latter matters because the hybrid
cooling callback writes `outside_molecular_domain`. Test-only bridge helpers
are guarded by `SNRT_CHIMES_JACOBIAN_TESTING`; normal builds are unchanged.
Any RHS failure or domain flag causes a reference serial reconstruction, so
diagnostic side effects and failure handling are not silently dropped. On
successful non-verifying calls the owner's last-column RHS is repeated to
preserve upstream mutable state. This adds one RHS evaluation; it does not
claim to reduce the total number of evaluations.

Active outer OpenMP teams retain upstream serial Jacobian construction, avoiding
nested oversubscription and preserving the existing cell-parallel path. Spectra,
gas-owned opaque hybrid state, and GPU backend selectors fall back as well.
Zero-spectra dark systems retain their original RT/case-B configuration; the
test does NOT equate N_spectra=0 with rt_update_flux=0.

The final bounded run on authorized LagEunha is:
`.cvode-gpu/parallel-jacobian-lageunha-20261001-r6/`, with parent log
`.cvode-gpu/parallel-jacobian-lageunha-20261001-r6.log`.
Reproduce with `bash simulation/snrt/tools/chimes_parallel_jacobian_test.sh NEW_DIR`
on LagEunha (not the login node). The runner links the existing native fixtures,
uses a new build directory, and does not submit a RAMSES simulation or emit dumps.

Validation results:

- The existing default/mixed 24-cell fixtures at dt=1e10 and 1e9 s passed
  direct bitwise Jacobian comparisons at 1, 2, 4 and 8 workers. Final species,
  temperature and rejection statuses matched; the two known status-41 failures
  at default dt=1e10 remain failures, not claimed fixes.
- The original Fortran captured-cell replay also passed bitwise Jacobian and
  complete endpoint-state comparison at all four worker counts. It retains
  the original group initialization, case-B switches and physical timestep.
- Every warm replay retained all four Jacobian refreshes; all 100 repetitions
  within each process reproduced the same full endpoint state. Worker counts
  were measured, not inferred solely from OMP environment settings.

Five alternating-order runs, each with one cold solve followed by 100 warm
solves, gave the following median of per-run mean warm times:

| Inner Jacobian workers | Warm cell time (ms) | Reduction vs reference |
|---|---:|---:|
| Reference serial | 10.7548 | — |
| 2 | 9.8485 | 8.43% |
| 4 | 10.0761 | 6.31% |
| 8 | 10.3858 | 3.43% |

These are wall times for one captured cell, not CPU-time savings or full
simulation throughput. More workers are not faster. The synthetic short-step
fixture likewise showed no substantial benefit; its existing outer-cell OMP8
path remained about 0.082 s per 20x24-cell batch. Verification runs include
extra serial reconstruction and are excluded from timing claims.

Decision: retain this correct parallel prototype as experimental evidence;
do not promote it, do not spend GPU resources on a direct port, and do not add
a namelist selector. The best measured 8.4% wall-time reduction consumes extra
cores and does not establish a throughput gain over cell-level parallelism.
The substantive unresolved target remains fewer expensive operations during
Jacobian construction, including all temperature/electron-density coupling,
not more workers or less frequent refreshes. No conservation/error tolerance
was relaxed and no physical model was changed.

Runner fixes during development: r1 rejected an invalid fixture label; r2
incorrectly preloaded the adapter into `timeout` instead of the test process.
r4 tried a zero-spectrum C++ initialization for an original RT-configured
capture and correctly hit the existing chat/c admission check. That shortcut
was removed; r5/r6 use the established native Fortran replay and its original
inputs. These failed harness attempts are not chemistry failures or performance
measurements. No test raw simulation outputs were generated or removed.

### Dependency-exact Jacobian work reduction

Implemented `simulation/snrt/tests/fixtures/phase0/chimes_jacobian_cooling_reuse.c`
and reused the existing native test runner, without adding a new benchmark
suite, physics switch, or accuracy tolerance. Two complementary reductions
operate only inside one dense difference-quotient Jacobian:

1. Reuse the live cooling-coefficient arrays only when temperature, electron,
   HI and HII number densities and the owning buffer are exactly unchanged.
   These are the evolving inputs of the upstream 1D/2D/4D cooling coefficient
   evaluator. Actual cooling/heating summation remains executed with every
   perturbed abundance; this does NOT replace NLTE evolution with a T-only table.
2. Build a reaction-to-species contribution list for the admitted dark network.
   After rates are updated normally, compare their exact stored values and
   recompute only creation/destruction sums touched by a changed reaction.
   Preserve the original per-species addition order, duplicate reactants and
   inactive destinations. Never update a sum by subtracting an old contribution
   and adding a new one, which would change rounding. All temperature and
   electron-density dependencies are retained through the original rate update.

All scratch state is thread-local and freed/reset at each Jacobian boundary.
The adapter remains pinned to SUNDIALS 5.8 private ABI and is experimental.
GPU-selected or illuminated paths remain upstream. Unknown secondary hooks
disable reaction-row reuse; a test-only bridge capability check recognizes
the current dark photoelectron hook as a no-op. Nothing is enabled in normal
production builds, and the standalone adapter must not be stacked with the
other Jacobian/cost interposers.

Evidence on LagEunha:

- `.cvode-gpu/jacobian-cooling-reuse-20261001/` (coefficient reuse only).
- `.cvode-gpu/jacobian-reaction-reuse-20261001-r2/` (combined, final hook guard).
- Corresponding `.log` files beside those directories ended with
  `JACOBIAN_COOLING_REUSE_TEST_PASS`.
- Reproduce with `bash simulation/snrt/tools/chimes_parallel_jacobian_test.sh
  NEW_DIR reaction_reuse` (or `cooling_reuse`) on authorized LagEunha.

Both variants passed direct bitwise matrix comparison and the existing
default/mixed-composition endpoint comparisons at dt=1e10 and 1e9 s, including
outer-cell OMP1/OMP8. Known default long-interval status-41 failures remained
unchanged. The original native captured-cell replay retained all four Jacobian
refreshes and exactly the same final temperature and 157 species, including
100 warm repetitions per process.

For this captured cell, cooling coefficient evaluations during its four
Jacobians fell from 632 to 28 (604 exact reuses, 95.57%). Of 98,596 subsequent
species-row updates, 92,660 (93.98%) reused unchanged sums. This reduces RHS
internal work, NOT the number of finite-difference columns or RHS callbacks.

Five alternating-order runs with 100 warm solves each gave median per-run mean
times of 10.7902 -> 10.1530 ms for coefficient reuse alone (5.91%), and
10.6590 -> 9.74085 ms for combined reuse (8.61%). Unlike column parallelism,
these improvements use the same single-core resource. The synthetic fixture's
cell-parallel timing also remained comparable or modestly improved; there is
no demonstrated full-simulation speedup yet.

Decision: correctness of this bounded implementation is established on the
existing comparisons, but the gain is not a major acceleration and it is not
promoted to production through a private-ABI preload. Previous captured-cell
profiling attributes about 3.83 ms to dense LU/setup and solves alone; reducing
only coefficient and row assembly cannot remove that cost. The remaining
algorithmic target is a structured Jacobian/linear solve, with its thermal
coupling and trace-species endpoints explicitly preserved. Do not represent
this exact-work-reuse prototype as an analytic or colored sparse Jacobian.
No simulation dumps were produced or deleted and no jobs remain running.

### Linear-solve follow-up: exact zeros and optimized direct LU

Implemented and evaluated the remaining bounded linear-algebra candidates on
LagEunha, using the existing native fixtures and captured production cell:

- `chimes_sparse_triangular.c`: compress exactly-zero entries after the
  existing LU, first as column scatter and then as row gather. Row gathering
  preserves the original subtraction order for each solution component.
  Nonfinite factors/RHS use the upstream solve; overflow during a trial solve
  also re-evaluates upstream from the saved RHS. No nonzero entry is dropped.
- The same adapter optionally performs partial-pivot LU with exactly-zero
  multipliers skipped. It retains pivot tie-breaking and arithmetic order;
  dense updates remain for nonfinite operands and densely populated columns.
  The density switch is a computational cost choice, not a drop tolerance.
- `chimes_lapack_dense.c`: compare oneMKL 2025.3 LAPACK direct LU/solve with
  one inner thread and LP64 interface. It keeps the existing BDF tolerances,
  Jacobian construction and refreshes, but may change floating-point rounding.

The runner accepts `sparse_solve`, `sparse_lu`, and `lapack` variants. These
are test interposers, not production namelist options. All use the previous
exact cooling/reaction reuse as part of the combined candidate.

All native default/mixed dt=1e10/1e9 comparisons and OMP1/OMP8 endpoint tests
passed under the existing criteria. The known two status-41 failures were
preserved. For the captured cell, the exact-zero implementation matched all
29 LU factorizations (numeric factor entries and pivot vectors), all 203
triangular solves, the four Jacobians and the complete final state. Zero signs
are not required to match in the LU/solve numerical checks; the captured
endpoint does match bitwise. The triangular solve operated on 852,375 entries
instead of 5,035,618 (83.07% fewer multiply/subtract terms). These checks do not
constitute general-purpose solver certification outside the tested CHIMES path.

LAPACK also passed; its largest existing fixture species-scaled discrepancy
was 2.55e-8 against a limit of 1, with temperature discrepancy 6.66e-16 relative.
The captured-cell endpoint matched exactly. Thus its rejection below is on
performance grounds, not a failed physics criterion.

Five alternating-order runs of 100 warm solves each gave these median mean
cell times; each row has its own paired upstream reference:

| Combined candidate | Reference (ms) | Candidate (ms) |
|---|---:|---:|
| Exact reuse + column-sparse triangular solve | 10.7245 | 10.0510 |
| Exact reuse + row-sparse triangular solve | 10.7100 | 9.87046 |
| Exact reuse + zero-aware LU + row-sparse solve | 10.8047 | 9.74415 |
| Exact reuse + oneMKL direct LU/solve | 10.7033 | 12.4013 |

An additional direct comparison of the already implemented exact reuse against
the best sparse-LU combination, using the same native replay executable and
alternating order, gave medians of **9.74227 versus 9.72097 ms** (only 0.22%
additional reduction). Every endpoint was identical. This separates the
previous reuse improvement from the effect of this turn's linear algebra.

Decision: do not promote these linear-algebra alternatives or infer a large
simulation speedup from their lower arithmetic counts. Sparse indexing,
packing and pattern construction offset the saved operations for these small
matrices. LAPACK was about 15.9% slower than its upstream reference on the
captured cell and also slower on the synthetic cell batches. Existing dense
production behavior remains unchanged. A substantial acceleration is still
unresolved; these results do not establish that every possible structured or
analytic solver has been exhausted.

Evidence directories under `.cvode-gpu/`: `sparse-triangular-20261001`,
`sparse-triangular-row-20261001`, `sparse-lu-20261001`, `lapack-dense-20261001`,
and `linear-paired-20261001`; parent `.log` files exist for the first four.
No full simulation or GPU job was launched, no raw dumps were generated or
removed, and all invoked processes finished. No tolerances, physical models,
production switches, shared context, commits or remotes were changed.

## Approved semi-analytic Jacobian experiment

Implemented `chimes_semianalytic_jacobian.c` as a bounded SUNDIALS-5.8
test interposer, not a production replacement. Atomic-metal chemical columns
use mass-action derivatives (including repeated reactants and three-body
density factors) only when the actual upstream DQ perturbation leaves the
computed temperature exactly unchanged. H/He, electrons, molecules, energy,
and temperature-changing columns retain numerical differentiation. This is
therefore a partial semi-analytic Jacobian, not an implementation of all
thermal/compositional coefficient derivatives. The full energy row remains
numerical, with fresh rates and full cooling/heating hooks. No species or
small nonzero rates are discarded. Unsupported paths fall back to upstream.

The adapter also checks identity species mapping and excludes columns which
appear in exceptional dust-H2, CO-CR or secondary-CR table dependencies.
Cooling interpolation reuse is scoped to one Jacobian and exact temperature,
electron, HI, HII densities and coefficient-buffer identity. All cooling sums
are still evaluated. No rounded cache keys or increased Jacobian age.

Tests on Lageunha, default/mixed 24-cell sets at the two existing timesteps,
passed the existing endpoint criteria. Largest species-scaled difference was
0.350906 (limit 1); maximum temperature relative difference was 5.33e-14.
The two pre-existing status-41 failures at the longer default timestep were
preserved, not fixed or counted as successful integrations. OMP 1/8 checks
were included. The real captured cell retained exactly the final temperature;
maximum species-scaled endpoint difference was 4.08e-39. The initial runner
stopped at a bitwise native-state comparison: analytic and DQ matrices need
not give bitwise equal endpoints. For these two explicitly semi-analytic modes
only, native endpoint comparison now uses the existing C++ fixture bounds
(species 1e-17 + 1e-6*abs(reference), temperature relative 1e-6). Solver and
conservation controls were not changed; older exact optimizations still
require byte-for-byte native state equality. Native replay's own admission
checks remain enabled.

For the real cell's four Jacobians, 520 of 632 columns used analytic chemistry,
112 used full numerical RHS, plus four baseline refresh calls: 116 full RHS
evaluations instead of 632 (81.65% reduction). Cooling/heating was still
evaluated for the analytic columns, so this is not an 81.65% work reduction.
Finite Jacobian differences were recorded diagnostically, not asserted zero;
the real-cell maximum absolute analytic-versus-DQ matrix difference was
1.99e-9. This dimensional maximum is not a relative accuracy certificate.

Five alternating-order samples of 100 warm native solves gave:

| Variant | Paired reference (ms/cell) | Candidate (ms/cell) |
|---|---:|---:|
| Initial partial semi-analytic Jacobian | 10.6817 | 11.9338 |
| With exact cooling interpolation reuse and zero-product skips | 10.7126 | 11.1772 |

An existing cost probe, reused for 101 native solves per arm, located the
offsetting cost: reaction-vector time fell 0.12925 -> 0.04981 s and cooling
0.23477 -> 0.16752 s, but LU rose 0.19046 -> 0.37080 s, with the same 2,929 LU
calls, 20,503 triangular solves, 18,887 integration steps and 404 Jacobians.
Triangular solve time also rose 0.18701 -> 0.21022 s. These profiled timings
include instrumentation and are diagnostic, not the unprofiled benchmark.
The LU increase is measured; a subnormal-arithmetic explanation is not yet
proved. Full RHS counts across those solves fell 84,638 -> 32,522.

Evidence: `.cvode-gpu/semianalytic-20261001` (initial bitwise comparison stop),
`semianalytic-20261001-r2` and `semianalytic-20261001-r3` (complete tests).
The r3 directory also contains `profile-reference.log` and
`profile-candidate.log`. No production promotion is justified by these
results alone.

The already available exact-zero LU/triangular implementation was then paired
with the semi-analytic candidate (`semianalytic_sparse` runner mode). It passed
the same endpoint checks and, on the captured cell, numerically matched all
29 upstream LU factorizations/pivots and all 203 triangular solves of the
candidate matrices. Median paired native times were 10.7628 ms reference and
11.4188 ms candidate: still slower. Evidence directory:
`.cvode-gpu/semianalytic-sparse-20261001`. No new sparse algorithm was added.

Disposition: retain this partial semi-analytic implementation as experimental
evidence only; do not enable it in production. The best measured version of
this turn is still about 4.34% slower than its own paired upstream baseline,
despite substantially fewer full RHS calls. A general analytic Jacobian,
including temperature/energy coupling and an effective linear solve, remains
unresolved; these measurements do not reject that broader approach. Do not
claim GPU/whole-simulation acceleration from this CPU experiment. All invoked
tests completed on Lageunha, with no simulation jobs/dumps or production
switch changes. No commit/push was performed in this turn.

## Requested remaining thermal-coupling and linear-solver tests

Added opt-in `SNRT_CHIMES_JAC_THERMAL_CHAIN` to the experimental adapter.
For atomic-metal columns this combines fixed-temperature analytic chemistry
with the ideal-gas derivative dT/dx = -T*nH/n_total at fixed energy. A single
numerical energy perturbation supplies the common temperature response of
both chemistry and cooling. Fixed-T energy dependence still uses the full
cooling evaluator. Temperature-clamped columns retain the earlier admission
rule. H/He/electron/molecular columns still use full DQ. This is a chain-rule
hybrid, NOT completion of a fully analytic thermal/chemical Jacobian.

The first experiment used the rounded temperature change under the species
perturbation; the second differentiates the EOS itself. Both failed the
existing default-composition, dt=1e9 s endpoint comparison. For the EOS version,
cell 5 CO was 1.6530074148e-14 versus upstream 2.5296233321e-14 (abundance per H),
scaled difference 874.404 against the unchanged limit 1. Temperature/energy,
charge and nuclei checks passed, and all statuses were zero at this timestep.
The longer dt=1e10 s comparison preserved the two previously known status-41
failures. No failed chemistry comparison was converted to a pass.

Compared the same thermal candidate with upstream dense LU, the previous
exact-zero sparse LU/solve, and sequential oneMKL LAPACK. All three produced
the same failing CO endpoint at displayed precision. Each full runner stopped
on that accuracy check, before its performance phase; no speedup claim is made
for these variants. This isolates the observed discrepancy from the choice
among these three linear algebra implementations, not from all possible
linear-solver effects.

To distinguish candidate error from an unconverged reference, the existing
cost probe gained test-only wrappers for public CVodeSS/SVtolerances. An
explicit scale in (0,1] tightens both CVODE relative and absolute tolerances;
the vector is cloned, not mutated. With the variable absent the setters pass
through unchanged. Production settings, explicit-step criteria, physical
admission and endpoint comparison bounds remain unchanged. Separate 24-cell
runs at dt=1e9 s yielded this cell-5 result:

| CVODE tolerance multiplier | Upstream CO/H | Thermal-chain CO/H |
|---|---:|---:|
| 1 | 2.5296233321e-14 | 1.6530074148e-14 |
| 0.1 | 1.7453169803e-14 | 1.7453158857e-14 |
| 0.01 | 1.7617316789e-14 | status 41 (rejected) |
| 0.001 | 1.7614076903e-14 | 1.7612977373e-14 |

The 0.01 candidate run rejected two cells consistently at OMP 1 and 8;
rejected states were not used as physical endpoints. At 0.1 and 0.001 both
versions accepted all 24 cells. Pairwise maximum species-scaled differences
were respectively 0.0322833 and 0.109760, within the existing bound of 1.
The default reference CO value is therefore tolerance-sensitive, not an
established physical truth. The tightening experiment supports a numerical
convergence issue but does not by itself prove the root cause or establish
global convergence across all species/states. Do not replace production with
the candidate or globally tighten tolerances based on this small fixture.

A single native replay of the captured real cell, using the thermal-chain
candidate at the unchanged production tolerances, passed native admission.
Its final temperature matched the reference exactly; largest species-scaled
difference was 4.08e-39. It used 120 Jacobian RHS evaluations instead of 632,
including four extra temperature-response calls. This one-cell result does
not override the mixed-state accuracy/stability failure.

Evidence under `.cvode-gpu/`: `thermal-chain-20261001` (first rounded-response
experiment), `thermal-eos-semianalytic_thermal-20261001`, the corresponding
`_sparse` and `_lapack` directories, and `thermal-convergence-20261001`.
The convergence directory retains all input-state results/logs, diagnostic
library, pairwise comparisons and `native-thermal.log`. Tests ran on Lageunha;
no GPU/full-simulation job, dump, tolerance change in production, commit or
push was performed. Result: the requested candidates were exercised, but
general thermal-Jacobian readiness and a net acceleration remain unresolved.

## Trace-species failure recovery bundle

Identified the previously ambiguous failure as **CV_ERR_FAILURE (-3)**,
not a diagnosed RHS callback failure: the thermal-chain run at tolerance
multiplier 0.01 stopped at t=69,549,048.3535123 s of 1e9 s, current step
12.305239 s, with CO (species 148) the largest weighted local-error component
(62.446). CHIMES had printed the same misleading RHS warning for every
negative CVODE result. Both maintained receiver baseline patches now print
the actual flag, reached time and requested time. This diagnostic change
does not alter solver acceptance.

CO's CR reaction contains sqrt(max(x_CO,0)), a steep/nonsmooth near-zero
dependency. It is a plausible contributor, not an isolated proof that this
reaction alone caused every failure. Scoped experiments, using the existing
cost probe and fixtures, rejected two tempting shortcuts:

- Nonnegative CVODE constraints cured the first case but introduced
  CV_CONSTR_FAIL (-15) in mixed-composition long steps (cells 9/21).
- A local CO DQ perturbation cap cured default-tolerance failures, but did not
  establish trace accuracy and failed again with tighter molecular tolerances.

Neither shortcut is in the selected native-library patch. No reaction was
removed, floored, changed, or solved in equilibrium; no abundance clipping was
added. Diagnostic switches remain confined to test interposers.

Implemented the selected candidate in `data/chimes_trace_solver.patch`,
applied only by the explicit `trace` mode of `tools/build_chimes_cpu_library.sh`:

1. For dark molecular thermal solves, tighten molecular **absolute** tolerances
   by 1e-3 only. Atomic species, thermal-energy/photon absolute tolerances,
   relative tolerance, explicit acceptance and conservation limits are unchanged.
2. On CV_ERR_FAILURE only, restart CVODE multistep history from its returned
   last accepted y,t. CVodeReInit preserves the configured tolerances, user data,
   root function/direction and stop time in the pinned SUNDIALS 5.8 source.
   No chemistry state or target time is changed. Maximum four restarts and
   stop immediately on no progress; other errors propagate unchanged.

The prototype at 1e-2 molecular tolerance passed integrations but did not pass
the tighter endpoint comparison (CO scaled error 3.69). It is NOT the selected
patch. The 1e-3 version passed a further factor-ten tightening comparison.
Recovery examples in the diagnostic experiment completed after one history
restart; this is not evidence that every possible future failure is recoverable.

Native-library results (no LD_PRELOAD adapter for the candidate solve):

| 24-cell configuration | Max species-scaled difference vs 10x tighter |
|---|---:|
| default, dt=1e9 s | 0.0675812 |
| default, dt=1e10 s | 0.0244940 |
| mixed, dt=1e9 s | 0.7782203 |
| mixed, dt=1e10 s | 0.0414253 |

All 96 configuration/cell cases and their tighter comparisons had zero
rejections, including the formerly rejected default cells 5/17. Existing
temperature, energy, charge/nuclei criteria passed; OMP 1/8 endpoint parity
passed. Limit on the reported species-scaled difference remains 1.
This is convergence evidence for the tested matrix, not universal solver
certification. Native Fortran dark mixed density/temperature checks also
passed `DARK_MIXED_THREAD_PARITY_CONSERVATION_PASS`.

Re-evaluated the thermal-chain Jacobian against this corrected reference:
all four configurations passed, maximum species-scaled difference 0.003078,
no status mismatch or rejections. This resolves the earlier tested parity
failure without declaring the nominal old reference to be the truth.

Five paired runs of 100 warm captured-cell solves: original library median
10.6917 ms, corrected trace library 11.8039 ms (+10.4%). This is an accuracy/
robustness repair, NOT a speedup. A separate three-sample alternating test
on the corrected library gave median times:

| Jacobian path, corrected library | ms/cell |
|---|---:|
| Default numerical | 11.9127 |
| Prior exact cooling/reaction-row reuse | 10.9894 |
| Experimental semi-analytic thermal chain | 14.6815 |

Native reuse endpoint matched exactly; thermal-chain final T matched and
species-scaled discrepancy was 4.08e-39. Thus the semi-analytic path is still
not promoted for performance. The reuse path remains an experimental private-
ABI interposer, not silently installed in production.
The corrected-library reuse adapter also passed all four 24-cell comparisons
with zero endpoint differences and bitwise Jacobian verification enabled.

Evidence directories: `.cvode-gpu/trace-failure-20261001`,
`trace-positive-20261001`, `trace-molecular-20261001`, `trace-local-co-20261001`,
`trace-local-tight-20261001`, `trace-restart-20261001`,
`trace-native-20261001` (rejected 1e-2 candidate), and
`trace-native-20261001-r2` (selected 1e-3 candidate, build/state/parity/timing logs).
Native candidate libchimes.so SHA256:
`c2c30ec08e8d4542a710a6106ff98ef5bb9601ffb0b8f15360b24ec3aaa93b63`.
The first build attempt `trace-library-20261001` failed on patch context before
compilation; the patch was corrected and successfully built in the native
directories. No failed build was used for testing.

All computation ran on Lageunha. No production executable/library was replaced,
no full simulation/GPU job or raw dump was generated, and no commit/push was
performed. Default build mode remains unchanged; the selected solver repair
is available as an explicit candidate build, with the limitations above.

## GPU requalification and cost decomposition bundle

Operator approved steps 1 and 2 together: requalify the Astra-era GPU LU and
GPU dark-reaction RHS against the corrected CPU library, then profile only
qualified paths. No device-resident integration redesign was authorized by
this bundle. Implemented additional modes in the existing native fixture and
`tools/chimes_gpu_requalify.sbatch`. No GPU physics/backend code was changed.

Slurm **409744**, partition h100, node syn09, **H100 NVL**, completed **0:0** in
37 seconds; request was one GPU, eight CPUs, 16 GiB, fifteen minutes maximum.
Build and native tests ran inside this allocation. It was not a RAMSES run:
no IC/evolution/output schedule or simulation dumps. Other users/jobs and the
existing simulation binaries were untouched. H200 was eligible but not used.

The job built the selected `trace` CHIMES library and both existing GPU
backends with FP64, SUNDIALS 5.8 / 64-bit indices and CUDA 12.8.1. Library SHA256
matched the prior CPU correction exactly:
`c2c30ec08e8d4542a710a6106ff98ef5bb9601ffb0b8f15360b24ec3aaa93b63`.
Native test executable SHA256:
`8a09a8159ec79af098dafe923403ae1c76c930292748330e5bf0b09848e3c4d6`.

### Correctness

Each backend separately ran default/mixed composition at dt=1e9 and 1e10 s,
24 cells per configuration: **96 cases per backend**, including formerly
rejected cells 5/17. The new modes require every CPU and GPU cell to return
success; matching rejection alone cannot pass. Diagnostic bypass is rejected.
All states passed the unchanged species, temperature, energy, nuclei and
charge criteria, with zero status mismatches and zero GPU errors.

- GPU LU: zero endpoint differences in all four configurations.
- GPU RHS: maximum species-scaled discrepancy 0.00565621 (limit 1), maximum
  temperature discrepancy 1.9984e-15 and energy discrepancy 2.1094e-15.
- Actual warm GPU work was required in every configuration. Warm table
  uploads were zero. LU calls still retain their normal per-call CPU fallback;
  passing is not a claim that every factorization ran on the GPU.

The previous CO disagreement did not recur with the corrected reference.
This establishes parity for these tested dark states, not illuminated chemistry,
all event-root transitions, MPI transactions or a production cosmological run.

### Unprofiled timings

Default 24-cell fixture, dt=1e9 s, OMP8 on both arms, three warm samples with
alternating CPU/GPU order. Medians below are **per 24-cell batch**, NOT per-cell
times for the separate captured cosmological cell. Initialization/table upload
and profiling overhead are excluded.

| Backend | CPU OMP8 (ms) | GPU + OMP8 (ms) | GPU/CPU time |
|---|---:|---:|---:|
| LU offload | 10.08054 | 44.11002 | 4.376x |
| Dark-reaction RHS offload | 10.04839 | 130.15344 | 12.952x |

These are shared-node bounded measurements, not isolated full-simulation
scaling results. No favorable speedup was inferred from successful correctness.

### Separate warm Nsight measurement

Nsight Systems 2024.6.2 captured one warmed GPU solve via cudaProfilerStart/Stop,
after initialization and cold qualification; CPU phases were outside capture.
CUDA and OS-runtime traces used no CPU sampling. Both profiled runs also
passed the endpoint checks. Profiling materially changes runtime, so its
component times below must not be added to or normalized by the unprofiled
timing table.

| Warm captured component | LU | RHS |
|---|---:|---:|
| Profiled batch wall time (ms) | 70.6884 | 223.9315 |
| GPU kernel launches | 161 | 1,147 |
| Total GPU kernel time (ms) | 15.1703 | 30.0240 |
| H2D copy activity (ms) | 1.3299 | 1.4517 |
| D2H copy activity (ms) | 0.9805 | 5.4097 |
| cudaStreamSynchronize host time (ms) | 16.4244 | 44.1882 |
| cudaMemcpyAsync host API time (ms) | 4.6550 | 14.5225 |
| cudaLaunchKernel host API time (ms) | 1.9504 | 10.8907 |
| Worker pthread_cond_clockwait, summed (ms) | 23.4633 | 86.0435 |
| Worker pthread_cond_wait, summed (ms) | 311.5673 | 1,167.5143 |
| Worker pthread_mutex_lock, summed (ms) | 44.5579 | 197.7703 |

The OS-runtime rows were filtered to the eight application worker thread IDs
that issued cudaLaunchKernel_v7000, excluding profiler/driver poll threads.
Timed condition waits are the existing coordinator's batch gathering; ordinary
condition waits include waiting for processing/completion. These are cumulative
**thread times**, can overlap across workers, and must not be summed as wall
time. CUDA synchronization overlaps GPU kernel/copy time too. This is a
component breakdown, not a disjoint wall-time accounting.

Kernel grids were tiny: LU 4–8 blocks (mean 5.422), 1,024 threads/block;
RHS 1–8 blocks (mean 5.915), 256 threads/block. Mean kernel durations were
94.23 us and 26.18 us respectively. This demonstrates limited work per launch;
hardware occupancy was not directly measured. Transfer-only optimization is
not enough: GPU kernel time alone in each profiled batch already exceeded
the corresponding unprofiled complete CPU batch, and host gathering/locking
plus per-batch synchronization adds substantial overhead. That comparison
motivates redesign, but profiler overhead prevents a quantitative speedup bound.

### Disposition

Steps 1 and 2 are complete on H100. Both backends remain experimental/opt-in:
keep CPU as the performance default. The next implementation target is
tile/leaf-level batching with evolving state retained on-device and a coherent
RHS/cooling/Jacobian/linear-solve path, preserving per-cell adaptive error control
and transaction failure handling. Simply increasing stream count or moving
only LU does not address the measured launch granularity. No such architectural
change, tolerance relaxation, commit, push or production deployment was made
in this bundle.

Evidence: `.cvode-gpu/requalify-409744.log` and
`.cvode-gpu/requalify-409744/`: all eight qualification logs, unprofiled
`lu-bench.log`/`rhs-bench.log`, the two profile logs, Nsight reports/SQLite
databases and `*-profile-stats.csv`. These are bounded profiling evidence,
not disposable raw simulation snapshots.

## Unified full-cell CUDA selection (2026-10-01)

Added `SNRT_CHIMES_COMPUTE_BACKEND=cpu|cuda_integrated` to the native receiver.
The new selector configures both existing CUDA components together for one
complete CHIMES cell solve: batched dense-LU setup and dark-reaction RHS. It
resolves the backend before CHIMES initialization and rejects a conflicting
legacy LU or RHS selector. With the unified selector unset, the two legacy
controls retain their existing independent behavior. CPU remains the default.

The CUDA selector still uses the existing per-cell host CVODE state/adaptive
control, CPU cooling and unsupported chemistry paths, and the tested fallback
and full-cell retry boundaries. Per-callback reaction inputs/results still
cross the host/device boundary. This is an integrated execution selector for
the existing complete cell call, not the proposed tile-resident solver. No
solver tolerance, chemistry, transaction or production default changed.

Implementation edits: `patch/lagRamses/snrt_chimes_bridge.c`,
`snrt_chimes_cvode_cuda.{cu,h}`, `snrt_chimes_rhs_cuda.{cu,h}` and the runtime
description in `simulation/snrt/NATIVE_RUNTIME.md`. The bridge and both CUDA
translation units compiled with the cluster GNU/MPI/CUDA toolchain.

### H100 end-to-end qualification

Job `409782` ran on an NVIDIA H100 NVL (SM 9.0), with one GPU, eight OpenMP
threads and the fixed CHIMES trace library/table. It exercised 24 live CHIMES
cells for default and mixed elemental compositions at `dt` factors 0.1 and 1.
All four CPU/GPU comparisons accepted every cell, had zero status mismatches
and zero CUDA errors. The largest scaled species difference was 5.66e-3;
relative temperature and energy differences were at most 2.11e-15, relative
nuclei error at most 2.12e-16, and absolute charge error at most 3.01e-20.
The dark-reaction GPU RHS ran in all cases. Batched GPU LU was used when the
batch admitted it; zero GPU LU factors in a cold mixed-composition case are
an allowed CPU fallback, not a test failure.

Three alternating warm CPU/GPU samples per case gave these median wall times:

| Fixture | CPU (s) | Integrated GPU (s) | GPU/CPU |
| --- | ---: | ---: | ---: |
| Default, `dt` 0.1 | 0.009954 | 0.186160 | 18.70x slower |
| Default, `dt` 1 | 0.016749 | 0.283169 | 16.91x slower |
| Mixed, `dt` 0.1 | 0.030837 | 0.336391 | 10.91x slower |
| Mixed, `dt` 1 | 0.037856 | 0.474188 | 12.53x slower |

Conclusion: the unified selector is functionally wired and numerically agrees
with CPU on this bounded fixture, but this measured path is not an acceleration.
Keep CPU as the default and do not recommend `cuda_integrated` for production
performance. These 24-cell results do not establish full-domain scaling; the
per-cell host CVODE state and repeated host/device RHS transfers remain likely
costs to address before a larger GPU trial. The first qualification attempt
(`409773`) stopped on overly strict fixture assertions for permitted LU CPU
fallback and a lazy table upload; those assertions were corrected without
relaxing physics checks before the passing rerun.

Evidence: `.cvode-gpu/integrated-409782.log` and the four `check-*.log` plus
four `bench-*.log` files under `.cvode-gpu/integrated-409782/`. The qualification
is bounded native solver evidence, not a RAMSES evolution or production GPU
qualification.

## OpenMP/rank-to-GPU multi-stream dispatch follow-up (2026-10-01)

To remove the single global worker bottlenecks in both CUDA adapters,
concurrent OpenMP cell requests now reserve disjoint batches before stream
acquisition. The RHS batch leader leases a shared stream and submits H2D,
reaction kernel and D2H work; the LU leader submits a same-size matrix batch
to cuBLAS. Different leaders may run on separate streams. The RHS
table/topology payload and per-stream staging/device buffers remain resident
after lazy initialization. If no RHS stream is available, a cell returns to
the original CPU RHS; if it had already committed to GPU arithmetic, it uses
the existing whole-cell CPU restart to avoid mixing implementations inside
one numerical Jacobian. LU setup independently falls back to native CPU LU.

The current driver allocates a GPU to each MPI local rank and OMP workers share
that rank's stream pool. Thus multi-GPU scaling is rank-to-device across the
node, while multiple streams overlap RHS and LU batches within each rank.
This does not make the host CVODE callback asynchronous: each caller still
waits for its own RHS outputs, and dynamic cell state/results still cross
host/device at each callback. It is a scheduling/overlap step, not a
tile-resident CHIMES integrator. The fixture records peak concurrent RHS and
LU batches plus distinct stream slots; the two-task qualification uses Slurm's
per-task GPU binding.

Local CUDA translation-unit compilation and fixture C++ syntax checks pass.

### A100 integrated qualification

Job `409844` completed on an NVIDIA A100-SXM4-80GB with zero qualification
failures. Default and mixed compositions at both timestep factors accepted all
cells, with zero status mismatches and CUDA errors. The largest species-scaled
difference was 5.66e-3 (limit 1); temperature and energy differences were at
most 2.11e-15, nuclei error 2.12e-16, and absolute charge error 3.85e-20.
The fixture observed up to three simultaneous RHS batches/streams and two LU
batches. This qualifies correctness and multistream overlap on the bounded
fixture, not a speedup: median integrated-GPU time was about 9.4–13.2 times
CPU OMP8 across the four cases.

Job `409855` is the separate two-task/two-A100 qualification; on 2026-10-02
accounting reports `TIMEOUT` on syn102, without a batch output log. It yields
no multi-GPU qualification or scaling evidence and remains an unverified
requirement. It had no scheduler dependency on `409844`. Job `408505`
was cancelled by its owner, but Slurm still reports it `COMPLETING` with
`Reason=Prolog` and its allocation listed, so node cleanup has not finished.
The earlier H100 job `409843` was cancelled before start after Slurm projected
a next-day start. CPU remains the production default; no stale RHS state or
relaxed convergence criterion was introduced.

## Full RAMSES comparison launch (2026-10-02)

The operator requested a larger actual simulation comparison because the
24-cell fixture cannot establish application performance. The isolated runner
`.cvode-gpu/chimes_real_ramses_128_compare.sbatch` builds the current working
sources and uses one executable and one effective namelist for the CPU and
`cuda_integrated` arms. It evolves a uniform 128^3 mesh for eight hydro steps,
MPI=1/OMP=8, with live two-size dust and non-equilibrium CHIMES dark chemistry.
The initial attempt incorrectly described `SNRT_RT_LEVEL=0` as disabling RT.
Source inspection and job 410481 prove it is an invalid level filter that is
ignored: the template's `radiation_transport='snrt_sn'` activates S_N transport.
This uniform template has therefore been withdrawn as evidence for dark
material performance. The corrected actual-IC profile is documented below.
Stream count is fixed at one, not tuned. CPU remains the production default.

Each arm requests one final HDF5 dump, estimated 2.93 GiB of hydro fields
before metadata. The runner retains inputs, logs, timings, backend counters
and a per-dataset difference report. Raw outputs are removed only after the
execution/comparison evaluation succeeds; failed runs retain their outputs
for diagnosis. The report measures floating differences without declaring a
new numerical acceptance tolerance or production qualification.

Launch failures are recorded rather than counted as simulations: job 409971
failed before build due to missing Intel module prerequisites; 410423 failed
in 17 seconds because parallel make compiled consumers before
`amr_parameters.mod` existed. The test build is now serial; Makefile VPATH is
unchanged. Job 410456 allocated syn102 but remained in Slurm Prolog with no
batch log or run directory and was cancelled before simulation launch.

The runner also had an omitted Python loop body in its post-run HDF5 report;
that is corrected. Both embedded Python blocks now parse, and a tiny HDF5
report check writes the expected two data rows and reports the known injected
floating difference. This verifies the reporter only, not the simulation.

Replacement job **410469** started on A100 node **syn101** at
2026-10-02 15:02:43 KST, requesting one GPU, eight CPUs, 64 GiB and 24 hours.
Its submission explicitly sets the active checkout as the working directory
and excludes syn102; A100/H100 are the allowed partitions. At launch there is
no application timing or numerical comparison result. Evidence is the
job-specific Slurm log and `.cvode-gpu/chimes-real-ramses-128-410469/`.

Job 410469 completed the full CUDA-linked build (binary SHA256
`410873e9810a29d685bd6240c6706f32b61573f6d57761e2dba4579ace703928`)
but both simulation arms exited 127 because `/usr/bin/time` is also absent on
syn101. No simulation result was produced. The runner now uses Bash's reserved
`time`, separating child stderr from wall/user/system timings and recording the
child exit status. A small failing-command check verified stderr routing and
exit-status retention. Slurm supplies aggregate peak memory separately.

Job **410475** reuses that completed build and its matched CHIMES library;
no solver source, tolerance, mesh or timestep was changed. It started on
syn101 and entered the CPU RAMSES arm at 2026-10-02 15:11:17 KST. The log
confirms creation of 262144 level-7 octs (2097152 cells) and entry into gas
initialization. This is actual application launch evidence, not yet a completed
CPU/GPU timing or accuracy comparison. The two arms share the same newly
generated effective namelist. Evidence is
`.cvode-gpu/chimes-ramses-128-410475.log` and
`.cvode-gpu/chimes-real-ramses-128-410475/`.

Job 410475 subsequently failed before chemistry: both arms hit fine-MG
nonconvergence after 100 iterations (reported residual 0.9365), aborted with
MPI code 914 (shell exit 146), and reported zero CPU CVODE setups and zero GPU
chemistry calls. Its launcher wall times therefore do not measure chemistry.
The inherited small comparison template had disabled the project-default
FFTW base-grid path. The enlarged uniform periodic comparison now explicitly
builds with `USE_FFTW=1` and sets `use_fftw=.true.` in both arms; Poisson and
hydro remain enabled, and the MG abort policy remains strict. The existing
FFTW implementation solves the same periodic discrete Poisson operator and
sets the zero Fourier mode to zero. This chooses the supported uniform-grid
reference solver rather than treating the failed MG result as acceptable.

The preload reporter now also obtains `getrusage(RUSAGE_SELF)` inside each
RAMSES rank, recording CPU user/system time (including OpenMP workers) and
peak RSS. Bash timing fields are labelled as launcher CPU times to avoid
mistaking MPI launch-process time for simulation CPU consumption.

Replacement job **410481** rebuilds with FFTW enabled; it does not reuse the
FFTW-disabled binary. Both embedded Python blocks, shell syntax and reporter
C syntax passed checks before submission. Actual FFTW execution, chemistry
calls, completion, numerical parity and speed remain to be measured.

### Corrected actual-IC dark material comparison

Job 410481 successfully built with FFTW and logged its direct-solve marker;
the earlier MG startup failure was avoided using the existing uniform-grid
reference. Its subsequent `SNRT RT enabled`/80-direction log revealed that the
uniform generator's `SNRT_RT_LEVEL=0` was not an RT-off switch. The job was
cancelled for that confirmed workload mismatch; it supplies no dark-path
performance result. The completed FFTW binary SHA256 is
`376dda634becbc32c9c889a7ad42056f74888d63a7b97d125cef8973dbad95e0`.

The runner now copies the already admitted actual-IC dark profile from
`/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/dark-profile/run.nml`, using
the registered pre-enriched 12.5-cMpc 128^3 gas/DM IC. It replaces the legacy
`rt=.false.` with `radiation_transport='none'`, enables FFTW, sets one CUDA
stream explicitly, and uses two complete coarse steps and one final dump per
arm. This tests spatially varying real cosmological cell states and the
independent CMB/dust/CHIMES path, rather than the uniform non-cosmological
SNRT comparison. Two steps bound resource cost before scaling a GPU campaign;
this is a controlled performance/endpoint comparison, not galaxy calibration.

Job **410506** requests one A100/H100, MPI1/OMP8, 64 GiB, 24 hours, excludes
syn102, and reuses the completed FFTW binary and matching CHIMES library.
The effective namelist resides in the job's `setup/chimes_real128.nml` and
is copied unchanged into both arm directories. The generator dry run and
both embedded Python syntax checks passed; deprecated RT environment enable
and level-filter variables are explicitly unset. Successful arms must show
both the FFTW marker and the dark-material entry marker, and must not show
`SNRT RT enabled`. Per-arm HDF5 output is estimated at 3–4 GiB including DM
particles/metadata (6–8 GiB total); `/gpfs` had 211 TiB free before submission.
Raw outputs remain retained until their comparison succeeds.

Job 410506 exited at namelist admission (9 seconds), before reading ICs or
chemistry: the reused non-cosmological comparison binary had `NENER=1`, while
`read_hydro_params.f90` requires `NENER=0` for the admitted cosmological
kind7 C/silicate DL01/D03/CHIMES path. The strict check was correct and is
unchanged. The runner now builds `NENER=0` (the Makefile still reserves
`NVAR=187` for CHIMES/live dust); replacement job **410509** performs that
fresh build with FFTW enabled and uses the actual-IC inputs above. Unsupported
CR/SGS/Fe/PAH/drift/sublimation and sinks remain disabled in that profile.
No comparison or production acceleration claim follows from these failed
launches. Slurm resources remain one allowed GPU/eight CPUs/64 GiB/24 hours.

Job 410509 completed the correct-profile build and started its CPU arm at
2026-10-02 15:42:41 KST on syn101. Binary SHA256 is
`fee1bdb251971f3e83e631a79b4fb044a64ad7c16ddf66223fea8e4c458ab858`;
effective namelist SHA256 is
`ae07f3969474caee80c2cf2f630b1a012b6f191e69d6c89196b38928ac86b9d1`.
The log confirms FFTW execution with eight threads, zero reported post-force
and post-sync non-finite counts, and actual dark-material entry with 2097152
local leaf cells, tile width 128 and 16384 tiles. The first 128-cell tile
reported CHIMES wall 0.053151 s, dust wall 0.00064397 s and status 0. This is
only one tile, not a full-level cost or a CPU/GPU speedup. At this observation
the CPU arm is still running and the CUDA arm has not begun. Inputs/logs are
under `.cvode-gpu/chimes-real-ramses-128-410509/`.

By 2026-10-02 15:58:24 KST, the CPU arm completed its first complete
2097152-cell material level: CHIMES wall **793.19 s**, dust wall **4.078 s**,
approximately **2644 cell solves/s** across eight OMP workers. This is a
measured full-level cost, not the earlier first-tile extrapolation. The
material ledger reports max_balance=0; main-step printed mass/energy
conservation errors are zero and the subsequent non-finite checks report
zero counts. Initial dust in this registered IC is zero, and the ledger's
IR/CMB/gas-dust energy exchanges are zero; these timings therefore chiefly
exercise gas chemistry, not nonzero dust heating or growth qualification.
The CPU arm is advancing its second step. CUDA execution and the final
endpoint comparison remain pending. Live Slurm step peak RSS was about
11.30 GiB; this is not the final per-rank peak or whole-job build peak.

At 2026-10-02 16:12:42 KST, job 410509 and its MPI step remained RUNNING.
The completed first coarse-step timer reports 863.98 s total, 823.228 s
in `snrt_advance`, and 793.19 s in CHIMES: approximately 95.3% of the
coarse-step wall time is in the material advance and 91.8% in chemistry.
This identifies chemistry, not Poisson or dust, as the dominant measured
cost for this IC. The second CPU material level is still active; step
accumulated CPU time increased from 03:17:08 at 16:09:03 to 03:43:59 at
16:12:42, so unchanged tile-level logging alone is not evidence of a stall.
The CUDA arm has not started, and no CPU/GPU performance or endpoint-parity
verdict is available yet. No duplicate or replacement job was submitted.

The CPU arm finished successfully at 2026-10-02 16:15:08 KST, EXIT=0,
with two complete material advances and one final HDF5 output. Wall time
was **1947.259 s**; per-rank user CPU time was 14435.101552 s, system CPU
10.317004 s and peak RSS 12054004 KiB (about 11.50 GiB). The second material
level took CHIMES 970.74 s and dust 4.1857 s. The final printed main-step
mass/energy errors were 7.14e-16 and -1.99e-9; no NaN was reported. The
CUDA arm started automatically at the same timestamp with the intended
`cuda_integrated` selector and FFTW marker.

One runner-verdict defect was found during this transition: it requires
CPU `cpu_setups>0`, but `snrt_chimes_cvode_linear_solver` only installs the
counted custom setup callback when CUDA is enabled. Native CPU SUNDIALS
therefore correctly leaves the adapter's counters at zero. This is NOT
evidence that CPU chemistry was skipped. The running script is left
untouched; its eventual CPU-counter rejection must not be reported as a
physics failure or used to restart the completed CPU baseline. Correct that
condition after the live script finishes, using actual CPU completion and
material-advance evidence instead. CUDA device-work evidence and endpoint
comparison are still required. Because this rejection precedes successful
cleanup, raw comparison outputs will remain available for end evaluation.

Read-only inspection of the completed CPU HDF5 metadata confirmed 187
hydro datasets at `/hydro/level_7`, each containing 2097152 leaf values;
dataset names are `uold_1` (not zero-padded), etc. Runtime indices are
ichem=7, idust=18 and ichimes=idust+11=29, so the 157 transported chemical
carriers occupy uold_29 through uold_185. These are m_H*n_species carrier
densities, not the per-H abundances used by the existing native parity
fixture. Its species tolerance (1e-17 + 1e-6*abs(reference abundance)) must
therefore not be applied directly to raw conservative field values. The
final full-field report remains descriptive until evaluated in matching
physical units. At 16:19:50 KST the CUDA MPI step was confirmed RUNNING;
its first tile returned status 0 after 0.070305 s of CHIMES work. Neither
this first-tile timing nor successful selector admission proves device
work; the exit counters and full-level timings are still pending.

### Independent-CHIMES CUDA pool wiring correction

A single read-only nvidia-smi sample inside job 410509 showed 0 MiB device
memory and 0% utilization. Source inspection then established a wiring
defect independently of that sample: `adaptive_loop.jaehyun.f90` initialized
the shared stream pool only for hydro/Poisson/FFT/sink/scalar/particle/MHD
GPU flags. All those flags are deliberately false in this isolated
RT-none material benchmark. CHIMES's selector configured its adapters but
never initialized the pool; `cuda_acquire_stream` therefore returned -1.
The current run cannot establish genuine CUDA performance. Job 410509 was
cancelled after its successful CPU arm; its logs, inputs and CPU snapshot
are preserved, and its cancelled CUDA timing is not an admitted comparison.

The bridge now exports a side-effect-free startup query that reuses the
existing unified selector admission (and the legacy experiment selectors).
The adaptive loop uses that query only when CHIMES dust chemistry is active,
and includes CHIMES in collective node-local-rank stream-pool initialization
before any OpenMP cell callback. A requested CHIMES CUDA mode with no pool
now prints an error and stops, rather than silently acting as a CPU-only
performance comparison. CPU mode does not request a CUDA context. No
physics tolerance, RHS, timestep, radiation selector or VPATH was changed.

The comparison runner replaces the invalid native-CPU adapter-count check
with evidence of both completed material levels, requires the CUDA pool
startup marker, and retains snapshots until numerical stage evaluation.
Its shell syntax check passed. Replacement job **410554** was submitted
with the same one-A100/H100, MPI1/OMP8, 64-GiB, 24-hour configuration and a
fresh isolated build. Both arms are rerun with the same corrected binary;
the earlier successful CPU result remains a baseline reference, not a
substitute arm under a different binary. Runtime/build validation is still
pending; no speedup or production-GPU admission is claimed.

Job 410554 completed the fresh build without compiler/linker failure and
started its CPU arm at 2026-10-02 16:29:53 KST on syn101. Corrected binary
SHA256 is `5d8823da1376c061cd84c433687f012ad2fa87222554d0262fe6c1db34cbdde7`.
The effective namelist retains `radiation_transport='none'` and nstepmax=2.
CPU startup admits `CHIMES integrated compute backend: cpu`, does not print
a CUDA-pool initialization marker, and reaches the first cooling operator
with zero preceding non-finite checks. This verifies the corrected sources
build and preserves CPU-mode startup without a GPU context; independent
CHIMES CUDA-pool initialization still requires the forthcoming CUDA arm.

By 16:45:37 KST, the corrected CPU arm completed its first full level:
CHIMES 799.48 s, dust 3.8979 s and coarse-step wall 868.58 s. Its printed
chemical energy debit (-1.9649e43 erg), zero exchange/balance diagnostics
and zero main-step conservation errors match the earlier CPU log at the
printed precision. Chemistry wall differs by about +0.79% from 793.19 s;
one measurement does not establish a performance regression or improvement.
The second CPU step is active. Full endpoint parity and CUDA work remain
unverified until the comparison finishes.

The corrected CPU arm finished EXIT=0 at 17:02:26 KST: wall 1952.894 s,
rank user CPU 14470.840832 s, system CPU 9.980400 s, peak RSS 12052844 KiB.
Its second full chemical advance took 973.99 s (dust 4.0495 s). The CUDA
arm began automatically at 17:02:26 with the same binary and namelist.
It now prints the A100 local-rank device mapping with streams=1 and
`Adaptive loop: CUDA pool early-init, available=T`. This verifies the pool
wiring correction in the intended RT-none production material path.
A single in-allocation device sample showed 59% GPU utilization and
445 MiB allocated, unlike the cancelled uninitialized run's 0 MiB. That
sample is not an averaged utilization or speedup. The first 128-cell CUDA
tile returned status=0 with CHIMES wall 0.15899 s; the corresponding CPU
first tile was 0.050204 s. Complete-level/device-counter/endpoint evidence
is still pending, and the first tile must not be extrapolated into a final
performance claim.

Read-only callback inspection during the real CUDA arm confirms that larger
domain size does not itself increase device batch width: both coalescers
target ceil(OMP team size / stream count), with a minimum of four. For
OMP8/one stream the target is eight pending cell callbacks, not a 128-cell
tile or the full leaf grid. Both RHS and LU paths enqueue asynchronous
copies/kernels but synchronize their own leased stream before publishing
results to host CVODE. LU returns the factorized matrices and pivots to
host memory; this is not a device-resident integrator. These are candidate
latency/occupancy limitations to interpret alongside final counts/timings,
not measured attribution or permission to use stale CVODE RHS values.
No solver arithmetic or live execution was changed. At 17:09:39 KST,
410554.1 was confirmed RUNNING with increasing CPU usage and no new error.

At 2026-10-02 17:37:26 KST, the real CUDA arm was confirmed still RUNNING,
without a completed first-level ledger or a new error. Since RUN_BEGIN was
17:02:26, elapsed arm time was already at least 2100 s, exceeding the same
binary's complete two-step CPU time of 1952.894 s (lower-bound runtime ratio
about 1.075). Thus this actual 128^3 configuration cannot demonstrate a
wall-time speedup over its CPU baseline, even before CUDA completion. This
is a lower bound, not a final CUDA duration or a numerical-parity verdict.
The job is retained to obtain the final endpoint and device counters.

### 2026-10-02: warning attribution and aborted-trajectory guard

At 17:47 KST, 410554 was still RUNNING but its CUDA log had grown to
approximately 1.19 GB with repeated CHIMES CVODE failure warnings and full
species dumps. These warnings are not completion evidence. Source inspection
shows that the same generic warning is printed for any negative CVODE result,
not exclusively recoverable RHS failures. In particular, an unavailable RHS
stream after a cell has used GPU arithmetic deliberately returns -1 from f;
the bridge subsequently retries the whole cell on CPU from original inputs.
Therefore the warning alone does not prove that the accepted cell failed or
that its state is invalid. Final status/parity and device counts remain needed.

A bounded adapter guard now keeps a requested whole-cell retry sticky: further
RHS calls on that discarded trajectory return -1 until the existing bridge
switches to its fresh CPU retry. Previously the flag was set but not checked
at RHS entry, allowing a later callback to acquire a device stream again.
This changes neither tolerances nor accepted-state conservation rules. The
running job uses the previous immutable executable; this guard does not alter
410554, and does not by itself establish the cause of its elapsed time.

The modified CUDA adapter compiled successfully on syn101 using CUDA 13.0.2,
SM80, FP64 CHIMES and SNRT_CHIMES_RHS_TESTING (srun diagnostic step 410554.4,
RHS_COMPILE_STATUS=0). The first compile-only diagnostic step 410554.3 failed
because its command omitted the MPI include directory required by parallel
HDF5; it was corrected without changing the simulation. This is compile
verification only: runtime regression and endpoint parity are still pending.
At 17:51:47 the original CUDA arm remained RUNNING, with no completed first
level and approximately 1.41 GB of warning/state-dump log. No cancellation,
restart, tolerance relaxation or silent suppression of true solver errors
was performed.

The existing native RHS regression was submitted as job 410614 using
chimes_integrated_gpu_qualification.sbatch with
SNRT_CHIMES_RETRY_REGRESSION_ONLY=1, A100/OMP8/16 GiB and syn102 excluded.
This bounded mode compiles the current adapter and runs the existing `rhs`
fixture, including busy-pool fallback, serial fallback and forced midcell
drop/restart parity. Default integrated qualification coverage is unchanged.
It produces no RAMSES snapshots and cannot establish large-run speedup. At
17:55:41 it was PENDING (Priority); 410554.1 remained RUNNING. No duplicate
real simulation was submitted.

While 410614 remained PENDING (verified again at 18:03:07), the existing
RHS fixture gained a test-only direct assertion of the pending-retry guard.
It presents no usable UserData after setting a pending GPU-trajectory abort,
requires rejection before any trajectory/table/stream access, and restores
the thread-local state. This catches re-entry independently of whether the
subsequent CPU solve happens to give a matching endpoint. It adds no
production counters, numerical gates or physical-model changes. The marker
is RHS_PENDING_RETRY_STICKY_PASS; its runtime result is still pending.

### Completed real 128^3 comparison and retry regression

At 18:23:18 KST both jobs were authoritatively COMPLETED, exit 0. Job
410554 completed both actual gas+DM cosmological arms, and 410614 completed
the modified adapter's existing RHS fixture. The sticky guard marker passed;
eight forced midcell CPU retries passed with zero endpoint discrepancy.
Cold/warm, busy-pool and serial RHS fixtures had zero rejected cells, zero
device errors and existing species/thermal/conservation criteria passed.
These are correctness results, not acceleration evidence.

Real-run wall times (same binary, IC, namelist and OMP8 placement):

| Quantity | CPU | CUDA integrated |
|---|---:|---:|
| Complete two-step arm, seconds | 1952.894 | 4795.146 |
| First level CHIMES, seconds | 799.48 | 2450.4 |
| Second level CHIMES, seconds | 973.99 | 2162.6 |
| Rank user CPU seconds | 14470.840832 | 17909.647098 |
| Rank system CPU seconds | 9.980400 | 1624.664045 |
| Peak RSS, KiB | 12052844 | 12177996 |

CUDA wall time is 2.4554 times CPU; CHIMES level time is 2.6011 times CPU.
The production default remains CPU. CUDA is an explicit experimental
correctness-compatible option in this tested domain, not a production speedup.
The large-run binary predates the sticky guard; only its separate regression
has tested the guard, so no full-run timing improvement is attributed to it.

The completed Slurm-side HDF5 comparison inspected 1395 datasets with equal
layouts, zero nonfinite pairs and zero discrete mismatches. All finest-level
hydro/carrier fields 1..28 were exactly equal, including total gas hydrogen,
energy and dust carriers. The largest pointwise relative difference was
8.566037894953708e-10 (a trace species); largest field-peak-normalized
difference was 4.361764659982712e-11, absolute difference 2.06887006e-36.
With the identical gas-hydrogen normalization, this bounds the existing
per-H species error norm using 1e-17 + 1e-6*abs(CPU) below 0.000857; the
positive absolute floor only reduces it. Hydro energy and density/momentum
were exactly equal. Both level chemistry-energy ledgers agree at printed
precision and max_balance is zero. Endpoint parity is therefore accepted
for this bounded cold, near-primordial, dark two-step test only. Initial dust
is zero; this does not qualify illuminated, metal-rich or nonzero-dust states.

Actual device work: 12197058 GPU LU factors, 2644323 LU batches, 24284909
GPU RHS calls, 6298256 RHS batches, one resident table upload, zero device
errors, peak concurrency one and one used stream. LU reported 1379458 busy
fallback cells and 37114308 CPU setups. RHS reported 773846202 CPU fallback
calls (callbacks, NOT distinct cells or full-cell retries). Millions of small
callback batches, host/device synchronization and restarts are supported
bottleneck candidates; their separate time attribution is not measured by
these counters. Generic CVODE warning dumps are not final failed-cell counts.
No tolerances were relaxed. Next optimization must address callback
granularity/residency and avoid discarded GPU work, not simply rerun this
configuration on a faster card or use stale CVODE states.

Stage evaluation is complete. Retain the effective namelist, table/build
identities, timing files, logs, field comparison TSV and output metadata;
remove only these completed, no-longer-needed raw comparison snapshots:

- .cvode-gpu/chimes-real-ramses-128-410554/cpu/output_00001/data_00001.h5
- .cvode-gpu/chimes-real-ramses-128-410554/cuda/output_00001/data_00001.h5

Each raw file is 3883831760 bytes; total authorized removal 7767663520 bytes
(approximately 7.23 GiB). Earlier failed/cancelled-job outputs are not targets.

Cleanup executed successfully and both exact raw-file targets were verified
absent. Metadata and compact evidence remain. These raw snapshots cannot be
recovered from this checkout; reproducing them requires rerunning the retained
inputs with the recorded binary. No earlier failed/cancelled output was removed.

### Contention follow-up after the completed real comparison

The experimental RHS coalescer now retains an already-GPU trajectory through
transient stream contention: its batch leader waits without a queue mutex or
stream lease, sleeping between acquisition attempts rather than spinning.
Executing RHS and LU batches drain their own work and release their leases
before returning to CVODE, so there is no lease-to-callback dependency cycle.
New trajectories still use CPU when leases are unavailable. Pool disappearance,
unsupported state, forced test drop and device failure still use the original
whole-cell CPU retry; no stale RHS or tolerance changes are introduced. The
MG/global pool implementation and Makefile VPATH order are untouched.

Job 410626 was submitted to rebuild and run the existing RHS retry fixture
and integrated RHS+LU fixture together, on A100/OMP8/16 GiB, without RAMSES
snapshots. At 18:30:50 it was PENDING. The change is not promoted and no timing
benefit is claimed pending tests. This small regression is not a substitute
for a subsequent actual-workload performance comparison.

### Staged cell-loop follow-up (2026-10-02)

410626 completed with exit 0: isolated RHS busy/serial fallback and eight
forced midcell retries passed; integrated RHS+LU passed with zero device
errors. The transient-lease waiting experiment was subsequently withdrawn
following the user's objection: it preserves arithmetic but blocks a host
worker rather than scheduling independent cells. Legacy callback contention
again uses the original bounded fresh CPU retry.

New `snrt_chimes_rhs_tile_run` retains independent CVODE host continuations,
collects their current RHS candidates, computes batches into temporary rate
arrays and resumes the dependent CPU stage. No extra OS threads, stale RHS,
tolerance relaxation, new chemistry model or Makefile VPATH changes. One
stream lease belongs to the complete cohort; busy-at-entry cohorts run CPU.
Host dense LU is CPU-only while a cohort owns the lease, avoiding a nested
lease dependency. Stack allocation is bounded (256 KiB per continuation);
exact rate cache and adapter TLS controls are switched between cells.

410633 completed exit 0 in 20 s, A100/OMP8/16 GiB. The 24-cell staged path
passed endpoint status/conservation checks: rejected cells 0, status mismatch
0, species scaled error 0.00128314 (existing acceptance <=1), temperature
relative difference 1.89e-15, nuclei relative difference 1.11e-16. Busy-pool
whole-CPU fallback was exactly identical to CPU. GPU RHS calls 9964, batches
746, device errors 0. Timings: CPU OMP8 0.0188202 s; staged single-worker GPU
0.149146 s. This proves the mechanism, not production acceleration.

The RT-off cold/transition Fortran operator now stages input preparation,
disjoint 16-cell OpenMP cohort callbacks, result validation, and existing
dust/level commit. It reuses the full native cold/thermal-root adapter; no
AMR writes or MPI occur in callbacks. Grey/band and CPU-selected paths remain
unchanged. CHIMES_CUDA now defines its Fortran build macro, without VPATH
reordering. Cohort stack capacity is 4 MiB/worker, 32 MiB at OMP8. Callbacks
write disjoint private result slots; other workers can compute CPU cohorts
while GPU cohorts run. Host CVODE/cooling/LU still run on CPU, and a worker
does not yet steal a new cohort during its own device synchronization.

410638 submits an isolated fresh native RAMSES build only (A100/H100,
OMP8/64 GiB, 40 min); it does not start a simulation. 410639 submits the
extended existing fixture with concurrent cohorts sharing one available
lease and forced midcell restart checks. Both are pending validation; no
actual-workload speedup or production promotion is claimed.

410640 is submitted with `afterok:410638:410639`, reusing exactly the fresh
410638 binary for the existing two-step 12.5-cMpc/128-cubed CPU/CUDA actual-IC
comparison. Thus it starts automatically only after both prerequisites pass.
MPI1/OMP8/one stream and all physical inputs remain those of 410554; this is
not stream-count tuning. If a prerequisite fails, this dependent comparison
must not be mistaken for a running simulation. Queue resource waits are not
numerical failures. No new external audit or larger resource campaign is
introduced for this ordinary scheduling implementation.

410638 completed exit 0 in 4 min 47 s: fresh native Intel Fortran/CUDA build
including the BIND(C) staged callbacks linked successfully. Binary SHA256:
`6d6a3a249afd328480fad62d3b2bb4bc9331fb3e57140356d8ec52a88267db26`.
410639 completed exit 0 in 20 s. Concurrent-cohort parity passed with status
mismatches/rejections 0 and species scaled error 0.000190169; forced midcell
retry was exactly CPU-identical. Staged GPU timing 0.146496 s vs CPU OMP8
0.0194237 s remains a small-fixture regression, not acceleration evidence.

410640 began its CPU arm at 19:01:16 KST, but was cancelled early by the
operator because the retained comparison harness incorrectly required
nonzero GPU LU work for the new CPU-LU staged design. This is a harness
assumption, not a numerical failure. The script now uses an explicit staged
qualification mode requiring positive GPU RHS work and zero GPU LU work;
the original integrated mode still requires positive work in both. No
physics tolerance was relaxed and final numerical comparison remains
required. Replacement 410642 reuses the same verified binary and namelist,
with `CHIMES_REAL_STAGED=1` and the successful prerequisite dependencies.
Raw files from the interrupted job are not deleted: no stage evaluation has
been completed. The replacement's two full coarse steps and final comparison
are outstanding.

Follow-up monitoring: at elapsed 26 min 28 s, 410642 was still in its second
CPU material step, not the CUDA arm. The first CPU CHIMES level took
800.36 s versus the prior 799.48 s reference, with identical printed thermal
exchange (-1.9649e43 erg) and zero material balance residual. Step accounting
reported about 3 h 17 min CPU time at this wall time, consistent with the
eight-core allocation; peak RSS approximately 11.3 GiB. These are partial
progress evidence, not a completed comparison or a proof of no future stall.

The new tile test's forced-drop hook calls are now guarded by
`SNRT_CHIMES_RHS_TESTING`, matching the older RHS fixture. Both ordinary GPU
fixture and CPU-only fixture configurations passed C++ syntax checks without
the test-hook macro. This fixture-only fix does not change the running native
binary, its source identity at build, or the science inputs. No extra runtime
campaign was introduced to check this compile-time guard.

410642 CPU arm completed exit 0 at 2026-10-02 19:35:27 KST, elapsed
1956.826 s (32 min 36.826 s), versus prior 1952.894 s. CHIMES level times
800.36 and 975.34 s; dust times 3.5371 and 3.7016 s. Printed energy changes
match the prior reference (-1.9649e43 and -1.9559e43 erg), balance residual 0.
Rank user time 14498.619765 s, system 9.856255 s, peak RSS 12046780 KiB.

The CUDA arm started automatically at 19:35:27 KST with the same binary and
effective namelist, and reported `cuda_integrated`. Its first real 128-cell
tile completed with status 0: CHIMES 0.076450 s and dust 0.00067711 s,
compared with this run's first CPU tile 0.051515 s. Ratio about 1.48 for this
single startup tile, not the total workload. This confirms the staged
Fortran/C continuation path can complete an actual-IC tile; GPU counter and
full-field parity evidence remain outstanding until the CUDA arm finishes.
Do not infer acceleration, general chemistry qualification or completion
from this first-tile timing. The job continues through its two full coarse
steps and automatic compact output comparison; raw data remain retained
until numerical stage evaluation.

The first complete CUDA chemistry level is measured at 980.90 s versus
800.36 s for CPU, a 22.6% increase in that level's CHIMES time. At total job
elapsed 01:01:46, the CPU-only arm had finished in 1956.826 s (32:36.826),
and the sequential GPU/hybrid arm had run about 1749.174 s (29:09.174) and
was still running. Compare those arm times, not the total Slurm job time:
the hybrid arm had used about 89.4% of the CPU-only wall time, leaving about
3 min 28 s before it would lose on elapsed time if it did not finish. This is
only a live bound, not the final speed result. The second CUDA coarse-step
material solve, exact arm wall time, counters and HDF5 state comparison remain
outstanding. Keep the run for endpoint parity and report the completed
per-arm timing files as the benchmark evidence.

Operator clarification (2026-10-02): the hybrid design objective is additive
throughput: CPU workers continue independent cells while the GPU processes
its assigned cells. If GPU throughput matches the CPU-only aggregate
throughput and both lanes remain fully occupied on disjoint work, ideal
wall-clock speedup is 2x. Per-kernel CUDA submission being asynchronous is
insufficient if the host worker blocks and CPU capacity is stranded. The
current `execute` path synchronizes each RHS batch before returning to the
tile scheduler; other OMP workers can do CPU cohorts during that wait, but
the GPU-owning worker cannot take another cohort. Evaluate this run using
whole-arm wall clocks. The next performance decision must determine whether
that worker wait/limited scheduling overlap explains the observed throughput;
do not describe this architecture as fully asynchronous.

### 2026-10-02: asynchronous dynamic R/B stream dispatch

The completed CHIMES comparison 410642 recorded CPU-only 1956.826 s and
the previous staged CUDA/CPU route 2305.353 s, both exit 0. The staged route
was 17.8% slower. These timings describe that chemistry implementation;
they are not timings for the new MG smoother.

The experimental `mg_dynamic_hybrid` path now dynamically claims independent
4096-grid batches in each red/black pass. A worker atomically leases an idle
stream from the existing rank-local pool, submits H2D identifiers, smoothing
and D2H results asynchronously, then continues claiming work. While its
submitted batch is pending, it computes subsequent batches on the CPU.
Workers without a stream also compute on the CPU. Nonblocking queries reap
finished results; each worker drains its outstanding batch at the color
boundary. Separate device and pinned-host buffers belong to each stream;
allocations happen before the parallel region. CPU work-list reservations
use one atomic operation per batch, not per grid.

Same-color batches write disjoint cells and read the opposite color. Completed
GPU values are merged to host phi; CPU color updates are scattered to device
phi before the next color. Existing MPI halo exchange updates both copies.
GPU coarse interpolation now refreshes host phi before hybrid post-smoothing,
so CPU batches see the coarse correction. The existing exclusive GPU and
CPU smoothers remain selectable. Namelist generator, mkrun warning and
namelist reference describe the experimental mode.

The existing stream pool selects one device per MPI rank. Several ranks can
map to several GPUs, and workers select their rank's device when acquiring a
stream. Arbitrary use of several GPUs by one rank is NOT implemented here:
MG state is presently resident on one device. That requires per-device state
and synchronization of color updates. Stream leases are process-local, not
a card-wide lock shared by MPI processes.

CUDA and ISO_C_BINDING translation units compiled successfully in isolated
temporary directories; diff whitespace and runner shell syntax checks pass.
Full production Fortran build and numerical/timing verification are pending
Slurm job 410745, using `.cvode-gpu/mg_rb_dynamic_compare.sbatch` and
`.cvode-gpu/mg_rb_base.nml`. This compares CPU-only, GPU-only and hybrid with
identical 128^3 gas/DM ICs, one rank, eight threads and three pool streams.
Cooling/feedback are disabled to isolate gravity; this is not the coupled
CHIMES performance result. Three streams are a test configuration, not an
optimized stream count. Output policy: noutput=1, aout=1.1 (outside the
five-step interval), tout=1e100, foutput=5, fbackup=1000000; one expected
endpoint dump per arm, roughly order 1 GiB each, retained through comparison.
The checked GPFS mount had 223 TiB available. Effective namelists and timing
files will reside in `.cvode-gpu/mg-rb-dynamic-410745/{cpu,gpu,hybrid}/`.
Do not claim speedup or numerical acceptance before this comparison finishes.

### 2026-10-02: CHIMES fixed GPU worker comparison

Operator requested OpenMP dynamic scheduling with worker 0 exclusively
owning GPU submissions. `SNRT_CHIMES_GPU_WORKER=thread0` now restricts both
staged cohort stream acquisition and the RHS entry point to OpenMP worker 0.
Other workers use the original CPU chemistry. The default `any` preserves
the previous lease-based scheduler. Unknown values fail configuration.
GPU batch execution checks ownership, and the existing exit reporter now
prints `CHIMES_GPU_WORKER_MASK`: value 1 means only worker 0 submitted GPU
batches. Nonzero device work must also be demonstrated by existing RHS
counters. Dynamic scheduling remains `schedule(dynamic,1)` in the Fortran
cohort loop. The GPU worker still performs host CVODE continuation and waits
for its own current-state RHS batch; this change fixes worker assignment,
not the GPU/host CVODE dependency or batching size.

The modified CUDA unit compiled successfully with FP64/64-bit SUNDIALS
headers. Shell syntax and diff whitespace checks passed. Submitted job
410870 (`chimes-thread0`) runs a fresh full production build on a compute
node, then CPU-only, previous `any` stream assignment and `thread0`
assignment using one common executable, GPU allocation and input. It uses
the previously verified 128^3 actual gas/DM IC profile, two coarse steps,
MPI=1, OMP=8, one stream, radiation_transport=none, cosmological CMB/dust
and CHIMES enabled. All three arm wall times and rank CPU resource counters
are retained. Compare thread0 to both same-job CPU and any baselines; the
previous A100 measurement is historical context only.

Runner: `.cvode-gpu/chimes_thread0_compare.sbatch`; partitions h200/h100/a100,
one GPU, eight CPUs, 64 GiB, maximum one day. Effective inputs will be at
`.cvode-gpu/chimes-thread0-410870/{cpu,stream,thread0}/effective.nml`.
Output policy remains noutput=1, aout=1.1, tout=1e100 (unreached), foutput=2,
fbackup=1000000: one endpoint dump per arm, estimated 2.93 GiB hydro fields
plus particles/metadata each. Checked GPFS mount free space was 224 TiB.
Regression compares all HDF5 fields, requires finite floats and exact
discrete fields, and declares field-peak relative drift <=1e-8 before launch;
this is a numerical regression criterion, not a galaxy-physics calibration.
Raw test dumps are removed only after the stage comparison and performance
report finish successfully, with a cleanup manifest; failures retain them.
Job was PENDING/Priority at the submission check; no new speed measurement
or numerical verdict is yet available.

2026-10-03 follow-up: 410870 FAILED with exit 127 after 4m43s. The full
production build and link completed, but the first CPU arm never launched:
the compute node lacked `/usr/bin/time`. This was a benchmark-runner error,
not a CHIMES numerical failure. Replaced the external timer with Bash's
reserved `time` and TIMEFORMAT; launcher CPU times are labeled explicitly,
while rank CPU time and peak RSS remain in the exit reporter. Syntax/diff
checks pass. Verified both the binary and CHIMES library against the original
build SHA256 manifest. Replacement job 411173 reuses that exact build via
CHIMES_THREAD0_REUSE_BUILD and writes fresh arm inputs/results under
`.cvode-gpu/chimes-thread0-411173/`. No arm timing or numerical comparison
was produced by failed job 410870; historical timings must not be reported
as the new thread0 result.

### 2026-10-03: Astra level-queue/event-broker implementation

Completed comparison 411173 on H100 NVL: CPU8 wall1954.837 s, any-worker
lease2353.285 s, thread0 fixed2361.603 s. Both GPU arms matched1395 datasets
with nonfinite0, discrete mismatch0 and field-peak relative6.841759761981513e-11.
Thread0 ownership mask1 passed, but assigning the worker alone did not improve
throughput. Raw outputs were removed after the regression/performance check.

Operator approved the Astra proposal. Added opt-in
`SNRT_CHIMES_CELL_SCHEDULER=level_queue` to the RT-off cold CHIMES path.
One atomic cell-claim queue spans the entire rank-local level. Native CPU
workers execute independent whole-cell solves; with CUDA enabled, worker0
owns up to128 fixed-thread CVODE continuations and two separate64-request
staging buffers on one leased stream. Events identify completion of each
buffer's H2D/kernel/D2H sequence. The broker resumes only the suspended
continuations whose own current-candidate results have arrived; other ready
continuations can progress while device work runs. Requests live on stable
fiber stacks and are checked against their owning fiber before publication.
Partial batches are submitted, and waiting happens only when no host
continuation can advance. CPU-only uses the same level queue with all workers.

CPU workers bypass the GPU RHS counter and custom LU adapter, preserving
native dense CVODE arithmetic without per-callback shared atomics. The GPU
broker retains host CVODE/LU/cooling. Existing whole-cell CPU restart from
original inputs handles rejected GPU candidates; runtime device failure
drains pending work and disables device submission for that level. No stale
RHS is substituted and no continuation is migrated between OS threads.

The Fortran caller keeps cell inputs read-only, computes chemistry into
level-lifetime private staging, then runs the existing dust tiles and MPI
admission/commit. Dust temperature is derived from original cell material
before chemistry, and preparation is repeated for the subsequent dust tiles.
The duplicated preparation cost is included in both matched arms. There is
one queue startup rendezvous and one chemistry level join, with no per-tile
chemistry rendezvous in this mode. Default cohort scheduling is preserved.

CUDA RHS/broker and CVODE units compiled in isolated temporary directories.
The modified Fortran runtime compiled against the completed production build
modules with matching NDIM3/NVAR187/NVECTOR32 flags. Diff whitespace and
runner syntax checks passed. Existing bounded tile regression now also
checks level CPU parity, broker parity, pool-busy CPU fallback and forced
midcell CPU restart; these GPU runtime checks await the submitted allocation.

Submitted411408 (`chimes-broker`) combines that bounded regression, a fresh
production build, and a matched actual128^3/two-step CPU8 vs CPU7+broker
comparison. It requests one GPU, eight CPUs,64 GiB on h200/h100/a100.
Both arms use the new level queue, identical binary/IC/tolerances and one
stream. Numerical criteria remain finite fields, exact discrete fields and
field-peak relative<=1e-8; performance is judged by total arm wall, not
callback counts. Preserve inputs/logs/counters/compact comparison; remove
raw output only after evaluation succeeds. Output policy remains one final
dump per arm (noutput1 with unreachable aout1.1/tout1e100, foutput2,
fbackup1000000), estimated2.93 GiB hydro fields plus particles/metadata per
dump. Checked free GPFS space220 TiB. Evidence root will be
`.cvode-gpu/chimes-broker-thread0-411408/`, runner
`.cvode-gpu/chimes_thread0_compare.sbatch` with
`CHIMES_LEVEL_BROKER_COMPARE=1`; job log `.cvode-gpu/chimes-broker-411408.log`.
Initial job status PENDING/Resources. No runtime parity or speedup claim yet.
# Four GPU brokers: operator-requested comparison (2026-10-03)

## Compact dark RHS results (2026-10-04)

Initial three-arm attempt 411973 reached the A100 regression fixture but
stopped at the explicit assertion that compact batches had executed. The
reason was an overbroad eligibility condition: the zero-photon cold solve
retains `rt_update_flux=1` so CHIMES uses case-B recombination. The registered
host molecular callback under this state writes exact-zero photo rates and
does not consume omitted reaction arrays; `N_spectra==0` is the correct
compact-path boundary. The guard now checks actual photon variables only.
No measured run was made in 411973. CUDA and fixture syntax checks passed
after the fix; matched retry **412380** is queued for regression followed by
the three-arm full workload comparison. Results remain pending.

Timeline 411868 completed on A100-SXM4-80GB in 7 minutes including build.
The 20.0745-second observed CUDA activity window contained 24,565 reaction
kernels. Kernel interval union was 0.67414 s; memory-copy interval union
14.12297 s; combined device-activity interval union 14.58882 s. H2D payload
was 2,065,315,560 bytes and D2H 29,600,434,600 bytes. Broker NVTX wait ranges
were 11.71–11.93 s per worker, host sweeps 6.78–6.93 s per worker. These
overlapping, instrumented intervals are not additive wall-clock phases or
an uninstrumented whole-run performance estimate. The submission-failure
fresh-CPU-retry regression passed, and no raw simulation dump was created.

The operator approved proceeding with transfer reduction. Full reaction
scratch stays on device; compact output retains separate creation and
destruction rates for all 157 species, critical densities, five molecular
reaction values read by CPU cooling, and cosmic-ray rates. Unused fixed-size
reaction/coefficient arrays are no longer copied to CPU on the admitted
dark path. Device finite checks retain coverage of the full computed result.
FP64 arithmetic, species sums, cooling formulas, current-candidate ownership
and CVODE tolerances are unchanged. Payload is 2840 vs 18920 bytes/callback
(84.99% smaller); this is a byte-count improvement, NOT a speedup claim.

`SNRT_CHIMES_GPU_RESULT=full` preserves the previous transfer; paired rate
verification, flux callbacks and unsupported table-index layouts use full
results. Existing fixtures now compare compact/full trajectories to the CPU
reference and require actual compact batch execution. Completed-result byte
and layout counters are exported through the existing exit reporter.

Matched comparison **411973** was submitted: same allocation/GPU/binary,
CPU8 versus CPU4+four brokers with full results versus compact results,
actual 128-cubed IC, two coarse steps, one final dump per arm. Estimated
hydro output 2.93 GiB/arm plus particles/metadata (three arms); free GPFS
space at submission approximately 182 TiB. No tracing is enabled in this
wall-clock comparison. CUDA object compile, C++/C fixture syntax, shell
syntax and diff checks passed before submission; device regression and
real-workload timing results remain pending.

Log: `.cvode-gpu/chimes-compact-411973.log`.
Effective inputs:
`.cvode-gpu/chimes-broker-thread0-411973/{cpu,stream,thread0}/effective.nml`.
Historical arm names `stream` and `thread0` here mean full and compact
results respectively; both hybrid arms use workers 0–3 and four streams.
This experiment reduces payload only; host CVODE and repeated RHS
round trips remain. A full device solver is not claimed implemented.

## Astra-followup bounded timeline

The driver received the Astra follow-up analysis: successful GPU RHS result
consumption is established, but broker cell throughput is below CPU-worker
throughput. Per-RHS host CVODE/LU, context switching and output copies remain;
the contribution of host supply gaps, event ordering and device runtime is
not yet measured. Whole-team repeated chemistry barriers were not found.

Operator-approved follow-up **411868** uses one real 128-cubed input and
four brokers on one GPU (OMP8, streams4), with a 20-second Nsight Systems
CUDA/NVTX capture starting after broker buffer preparation. Ranges identify
host sweeps, output consumption, submission sequence/flight and waited event
sequence/flight. The profiler intentionally terminates the application at
capture end; this is NOT a completed evolution or a speedup measurement.
No simulation dumps are scheduled; allow a 5-GiB diagnostic trace budget.
Free GPFS space at submission was approximately 188 TiB.

The repeated-batch failure path was corrected to skip querying/consuming
reusable events after an unsuccessful submission. An isolated regression
injects a failed submission after prior event records and checks fresh CPU
retry against the existing reference. No transport, physical model,
tolerance, or normal broker scheduling was changed. CUDA object compilation,
fixture C++ syntax and shell/diff checks passed before submission. Device
regression and timeline results await allocation.

Runner: `simulation/snrt/tools/chimes_broker_timeline.sbatch`.
Log: `.cvode-gpu/chimes-timeline-411868.log`.
Effective input: `.cvode-gpu/chimes-timeline-411868/effective.nml`.
Expected results: `broker.nsys-rep`, `broker.sqlite`, `timeline-summary.txt`
in that run directory. Trace-on measurements must not be substituted for
the matched CPU-only/hybrid wall-clock benchmark.

Job 411744 failed after 22 seconds on syn09 (H100 NVL) during fixture
compilation: the four-worker participation assertion called
`snrt_chimes_rhs_worker_mask` without its public header declaration.
No production evolution or timing comparison ran. The driver added the
existing implementation's declaration to `snrt_chimes_rhs_cuda.h` and
resubmitted the same comparison as **411811**. Diff whitespace and shell
syntax checks passed; compute-node compilation/runtime are still pending.
Replacement log: `.cvode-gpu/chimes-broker4-411811.log`.
Replacement effective inputs:
`.cvode-gpu/chimes-broker-thread0-411811/{cpu,thread0}/effective.nml`.

Job 411408 completed on H100 NVL: CPU8 wall 1946.113 s versus
CPU7 plus one GPU broker 1954.225 s (0.42% slower). All 1395 datasets
passed the existing regression criterion; maximum field-peak-relative
difference was 6.651014300564651e-11, with no device errors.

The operator requested four GPU-compute threads. Job **411744** was
submitted with MPI=1, OMP=8, one GPU, four CUDA streams, and
`SNRT_CHIMES_GPU_BROKERS=4`, `SNRT_CHIMES_GPU_WORKER=any`,
`SNRT_CHIMES_CELL_SCHEDULER=level_queue`. Workers 0–3 each lease a
distinct stream and private staging buffers/continuations; workers 4–7
solve native CPU cells from the shared level queue. Host CVODE remains
on each owning worker; this is not full-device cell integration.

The allocation first executes bounded regression including four-worker
GPU participation and CPU-reference parity, then rebuilds and runs the
same actual 128-cubed gas/DM IC for two steps in CPU8 and four-broker
arms. Both effective namelists use `n_cuda_streams=4`; all physics and
tolerances are unchanged. One final dump per arm is expected (about
2.93 GiB hydro fields per arm plus particles/metadata). Free GPFS space
before submission: approximately 204 TiB. Accepted raw outputs are
removed only after the matched comparison, retaining logs and results.

Log: `.cvode-gpu/chimes-broker4-411744.log`.
Effective inputs: `.cvode-gpu/chimes-broker-thread0-411744/{cpu,thread0}/effective.nml`.
The inherited `thread0` directory label denotes the hybrid arm here;
its explicit policy is `any` with four brokers, requiring worker mask 15.
No four-broker runtime or speedup result is available at submission.
