# CHIMES solver cost and optimized SUNDIALS

Follow-up to the approved P0/P1 bundle; no new physics or solver tolerances.
Evidence lives in `.snrt-performance.jDx9Iz/`.

## Finding

Existing `.dust-extension.AOz7mU/sundials-build/CMakeCache.txt` has empty
`CMAKE_BUILD_TYPE` and `CMAKE_C_FLAGS`. Actual CVODE object flags contain
only `-fPIC`. The existing SUNDIALS library was compiled without optimization.
This is a concrete build issue, not a supposition about OpenMP efficiency.

Diagnostic-only `chimes_cost_probe.c` interposes on existing native calls via
LD_PRELOAD. It does not modify production code, solver tolerances, Jacobian
or callbacks' calculations. Use a serial probe: it prints main-thread totals
only. RHS timing includes Jacobian RHS evaluations; categories are nested.
Counters are sums per CVode invocation, appropriate to this probe's one call
per fresh solver. Do not apply the sums unchanged to resumed CVode sessions.

288 neutral matrix integrations, diagnostic timings in seconds:

| Item | Old unoptimized | Release |
|---|---:|---:|
| Whole chimes_network | 1.281338 | 0.818292 |
| CVode | 1.272354 | 0.810242 |
| All RHS | 0.694414 | 0.691418 |
| Dense LU setup | 0.169383 | 0.027859 |
| Dense solve | 0.158649 | 0.039613 |

Each run: 20992 accepted steps, 24800 ordinary RHS evaluations, 24832
linear-solver RHS evaluations, 544 Jacobians, 5056 LU setups, 23968 linear
solves, 384 error-test failures, zero nonlinear convergence failures.
288 network calls = 288 solver integrations; no extra bridge integrations
were observed in this matrix. No inference about live retry frequency.

Non-CVode work is 0.0090s (0.70%) BEFORE optimization, including more than
allocation alone. This workload does not justify workspace reuse as the
first large performance investment. Release whole-network time falls 36.1%;
same RHS time and counters support attribution to the numerical library.
Remaining optimized cost is predominantly RHS (about 84.5%). Counts alone
do not assign equal cost to ordinary versus numerical-Jacobian RHS calls.

## Build and precision

Same SUNDIALS5.8.0 source (e8a3e67e3883bc316c48bc534ee08319a5e8c620),
GCC13.2, DOUBLE, 64-bit indices, CVODE BDF/dense.
New prefix `.snrt-performance.jDx9Iz/sundials-release`; old prefix untouched.
Release flags explicitly `-O3 -DNDEBUG -ffp-contract=off`; no fast-math,
native-architecture requirement, mixed precision, OpenMP or CUDA change.
Source has a one-line CMake4 compatibility repair in the generated POSIX
timer test; retained as `simulation/snrt/data/sundials_cmake4_timers.patch`.
No mathematical source changes to SUNDIALS.

Reproduce configuration with `cmake -S SOURCE -B NEW_BUILD` and:

```
-DCMAKE_INSTALL_PREFIX=NEW_PREFIX -DCMAKE_BUILD_TYPE=Release
-DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG -ffp-contract=off"
-DCMAKE_POLICY_VERSION_MINIMUM=3.5
-DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF
-DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF
-DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF
-DSUNDIALS_PRECISION=DOUBLE -DSUNDIALS_INDEX_SIZE=64
```

Then `cmake --build NEW_BUILD -j8` and `cmake --install NEW_BUILD`.
On CMake4 apply the timer patch first. Confirm actual generated flags.
Use `SUNDIALS_DIR=NEW_PREFIX` for subsequent lagRamses builds. For identical
binary comparisons prepend NEW_PREFIX/lib64 and NEW_PREFIX/lib to
LD_LIBRARY_PATH and verify all five SUNDIALS libraries with ldd; never
assume a build flag changed the library selected at runtime.

CVODE SHA256: bb1b9236c5fb573dc895681ee131f7d0f505003d544010a534658cc117970e19.
Dense solver SHA256: dc8509557035d5644184ff83f967ad41ed4932edae8c956d6fd73128d90b3764.

## Validation

Neutral 9-state and mixed 6-state before/after printed states match exactly.
Mixed OMP4/serial output comparison passes; native thermochemistry root,
atomic remainder and long photo/CMB regressions pass. Diagnostic probe
reports identical solver counters. This is bounded parity, not proof for
all possible galaxy states.

Uninstrumented interleaved repeats (sum of nine native TIME entries, seconds):
old 1.280557/1.280454/1.278519; release 0.820735/0.817860/0.819667.
Median reduction 36.0%; this confirms the diagnostic wrapper did not create
the observed gain. These are native probes, not independent galaxy runs.

Integration job542047: same binary, effective NML, MPI32xOMP4, and requested
grammar[109-112] as job542031; only SUNDIALS runtime selection changes.
`.snrt-performance.jDx9Iz/release-profile/run.sbatch` verifies loaded library
paths. nstepmax=1, noutput=1, aout=1.1, tout=1e100,
foutput=fbackup=1000000; expected dumps zero; free at launch91TiB.
Job542047 COMPLETED 0:0, allocated173s. RAMSES wall162.228304s;
TotalCPU03:28:05 (3.4681core-hours). Cold wall63.739s; dark summed worker
246.978s; primary72.152s, IR57.811s, tail4.694s; decision collective0.000s
at printed precision. All32 CMB commits, one IR commit, reported
mcons=econs=0, no rejection, no output directories. Cleanup: raw outputs0,
deleted0; logs/inputs/binaries retained.

Same-node comparison against542031:

| Quantity | Old | Release | Reduction |
|---|---:|---:|---:|
| Cold chemistry wall | 193.346s | 63.739s | 67.0% |
| Whole RAMSES wall | 325.976843s | 162.228304s | 50.2% |
| Actual CPU | 08:20:25 | 03:28:05 | 58.4% |

The decision collective also fell31.863s ->0.000s. Whole-wall improvement
must not be attributed solely to optimized arithmetic. Same-node sequential
comparison is stronger than different-node tests but still not an exclusive
node/repeated full-run experiment. Native repeated comparisons and unchanged
state/counters independently support the numerical-library improvement.

## Runtime adoption and next priority

Retained linked binary `.snrt-performance.jDx9Iz/ramses_chimes_release_3d`
uses SUNDIALS_DIR pointing to the new prefix; ldd verifies all five libraries
without requiring an LD_LIBRARY_PATH override. SHA256:
`9a8feb2ac22152cd42a60a68e617eb0c4b76a1e922522ab1ed01595cae87c14c`.
It is built from the same objects/physics configuration with updated build
metadata and library paths. The integrated A/B run deliberately used the
previous identical binary plus explicit runtime selection, not this relink.
Use this release prefix for subsequent CHIMES builds/runs. Old dependency
prefixes are preserved for reproducibility; no shared/global environment was
rewritten. No new namelist, Makefile VPATH, or stream-count changes.

Workspace reuse is not the next priority on measured evidence. If further
CHIMES optimization is needed, investigate repeated RHS/rate-table work and
numerical Jacobian evaluation on the optimized baseline, preserving thermal,
electron and shielding derivatives. Do not launch a full CUDA solver rewrite
or automatically turn this into another audit gate. No commit/push this turn.
