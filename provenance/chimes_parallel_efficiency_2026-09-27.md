# CHIMES parallel efficiency — 2026-09-27

Scope: improve the approved gas chemistry implementation without changing
species, reaction data, timestep, solver tolerances or physical admission.
Work in `/gpfs/kjhan/LRD_JWST` (origin kjhan0606/LagRamses).
No direct GPU-node SSH/debugger access; any debugger work belongs on
lageunha by operator instruction. Slurm and shared logs remain available.

## Measured bottleneck

Real-IC job 406263 finished its intentional 1200-second sampling window.
994 tiles of 256 cells: mean preparation 0.962014 s, material 0.006550 s,
publication 0.002238 s. Dark integrator summed worker wall / preparation
wall = 3.88791 on four workers. This is concurrent overlap, not serial
speedup. Additional nested OpenMP teams inside each cell would oversubscribe
this already occupied CPU allocation. GPU M5 speedup is not chemistry
speedup; the implicit chemistry integrator is still CPU CVODE.

## Implemented candidate

Supplemental patch `simulation/snrt/data/chimes_parallel_rates.patch`,
applied AFTER `chimes_native_receiver.patch`:

1. Independent reaction-rate products use OpenMP SIMD. GCC x86_64 runtime
   target clones select AVX2 or the portable default; no binary-wide AVX2
   requirement. Species scatter accumulation is not parallelized unsafely.
2. Skip only exactly-zero reaction rates in the four principal accumulator
   loops. Preserve the order of every nonzero addition; no tiny-rate cutoff.
3. Reuse temperature-only cooling coefficients in a thread-local cache,
   only for a dark system, the same rate buffer and exactly equal T. Reset
   at every initial-rate/cell entry, including allocator-address reuse.
   All abundance weighting and 2D/4D density/electron/HI/HII tables are
   still recalculated. No temperature binning, shared writable cache,
   whole-cooling-value cache, or lost Jacobian dependence.

Compile update_rates.c with `-O3 -ffp-contract=off -fopenmp-simd`, other
cooling code with the existing O2 policy. Keep double precision, no
fast-math, the same release SUNDIALS and the same remaining CHIMES objects.
No ABI change. Existing live/frozen library is NOT overwritten.

Candidate workspace: `.chimes-parallel-20260927.5qb6Qv`.
Baseline library SHA256:
`ae0d11e649780d01c133d5c7b4ce851eec748eee39224bb7d4747860c254bc49`.
Candidate including cooling cache:
`4e81761da2cdccb136f4788460bed827f4053e14083b5f40c9d177bf9932fa7b`.
The compiler report confirms both 16-byte and 32-byte vector loops.

## Native results

Reuse the existing six temperature/density cases, each 32 independent
cells, with existing nucleus/charge checks and a serial reference in every
case. The fixture now obeys OMP_NUM_THREADS rather than forcing four workers.
Each timing includes 32 parallel calls and the one serial reference;
therefore it is not a pure parallel-region scaling measurement.

Job 406266 (SIMD/zero work only) passed all state comparisons and gave
about 6% lower four-worker matrix time. Job 406267 additionally enabled
cooling coefficient reuse. Three-repeat median whole-matrix seconds:

| Workers | Baseline | Candidate | Time reduction |
|---|---:|---:|---:|
| 1 | 3.081390 | 2.803050 | 9.03% |
| 4 | 0.870857 | 0.786691 | 9.67% |

All printed 157-species and temperature states are byte-for-byte identical
between libraries in all cases and repetitions. In-process OMP/serial
state comparisons are exact. This is a bounded regression result, not a
claim for all chemistry states or whole-simulation speedup. The change
reduces per-worker arithmetic work; it does not demonstrate an increase
in the already-high thread utilization of the source-dark pilot.

## Remaining verification for adoption

406269 reuses the rate-reset and long irradiated/CMB interval fixtures,
plus the actual M5 substep-duration neutral probe in both libraries.
406271 is dependent on its successful exit, using the identical frozen
RAMSES binary and namelist as 406263, with only the CHIMES library changed.
Effective input `.chimes-parallel-real128-20260927.z5lH2G/run.nml`:
levelmin=levelmax=7, nstepmax=1, original 128^3 pre-enriched IC.
One A100 on syn101, MPI1 x OMP4, 320 GiB (baseline RSS 186.11 GiB),
600-second launcher limit in a 12-minute allocation. Zero full dumps:
noutput=1, aout=1.1, tout=1e100, foutput=fbackup=1000000,
walltime_hrs=-1. GPFS free space 274 TiB. A timeout is a sample, not a
completed coarse step or completed science run. At initial check 406269
was waiting on QOSMaxGRESPerUser; no extra allocation was launched to
bypass that restriction.

### CPU regressions completed on lageunha

The operator-designated debugging server sees the same verified GPFS
repository and Intel/GCC dependencies. The bounded native script
`verify-lageunha.sh` (120-second ceiling, at most four workers) passed rate
reset and the original 68.65-Myr photo/CMB fixture for both libraries.
Nuclear/charge/photon/energy errors are identical at printed precision:
1.7387478265e-10 / 6.8827512767e-14 / 5.5559888611e-12 / 1.3289452254e-11.
The long photo case takes about 11.5 CPU seconds in either library: no
speedup claimed for this unrelated photo-dominated workload.

Neutral target-duration probe (64 calls, dt=8.158275148690823e11 s),
three-repeat median: baseline 0.720419 s, candidate 0.646861 s (-10.21%).
Both return T=165.399999991399 K and the full requested interval; native
charge/nuclei checks pass. These lageunha timings are compared only with
the same-server baseline, not A100-host timings.

Thus duplicate pending CPU-test job 406269 was cancelled. After actual
regression success, 406271's afterok dependency was cleared; the integrated
A100 comparison remains subject to the ordinary GPU quota/resource queue.
Do not report its result before it runs. Initial script portability faults
(missing explicit HDF5 include path, unavailable rg on lageunha) were
corrected; the final script exits 0 with CHIMES_PARALLEL_REGRESSION_PASS.

## Actual-IC fixed-work OpenMP scaling on lageunha

Operator approved 1/2/4-worker comparison after the integrated memory sample.
Completed sequentially on LagEunha on 2026-09-27, starting 13:18 KST, using
`.chimes-scaling-20260927.hG8d1t/run.sh` and `measure.py`. Reused the frozen
combined executable and optimized CHIMES from
`.integrated-opt-real128-20260927.LhVxyl/inputs.sha256`; no production source
or physics changes for this measurement. Each t1/t2/t4/run.nml has SHA256
3a2017bbc36a46b2f7522ebe02573b42311616b12d93cc23d42efc556698ec79.
Same real 128^3 IC, L7, MPI1, hybrid, RTX 5000 Ada, OMP_PROC_BIND=FALSE,
OMP_DYNAMIC=FALSE; only OMP_NUM_THREADS differs. No CPU affinity tuning.

Comparison uses exactly the first 64 diagnostic tiles, 256 cells each,
with first/last cell indices checked: 16,384 identical initial cell inputs
and the same integration interval. A bounded external launcher reads existing
diagnostics and signals only its own launched process group after collecting
the prefix. Up to a small tail may run before signal delivery; it is excluded
from tile timing. No complete step/state parity claim is made. Individual
cell-result parity was covered by earlier native regressions, not measured
by these timing-only logs. A 900s ceiling per run was not reached.

| OpenMP workers | Sum of 64 tile wall times (s) | Speedup | Efficiency |
|---|---:|---:|---:|
| 1 | 164.474 | 1.0000 | 100% |
| 2 | 82.567 | 1.9920 | 99.60% |
| 4 | 41.751 | 3.9394 | 98.49% |

Preparation (chemistry plus per-cell setup) totals are
164.044 / 82.164 / 41.287 s: four-worker efficiency 99.33%.
Summed dark chemistry worker wall times are 160.204 / 160.116 / 160.211 s,
consistent with negligible work-time inflation over this 1–4-worker range.
These worker wall sums are not CPU counters. The benchmark is one pass
per worker count, using existing timers rounded to milliseconds, not a
multi-repeat confidence interval or evidence for large-rank scaling.

Launch-to-termination wall times are 387.859 / 247.653 / 179.580 s;
child user+system CPU seconds are 383.360 / 386.061 / 395.352.
The short whole-process comparison is only 2.16x faster at four workers
because initialization and whole-grid setup dominate this short prefix;
it is not a whole-coarse-step speedup measurement. First completed tile
arrived at 226.258 / 166.196 / 138.166 s, respectively.

All three stop traces report SIGTERM from the intentional prefix stop.
No prior runtime/physics rejection was found. No output_* dumps were produced;
there are no raw simulation snapshots to clean up. Preserve inputs, launchers
and logs. Evidence: scaling.log and t{1,2,4}/live.log in the above directory.

Conclusion: cell preparation has near-ideal 1–4-thread scaling here. Low
OpenMP utilization is not the dominant cause of this source-dark workload's
cost. Initialization/full-field setup and per-cell integration cost remain
separate optimization targets. This does not establish GPU speedup (zero
radiation short-circuits angular work), later source-active behavior, or
efficiency beyond four workers; no extra benchmark bundle is implied.

## Next efficiency focus: level chemistry preparation

Follow-up on 2026-09-27. The initialization log itself reaches `Starting
time integration` after about 17 seconds on LagEunha; the much longer delay
to the first M5 tile also includes hydro/gravity work and level chemistry
preparation. Do not label that whole interval as input initialization.

`chimes_prepare_level` was serial over independent active leaf cells. It is
called at the hydro boundary and dust level update. The new implementation
collects the same ordered leaf list, evaluates each cell into a separate
staging column with OpenMP static work sharing, reduces failures, and keeps
the MPI rejection collective before any publication. After admission, it
copies the independent staged columns back in parallel. Each cell appears
once in the leaf list. The C reconciliation uses local work arrays and reads
the already initialized shared tables; `chimes_live_initialize` and the CMB
clock update still occur before the team. No chemistry equations, solver
tolerances, element/charge tests, MPI synchronization, physical controls or
output contract changed.

The existing fixed-first-64-tile (16,384-cell) real-IC runner was reused
on the same LagEunha host, RTX 5000 Ada, MPI1/OMP4, identical namelist,
optimized CHIMES library and physical inputs. Three frozen binaries, one
run each: baseline `aedda69ea0cb7f60c51df6a30378d465c4b5036bb3ce1bd7510d0260b9461165`,
cell-parallel only `30c31f6641022289cd6e142484ac8ec586f2b15b32c53be3e2f6c87df6d7c833`,
cell-parallel plus parallel publication
`98b50295470730eb7a30afede292367cad31d3ae295822de0e3537018a2a0c56`.

| Variant | First tile (s) | 64 tiles from launch (s) | Sum of tile wall (s) |
|---|---:|---:|---:|
| Serial level preparation | 138.166 | 179.580 | 41.751 |
| Parallel cell calculation | 125.165 | 166.579 | 41.730 |
| Parallel cell calculation + publication | 123.157 | 164.569 | 42.100 |

The final short prefix is 15.011 s or 8.36% shorter than the baseline. Tile
chemistry work is effectively unchanged; its summed worker wall was
160.211 vs 161.685 s. These are single runs and the external launch timing
contains all pre-tile work, so the 15.011 s cannot be assigned exclusively
to one subroutine or extrapolated to a whole coarse step. The final
executable is a frozen candidate, not an automatically deployed production
binary. Source SHA256:
`bae940c9e0a2e4662e850e7eb7c9984a2d4ac88ebc4d3fb7747037191aaef892`.

All three matched-prefix runs reached M5 without a reported admission or
physics error; the pre-M5 printed mass/energy totals and zero NaN checks
matched. The runner intentionally sends SIGTERM after the 64th diagnostic;
that traceback is expected. This comparison did not complete a coarse step
or perform a bitwise chemical-state comparison. No `output_*` raw dump was
generated. Evidence and effective inputs are retained at
`.chimes-scaling-20260927.hG8d1t/{scaling.log,candidate4.log,candidate2-4.log}`
and the corresponding `t4/live.log` in each run directory. Build source
selection can be checked in `build-next2.log`; the build used the
`patch/lagRamses` VPATH. Existing native CHIMES serial/OpenMP state parity
tests apply to the C receiver; the level-wide rearrangement relies on
independent staged cells and collective admission. Source-active and
completed-step verification remain to be established before a production
speedup claim.

### Completed-step state test of level-preparation parallelism

The same frozen serial/parallel binaries and CHIMES library were run on
LagEunha in a bounded 4^3, MPI2 x OMP2, four-coarse-step, source-dark
CHIMES+dust+RT control. The old AGN fixture was not silently reused: its
first attempt failed current dust/AGN admission, and disabling stellar
physics alone failed the stellar-source admission. The final input keeps
stellar physics enabled but sets `n_star=1e30`, with no sink/AGN and thus
no live source. This is a numerical regression fixture, not physical
calibration. Input SHA256:
`8e9bd47b72d5ee38ce7eed8f60c406bcf2c7a350a0b310d9e721053c55b983a3`.

Both binaries reported `Run completed`, four coarse steps, four SNRT
transaction/closure passes, and no reported runtime/physics errors. HDF5
snapshots at steps 2 and 4 contain the same 438 physical datasets in each
binary; every dataset under amr/gravity/hydro/particles/sinks/snrt is
bitwise equal, with no missing or extra names. The comparison used
`.chimes-level-parity-20260927.zn9Nyj/compare.py`; inputs and logs are in
that directory. Runtime was 11.395 s serial versus 9.406 s parallel, a
single tiny run dominated by overhead, so it is not a production speedup
claim. This regression does not validate source-active radiation/chemistry.

The source-active follow-up reused the historical 4^3 existing-sink/Bondi
fixture with the now-required explicit `SNRT_AGN_MODEL=partition_reference_v1`.
Both frozen binaries completed four coarse steps, with 32/64/96 active
sources after the first source-free step, nonzero CHIMES gas and grain
absorption, and four SNRT closure passes. At steps 2 and 4, all 485 physical
dataset names matched. Of these, 482 datasets including *every* hydro,
chemistry, dust, radiation, gravity, AMR and particle field were bitwise
identical. Only `sinks/jsink_{1,2,3}` differed, by at most 2.17e-34 at
step 2 and 3.94e-31 at step 4, roughly roundoff-relative to their
1e-20--1e-15 values. An unchanged-baseline repeat also varied in the same
three sink-angular-momentum fields (up to 7.89e-31); therefore full-state
bitwise equality is not a valid admission criterion for this fixture.
This does not prove every sink-ordering effect is unrelated to the patch,
but the modified hydro/CHIMES/RT state itself matched exactly. Whole
tiny-run wall times 27.066 s serial, 27.402 s candidate, and 27.513 s
unchanged-baseline repeat show no measurable end-to-end speedup in this
source-active case. No production performance claim follows. The logs,
comparison helper and exact effective inputs are retained in
`.chimes-level-parity-20260927.zn9Nyj`. After the evaluation, the ten
exact HDF5 snapshot files (about 1.01 GB total) were deleted; run logs,
inputs, binaries, scripts and small output metadata remain. No other
simulation data were removed.

### Exact-zero molecular line-cooling candidate (not adopted)

`simulation/snrt/data/chimes_zero_molecular_cooling.patch` is an optional
incremental patch **after** `chimes_parallel_rates.patch`. It skips H2
rovibrational interpolation only for exactly zero H2, CO/H2O line tables
only when their H2 or emitter abundance is exactly zero, and OH rotational
cooling only for exactly zero OH. Nonzero arithmetic, H2 formation and
dissociation, other heating/cooling, and solver controls are unchanged.
The frozen optimized library remains the default; this patch was tested
only in `.chimes-zero-cooling-20260927.wC7hvr` on LagEunha.

Baseline/candidate library SHA256 respectively:
`4e81761da2cdccb136f4788460bed827f4053e14083b5f40c9d177bf9932fa7b`
and `038fcb591d7cc6ed397800c75519b976c50f95fca6fed8e0f170b33e1948e642`.
The identical native binary exercised the standard six-state dark matrix,
an exact captured live cell, a molecular-rich state with H2/CO/H2O/OH all
nonzero, and OH with H2O/CO absent. All returned temperatures and 157-species
states matched byte-for-byte, including serial/OpenMP comparisons. The
three-repeat median matrix times on the same server were 2.14805 vs
2.11662 s (one worker, 1.46% lower) and 0.588995 vs 0.587273 s (four
workers, 0.29% lower). These are small native-only effects.

The same frozen 4^3 MPI2/OMP2 source-active binary/input then completed
four coarse steps with either library and all four SNRT closure checks.
At steps 2 and 4, 485 physical dataset names matched. All 482 non-sink-
angular-momentum datasets, including hydro, chemistry, dust and radiation,
were bitwise identical. Only `sinks/jsink_{1,2,3}` varied by at most
3.70e-32 (step 4); the unchanged-baseline repeat above already varied
in these fields by up to 7.89e-31. The baseline/candidate run walls were
27.255/27.058 s, within the existing same-binary run variation. The sum
of printed cold-chemistry wall times was 13.434/13.439 s: no integrated
chemistry speedup was observed. The four evaluated 97 MB raw HDF5 files
were removed after comparison; logs, effective inputs, patch and candidate
library remain. Do not promote this candidate on performance grounds.

## Routine-cost triage for the next optimization (2026-09-27)

Use three distinct denominators; do not add nested rank-local or summed-worker
timers into whole-job elapsed time.

1. Real 128^3 IC, first 64 rank-1 M5 tiles (16,384 cells), MPI1/OMP4,
   source-dark, frozen pre-level-parallel binary: 41.751 s tile wall =
   41.287 s cell preparation, 0.320 s material/IR stage, 0.144 s commit.
   The CHIMES dark *summed worker wall* was 160.211 s; divided by four
   workers, this is about 40.053 s, consistent with ~97% of preparation.
   OMP1/2/4 tile scaling was 1/1.992/3.939, so adding another nested
   OpenMP team is not the near-term answer. The newer level-preparation
   binary reduced the launch-to-64-tile prefix by 15.011 s, but its
   41.696 s preparation/0.275 s stage/0.129 s commit remain essentially
   the same per-tile work. The ~123 s to its first tile includes startup,
   gravity/hydro and level preparation, not one identified routine.
2. Bounded 4^3 source-active control, rank 1 over four completed coarse
   steps, MPI2/OMP2: SNRT coupling 16.906 s, transport 2.763 s, source
   0.038 s, NLTE 0.004 s. The *nested* cold chemistry wall within coupling
   was 13.598 s; photo/dark summed worker walls were 25.910/0.913 s.
   The *nested* IR halo/solve/scatter totals were 0.0904/0.3163/0.0153 s.
   The ordinary RAMSES phase timer reports broad cooling at 16.995 s
   (86.9% of its 19.566 s rank-averaged timer total), consistent with
   chemistry/coupling dominance but not directly additive to rank-1 SNRT
   intervals. This tiny source-active case is a path check, not 128^3 scaling.
3. One exact live source-dark cell (captured cell 18447) replayed on LagEunha
   with the *current optimized* CHIMES SHA256
   `4e81761da2cdccb136f4788460bed827f4053e14083b5f40c9d177bf9932fa7b`
   and Release SUNDIALS. Network 10.268 ms, CVODE 10.108 ms, RHS 4.972 ms,
   dense LU 1.913 ms, dense solve 1.877 ms. It required 187 accepted steps,
   838 RHS calls (206 ordinary, 632 numerical-Jacobian), four Jacobians,
   29 LU setups, 203 solves, one error-test failure, zero convergence
   failures. A diagnostic wrapper measured coefficient/rate/vector/cooling
   calls at 0.776/0.964/0.781/1.720 ms respectively. These subroutine
   intervals can include initialization or final calls outside CVODE RHS;
   they are not an exact disjoint decomposition of 4.972 ms. Instrumented
   and plain replay walls were 10.335 and 10.475 ms, with the same accepted
   state. Evidence: `.chimes-scaling-20260927.hG8d1t/current-replay*.log`.

Cost/effort judgment: first test exact-zero molecular-emitter skips and
other local, algebraically exact work elimination inside
`calculate_total_cooling_rate`/`update_cooling_rates`, then the remaining
`update_rates`/`update_rate_vector` loops. Temperature-only cooling reuse,
temperature-only coefficient reuse, exact-zero reaction skips, SIMD and
optimized SUNDIALS are already installed: do not count them as future work.
Do not skip a nonzero species, loosen tolerances, or cache a density/species
dependent cooling value. The cooling routine is the largest measured RHS
subroutine in this one real cell, but even removing it entirely would not
remove CVODE's dense linear cost; require repeated native and one integrated
same-input A/B before claiming a gain. The higher-upside, higher-effort
option is to reduce the 632 numerical-Jacobian RHS evaluations per cell
(sparsity/coloring or a validated alternative linear solve); this needs a
separate bounded correctness/performance pilot, not an immediate solver
replacement. M5 transport, IR kernels and stream-count tuning are low ROI
for the measured source-dark 128^3 tile; source-active transport is a
secondary target. A single minimal phase timing of the still-unassigned
pre-first-tile 128^3 interval should precede optimization of hydro or
level setup. No production physics/code was changed by this triage and
no raw simulation output was generated.
