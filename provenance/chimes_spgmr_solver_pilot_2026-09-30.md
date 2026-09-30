# CHIMES dense-vs-SPGMR solver pilot — 2026-09-30

## Scope

Investigate the measured CHIMES cost without changing production chemistry,
species, reaction data, integration interval, or tolerances. The existing
solver remains the dense CVODE BDF path. A separate candidate uses SUNDIALS
SPGMR (no preconditioner, Krylov dimension 20, two restarts) only when the
cell has zero photon ODE groups and `nH_tot <= 0.2 cm^-3`; higher-density and
photon-coupled cases retain dense LU. Candidate source/build/logs are isolated
under `.chimes-spgmr-pilot.jQ3hgB/`.

## Why test an iterative linear solver

The dark 128^3 smoke profile (job 409088) used 157 chemical abundances plus
thermal energy (`nstate=158`). A sampled cell required 232–233 BDF steps,
261–262 ordinary RHS evaluations, 790 finite-difference linear RHS calls,
five Jacobian evaluations, 35–36 linear-solver setups, and about 258–259
nonlinear iterations. The dense finite-difference Jacobian therefore
re-evaluated the expensive chemistry RHS many times; the GPU does not
accelerate this CPU CVODE solve.

On the captured cell, SPGMR reduced linear RHS evaluations from 632 to 281
(55.5%) while BDF steps remained 187 vs. 185 and convergence failures stayed
zero. This is the measured algorithmic saving behind the roughly 48%
single-cell wall-time reduction; ordinary RHS work increased slightly, so
the benefit comes primarily from avoiding full dense Jacobian construction.

## Native same-input evidence on LagEunha

The exact captured dark cell 18447 (nH=0.165361 cm^-3, T=158.854 K,
dt=8.3070e11 s) was replayed five times per solver from the same preserved
input. Dense median integration time was 9.972 ms; the density-limited SPGMR
candidate was 5.153 ms (48.3% lower). One profiled replay reported:

| Metric | Dense | SPGMR |
|---|---:|---:|
| CVODE steps | 187 | 185 |
| ordinary RHS evaluations | 206 | 253 |
| linear RHS evaluations | 632 | 281 |
| Jacobian evaluations | 4 | 0 |
| error-test failures | 1 | 5 |
| convergence failures | 0 | 0 |

The complete 158-value output state was compared: temperature differed by
8.53e-14 K (5.37e-16 relative); all species above 1e-15 abundance matched at
printed precision. The native replay's element and charge acceptance checks
passed for both solvers. A later replay with the candidate's reconciliation
guard enabled measured 10.291 ms dense versus 5.399 ms SPGMR; the temperature
difference remained 8.53e-14 K and the maximum absolute species difference
was 6.97e-20 (zero among species above 1e-15).

The dark mixed regression was extended to a 12-case density/temperature grid:
T=100/1e4 K and nH=1e-4, 0.2, 0.5, 1, 10, and 100 cm^-3. Each case passed
same-solver thread/state and element/charge checks. Dense, adaptive, and
all-SPGMR were each run three times. The median per-case all-SPGMR/dense wall
ratios were 0.744–0.746 at 1e-4, 0.892–0.928 at 0.2, 0.950/1.040 at 0.5,
0.998/1.072 at 1, 1.230–1.286 at 10, and 1.609–2.113 at 100 cm^-3 (the
two values at 0.5 and 1 are T=100/1e4 K). Thus SPGMR's advantage declines
with density and becomes a loss above the crossover; the conservative
nH<=0.2 route is supported by this grid, but behavior between sampled points
is not established.

The aggregate 108-integration profile explains why the low-density cell
speedup does not translate into a broad tile speedup. Dense took 1.267 s in
CVODE, of which RHS evaluations accounted for 0.975 s, linear setup for
0.058 s, and linear solves for 0.112 s. Adaptive took 1.293 s (RHS 1.011 s;
setup 0.054 s; solves 0.105 s). All-SPGMR took 1.812 s and raised RHS work
from 102,501 to 157,509 calls (+54%), despite removing dense setup/solve
work. The instrumented RHS subcosts were primarily coefficient updates and
cooling evaluation (dense: 0.313 s and 0.325 s; all-SPGMR: 0.548 s and
0.515 s). In this fixture, therefore, repeated chemistry RHS evaluation—not
dense factorization alone—is the dominant optimization target. This is a
synthetic 12-case profile, not a cosmological cell-population measurement.

## Algorithm-level diagnosis

Turning RT off removes photon-coupled groups; it does not turn off the
157-species chemistry network or its thermal equation. Every eligible cell
still advances that stiff, coupled system with CVODE BDF. The solver and its
rate/cooling RHS run on CPU; having a GPU allocated does not accelerate this
CVODE work. In the RAMSES dark-material path, OpenMP distributes cells in a
tile, and each cell invokes the CHIMES stage separately. The logged
`CHIMES_s` timer covers that cell loop plus its per-cell dust-composition
preparation; the dust material solve is separately timed. On the exact-patch
first tile (128 cells), those timers were 3.4632 s and 0.000356 s
respectively, so the dust solve itself is not the observed bottleneck.
If the dense-control first-tile time (3.2893 s) repeated across roughly
1,000 tiles assigned to one rank, the arithmetic is about 3,290 s (55 min);
the exact-patch candidate gives about 3,460 s (58 min). This is consistent
with the user's order-of-magnitude concern, but it is only a linear
extrapolation from one tile—not a measured level runtime. Cell density,
fallback frequency, first-tile effects, and rank load balance can change the
mean substantially.

The source path is `CVode BDF -> RHS -> update_rate_coefficients/update_rates
-> update_rate_vector -> calculate_total_cooling_rate/update_cooling_rates`.
The dark path already reuses temperature-only coefficients when the cell's
temperature is unchanged and caches the 1-D cooling-table interpolation for
an unchanged temperature. It cannot reuse abundance-dependent reaction
rates or density-dependent 2-D/4-D cooling interpolations when CVODE perturbs
the species state. The profile agrees: coefficient and cooling work dominate
the measured RHS subcosts, whereas dense matrix setup/solve is a smaller
share. Therefore a useful next optimization should target repeated RHS/table
work or exploit the reaction Jacobian's sparsity with an effective
preconditioner; merely replacing dense LU with unpreconditioned Krylov is not
enough. Any reuse across perturbed states must preserve the exact RHS and
conservation behavior—no interpolation approximation is justified by these
timings alone.

## Integrated candidate failure history and bounded A/B

Job 409105 was a same-input, output-suppressed prefix comparison of the
128^3 no-RT/dust+CHIMES run, using the existing RAMSES executable and only
the isolated candidate `libchimes.so`. It requested one H200 node, 2 GPUs,
4 MPI ranks x 3 OpenMP threads, 64 GiB, and was capped at 360 seconds. The
effective namelist has `nstepmax=1`, `noutput=1`, `aout=1.1`, `tout=1e100`,
and large periodic output intervals. During its first dark material level,
four rank-local cells were rejected by the unchanged strict charged
reconciliation gate with negative electron balances from -3.81e-23 to
-1.08e-19. The candidate prefix exited nonzero and was stopped; no
`output_*` dumps were produced. This is a real integrated admissibility
failure despite the native replay passing, so the SPGMR candidate is not
acceptable without a safe fallback.

The isolated candidate snapshots each cell's pre-solve abundances and
temperature, evaluates the expected net electron abundance using CHIMES'
existing charge bookkeeping after constraint enforcement, and—only if that
budget is negative or CVODE reports failure—retries from the untouched input
using the dense reference path. The bridge's full charged reconciliation
remains unchanged and still controls final acceptance. No tolerance,
conservation threshold, or production source was changed. Native replay and
the dark/mixed conservation suite pass with this guard enabled.
The dense-retry branch was also exercised once deliberately on the captured
live-cell replay (test-only compile flag, absent from the integrated candidate
library). It emitted the retry trace, completed successfully, and its full
158-value state was byte-for-byte identical to the dense baseline. The
forced SPGMR-then-dense replay took 15.88 ms versus 10.29 ms for dense alone,
consistent with paying for both attempts on a rejected cell.

The first A100 retest request (409111) was cancelled before it ran because the
A100 nodes had no fitting GPU slot. A one-GPU H100 request (409131) reached the
first dark-material call, but the candidate library then failed with loader
exit 127: it imported `snrt_chimes_reconcile_charged`, a host executable symbol
not exported for dynamic linking. No chemistry result or output dump was
produced. This was a pilot wiring defect, not a solver or physics rejection.
The candidate was corrected to use the in-library charge balance above, and a
preflight now rejects that unresolved host-only import before launching.
Job 409143 reran the bounded prefix on H100 NVL `syn08` (4 MPI x 3 OpenMP,
64 GiB), using the same namelist SHA256
`7b3bf903b526bb7876e8b44d67fdcf758588e0eab2b07314beaf4009cca534b0` and
RAMSES binary SHA256
`44b9195a1a1b1791d0357212cd222676d529adc14cf21259198667b8b455b39e` as the
dense control. It reached the first dark-material tile and logged a negative
expected-electron balance on each of the 12 worker threads; these are lower
bound retry observations, not a total retry count. The untouched-input dense
fallback cleared them, the first tile reported `status=0`, and no strict
reconciliation rejection was logged. Its tile timing was CHIMES 2.2929 s and
dust 0.000545 s. Dense job 409149 then ran the same input and binary on the
same `syn08` H100 NVL with the same 1-GPU, 4 MPI x 3 OMP allocation. It loaded
the frozen dense library (SHA256
`4e81761da2cdccb136f4788460bed827f4053e14083b5f40c9d177bf9932fa7b`) instead
of candidate library `33df2cc9cd3b5d2564d8b17c715eadf413f2dc72acc8e0dbaaef634f7c796191`.
The dense tile took 3.2893 s CHIMES and 0.000343 s dust, `status=0` in both
cases. The candidate reduced first-tile CHIMES time by 30.3%. This is a
matched first-tile result, not a completed level/coarse-step speedup. Both
360-second caps fired as intended (run status 124, successful batch wrapper),
and neither prefix created `output_*` raw dumps.

The exact final patch build was then tested in job 409159 on the same node,
binary, namelist, and allocation. Its candidate library SHA256 was
`362d88e93262527512efdbdb6a91f40bbb4c21c98d376462cffd747e56bf66cd`; it
resolved cleanly, and the first 128-cell tile completed with `status=0`, no
reconciliation rejection, CHIMES 3.4632 s, and dust 0.000356 s. The matching
dense control was 3.2893 s, so this exact-patch observation is 5.3% slower,
not faster. The run again stopped at its planned 360-second prefix and wrote
no raw output dumps. This contradicts the earlier 30.3% pilot result; the
integrated candidate speed remains unresolved and neither figure is a
multi-tile estimate. Do not treat the earlier pilot timing as validation of
the final patch.

## Reproducible patch and post-patch verification

The isolated implementation is now captured in
`simulation/snrt/data/chimes_spgmr_adaptive_solver.patch`. It is a follow-on
patch for the staged CHIMES source used in this experiment—not pristine
upstream CHIMES. Its required input `src/chimes.c` SHA256 is
`00cc645f24d0f334a3c5b7c11341b404e73cde851c4ea498fd2bc789dd862bc3`; the
patched source SHA256 is
`468fdd0dd953c09497b3d2f13a82bd1867f9b951da8c622ea587be3bc50843dd`.
All nine hunks dry-run and apply, and the applied file matches the measured
candidate source byte-for-byte. Compiling without
`SNRT_CHIMES_ADAPTIVE_SPGMR` retains the dense code path; enabling the macro
requires linking `-lsundials_sunlinsolspgmr` in the separately built CHIMES
library. The project's RAMSES Makefile already links that SUNDIALS solver.
No frozen library or default runtime was replaced.

Both macro-off and macro-on source variants pass GCC 13.2 syntax checks. The
isolated macro-on CHIMES library links against the pinned SUNDIALS 5 and HDF5
stack (candidate library SHA256
`362d88e93262527512efdbdb6a91f40bbb4c21c98d376462cffd747e56bf66cd`) and
passes `ldd -r` without unresolved symbols. An
initial packaging draft also required `rt_update_flux==0`; a live replay
showed that this unnecessarily routed a tested zero-photon state back to
dense. That extra condition was removed, leaving the measured routing rule
`N_spectra==0 && nH_tot<=0.2`. Its 12-case T=100/10^4 K,
nH=10^-4/0.2/0.5/1/10/100 cm^-3 dark mixed regression passes thread-parity
and element/charge-conservation checks at four OpenMP threads. The density
crossover was measured in three repeats per solver. At 0.5 and 1 cm^-3 the
result is temperature-dependent and near parity; at 10 and 100 cm^-3 SPGMR
is slower. This bounds the conservative low-density route but does not claim
that the threshold is globally optimal.
In one paired four-thread pass of the same 32-cell-per-case fixture, adaptive
versus dense times were 7.099 vs 8.721 ms and 6.800 vs 8.717 ms at
nH=10^-4; 38.004 vs 40.993 ms and 39.418 vs 41.058 ms at nH=0.2 (T=100 and
10^4 K respectively). At nH=100 it routes to dense; small single-pass timing
differences there are not significant.
On the captured live-cell replay (cell 18447), the patched library and frozen
dense library both pass. Under the same RHS-cost probe, wrapper time was
7.335 ms for the adaptive solver and 13.423 ms for dense. Their 158-value
states differed by 8.53e-14 K in temperature and at most 6.97e-20 in species;
none of the 11 species above 1e-15 changed at printed precision. This single
replay corroborates the earlier independent profiled result, but is not a
population-level speed estimate. The exact-patch H100 integrated run 409159
also completed its bounded prefix and passed the first tile, but measured
5.3% slower than the dense control; see the integrated comparison above.
Functional admissibility is demonstrated for the sampled prefix, while an
integrated speed benefit is not.

## Decision

Keep this as an isolated, opt-in candidate; do not overwrite the
frozen/production CHIMES library. The native live-cell result shows that
SPGMR can reduce cost in a low-density state, and the 12-case grid supports a
conservative low-density-only route. However, the final-patch integrated
first-tile A/B was 5.3% slower while the earlier logically equivalent pilot
was 30.3% faster. That disagreement is unresolved. More importantly, the
mixed-grid profile identifies repeated chemistry RHS work—especially
coefficient and cooling evaluation—as a larger aggregate target than dense
linear algebra; all-SPGMR increases total RHS calls. The next efficient step
is a bounded same-node, paired multi-tile timing with the exact patch and
frozen dense library, alongside a small attribution check for fallback/RHS
cost. Do not run a full level or change the production default until that
comparison shows a repeatable benefit. No production solver, tolerance, or
physical admission rule was changed.

## Follow-up on LagEunha: exact-patch first-tile cross-check

On 2026-09-30, I attempted the proposed 64-tile collector on LagEunha using
the existing `measure.py`. That collector expects the radiation-transport
`SNRT PERF TILE` record, while this no-RT dark-material path emits only one
`SNRT dark material first_tile` record. I stopped both misconfigured
collection attempts rather than allowing a long one-step run; neither was
used as multi-tile evidence. They produced no `output_*` dumps. One preliminary
run also selected the older SPGMR library by a mislabeled path; although it
emitted a tile timing, it was the wrong arm and is excluded from all
comparisons.

A corrected same-input first-tile cross-check used the same RAMSES binary
(SHA256 `44b9195a1a1b1791d0357212cd222676d529adc14cf21259198667b8b455b39e`),
the same no-RT/CPU namelist (SHA256
`3bfe87bd0dfa8ece55b023879751b46bdaca93b5e44629c75da201a5e21d74e3`), and
four OpenMP threads on LagEunha. `ldd` confirmed the frozen dense library
(`4e81761da2cdccb136f4788460bed827f4053e14083b5f40c9d177bf9932fa7b`) for
the control and the exact SPGMR patch library
(`362d88e93262527512efdbdb6a91f40bbb4c21c98d376462cffd747e56bf66cd`) for
the candidate. First 128-cell tile CHIMES times were 0.40940 s dense and
0.42305 s adaptive, i.e. the candidate was 3.33% slower. Both first tiles
returned `status=0`; the adaptive prefix ran to its planned 180 s cap, and
neither arm produced output dumps. This remains a one-tile observation, not
a tile-population or level-speed estimate. Its direction agrees with the
exact-patch H100 result (+5.3% slower), so there is no evidence to promote
SPGMR beyond its isolated opt-in status.

The next efficiency candidate should target ordinary chemistry RHS cost,
not add more tile instrumentation or widen the SPGMR trial. In particular,
the existing exact temperature-only coefficient and 1-D cooling caches must
be preserved; any Jacobian/preconditioner optimization must be checked
against the dense CVODE result for state, thermal energy, elemental/charge
conservation, and CVODE convergence. No production source or solver default
was changed by this follow-up.

## Exact temperature-only molecular-rate cache pilot

Because RHS coefficient/cooling work dominated the prior profile, an isolated
exact-reuse candidate was also built from the existing CHIMES source snapshot.
It caches the five H2 low-density coefficients, H2 LTE coefficient, and
gas-grain transfer coefficient only while both the per-thread current-rate
buffer identity and gas temperature are exactly unchanged. Abundances,
density-dependent tables, and all state-dependent cooling sums remain live.
The cache resets at the existing per-cell initialization point. No rounded
keys, interpolation approximations, or production defaults were introduced.
The isolated source is `rhs-temperature-cache/src/chimes_cooling.c`
(SHA256 `1ab6a90078d1b267a5f77ac884f553f20b68732ffac6faf6384cffae299e739d`);
its dense-path library is `rhs-temperature-cache/build/libchimes.so`
(SHA256 `6efed4b6884a528a1fa8f3893fb2f3d186370787c97d656bec5c501c2cb14d92`).

On LagEunha, the dense library and candidate were paired five times on the
same captured live cell with one OpenMP thread. All runs reported identical
CVODE work (838 RHS calls, 187 steps, 4 Jacobians) and the same final
temperature. A fresh replay driver compared the full 159-value `REPLAY_STATE`
record byte-for-byte; both hashes were
`16b4ec2ccf2a1090cc9d3baf008a18fef0262875e45d179cadcf604f0761c09`.
Median network time was 10.384876 ms dense versus 10.397105 ms cached
(candidate +0.12%); median cooling sub-time was 1.700640 ms versus 1.695871
ms (candidate -0.28%). These are negligible/noise-level changes for one
captured cell, not a population speedup. Keep the cache as an isolated
experiment only; do not promote it.

Taken together, the exact-patch SPGMR tile cross-check and this exact cache
pilot do not justify a production change. The next high-value candidate is
reducing the expensive finite-difference Jacobian/RHS work without changing
the governing chemistry or cooling equations. Before implementation, inspect
the dark-network dependency graph and verify whether conservative Jacobian
coloring or an effective preconditioner can reduce CVODE RHS calls; retain
the dense solver as the numerical reference. No production source or runtime
default was modified in these pilots.
