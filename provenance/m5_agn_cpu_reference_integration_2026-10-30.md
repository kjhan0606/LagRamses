# M5 AGN-only CPU integration evidence — 2026-10-30

## Result

Slurm job `1244059` completed successfully on `grammar029` (`normal`): 1 MPI
rank, 8 OpenMP threads, CPU CHIMES, 16 GiB requested, 1,368,528 KiB peak RSS,
55:04 elapsed. The runner reported
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS` after three coarse steps, with
M5 order 5, live two-size dust, CHIMES, an actual accreting BH sink, and stellar
SED disabled to isolate the AGN radiation source. The profile uses reference
controls and is explicitly nonproduction.

The effective namelist is retained at
`.m5-source-dust-chimes-1244059/effective.nml`; input and executable digests
are in `.m5-source-dust-chimes-1244059/inputs.sha256`. The CPU executable was
built by job `1243944` from source commit
`4e4c323450ef01baeb8848011d286b58ae2a7628`; its build manifest is in
`.m5-ramses-momentum-build-1243944/build-result.txt`.

## Coupling evidence

- The AGN sidecar records positive BH accretion and positive thermal luminosity
  at coarse step 2, satisfying the source-attribution condition. No stellar
  radiation population was loaded.
- The runtime reached `Run completed`; the wrapper result records 16,896
  momentum-assertion cells and `raw_outputs=0`.
- At the source-active M5 commit, photon-number residual was
  `1.79e-17` and radiation-energy residual was `9.98e-18` in the printed
  ledger convention. At the following commit, the energy residual was
  `-5.62e-18`. No M5 accounting rejection occurred.
- The run also records `Main step=3`; the output schedule remained in the
  future and no `output_NNNNN` directories were created.

## Runtime and limits

Total reported runtime was 3,291.5 s. RAMSES attributed 3,286.2 s (99.8%) to
its `cooling` bucket. This bucket includes the coupled M5/CHIMES path and is
not a standalone CHIMES timing. It identifies that coupled path as the next
runtime comparison target.

The three level-3, order-5 M5 commits took 449.1 s (16 substeps), 2,361.0 s
(16 substeps, source-active), and 437.2 s (1 substep), totaling 3,247.3 s or
98.7% of the run. For the source-active 16-substep call, existing tile timers
sum to 848.1 s of material preparation and 499.0 s of staged commit work;
about 1,013.8 s remains outside those tile intervals and is not split further
by current timers. The summed per-cell photo durations (5,942.2 s) overlap
across OpenMP workers and are not wall time. These measurements identify the
M5+material path as the target but do not justify attributing all its cost to
CHIMES.

The printed RAMSES `econs=9.98e-1` is computed from kinetic energy, potential
energy, and the integrated potential-work term in `amr/update_time.f90`; it
does not include the gas internal, dust, radiation, or explicit feedback-source
terms. It is therefore not a closed-system energy audit for this run. The M5
photon/radiation ledger closes to the residuals above, but a total
gas+radiation+dust+AGN energy budget was not measured here.

This establishes CPU-only coupled execution and AGN source attribution under
the declared reference profile. It does not establish production physical
calibration or GPU/hybrid runtime behavior. H100 CUDA build job `441960`
passed at source commit `ca40f0e`; its dependent runtime job `442061` was
superseded before allocation. A build of current commit `598a302` and its
dependent hybrid runtime job are now queued as `442109` and `442110`;
`442109` is pending `QOSGrpGRES` capacity.

## M5 AMR slope OpenMP trial — 2026-10-30

The M5 slope/sample cell loops were parallelized in commit
`598a302b0303c0d763d4daa597b507b175f3df5e`. CPU build job `1244335`
passed with `CHIMES_CUDA=0`. The same three-step AGN-reference profile then
completed as job `1244378` on `grammar017`: 1 MPI rank, 8 OpenMP threads,
16 GiB requested, 1,344,664 KiB peak RSS, and 53:44 Slurm elapsed (the RAMSES
timer reports 3,213.772 s). Its result sentinel is
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS`, with `raw_outputs=0`.

Against job `1244059` on `grammar029` (1 MPI rank, 8 OpenMP threads, 16 GiB,
55:04 elapsed), the three M5 commits were:

| M5 operation | Baseline (s) | OpenMP trial (s) | Change |
|---|---:|---:|---:|
| 16 substeps, source-free | 449.093 | 416.516 | -7.3% |
| 16 substeps, source-active | 2,360.983 | 2,310.713 | -2.1% |
| 1 substep | 437.222 | 442.189 | +1.1% |
| Total M5 | 3,247.298 | 3,169.418 | -2.4% |

Full RAMSES elapsed changed from 3,291.5 s to 3,213.772 s (-2.36%), while
accounted CPU time increased from 4:22:04 to 4:24:42 (+1.0%). The runs used
different compute nodes, so this is not a controlled same-node speedup; the
measured full-run difference is below the project's useful-improvement bar
and is not claimed as a confirmed performance gain. The source-active photon,
radiation-energy and M5 dispatch counts match exactly between runs. In the
trial, the dispatch counter records 51,019,776 CPU face operations and zero
GPU face operations during the active 16-substep commit.

The existing `SNRT PERF TILE stage_s` timer measures
`snrt_dust_live_moment_tile`, not `sample_slopes`; it cannot isolate this code
change. This trial verifies compilation and coupled CPU execution only. It
does not qualify GPU/hybrid behavior, production physical inputs, or
production-scale scaling. Keep the small timing change provisional and direct
the next performance work toward the dominant M5 transport path, using a
representative CPU/GPU allocation when available.

## Reconciliation of the earlier M_N band-split failures

CPU integrated job `1222367` and H100 integrated job `436882` both used
executables built from source commit `00f7a4f`, before commits `0b69763`
(`Fix subnormal M5 spectral endpoint split`) and `e63fc6d`
(`Project subnormal M5 band packets`). Their group-4 rejection logged a
negative endpoint weight of `-1.976262583364986E-323` and a zero opposite
weight, consistent with the subnormal endpoint case covered by those changes.

The later CPU integration job `1244378` built from `598a302` and completed the
same three-step AGN-reference profile with the pass sentinel above. Its
effective namelist, sink IC, and runner hashes match exactly the attempted
revalidation job `1244569` (`6c8e0a0be400110576835fdc218c0124f8394f7670b93629590d56db4c5772fa`,
`971d2f0f2b0c4731fd2b9801ebc5c363dd9cd51268660eb3cf801e7a34d88a54`, and
`df8a0f6dd6b7863bf67852c344e06f9cd6568443b74b511487a42b38bb4a2f0c`). The
revalidation was canceled after 39:21 to avoid repeating an identical
completed run; its partial log is retained. Thus the earlier failures are not
evidence of failure in the current `598a302` CPU path. Current-source hybrid
execution remains unverified while job `442109` waits for GPU capacity.

## Partial CHIMES photo-cost details from job 1244569

The partial `SNRT PERF` rows separate the slow material path. Across observed
256-cell tiles, `chimes_photo_s` is 201–227 summed per-cell elapsed seconds;
`chimes_dark_s` is 1.3–1.4 summed per-cell seconds. Photo-cell maxima are
2.0–4.5 seconds, with 42–69 cells over one second and none over five seconds.
These are summed cell-seconds, not tile wall time. The enclosing 512-cell
preparation loop takes about 58–60 wall seconds, while the M5 dust/IR tile
stage takes about 12–15 wall seconds per 256-cell tile. This remains partial
evidence from the canceled revalidation, not a completed run benchmark.

The photo timer surrounds `chimes_band_photo_molecular_*`; its driver uses
OpenMP `schedule(dynamic,1)`. The called `snrt_chimes_photo.cpp` solver creates
a per-cell CVODE BDF instance with a serial N_Vector and SPGMR (`maxl=30`),
then advances with `CV_ONE_STEP` and `rtol=1e-10`. Its state includes 157
chemical species, up to 1,152 active photon-node optical depths, and the
conservation ledgers (plus optional phase counters). The current timer does
not expose accepted-step, RHS, or linear-iteration counts. Commit `c116941`
adds those counts to the existing opt-in photo performance row. CPU-only
build job `1244657` completed on grammar's `normal` queue in 4:16 with
`CHIMES_CUDA=0`, 26,834,488 KiB peak RSS, and the
`M5_MOMENTUM_CLOSURE_BUILD_PASS` sentinel. Its executable is
`.m5-ramses-momentum-build-1244657/bin/ramses_m5_ramses_momentum_closure_v13d`
(SHA256 `a4c15b4da13c9e00b0a1e140e72479a69a77f38605d2235c297da05f55310140`);
`ldd` shows no CUDA runtime dependency. The separate syntax GPU build/run
chain `442109` → `442110` remains pending (`QOSGrpGRES` and dependency,
respectively); no replacement jobs were submitted.

## Current syntax hybrid execution chain — 2026-10-30 15:48 KST

The latest scheduler query confirms the submission contract for the current
source revalidation. Build job `442109` requests one H100, 4 CPUs, 32 GiB and
30 minutes; its submission does not override `M5_BUILD_CUDA`, so the runner
uses its default `1`. The build script snapshots `git rev-parse HEAD` when it
starts; the actual source commit must therefore be taken from its
`build-result.txt`, not inferred from the submission-time branch head.

Dependent job `442110` requests one H100, 8 CPUs, 64 GiB and two hours, with
`afterok:442109`. Its recorded `SubmitLine` explicitly exports
`M5_CHIMES_BACKEND=hybrid`, `M5_CHIMES_CELL_SCHEDULER=level_queue`, and the
job-specific binary path under `.m5-ramses-momentum-build-442109`. The runner
profile was not overridden, so it uses `agn_reference`: active accreting-BH
source, stellar SED disabled, nonproduction reference controls, and no raw
RAMSES outputs expected. At the 15:48 KST query, `442109` was still pending
`QOSGrpGRES` and `442110` was pending on its dependency. No duplicate jobs
were submitted.

Correction from the later QoS/partition inspection below: `normal` has
`MaxTRESPerUser=gres/gpu=8` and no configured `GrpTRES`; it is not an
account-level cap. The H100 partition separately applies QoS `part_h100` with
`GrpTRES=gres/gpu=8`.

Runtime acceptance remains the existing bounded runner contract: successful
M5/CHIMES/dust/source checks, a `CHIMES_LEVEL_QUEUE` report with all cells
completed and at least one GPU broker, the pass sentinel, and no raw output
dumps. Per-cell CVODE performance rows remain intentionally disabled for the
level-queue continuation path; this run verifies hybrid execution, not a
per-cell CVODE timing profile or production calibration.

## Current-source CPU revalidation and trace-recovery candidate — 2026-10-30

CPU reference job `1244668` completed on `grammar029` (`normal`) in 53:37:
1 MPI rank, 8 OpenMP threads, 16 GiB requested, 4:26:29 accounted CPU and
1,326,388 KiB peak RSS. It used the CPU executable built by `1244657` from
`c116941e7ac1774a61a8e7be80b0ad6c5e94a6f5`, SHA256
`a4c15b4da13c9e00b0a1e140e72479a69a77f38605d2235c297da05f55310140`.
Between that source commit and the current `de49876`, the only tracked changes
are provenance and the opt-in build script; the physical source is unchanged.

The runner emitted
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS`, with M5 order 5, live dust,
CPU CHIMES and a nonzero AGN source ledger. The AGN sidecar records the expected
sequence: inactive at coarse steps 0–1, thermally luminous accretion at step 2
(`dM_BH=1.0104259799605118e-8` code units,
`L_bol=3.0561937668561012e40 erg/s`), then jet mode at step 3. The active M5
photon-ledger residual is `1.7888835600560e-17`; all 16,896 reported
momentum-assertion cells passed. The runner found no raw outputs and no positive
`NaN_CHK` state counters.

The three level-3, order-5 M5 commits took 401.936, 2,315.817 and 447.083 s
(16, 16 and 1 substeps), totaling 3,164.836 s of the 3,217 s Slurm elapsed.
Across 66 tile records, the existing summed per-cell CHIMES timers reported
6,195.181 s for spectral photo solves and 57.501 s for dark solves; these cell
times overlap across OpenMP workers and are not wall time. Photo rows counted
444,629 CVODE steps, 590,818 RHS calls and 411,048 linear iterations; 9 cells
exceeded 5 s and the maximum was 10.996 s. This localizes CHIMES solver work
to the spectral photo path, while the full-run wall time remains dominated by
the coupled M5/material commits. The dark-trace recovery candidate is a
robustness change, not a speedup for this dominant photo path.

The non-cosmological smoke profile has two `SFRD=NaN` diagnostic lines when
`SFR=0`; this is outside the RAMSES state NaN checks and did not fail the
runner. Keep it classified as an undefined SFR-density diagnostic for this
profile, not as evidence of production SFRD validation.

Trace library build `1244720` completed on `grammar` in 4:17 with
`M5_BUILD_CUDA=0` and `M5_CHIMES_LIBRARY_MODE=trace`. Its build record pins
CHIMES base commit `a58e5c0311993b51abc63d84fff2958a0104f6d6`, the six-file
source-diff SHA256 `e33eac86687d16ed5894c55d4d3f0125771608df27cc5cece66a00a574d55dd9`,
and trace-patch SHA256
`dd625057b6f58549646c69a1ac78e400a8dab5e540ed079216efba6e867035c4`. Its
exported symbol list exactly matches the reference library and `ldd` resolves
all dependencies. Candidate `libchimes.so` SHA256 is
`e59e7d5365681631186bb0dad4e95594c29c4561bac668dcb5e6c716197cc3e2`; the
adjacent reference copy is preserved and hashes to
`64ad78dbd531621904a9fc00dea2c81d8068198fa245fff659bede9ca4d24ab2`.
The production/default library selection remains `reference`.

## Integrated trace run closeout — 2026-10-30

CPU job `1244750` completed on `grammar024`/`normal` with exit code 0, 1 MPI
rank, 8 OpenMP threads, 53:54 Slurm elapsed, and 1,487,500 KiB peak RSS. The
RAMSES log reached `Run completed`; `result.txt` records `raw_outputs=0`, and
no `output_00001` or `output_00002` directory was created. Its effective
namelist and sink IC hashes exactly match reference job `1244668`. Relevant
runtime source under `patch/lagRamses` and `simulation/snrt` is identical
between the paired builds; the only source-tree difference in those paths is
the later build script.

The run used trace library build `1244720`: its `build-result.txt` says
`CHIMES_LIBRARY_MODE=trace`, `ldd` resolves the executable to that build's
`chimes-install/libchimes.so`, and `inputs.sha256` pins the library to
`e59e7d5365681631186bb0dad4e95594c29c4561bac668dcb5e6c716197cc3e2`. This
matches the trace-library provenance above. The old runner emitted the fixed
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS` label even for this trace
library; a follow-up makes the runner read the build mode, include
`build-result.txt` in the input hashes, and label future results as
`..._TRACE_PASS` or `..._REFERENCE_PASS` accordingly.

The reference and trace runs each produced three level-3, order-5 M5 commits
with 16, 16, and 1 substeps. Their commit times were 401.936/2,315.817/447.083
s and 407.442/2,323.601/446.971 s, respectively. M5 time increased by 13.178 s
(0.42%); total RAMSES elapsed increased from 3,209.604 to 3,222.777 s
(0.41%). These runs used separate normal-queue nodes, so this is a runtime
compatibility comparison, not a controlled performance result or speedup.

At the source-active commit, total emitted photon number agrees to printed
precision (`8.5315194363883e-3`). The trace/reference differences in absorbed
number, gas heat, dust heat, and peak `xHII` are about `1.5e-10` or smaller
relative (peak `xHII` differs by `2.3e-13` absolute). Printed number and energy
closure residuals remain at or below about `1e-16`, with no M5 transaction
rejection. Both logs report the same RAMSES mass/energy conservation values at
the corresponding coarse step.

Across 66 photo rows, the trace run counted 443,641 CVODE steps, 589,844 RHS
calls, and 409,899 linear iterations, versus 444,629, 590,818, and 411,048 in
the reference run. Summed per-cell timers were 6,199.089 photo seconds and
61.105 dark seconds for trace, versus 6,195.181 and 57.501 seconds for
reference; these overlap across OpenMP workers and are not wall time. The
small count/timing differences do not qualify a speed improvement.

This confirms that the trace candidate completes the integrated source-active
M5 + live-dust + CHIMES reference-control profile without a numerical or
runtime regression at the measured precision. It does not prove the
`CVodeReInit` retry branch was exercised: the runtime reporter exposes RHS
counts and worker mask, not retry counts. Keep that recovery branch
unqualified until direct branch evidence exists. The profile remains
`NONPRODUCTION`, and the default CHIMES library remains `reference`.

## GPU queue correction and A100 reroute — 2026-10-30

The H100 pending reason was traced to the partition QoS, not the user's
`normal` per-user limit. `sacctmgr` reports `normal` with
`MaxTRESPerUser=gres/gpu=8`, while the H100 partition uses `part_h100` with
`GrpTRES=gres/gpu=8`. A read-only running-job query counted eight GPUs in the
H100 partition QoS group. Thus build `442109` correctly remained pending with
`QOSGrpGRES`; no H100 job had started.

The A100 partition has no partition-level QoS cap. Its two nodes each expose
eight A100 GPUs; the live Slurm allocation query showed seven in use on each
node, leaving one on `syn102`. That node also reported 15 idle CPUs and ample
memory for this job. The Makefile's CUDA targets include `sm_80`, so the A100
build does not require an architecture change.

Replacement chain: build `442271` (A100, 1 GPU, 4 CPUs, 32 GiB) and dependent
hybrid run `442272` (A100, 1 GPU, 8 CPUs, 64 GiB, `agn_reference`,
`level_queue`). At submission, the build was pending `Priority` and the run
pending on its build dependency. The previous H100 jobs `442109` and `442110`
were canceled before start to avoid a duplicate chain.

Lageunha's RTX 5000 Ada has 32 GiB total, with 26 GiB already allocated and
zero reported GPU utilization during the check; no manual run was launched
there because the remaining memory is not a safe basis for the integrated
test. CPU trace run `1244750` remained active at the last check; its log had
advanced through multiple M5 tiles and photo-CVODE records, but no result
sentinel was present yet.

## CPU trace and A100 hybrid integrated reruns completed — 2026-10-30

The current-source CPU trace rerun `1244750` completed on grammar `normal`
(1 MPI rank, 8 OpenMP threads) in 3,222.757 s (53:54), with 1,487,500 KiB peak
RSS and 4:27:12 accounted CPU. It emitted
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS`: M5 order 5, live dust, a
nonzero AGN source ledger, 16,896 passing momentum-assertion cells, and no raw
RAMSES outputs.

A100 hybrid job `442272` on syntax `a100` (1 MPI rank, 8 OpenMP threads, one
A100) completed with exit code 0 in 3,036.995 s (50:48), with 2,075,020 KiB
peak RSS. It emitted the same integration sentinel with
`CHIMES_backend=hybrid`. Its `CHIMES_RHS_COUNTS` report is
`gpu_calls=759329`, `cpu_fallback_calls=26501`, `gpu_batches=17117`,
`device_errors=0`, `warm_table_uploads=2`, `peak_inflight_batches=2`, and
`streams_used=1`. The per-tile level queue completed all 512 cells with one
GPU broker (`gpu_worker_mask=0x1`) and CPU workers. The run therefore exercised
both CPU and GPU lanes in the integrated material path.

The M5 dispatch counters, summed over the three level-3 order-5 commits, were:

| Run | Closure CPU | Closure GPU | Face CPU | Face GPU |
| --- | ---: | ---: | ---: | ---: |
| CPU reference `1244668` | 5,987,328 | 0 | 88,547,328 | 0 |
| A100 hybrid `442272` | 2,327,856 | 3,659,472 | 88,547,328 | 0 |

The photon/energy ledger residuals remained small: the A100 active-source
step reported `1.1246e-17` photon residual and `1.2270e-16` energy residual;
the final step reported `2.4627e-18` and `2.0618e-16`, respectively. CPU and
A100 AGN sidecars cover coarse steps 0–3 and follow the same state sequence:
inactive, inactive, thermal, jet. BH mass and radiative efficiency match to
roundoff. At the final jet state, inflow and bolometric luminosity differ by
`8.632e-11` relative; among the selected sidecar fields, the largest relative
difference is `2.433e-7` in the small gas-angular-momentum norm. Both
efficiency contracts are valid at all four states.

Using comparable RAMSES process `Total elapsed time` values, CPU reference job
`1244668` took 3,209.604 s and A100 hybrid took 3,036.995 s: a 5.38% reduction
(1.057x) on this bounded 512-cell, one-rank smoke. The three M5 commits summed
to 3,164.836 s on CPU and 3,001.129 s on A100 (5.17% lower). The Slurm elapsed
records were 53:37 and 50:48, respectively (5.25% lower). This is a modest
single-pair result, not a substantial or replicated performance gain. It is
not a scheduler-matched benchmark: the CPU reference used `cohort`, while
hybrid used `level_queue`. The effective namelist SHA256 is identical
(`6c8e0a0b…5772fa`), both runs use the same reference CHIMES library SHA256, and
the `patch/lagRamses` source is unchanged between their recorded source
commits. Backend and cell-scheduler differences remain part of the measured
configuration.

Slurm also lists a two-second auxiliary bash step `.9` with exit 1. The main
application step `.0` and batch step both completed with exit 0, and the
runner wrote its pass sentinel. The failed auxiliary command is not preserved
in the accounting record and is not attributed to the simulation.

Both results are **integration passes only**, not production qualification.
The profile deliberately disables stellar SED emission, uses the AGN/dust
reference controls, and declares `source_scope=AGN_only_point_source` and
`reference_controls=NONPRODUCTION`. Raw outputs were suppressed in both runs;
the run directories retain only compact inputs, logs, sidecars and result
sentinels. This supersedes the in-progress status recorded immediately above.

## Current-runtime-source effective-population hybrid integration — 2026-10-30

A40 job `442357` completed on `syn05` with exit code 0 (3:04 Slurm elapsed;
1 MPI rank, 8 OpenMP threads, one typed A40, 64 GiB requested; application
step MaxRSS 1,025,736 KiB). It used the `effective_population` profile with
stellar SED enabled, a real accreting BH sink, M5 order 5, live dust, CHIMES
hybrid level-queue execution, and 256 passing momentum assertions. Four
coarse steps completed; the process reported 170.060 s total elapsed. The
output policy suppressed RAMSES dumps (`raw_outputs=0`).

The executable is the CUDA-enabled build from job `442271`, source commit
`a572d1dbfd12ecf4b2ce68773ff53f13b103dfa5`, using the reference CHIMES library
SHA256 `64ad78dbd531621904a9fc00dea2c81d8068198fa245fff659bede9ca4d24ab2`.
The current branch changes after that source commit touch only the integration
run script under `simulation/snrt/runs`; the runtime Fortran/C++ source is
unchanged. The run's `inputs.sha256` also binds the executable, effective
namelist, sink IC, SED source, CHIMES tables, group tables and selected
contracts.

The level queue completed all 64 cells at each M5 commit with eight workers
and one GPU broker. CHIMES counters report 61,258 GPU calls, zero CPU
fallbacks, 2,150 batches, zero device errors, two warm-table uploads, and a
peak of two in-flight batches. The M5 source ledger is nonzero after startup;
the runner completed its M5/dust/CHIMES and momentum checks. The stellar
population binding is `common_imf_mixed_evolution`, and the configured SNIa
accounting is `effective_ssp`. The run formed 64 new stars at each coarse
step, reaching 256 total. The AGN sidecar records positive BH accretion and
bolometric luminosity with `JET` feedback mode at coarse steps 2–4; each
record's efficiency contract is valid.

The result line says `AGN_state=not_requested` because the runner's explicit
AGN-state pass condition is scoped to the separate `agn_reference` profile.
It does not mean the effective-population run had no AGN activity; the
sidecar evidence above records that activity. The four-step horizon is only
about 0.206 Myr, so this run does not establish a new SNII or SNIa explosion
and deposition event. It verifies the current-runtime-source combined stellar-radiation,
AGN, dust and CHIMES execution path under its declared controls.

This evidence is not production qualification. The admitted run explicitly
uses the reference stellar spectral, nine-group, AGN and dust contracts
(`reference_controls=NONPRODUCTION`). It must not be used to claim those
physical inputs are production-approved. No raw dumps were retained.

## Source-bound Kroupa package follow-up — 2026-10-30

CPU Slurm package-build job `1245335` completed successfully on grammar's
`normal` queue in 3:38 (peak RSS 642,500 KiB). It generated
`.m5-kroupa-source-package-1245335/input-v5` with `imf_id=1` and a separate
effective-SSP overlay. The manifest SHA256 is
`0d8fa29295f4bc9df6295c0965355340de20c511dcb318933a6a716fca8951e7`; the
binding records the SNIa runtime-contract hash and retains
`reference_control_not_production_approval` status.

The package was compared to the inputs used by A40 hybrid integration job
`442357`. `history.nml` is byte-identical (SHA256
`ac72b8022dcd7258926b2cd891d09983dc61c91c978747ca4122c619cf880874`), as are
`yields.dat` (`a442f8b5637d8cb80260b835f02efae93908b6989029b1dd511f132fe4c140c6`)
and `nodes.dat`
(`f319b2e55d64e61354fe9f8943a5ab17f5375e292bbf5a7b825203c15713d7c3`). The
effective-SSP SED has the same embedded node digest and population controls;
its substantive difference is the absolute node-file path into the newly
bound package (plus a trailing blank line). The package validator and
effective-namelist history-path rewrite both passed on the generated package.
Thus the earlier live hybrid result
provides runtime evidence for the same numerical stellar inputs: 61,258 GPU
calls, zero CPU fallbacks, and zero device errors.

A follow-up H100 job `442427` was submitted on syntax but remained pending
with Slurm reason `QOSGrpGRES`. It was cancelled while still pending after the
exact input comparison and package/path checks showed that another coupled
run would add no new numerical evidence. It consumed no GPU time and is not a
runtime failure.

This closes the Kroupa package/IMF metadata and runtime-path consistency
repair. It does not promote any yield, SED, SNIa, AGN, dust, or nine-group
reference contract to production status; the preceding nonproduction limits
remain in force.

## Three-step AGN reference revalidation — 2026-10-30

CPU Slurm job `1245972` completed on grammar `normal` (grammar078) with exit
code 0 after 54:09 elapsed. It used one MPI rank, eight OpenMP threads, and
1,392,612 KiB peak RSS. The build was job `1245437`, `CHIMES_CUDA=0`, source
commit `a0a557387f306249e7c9591a6bd311af8b771eb3`; executable SHA256 is
`47c2a56a12c2b6297b063541c9f5fbc59bd7a934ca236afb955e0d1f6a54834e`. The
effective namelist SHA256 is
`b63d24b1b9274f03151377b0e3145f2fbc7310e5ac0f05a0b3ba8554dfc1a349`, and the
sink IC SHA256 is
`971d2f0f2b0c4731fd2b9801ebc5c363dd9cd51268660eb3cf801e7a34d88a54`.

The result sentinel is
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS`: three completed coarse steps,
M5 order 5, live two-size dust, CPU CHIMES, a real accreting BH sink, and
stellar SED disabled to isolate the AGN source. The step-2 sidecar records
positive BH growth, inflow `5.39665e-6 Msun/yr`, Eddington ratio `83.9178`,
bolometric luminosity `3.05619e40 erg/s`, and active thermal feedback. Step 3
records jet mode. The 16,896 momentum-assertion cells passed. The active M5
photon-ledger residuals are `2.68312e-17` at step 2 and `3.79273e-18` at step
3; corresponding radiation-energy residuals are `9.62051e-17` and
`-5.52590e-17` in the printed ledger convention. No nonzero NaN checks, M5
rejections, dust-transaction failures, or fine-MG nonconvergence warnings
were found. `abort_on_mg_nonconvergence=.true.` was active. The output policy
produced no raw snapshots (`raw_outputs=0`).

This is an AGN-source and coupling smoke under the declared nonproduction
reference controls, not production physical calibration. RAMSES attributed
99.8% of its 3,241.2 s runtime to the combined `cooling` bucket; this does not
separate CHIMES from the coupled M5/material work. Printed `econs` is not a
closed gas+radiation+dust+AGN energy budget. The retained compact run record is
`.m5-source-dust-chimes-1245972/driver-evaluation.txt`.

## Effective-population CPU/GPU comparator retention update — 2026-10-30

The CPU effective-population run `1245438` completed all 128 coarse steps and
passed its physical postconditions, although its old runner selected the
hybrid-only `CHIMES_LEVEL_QUEUE` log check for a CPU build and therefore
returned batch exit code 1. Its final HDF5 dump remains the CPU comparison
reference. H100 job `442500` was cancelled before execution (zero runtime), so
it produced no comparison. The current-source A40 hybrid counterpart is job
`442563` on syntax, pending GPU resources with tentative start
2026-10-31 12:06 KST. Preserve the CPU dump until `442563` and the state
comparison are evaluated; remove only the raw dump after that evaluation.
