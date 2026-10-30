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

The read-only QoS query reports `MaxTRESPerAccount=gres/gpu=8` for `normal`;
the user association has no separate GPU cap reported. This establishes the
account-level ceiling associated with the pending reason, but does not identify
which account jobs currently consume those GPU slots.

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

Integrated trace job `1244750` now uses that exact candidate library and the
same AGN-reference runtime profile on the CPU-only `normal` queue. It was
submitted to start after the reference run terminates (regardless of its exit
status) and the trace build succeeds. At this record update it is running;
no trace-integrated verdict is claimed yet. GPU hybrid jobs `442109`/`442110`
remain a separate `syntax` H100 chain pending `QOSGrpGRES` and its build
dependency.
