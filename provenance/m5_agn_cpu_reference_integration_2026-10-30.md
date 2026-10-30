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
