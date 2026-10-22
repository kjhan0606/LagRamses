# Current bridge-state requalification after the hybrid speed handover

Source baseline: `1b527042810d3fae09b69f26611a88449883b435`, branch
`agent/chimes-cvode-gpu-20261001`, origin `kjhan0606/LagRamses`.

The Oct 6 hybrid speed handover remains valid for its separate fixture:
512-cell speedup about 1.02x, no production promotion. Its inputs bypass the
bridge's element/charge admission and do not resolve the earlier logical-cell-5
failure. The following is numerical requalification, not a speedup claim.

## Changes now applied to the /gpfs working tree

The formerly staged fixed-worker diagnostic is applied. It runs six copies
of one bridge-admitted state on workers 0 and 3, CPU/GPU RHS, two repetitions,
at tolerance factors 1, 0.1, 0.01. Test-only CVode interposition forwards the
original call without changing the solve and records terminal statistics.
It does not trace the first divergent step. Failed physical comparisons now
state the compared count and print NaN when no pair succeeded.

The runner has an explicit `SNRT_CHIMES_DIAGNOSTIC_LIBRARY=reference|trace`
choice, default `trace`. The trace option applies the existing tracked
`simulation/snrt/data/chimes_trace_solver.patch` to the job-local library
build, without changing the external CHIMES checkout or production library.

The existing bundle fixture now also supports `check-level` / `bench-level`:
same bridge-admitted cells, host CVODE and dense LU, CPU-only versus current
8-worker/8-stream GPU RHS/cooling level queue. Every cell must succeed before
timings are usable. Existing species, temperature, energy-comparison, element
and charge checks remain enforced. This mode passed job 414412 below. It is
not a complete RAMSES or full GPU integrator test.

## Job 414219: current GPU code, uncorrected reference library

H200 NVL, syn104, COMPLETED 0:0, 27 s. Log:
`.cvode-gpu/tile-rhs-414219.log`.

| Tolerance factor | CPU accepted | GPU accepted | Worker/repeat verdict dependence |
|---|---:|---:|---|
| 1 | 0/6 | 0/6 | none observed |
| 0.1 | 6/6 | 6/6 | none observed |
| 0.01 | 0/6 | 0/6 | none observed |

Each row was repeated twice on each of workers 0 and 3. Failure was actual
`CV_ERR_FAILURE=-3`, with zero reported nonlinear convergence failures.
At factor 1 the CPU failed at t=970708788.55756474 s, last accepted order 5;
the GPU failed at t=971959896.69638705 s. At factor 0.01 failures occurred
earlier (CPU 58874758.892786466 s, GPU 55898149.813500293 s).
At factor 0.1 all admitted states passed the fixture budgets; CPU/GPU maximum
species-scaled discrepancy was 0.1226535447 (existing bound 1), temperature
relative discrepancy 2.44e-15. This single successful tolerance is not a
converged reference and is not selected as a production fix.

Current GPU cooling changed the old CPU-fail/GPU-pass pattern to both failing
at the nominal tolerances. No thread-routing explanation was observed in this
bounded test. Success/failure alone does not establish accuracy.

Build hashes from this run:

- test: `889e14a1c9401c81d1d2fe4ca8acb60847143ac19aa7322d97ec0110a0caf3c8`
- libchimes: `56ab7b5e051c1a05b2d2899fb4b602cfb80ab6fe382917d650cb3a2b8e5df61c`

## Existing correction found in the repository

`provenance/chimes_gpu_followup_2026-10-01.md`, sections “Trace-species failure
recovery bundle” and “GPU requalification and cost decomposition bundle”,
already document a selected correction and prior CPU/H100 tests. The new
diagnostic runner initially omitted this correction from its library build.

The selected trace patch tightens only molecular absolute tolerances by 1e-3
for dark molecular thermal solves and allows at most four CVODE history
restarts on CV_ERR_FAILURE from the returned last accepted y,t. Relative and
thermal/atomic tolerances, target time and conservation checks are unchanged
by that patch. No-progress and other failures remain failures. Earlier
qualification does not automatically qualify the subsequently changed GPU
cooling and scheduler; that is the purpose of the current follow-up.

## Job 414224: corrected trace library

Submitted with one H200 (`--gres=gpu:H200:1`), 8 CPUs, 16 GiB, 15 min,
excluding syn102, fixed-route mode and six cells. Expected log:
`.cvode-gpu/tile-rhs-414224.log`.

Completed 0:0 on syn104 (H200 NVL) in 28 s, after initially waiting for
`QOSMaxGRESPerUser`. The normal QoS limit is eight GPUs per user across
partitions; no unrelated jobs were cancelled or duplicate submission made.

All 144 cell executions succeeded: six identical states, two workers, two
repetitions, two backends, and three tolerance factors. Within each backend
and tolerance the reported endpoint differences across workers/repeats were
zero. Every accepted state passed finite/nonnegative, element and charge
checks. The previous status-41 reproducer is resolved for this tested matrix
when using the existing corrected library.

| Tolerance factor | CPU/GPU species-scaled difference (bound 1) | T relative difference |
|---|---:|---:|
| 1 | 1.7643228895e-4 | 2.44e-15 |
| 0.1 | 1.6266592133e-6 | 3.33e-15 |
| 0.01 | 1.2964969299e-8 | 8.88e-16 |

Successive tolerance-tightening endpoint differences decreased from about
0.02440 (1 to 0.1) to 0.0001975 (0.1 to 0.01) in the same species-scaled
metric, for both backends. This is scoped convergence evidence, not universal
certification. The initial call completed for each of these corrected solves;
this run does not separately qualify the history-restart branch.

- test SHA256: `3e63075f40fece62ca54133b64c7b502e2b5f9bab4c1170613207f6f6833db00`
- libchimes SHA256: `2067fcd0c34b6a6443be220b620b00263d34736094be7b2203d7b13c449a6d99`

## Job 414412: bridge-admitted 512-cell level queue

Submitted after inspecting 414224, with the same H200/8 CPU/16 GiB request.
Runner mode `level-qualify`, library `trace`, 512 cells. Default and mixed
composition each require every CPU/GPU cell to succeed and pass the unchanged
endpoint criteria. Warm timings use three alternating CPU/GPU samples. Log:
`.cvode-gpu/tile-rhs-414412.log`. Completed 0:0 on syn104, H200 NVL, in 43 s.

Every 512-cell solve succeeded in both compositions and both arms, including
the cold comparisons and three warm alternating pairs. All 8,192 cell solves
completed; there were no mismatched rejections or excluded failed endpoints.
Species, temperature, energy-comparison, element and charge checks passed.

| Composition | CPU median, seconds | Hybrid median, seconds | CPU/hybrid speedup | Max species-scaled error (bound 1) |
|---|---:|---:|---:|---:|
| Default | 0.459555727 | 0.435920322 | 1.0542 | 0.0005945003 |
| Mixed | 1.38510728 | 1.75874853 | 0.7876 | 0.0032048625 |

Warm sample ranges: default CPU 0.4575–0.4596 s, hybrid 0.4277–0.4530 s;
mixed CPU 1.3827–1.3863 s, hybrid 1.7461–1.9553 s. These are three-sample,
shared-node bounded observations, not isolated scaling measurements. Mixed
composition makes the hybrid about 27.0% slower. Default's 5.4% improvement
does not justify changing the production default.

Warm GPU work was real: eight streams were used, no warm table uploads and
no device errors. Default had 639,342 RHS calls / 35,518 batches; mixed had
826,104 / 97,762, over the three warm GPU samples. GPU LU factors were zero.
The test did not use the `peak_concurrent_batches` counter as a performance
or correctness gate and did not modify that counter.

- test SHA256: `0b5f007048373bfff18bfea7f7af2faa6a66e37b063c3932faf818537e41d675`
- libchimes SHA256: `2067fcd0c34b6a6443be220b620b00263d34736094be7b2203d7b13c449a6d99`

The fixture is now bridge-admitted and uses the corrected solver; its values
must not be conflated with the separate speed-handover fixture (414170).
Both fixtures support retaining CPU as default, not a universal speedup.

## Next actions

1. The fixed-route and level-queue requalification bundle is complete. Do not
   repeat these tests without a new code change or a concrete uncovered case.
2. Use the corrected `trace` library explicitly for further chemistry
   comparisons; plain `reference` intentionally reproduces the old failure.
3. Keep CPU as the production default. This bundle qualifies the tested
   experimental GPU RHS/cooling path but does not install a production library
   or qualify a full RAMSES evolution. Full device-resident BDF/Newton control
   has not been implemented, and remains separate algorithm work rather than
   a prerequisite invented by this bounded test.
4. For further acceleration, address measured host integration/linear algebra
   and device residency or useful per-cell GPU parallelism. Repeating the
   rejected single-stream schedule or loosening tolerances is not justified.

No production executable replaced, no cosmological run or raw snapshots,
no commit/push in this follow-up. The old filesystem/network restriction was
lifted; it is no longer the blocker. Jobs 414219, 414224 and 414412 are all
terminal COMPLETED; this bundle has no remaining live submitted job. Retain
inputs, build identities and logs. No simulation raw snapshots were created.
