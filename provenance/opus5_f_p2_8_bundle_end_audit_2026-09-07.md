# Opus5 F-P2.8 bundle-end audit

Date: 2026-09-07 (Asia/Seoul)

## Initial verdict

**CONDITIONAL PASS.** Opus5 confirmed that `--mode gui` reaches the shared
generator and that the native harness reaches
`amr_step -> snrt_ramses_advance_level`. It also confirmed the baseline
transaction commit and forced receiver rollback control flow.

The initial audit identified three conditions before accepting the evidence:

1. The baseline runner accepted a zero-status `clean_stop` after the AGN
   reference path reported an optional-array size error. The cause was an
   implicit interface at the external `average_AGN`/`AGN_blast` calls.
2. `SNRT_DRIVER_TEST_SEED_SOURCE` was not gated by the existing diagnostic
   mode and was compiled into the normal binary without a complete array
   bound check.
3. The GUI real-widget test was skipped on the headless node and would have
   called the scalar `next()` method at the new level-pair screen. The
   provenance also cited an older run and binary hash.

Lower-severity findings were: the level replay did not verify that the
dependent prompt was actually `levelmax`, CLI help/error wording still
centered on `--gui`, and the rollback test was marker/control-flow evidence
rather than a post-rollback state observation. Opus5 found no
over-instrumentation and accepted the explicit limitations on physical sink
formation, physical AGN SED/MAD, MPI/restart persistence, and live dust.

## Closure performed in this bundle

- Added explicit caller interfaces with explicit-shape mandatory arguments and
  assumed-shape optional reference arguments in
  `patch/lagRamses/sink_particle.kjhan.f90`.
- Latched the test-seed environment check, required
  `SNRT_RT_TX_DIAGNOSTIC_MODE=1`, initialized the minimal sink state needed by
  the downstream AGN maintenance pass, and checked allocation and dimensions
  before writing.
- The runner now rejects Fortran fatal, SIGSEGV, invalid-AGN, SNRT preflight,
  and sink-merge failure markers. Its baseline is explicitly diagnostic-mode
  test execution with failure injection unset.
- The GUI verifies the replayed dependent prompt, exposes an explicit active
  level-pair state, tests the paired widget path when a display is available,
  accepts the redundant `--gui --mode gui` spelling, and uses updated CLI
  wording.
- Replaced the stale evidence run with the current run below.

## Post-closure execution

Runner:
`simulation/snrt/tests/run_snrt_agn_driver_faithful_smoke.sh`

Run root:
`simulation/snrt/runs/fp2_8_agn_driver_faithful_smoke/job_20260907T080854_2554271`

Results:

```text
F-P2.8_CASE baseline PASS return_code=0
F-P2.8_CASE injected PASS return_code=0 expected_mpi_abort=1
F-P2.8_PASS run_root=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/fp2_8_agn_driver_faithful_smoke/job_20260907T080854_2554271
```

The baseline contains `SNRT_RT_TRANSACTION_COMMIT_PASS`,
`SNRT_RT_CLOSURE_PASS`, `active sources: 1`, and a clean AGN maintenance
sequence through `average_AGN` and `AGN_blast`. The injected case contains
the receiver rollback and diagnostic fail-closed markers and no commit
marker. The runner found no fatal/error marker and no `output_*` directory.

Post-closure binary SHA256:
`4c2b9597f785901f3956a7c1039644d67f965bc10bf8a64e0d1e110559506d0d`

GUI result:
`17 passed, 1 skipped` (the only skip is the real-widget test because this
node has no usable display; its level-pair branch is now explicit in the
test).

## Final scope

The Opus5 verdict remains a **bundle-level CONDITIONAL PASS**: the closure
fixes the evidence and harness conditions, but no claim of full production
science readiness is made. The physical sink-formation blocker, MPI/restart
coverage, physical AGN SED/MAD, and live-dust production validation remain
separate follow-up work.
