# P0.4 fail-closed runtime gate — 2026-09-02

## Decision before independent audit

**PASS.**  The initial read-only Claude Opus 5 audit returned CONDITIONAL PASS;
its in-scope findings were remediated and the focused re-audit returned PASS
with no new blocker.

P0.4 is an admission-control gate.  It ensures that incomplete source physics
cannot silently enter the production executable; it does not approve a
physical yield asset and does not close the mandatory population, SNIa DTD, or
PISN/PPISN gates.

## Implemented contract

- Channel-resolved production startup requires a nonempty, existent external
  `PHASE0_YIELD_TABLE`; it cannot select the compiled synthetic fixture.
- Mass, birth-metallicity, and age queries outside the table domain, and all
  nonfinite queries/results, return hard errors rather than endpoint clamps.
- `imf_id`, `population_model`, and the five channel mass windows are parsed
  transactionally and validated; effective values and the external table path
  are emitted to the RAMSES info file.
- The production main program runs this preflight immediately after parameter
  parsing.  Every enabled channel must exist and span its configured mass
  window before IC allocation or time integration.
- The namelist group is mandatory in a Phase-0-enabled build.  The actual
  loaded path, row count, channel enables, active elements, population, IMF,
  and mass windows are recorded rather than re-reading the environment later.
- The currently implemented source model is limited to `single_star_ssp` with
  SNIa and PISN disabled.  Unsupported binary/SNIa/PISN configurations fail at
  production initialization until their mandatory active gates pass.

## Reproducible evidence

- `tests/run_stellar_feedback_policy_unit.sh`:
  `stellar feedback policy: PASS`.
- `tests/run_stellar_yield_fail_closed_unit.sh`:
  `stellar yield fail-closed policy: PASS`.
- `simulation/snrt/tests/run_g2_preflight.sh`: PASS.
- `P0_DIAGNOSTIC=1 simulation/snrt/tests/run_g1_native_contract.sh`: PASS with
  JAX 0.11.1 CPU differential tests.
- `simulation/snrt/tests/run_g1_native_contract.sh`:
  `STELLAR_SOURCE_PARITY_PASS blocked=none` and
  `G1_NATIVE_CONTRACT_RUN_OK`.
- `P0_BUILD=1 simulation/snrt/tests/run_p0_production_linked_contract.sh`:
  `P0_PRODUCTION_BUILD_EVIDENCE_RECORDED`,
  `STELLAR_SOURCE_PARITY_PASS blocked=none`, and
  `P0_PRODUCTION_LINKED_CONTRACT_OK`.
- The binary audited for P0.4 had SHA-256
  `46e408e60ab289eccffde6b00c3e7f93c4669e4ada25f455448bd2354deedb05`.
  The descendant F-P1 production build has SHA-256
  `188ddd9fc698730b58fe4b3736c395c124c42a234e466bd5dc436786a59c6665`
  and matches the current machine-readable build evidence below.
- Machine-readable build evidence:
  `simulation/snrt/data/p0_production_linked_build_evidence.json`.
- `simulation/snrt/tests/run_p04_production_negative.sh` executes the real
  production binary without entering time integration: an unset table exits 1
  and the intentionally incomplete fixture exits 121 at channel-coverage
  preflight.

The G1/G2 numerical suite runs the declared native mirror and is a differential
oracle, not a substitute for the production binary.  The production-linked
harness proves forced compilation/linkage; the separate P0.4 negative runner
executes the admission branches in that binary.  No positive live-feedback run
is claimed because a physically approved full-grid asset is a later P0.5 gate.

## Initial Opus 5 finding disposition

- C1 OpenMP/MPI error exit: fixed locally in the stellar-feedback call path.
  Worker errors are collected through the end of the OpenMP region and the
  parent thread then uses `MPI_ABORT`; no global MPI redesign was made.
- C2 deferred missing-channel failure: fixed by startup enabled-channel/window
  coverage preflight.
- C3 silently absent namelist: fixed by requiring the group and echoing IMF and
  population identifiers.  The legacy global `clean_stop` zero exit status is
  not a silent acceptance and is tracked as supporting infrastructure.
- C4 environment re-read/incomplete output identity: fixed by recording the
  loaded runtime identity, row count, channel enables, and active elements.
- C5 evidence overstatement: corrected above and covered by production-binary
  negative execution.
- C6 thin unit matrix: added NaN mass/metallicity/age, negative arguments,
  nonfinite result, transactional-state, missing-group, and PISN cases.
- C7 nonfinite input columns: fixed in the production loader before row commit.
- C8 production/native driver finiteness parity: fixed in both trees.

The independent audits judged only this P0.4 admission-control contract.
Physical DTD/population/yield approval is intentionally not claimed here and
is tracked as active mandatory work in
`feedback_population_dtd_active_roadmap.md`.
