# Claude Opus 5 audit prompt: P0.1 source parity implementation

You are the final independent gate auditor for `/gpfs/kjhan/LRD_JWST`.
Perform a read-only audit. Do not edit files, compile RAMSES, launch jobs, or
change running job `330195`.

Read completely:

1. `simulation/snrt/config/stellar_source_identity_v1.json`
2. `simulation/snrt/tools/validate_stellar_source_parity.py`
3. `simulation/snrt/tests/test_stellar_source_parity.py`
4. `simulation/snrt/tests/run_g1_native_contract.sh`
5. `simulation/snrt/tests/run_p0_production_linked_contract.sh`
6. `simulation/snrt/tools/record_p0_production_linked_build.py`
7. `simulation/snrt/P0_SOURCE_PARITY.md`
8. `provenance/p0_source_parity_gate_2026-09-02.md`
9. `provenance/feedback_implementation_plan.md`
10. `provenance/production_publication_readiness_plan.md`

Spot-check `bin/Makefile`,
`patch/lagRamses/stellar_enrichment_sources.mk`, all production stellar objects under
`patch/lagRamses`, the native mirror source list, and the shared source files
whose hashes are compared. Confirm that the validator itself is read-only and
does not treat stale build objects as evidence. The production-linked harness
is an opt-in build path; do not run it, compile RAMSES, launch jobs, or change
running job `330195`.

Audit questions:

- Does the gate correctly parse and resolve the production Makefile source
  path and all required stellar objects?
- Does the source-order sidecar remain bound to the exact production source
  set, including the backend module?
- Does the current BLOCK accurately identify the remaining P0.1 closure, or
  is any criterion too strong, too weak, or incorrectly implemented?
- Does inserting the gate into `run_g1_native_contract.sh` prevent a false
  native-only PASS without altering running calculations?
- Does the production-linked harness actually force a fresh `bin/Makefile`
  build, verify the binary linkage symbols and startup smoke, and does its
  evidence become invalid after any production source, config, harness,
  validator, recorder, binary, log, or compile-contract change?
- Does the harness/source sidecar relationship load the exact effective make
  object order, and does the build-log parser reject stale or overridden
  compile definitions?
- Does the smoke contract use the source-confirmed expected exit code (RAMSES
  `clean_stop` code 0), and is the free-subroutine `thermal_feedback_` symbol
  pattern valid for an ifx binary?
- Is the config's source-of-truth policy precise enough to support a future
  production-linked test and avoid a new mirror divergence?
- Are the plan and implementation record consistent with the actual gate
  output and with the Fable/independent reproduction record?

Return a concise report containing:

1. `PASS`, `CONDITIONAL`, or `BLOCK` for the P0.1 implementation;
2. evidence-based findings;
3. required corrections, if any;
4. whether P0.1 is closed or remains blocked, and the exact next action.

Treat a correctly failing fail-closed gate as a valid implementation result,
but do not call P0.1 complete while production source/runtime parity remains
unproven.

Re-audit the post-review corrections explicitly: verify that both logs carry
and contracts compare the exact binary SHA-256, that `forced_rebuild` is
derived from the logged `-B` marker, that evidence remains valid after a
descendant commit without weakening the relevant source/tool/binary/log hash
pins, and that shared-contract convergence no longer fails the gate while
production resolution and partition bounds remain enforced.
