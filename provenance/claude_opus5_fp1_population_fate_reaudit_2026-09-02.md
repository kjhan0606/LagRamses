# Claude Opus 5 F-P1 consolidated re-audit — 2026-09-02

Auditor: Claude Opus 5 CLI, high effort, read-only. This was one consolidated
re-audit of the initial D1--D5 repair bundle. No files were edited and no
RAMSES time integration was launched by the auditor.

## Reproduced evidence

The auditor reproduced the population contract (34/34 checks), production
yield-audit mutation test, feedback-policy test, interpolation fail-closed
test, JAX 0.11.1 CPU IMF differential, G2 population-ledger runner, G1
native/JAX diagnostic, P0.4 production negative startup, source-identity
classification, all 136 recorded build-input hashes, and binary SHA-256
`5c999dc5b34812a92d2e3a0bf4f538bad9e15ed2bb803b06e549ae0f2d91b815`.

## Finding disposition

- **D1 closed.** The auditor independently re-derived the four IMF mass
  antiderivatives and confirmed the Chabrier continuity amplitude. Scratch
  fault injection into the Fortran antiderivative and restoration of the old
  high-branch bug both caused the unchanged JAX differential to fail.
- **D2 closed on the live runtime, with one condition.** Production startup
  and inline RAMSES deposition correctly accept
  `sum(tracked)<=returned`; generic metallicity receives
  `returned-H-He`, independent of individual-element switches. However,
  three linked but currently uncalled adapters still required equality and
  the Phase-0 design document repeated that obsolete rule.
- **D3 structurally closed, with one condition.** The production timestep
  path calls the ledger unconditionally using already computed cumulative
  states, clears rejected sources, and checks particle mass before writes.
  The runtime cross-check nevertheless used `max(1,...)` in code-mass units,
  weakening a nominal relative bound to an absolute fraction of box mass.
  Similar absolute terms appeared in returned/SNII checks. A minor related
  observation was the ledger's absolute negative-mass tolerance.
- **D4 closed.** Production and native strict audit modules were
  byte-identical and covered duplicate coordinates, Cartesian completeness,
  age-zero anchors/nullity, cumulative material/residual and energy
  monotonicity, mass bounds, and remnant ownership.
- **D5 closed.** Legacy mode commits all runtime switches after a successful
  read, while every failed channel-resolved parse returns before commit.

Regression review found no API, zero-width increment, fail-closed clearing,
unit-conversion, field-map, or source/particle semantic regression beyond the
code-mass tolerance weakness above.

## Auditor verdicts

**D1--D5 ENGINEERING VERDICT: CONDITIONAL PASS.** Conditions:

1. Replace the dimensionally inappropriate unit code-mass tolerance floor by
   an initial-particle-mass-relative scale.
2. Correct or explicitly de-scope `stellar_cell_deposition`,
   `stellar_ramses_bridge`, `stellar_ramses_mapped_bridge`, and the obsolete
   equality statement in `PHASE0_STELLAR_ENRICHMENT.md`.

**OVERALL F-P1 SCIENTIFIC VERDICT: BLOCK.** The source package, IMF
convention, lifetimes, fallback/remnant/fate map, license, checksum, and
approval sidecar are not approved. The auditor independently obtained 6.75%
of initial Kroupa SSP mass in the 40--120 Msun winds-only/no-terminal-owner
gap. The unapproved source and this gap are reported but not runtime-enforced;
SNIa and PISN are genuinely fail-closed. PISN remnant semantics remain an
F-P3 decision.

## Independent post-audit disposition

Both engineering conditions were independently reproduced and repaired after
this audit:

- production and native runtime bounds now scale to `abs(mp0(ipart))`, with no
  box-mass-sized `1.0` floor; ledger negative-mass tolerances scale to initial
  population mass;
- all three linked adapters now accept non-negative untracked residual mass;
  the mapped and optional generic-metal outputs deposit tracked C--Fe plus the
  residual, and the design document states the same contract.

`tests/run_stellar_residual_deposition_unit.sh` passes, as do the full
production forced build/link/smoke, P0.4 negative startup, policy, yield-audit,
population-ledger, IMF differential, and source-parity checks. The rebuilt
binary SHA-256 is
`4e7bc3a4e8b2be6681cdd24e1acc2d40dac53e68110304929f164cb4bcfa555e`.
Per the reduced audit-frequency decision, these post-audit repairs will be
included with subsequent F-P1 work in the next bundled audit rather than
triggering an immediate micro-audit.
