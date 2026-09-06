# F-P2 blocker reconciliation — 2026-09-03

Status: **physical baseline approved; runtime activation remains gated**.

## Completed in this bundle

- HESMA `n300c` is classified as a source-data anomaly and quarantined.
- The adapter and admission sidecar reject all selected models carrying
  unresolved physical warnings.
- HESMA audit, comparison, selection, adapter, and Keegans audit payloads use
  repository-relative provenance paths.
- The sidecar has one canonical 19-field promotion list, mirrored by
  `required_for_promotion`; validator tests reject a stale mirror or malformed
  list.
- The runner resolves the project root from its own script path and executes
  both admission and contract negative paths.
- The DTD primitive uses a cancellation-safe representation near `alpha = -1`
  and passes an independent logarithmic-series reference test at
  `alpha = -0.99999999`.
- A review-only physical event-contract interface now validates explicit WD
  reservoir debit, returned/remnant mass closure, event energy, and either
  signed source-frame or isotropic-zero-vector momentum.  Its native and
  production mirrors are hash-checked by the F-P2 contract audit.
- The population ledger now accepts an explicit WD-reservoir assignment and
  applies a validated SNIa event budget transactionally, preserving aggregate
  initial = living + remnant + returned closure.  Radial momentum requires a
  unit deposition direction and is converted before ledger consumption.
- The cell-local SNIa deposition adapter now converts the validated physical
  budget into mass, bulk-plus-event momentum, event energy, bulk kinetic
  energy, and total-energy density increments using the same cgs convention as
  the existing generic runtime path.  It admits only an explicit
  all-to-total-energy thermal policy; unsupported fractional policies fail
  closed.  Native/production hashes and `-check all` unit tests pass.
- The production RAMSES bridge now consumes those increments with explicit
  length/density/velocity code-unit scales and normalized multi-cell weights.
  Its validation is pre-write and the bridge unit test covers conservation,
  scale conversion, and unchanged `uold` on an unapproved policy.
- The bridge now has a RAMSES-facing `unew(local_cell,variable)` entry point
  that validates target-cell uniqueness, local MPI ownership, and target
  bounds before depositing to variable-major scratch storage and scattering
  to the selected row-major cells.  The production bridge, runtime, and
  feedback caller objects compile under the actual `/gpfs` Makefile.
- The population realization contract now makes binary source identity, IMF
  conversion, DTD support, metallicity factor, expectation/Poisson choice,
  immutable source revision, and approval id explicit.  Its approved baseline
  is checked by the F-P2 audit.
- The full F-P2 contract runner passes, including native Fortran tests,
  source-format/profile/selection tests, admission negatives, and production
  mirror compilation.

## Still blocking production

- HESMA `yysd4-xap92/n100` is selected and approved as the physical baseline;
  Keegans remains a review-only comparison because its project-element rows are
  incomplete.
- The approved Maoz field DTD supplies the binary-population normalization and
  the approved HESMA event source supplies the WD debit and zero terminal
  remnant policy.  The runtime now loads the three-group handoff and applies
  the ledger transaction behind the production gate.
- The event-frame momentum convention is explicit and radial conversion is
  tested.  The cell-local conversion adapter and production-side RAMSES array
  bridge are called by the runtime's located AMR leaf-cell path; this is a
  one-cell NGP handoff, not a distributed neighbour stencil or MPI exchange.
- SNIa thermal coupling is approved as all-to-total-energy for this baseline but
  is not called by the runtime; existing generic stellar energy deposition and
  SNII delayed-cooling are separate paths.
- DTD normalization, IMF conversion, event realization, source energy, WD
  debit, momentum convention, and approval id are populated in the approved
  baseline.  Metallicity sensitivity and full net-yield semantics remain
  explicit follow-up work.
- Approval id `FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ` is bound to source
  staging commit `c6c8042b03406b9d69bc50434fe5d6af7f542be6` and the promotion
  tool checksum in the sidecar.

## Independent audit disposition

AGY (`gemini-3.8-flash-high`) assessed the review-only scaffolding as
coherent and production activation as blocked.  Fable assessed it as a
conditional review-only pass and also blocked promotion, with stronger focus
on WD debit, momentum semantics, portability, and evidence scope.  The
audits were performed before the final mirror/warning hardening in this
bundle; the post-hardening full runner is the current local evidence.

The profile discrepancy values are recorded from the local audit rather than
the AGY prose: `n300c` is 6.4104173893 relative discrepancy (641.04%) and
`n1600c` is 0.05112189345 (5.112%).

## Qualification bundle result

The connected caller qualification is complete: normal retry/restart
reconstruction, weighted bridge conservation, and the linked production
runtime-disabled negative paths all pass.  The actual production binary was
rebuilt from the current `/gpfs` tree and source parity remains closed.

The F-P1 terminal-fate gate and runtime activation remain closed.  A
hard-crash pending-event journal, distributed neighbour/MPI deposition, and
full net-yield/metallicity population extensions are separate follow-up
gates; none is silently claimed by this result.
