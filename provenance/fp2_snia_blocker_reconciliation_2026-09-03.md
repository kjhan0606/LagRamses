# F-P2 blocker reconciliation — 2026-09-03

Status: **implemented review-boundary hardening; physical promotion remains blocked**.

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
- The full F-P2 contract runner passes, including native Fortran tests,
  source-format/profile/selection tests, admission negatives, and production
  mirror compilation.

## Still blocking production

- No HESMA or Keegans event model is selected for production.
- The ledger interface is tested but no approved binary-population source yet
  supplies the WD reservoir or terminal-remnant policy for runtime.
- The event-frame momentum convention is explicit and radial conversion is
  tested.  The cell-local conversion adapter is implemented, but its output is
  now wired to a tested production-side RAMSES array bridge, but the bridge is
  not yet called from AMR/MPI neighbour selection.
- SNIa thermal coupling is represented by a guarded all-to-total-energy cell
  increment policy, but it is not approved or called by the runtime; existing
  generic stellar energy deposition and SNII delayed-cooling are separate
  paths.
- DTD normalization, IMF conversion, event realization, thermal coupling, and
  metallicity dependence remain requirements, not populated physics fields;
  only the numerical evaluation kernel is now complete.
- A source commit binding and named approval id do not exist because no
  production commit has been created for this contract.

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

## Next bundle

Select and bind the physical source/population realization, then wrap the
tested cell adapter with the actual AMR/MPI target-cell and RAMSES conserved
array bridge.  Add multi-cell weighted conservation, restart/idempotence, and
runtime-disabled negative tests before any activation change.  A source
commit binding and named approval id remain mandatory.
