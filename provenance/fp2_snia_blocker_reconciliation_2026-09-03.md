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
- The full F-P2 contract runner passes, including native Fortran tests,
  source-format/profile/selection tests, admission negatives, and production
  mirror compilation.

## Still blocking production

- No HESMA or Keegans event model is selected for production.
- The new contract interface is not yet wired into the population ledger;
  WD-reservoir debit and terminal-remnant ownership therefore remain
  unapproved for runtime.
- The signed event-frame momentum convention is represented as an explicit
  interface, but the scalar/radial versus vector cell-deposition policy is
  not approved or wired into AMR deposition.
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

Wire the review-only interface into the population ledger and a separately
tested deposition adapter, then bind the populated source/approval record to
an immutable commit.  Add conservation and restart tests before any runtime
activation change.
