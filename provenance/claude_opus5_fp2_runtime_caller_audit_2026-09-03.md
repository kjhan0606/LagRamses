# Claude Opus 5 F-P2 runtime-caller audit — 2026-09-03

Auditor: Claude Opus 5 CLI, high effort, read-only. Full external artifact:
`/home/kjhan/.claude/plans/claude-opus-5-virtual-yeti.md`.
The audit treated the dirty checkout as authoritative and did not modify files,
build, or launch a simulation.

## Verdict

Top-level: **CONDITIONAL PASS** for a coherent production-linked,
runtime-gated implementation. SNIa activation: **BLOCK**. Publication:
**BLOCK**. Opus found no path by which SNIa can activate in the current
production binary.

## Commonly verified behavior

Opus confirmed the fail-closed production gate, ordered contract structure,
DTD interval/restart mathematics, AGB-owned WD debit and ledger closure,
generic-driver SNIa exclusion, normal restart reconstruction, real AMR leaf
lookup, local row-major `unew` indexing, source parity, and the fresh linked
production build. It explicitly distinguished normal restart evidence from
unproven hard-crash exactly-once behavior.

## Additional Opus findings

These are additional to the AGY findings and must be treated as unresolved
until independently dispositioned:

1. **B1 — blocker, binding provenance.** The namelist binds to
   `c6c8042b03406b9d69bc50434fe5d6af7f542be6`, but that commit does not contain
   the current population, cell-deposition, or new runtime-accounting modules.
   The runtime currently checks only that three strings agree and that the
   string is hexadecimal; it does not prove ancestry against the build
   revision. The current working tree is also dirty.
2. **E1 — major, model contradiction.** The SNIa contract declares
   `population_binary_ssp`, while `production_source_model_supported()` still
   requires `population_single_star_ssp`. This is safely fail-closed now but
   must be reconciled before activation.
3. **E2/E4/E5 — major, evidence/build gaps.** The runtime loader itself is not
   directly executed by the contract test; linkage patterns do not explicitly
   assert SNIa symbols; and the new accounting object lacks its own explicit
   Makefile prerequisite on `stellar_enrichment_config.o`, relevant to parallel
   builds even though the forced serial build passed.
4. **E6/E7 — major, caller correctness.** The SNIa scatter occurs before later
   possible failure/abort points, so a future nonzero-momentum or late-failure
   path could duplicate deposition on retry. Element variables are built
   without the generic path's `active_element` mask.
5. **E8/E9/E12/E15/E16 — further correctness/qualification.** The combined
   closure check is partly tautological after reconstruction; re-callable
   initialization has latent shared-state race exposure; the caller does not
   reject `enable_snia` without AGB; the metallicity source id is decorative
   while the caller uses unity; and no IMF consistency check ties DTD
   normalization/binary fraction to the run.
6. **B2 — blocker, physical provenance.** The normalized HESMA record marks the
   event contract and canonical conversion as blocked, with five approvals
   still open. Adopted returned mass and energy are review-derived values; the
   profile mass estimate differs from the integrated stable-element mass by
   3.7683%, and isotope-to-element conversion, decay horizon, terminal-remnant
   ownership, signed momentum, and DTD population weighting remain unapproved.

## Activation gates requested by Opus

Commit and rebind all consuming modules, prove binding ancestry/build identity,
resolve the binary-vs-single SSP and IMF consistency, execute loader and
production-binary negative tests, add explicit SNIa linkage symbols and Makefile
dependencies, make scatter/update all-or-nothing, apply the active-element
mask, then close the hard-crash journal and multi-cell/MPI qualification. For
publication, close the five HESMA approvals and replace the review-only mass,
energy, isotope/decay, momentum, and population assumptions with authoritative
contracts.
