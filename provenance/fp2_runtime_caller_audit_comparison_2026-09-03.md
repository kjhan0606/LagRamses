# F-P2 runtime-caller audit comparison — 2026-09-03

The automatic evidence bundle was strengthened first, then the same checkout
was audited independently by AGY and Claude Opus 5. Audit records:

- [`agy_fp2_runtime_caller_audit_2026-09-03.md`](agy_fp2_runtime_caller_audit_2026-09-03.md)
- [`claude_opus5_fp2_runtime_caller_audit_2026-09-03.md`](claude_opus5_fp2_runtime_caller_audit_2026-09-03.md)

## Consensus

Both auditors agree that:

- the production binary is freshly linked and source parity passes;
- DTD interval mathematics, WD-reservoir ownership/debit, generic-driver
  exclusion, normal restart reconstruction, and the fail-closed gate are
  coherent and automatically evidenced;
- the connected runtime caller is one-cell NGP, not a qualified multi-cell/MPI
  deposition algorithm;
- normal restart does not prove hard-crash exactly-once;
- SNIa activation must remain blocked;
- the current baseline is not publication-ready because net yields,
  metallicity dependence, and broader progenitor/population approval are not
  complete.

## Difference in audit emphasis

AGY primarily validated the intended algorithm and current production-linked
fail-closed wiring. Its additional concrete issue was the unit/production
environment variable name mismatch and the ambiguity of `production_ready`.

Opus performed the more adversarial promotion audit. It additionally found the
stale/uncommitted source binding, the `population_binary_ssp` versus
`population_single_star_ssp` contradiction, lack of direct runtime-loader and
production-binary negative execution, missing explicit accounting-object
parallel-build prerequisite, possible late-failure duplicate scatter, missing
active-element masking, and the five still-open authoritative HESMA physics
approvals. Opus also quantified the profile-versus-integrated-mass discrepancy
as 3.7683%.

## Disposition

The safe combined disposition is **CONDITIONAL PASS for review-only,
production-linked scaffolding; BLOCK for production activation and publication**.
The Opus additions are retained as real blockers/qualification gates, not
discarded because AGY's narrower implementation verdict was favorable.

The next implementation batch should prioritize binding/provenance and direct
runtime evidence, then transactional caller correctness and MPI/multi-cell
qualification, while the HESMA physical-contract approvals and publication
yield/population work remain explicit scientific gates. No activation, commit,
push, or simulation launch was performed as part of this audit turn.
