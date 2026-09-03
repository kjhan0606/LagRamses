# Fable pre-implementation audit: next F-P1H bundle

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit files, commit, contact
authors, download or replace source data, or launch simulations. This is a
plan-legitimacy and sequencing audit before implementation.

## Final project objective

The final objective is a production-ready and publication-ready lagRamses
high-level hydrodynamics physics stack, focused on radiative transfer, stellar
and AGN feedback, and dust. The current track must move stellar feedback from
scientifically identified source data through population/fate selection,
canonical conversion, momentum/thermal deposition, and an audited runtime
consumer. Previously completed RAMSES topology, headers, HDF5/restart, and
generic distributed-payload work are out of scope unless a concrete new
feedback dependency is demonstrated.

## Current state

The uncommitted source-identity/rights gate bundle is described in
`provenance/fp1h_f_source_identity_rights_bundle_2026-09-03.md`. Its generated
gate pass is not accepted. Read these independent audit records:

- `provenance/agy_fp1h_f_source_identity_rights_audit_2026-09-03.md`
- `provenance/gpt56sol_fp1h_f_source_identity_rights_reaudit_2026-09-03.md`
- `provenance/fp1h_f_source_identity_rights_audit_comparison_2026-09-03.md`

AGY and Codex `gpt-5.6-sol` independently reproduced a candidate-substitution
pass. The Codex re-audit additionally demonstrated empty/self-consistently
rewritten inventories, circular local DOI/rights evidence, missing blocked
publication and path-confinement invariants, and malformed-type/exception
gaps. Production safety remains closed: no physical nodes exist and canonical
conversion, runtime deposition, production, and publication flags are false.
Grok is unavailable and intentionally excluded at the user's direction.

## Proposed next bundle

Read `provenance/fp1h_next_bundle_plan_2026-09-03.md` in full and inspect the
referenced source, admission, sidecar, candidate-grid, and runtime contracts as
needed. The proposal combines:

1. mandatory trust-root/adversarial remediation of the failed rights gate; and
2. a substantive source-byte-derived coordinate/population inventory plus a
   fail-closed `coordinate_hull_and_population` executable validator.

The production metallicity domain and rotation/binary population are currently
unselected, so Part B is expected to remain blocked unless legitimate approved
decisions already exist. The failed-model wind anomaly also remains unresolved
pending an author explanation or corrected release.

## Audit questions

1. Does this bundle directly and efficiently advance the final RT/feedback/
   dust objective, rather than becoming another schema/audit-infrastructure
   loop?
2. Is Part A sufficient and minimal to make source identity and redistribution
   rights non-circular, reproducible, and publication defensible?
3. Is Part B scientifically justified now? Decide whether to:
   - keep Part A and Part B together;
   - split them but execute both in a stated order; or
   - replace Part B with a checksum- and row-specific failed-wind
     author/corrected-release inquiry and source-acquisition packet.
4. Can a blocked coordinate validator still be worthwhile physical progress?
   If yes, state the exact evidence it must produce. If no, name the better
   concrete next deliverable.
5. Are the proposed pass/block conditions strong enough to prevent source hull
   coverage from being mistaken for production-domain/population approval?
6. Identify any missing physical prerequisite that should take priority over
   coordinate admission, especially age-resolved winds, fate/remnant
   structure, energy/momentum/deposition, PPISN/PISN, or runtime coupling.

Return exactly one plan verdict: **APPROVE**, **APPROVE WITH CHANGES**, or
**REJECT**. Then provide:

- required changes, ordered by priority;
- the exact recommended bundle boundary and execution order;
- machine-checkable acceptance criteria for the approved bundle;
- what must explicitly remain blocked;
- a concise purpose-fit judgment tied to the final project objective.

Do not audit unrelated completed infrastructure and do not implement anything.
