# Stage-3 closure audit

Read `/gpfs/kjhan/LRD_JWST` without editing it. Use `claude-opus-5` to decide
whether the only two conditions in
`provenance/claude_opus5_rsla_refinement_final_audit_2026-09-02.md` are now
closed by the current `simulation/snrt/RSLA_REFINEMENT_VALIDATION.md`.

Specifically verify:

1. the predeclared acceptance bullet now names the larger padded upper bound
   from both `1/(ĉ/c)` and photon-storage extrapolations, not inverse-ĉ alone;
2. the v3 coordinate set is explicitly closed to exactly those two physically
   motivated coordinates, and changing it requires a new schema/version and
   fresh predeclared run rather than post-hoc coordinate shopping;
3. those statements agree with the already-audited validator, v3 JSON, and
   narrow `0.019764600042460632 < 0.02` gate.

Check the actual diff and relevant files; do not trust this prompt. No solver
or artifact was changed after the prior full audit. Return exactly one verdict:
`PASS`, `CONDITIONAL PASS`, or `BLOCK`, with a concise reason. Do not reopen
later dynamic-source, hydro, helium, or dust work that is explicitly outside
stage 3.
