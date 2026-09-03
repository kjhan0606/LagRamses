# Claude Opus 5 final G2 source-package staging audit request

Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, commit, push,
launch a RAMSES time integration, or promote a candidate. This is the final
bundled audit after the residual findings in
`provenance/claude_opus5_g2_source_package_staging_reaudit_2026-09-02.md` were
addressed.

Inspect the current versions of:

- `simulation/snrt/tools/audit_g2_candidate_sources.py`
- `simulation/snrt/tools/audit_g2_sukhbold2016_candidate.py`
- `simulation/snrt/tools/build_g2_sukhbold_channel_projection.py`
- `simulation/snrt/tests/g2_candidate_sources.py`
- `simulation/snrt/tests/g2_sukhbold2016_candidate.py`
- `simulation/snrt/tests/g2_sukhbold_channel_projection.py`
- `simulation/snrt/data/g2_candidate_source_audit.json`
- `simulation/snrt/data/g2_sukhbold2016_candidate_audit.json`
- `simulation/snrt/data/g2_sukhbold_channel_projection_review.json`
- `simulation/snrt/config/g2_sukhbold_channel_projection_contract_v1.json`
- `provenance/fp1_source_package_selection_plan_2026-09-02.md`

Check these dispositions independently:

1. Aggregate acquisition-manifest and Limongi/NuGrid parser failures are
   fatal at the top level and return a non-zero CLI code; clean input still
   remains review-only.
2. The projection requires the complete fail-closed firewall map and approval
   prerequisites, the component map is complete, and each component's source
   column, proposed channel, and ownership are coupled to a fixed expected
   identity. Mutation tests must cover missing firewall/component and bad
   ownership/channel.
3. Implosion-wind parsing does not fabricate an ejecta value, derives its
   wind-only and nonnegative claims from parsed data, preserves radioactive
   overlap markers, and rejects invalid input.
4. High-mass W18/N20 yield tables and W18 implosion-wind records are exposed
   with source-branch provenance, per-engine counts, explicit missing masses,
   source fingerprints, and zero canonical/deposition permission. No
   interpolation or source values are invented.
5. Run bounded tests if possible. Treat stale production-linked evidence as
   stale. Confirm that all scientific blockers remain honest: solar-only and
   permission-limited Sukhbold source, missing age/decay/momentum semantics,
   unapproved transition seams, and no project physics approval.

Return exactly two verdicts:

- `G2 SOURCE-PACKAGE ENGINEERING FINAL VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.
- `G2 SOURCE-PACKAGE SCIENTIFIC FINAL VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.

List any remaining actionable defect with file/line evidence. Do not edit the
workspace. The scientific verdict must remain BLOCK while no approved physical
source package exists.
