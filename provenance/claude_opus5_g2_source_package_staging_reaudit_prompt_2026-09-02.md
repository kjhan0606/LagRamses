# Claude Opus 5 G2 source-package staging re-audit request

Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, commit, push,
launch a RAMSES time integration, or promote a candidate into canonical
production tables. This is the bundled re-audit after the dispositions from
`provenance/claude_opus5_g2_source_package_staging_audit_2026-09-02.md`.

Read the previous audit and the current implementation:

- `provenance/claude_opus5_g2_source_package_staging_audit_2026-09-02.md`
- `provenance/fp1_source_package_selection_plan_2026-09-02.md`
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

Independently check:

1. A mutated acquisition-manifest fingerprint or Limongi/NuGrid parser error
   propagates to the aggregate top-level status and a non-zero CLI exit code.
2. The projection rejects a swapped component/source-column contract, retains
   source-package file fingerprints, and does not claim exact mass closure.
3. W18/N20 high-mass yield tables and W18 implosion-wind records are parsed as
   review-only source components; missing branch masses are explicit, source
   components are not interpolated, and no canonical rows/runtime deposition
   can result.
4. The Sukhbold source remains solar-only, permission-limited, and scientifically
   incomplete for production: no age-resolved wind history, complete decay,
   canonical momentum, transition-seam ownership, or project approval.
5. Run bounded tests if possible. Do not treat a stale production-linked build
   artifact as current evidence.

Return exactly:

- `G2 SOURCE-PACKAGE ENGINEERING REAUDIT VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.
- `G2 SOURCE-PACKAGE SCIENTIFIC REAUDIT VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.

List any remaining actionable defect with file/line evidence. Distinguish code
defects from intentional scientific/licensing blockers. The scientific verdict
must remain BLOCK while no project-approved source package exists.
