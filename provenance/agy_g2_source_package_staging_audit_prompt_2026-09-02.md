# AGY independent G2 source-package staging audit

You are the independent AGY/Gemini reviewer. Work read-only in
`/gpfs/kjhan/LRD_JWST`. Do not edit files, commit, push, launch RAMSES, or
promote any source candidate. Review the current tree itself; do not rely on
another auditor's verdict.

This is an algorithm, wiring, and scientific-semantics audit of the G2
review-only source-package staging for the F-P1 40--120 M☉ program. Read:

- `provenance/fp1_source_package_selection_plan_2026-09-02.md`
- `simulation/snrt/tools/audit_g2_candidate_sources.py`
- `simulation/snrt/tools/audit_g2_sukhbold2016_candidate.py`
- `simulation/snrt/tools/build_g2_sukhbold_channel_projection.py`
- `simulation/snrt/tests/g2_candidate_sources.py`
- `simulation/snrt/tests/g2_sukhbold2016_candidate.py`
- `simulation/snrt/tests/g2_sukhbold_channel_projection.py`
- `simulation/snrt/config/g2_sukhbold2016_candidate_contract_v1.json`
- `simulation/snrt/config/g2_sukhbold_channel_projection_contract_v1.json`
- `simulation/snrt/data/g2_candidate_source_audit.json`
- `simulation/snrt/data/g2_sukhbold2016_candidate_audit.json`
- `simulation/snrt/data/g2_sukhbold_channel_projection_review.json`

Audit independently:

1. Whether acquisition-manifest and Limongi/NuGrid parser integrity failures
   propagate to an aggregate failure status and nonzero CLI exit, while clean
   candidate input remains review-only.
2. Whether the projection contract is fail-closed against missing safety keys,
   missing components, and source-column/channel/ownership swaps.
3. Whether wind-only implosion tables preserve only quantities actually
   present, derive nonnegative/wind-only claims from parsed input, preserve
   radioactive overlap markers, and never manufacture terminal ejecta.
4. Whether W18/N20 high-mass yield tables and W18 implosion-wind vectors are
   parsed and provenance-labelled without interpolation, with explicit missing
   masses and source fingerprints.
5. Whether any path can emit canonical rows, enable runtime deposition, or
   falsely imply a production source. Treat all missing age-resolved wind,
   decay, momentum, transition-seam, licensing, and project-approval fields as
   scientific blockers, not values to be guessed.
6. If possible, run only bounded tests. Do not use stale production-linked
   evidence as a current build proof.

Report exact file/line evidence and distinguish defects from intentional
scientific blockers. Return exactly these two lines:

`AGY G2 SOURCE-PACKAGE ENGINEERING VERDICT`: PASS, CONDITIONAL PASS, or BLOCK

`AGY G2 SOURCE-PACKAGE SCIENTIFIC VERDICT`: PASS, CONDITIONAL PASS, or BLOCK

The scientific verdict must remain BLOCK while there is no project-approved
physical source package and no complete 40--120 M☉ physical fate/yield,
remnant, lifetime, energy, or momentum prescription.
