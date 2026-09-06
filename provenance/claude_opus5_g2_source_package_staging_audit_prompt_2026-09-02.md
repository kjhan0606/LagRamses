# Claude Opus 5 G2 source-package staging bundle audit request

Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, commit, push,
launch a RAMSES time integration, or promote any candidate into canonical
production tables. This is one bundled audit of the G2 source-package
selection and review-only staging work for the F-P1 40--120 M☉ plan.

Read first:

- `provenance/fp1_source_package_selection_plan_2026-09-02.md`
- `provenance/fp1_mass40_120_application_record_2026-09-02.md`
- `provenance/feedback_population_dtd_active_roadmap.md`
- `simulation/snrt/config/g2_source_selection_matrix_v1.json`
- `external/g2_candidates/README.md`
- `external/g2_candidates/acquisition_manifest_v1.json`

Inspect the complete G2 candidate audit and its candidate-specific auditors:

- `simulation/snrt/tools/audit_g2_candidate_sources.py`
- `simulation/snrt/data/g2_candidate_source_audit.json`
- `simulation/snrt/tools/audit_g2_sukhbold2016_candidate.py`
- `simulation/snrt/tools/build_g2_sukhbold_channel_projection.py`
- `simulation/snrt/config/g2_sukhbold2016_candidate_contract_v1.json`
- `simulation/snrt/config/g2_sukhbold_channel_projection_contract_v1.json`
- `simulation/snrt/data/g2_sukhbold2016_candidate_audit.json`
- `simulation/snrt/data/g2_sukhbold_channel_projection_review.json`
- `simulation/snrt/tests/g2_candidate_sources.py`
- `simulation/snrt/tests/g2_sukhbold2016_candidate.py`
- `simulation/snrt/tests/g2_sukhbold_channel_projection.py`

Audit these claims independently:

1. The staged files are treated as immutable review inputs with exact source
   fingerprints, license/redistribution terms, and candidate identity. Check
   that the all-candidate audit is fail-closed and emits no canonical rows.
2. The recommended staging order (Sukhbold first as solar engine-specific
   validation, LC18 next as a multi-Z/rotation comparison) is scientifically
   defensible as a comparison plan, but is not falsely presented as a project
   source choice. Check that no public availability is used as physics
   approval.
3. The Sukhbold adapter/projection preserves source components and nulls:
   integrated presupernova wind versus terminal CCSN ejecta, stable tracked
   elements plus untracked mass, selected radioactive inventory, lifetime and
   release age, decay-complete mass, and canonical momentum. Check exact
   source-node identity, component ownership, closure arithmetic, and
   resistance to accidental canonical conversion or runtime deposition.
4. Check the physical semantics and seams: Sukhbold is solar-only, begins at
   9 M☉, uses one Z9.6 engine branch, and does not supply age-resolved wind,
   complete decay, canonical momentum, or a project-approved transition to
   Stockinger/F23. LC18 anomalies and missing energy/momentum/age histories
   must remain blockers; no interpolation or invented values may be allowed.
5. Run bounded candidate/projection tests if the session permits. Distinguish
   a code defect from the intentional scientific/permission blockers. Do not
   claim production-linked build evidence unless it is fresh and actually
   linked to this bundle.

Return exactly two verdict lines:

- `G2 SOURCE-PACKAGE ENGINEERING VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.
- `G2 SOURCE-PACKAGE SCIENTIFIC VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.

Give file/line evidence for every remaining defect, state whether it is
actionable in this bundle, and state the next in-scope step. The scientific
gate must remain blocked while no source package has project physics approval
and while 40--120 M☉ physical fate/yield/remnant/lifetime/energy/momentum
values are absent.
