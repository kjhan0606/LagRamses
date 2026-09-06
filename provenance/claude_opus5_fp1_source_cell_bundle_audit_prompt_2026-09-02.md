# Claude Opus 5 F-P1 source-cell bundle audit request

Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, commit, push,
launch a RAMSES time integration, or expand into generic HDF5/AMR/ksection/
CPU-box work. This is one bundled audit of the just-completed F-P1 runtime
source-cell/admission work.

Read first:

- `provenance/fp1_mass40_120_application_record_2026-09-02.md`
- `provenance/feedback_population_dtd_active_roadmap.md`
- `provenance/claude_opus5_fp1_population_fate_reaudit_2026-09-02.md`

Inspect the exact production files in `patch/lagRamses`, their native mirror in
`simulation/snrt/native/phase0`, the F-P1 tools/tests, and:

- `simulation/snrt/config/fp1_population_fate_map_v1.json`
- `simulation/snrt/config/fp1_fate_resolver_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/config/g2_physics_contract_v1.json`
- `simulation/snrt/config/stellar_feedback_contract_v1.json`

Audit these claims independently:

1. The table's piecewise-constant source-cell mode implements the declared
   half-open `[m_i,m_{i+1})` convention and exact upper-edge behavior, rejects
   invalid modes, and is actually selected by the production runtime after a
   successful table load. Confirm that this mode does not falsely authorize
   metallicity/rotation/engine/age interpolation.
2. The unresolved mass bucket is computed from the configured IMF and explicit
   non-overlapping F-P1 intervals, is reflected in the population ledger, does
   not enter returned/living/remnant closure, and is never deposited as
   feedback. Check malformed/overlapping intervals and IMF-support mismatch
   fail closed.
3. The admission sidecar verifies SHA256 for fate map, resolver contract,
   source contract, and physics contract; compares its unresolved intervals
   and approval id to the map; and keeps production/canonical conversion
   blocked while the source has zero physical nodes or unresolved intervals.
   Check that the sidecar test does not merely test its own constants.
4. Check production/native source parity and API compatibility for the new
   ledger and IMF routine. Run the bounded native/population tests and
   differential test if useful, but do not claim a production-linked PASS
   unless the real production-linked build evidence is fresh and valid.
5. Review the 40--120 Msun scientific status: no fate/yield/remnant/energy/
   momentum/lifetime values may be invented, and the next required step must
   remain licensed source-package selection and physical-node staging.

Explicitly distinguish a code defect from the acknowledged production
blockers (zero physical fate nodes, unresolved 0.8--1 and 40--120 Msun
intervals, and the pre-existing production-linked evidence requirement).

Report:

- `F-P1 SOURCE-CELL/ADMISSION ENGINEERING VERDICT`: PASS, CONDITIONAL PASS,
  or BLOCK.
- `OVERALL F-P1 SCIENTIFIC VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.

Give file/line evidence for every remaining defect and state whether it is
actionable in this bundle. Do not edit the workspace.
