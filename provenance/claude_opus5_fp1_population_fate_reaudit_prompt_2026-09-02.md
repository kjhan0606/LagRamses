# Claude Opus 5 F-P1 consolidated re-audit request

Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, commit, launch a
simulation, or expand into generic HDF5/AMR/ksection/CPU-box work.  This is one
consolidated re-audit of the complete D1--D5 repair bundle, not five separate
micro-audits.

Read first:

- `provenance/claude_opus5_fp1_population_fate_audit_2026-09-02.md`
- `provenance/fp1_population_fate_contract_2026-09-02.md`
- `provenance/feedback_population_dtd_active_roadmap.md`

Then independently inspect the relevant production files in
`patch/lagRamses`, exact/shared native mirrors, tests, and these contracts:

- `simulation/snrt/config/g2_physics_contract_v1.json`
- `simulation/snrt/config/stellar_source_identity_v1.json`
- `simulation/snrt/data/p0_production_linked_build_evidence.json`

Re-evaluate every original finding:

1. D1: Chabrier continuity and normalization, including whether the JAX
   differential is now independent enough to catch a Fortran antiderivative
   error rather than repeating it.
2. D2: `sum(tracked)<=returned` residual semantics and whether production
   generic metallicity receives tracked metals plus the residual even when an
   individual metal field is disabled.
3. D3: whether the population ledger is mandatory on the production timestep
   source path, uses already computed cumulative states without another IMF
   convolution, rejects bad ownership/over-return, and is cross-checked against
   RAMSES particle mass before writes.
4. D4: production/native strict audit equivalence and coverage of duplicate
   coordinates, Cartesian completeness, age-zero anchors, cumulative material
   and energy monotonicity, mass bounds, and remnant ownership.
5. D5: whether a successful legacy namelist read transactionally retains its
   element/channel switches and whether failed channel-resolved reads preserve
   prior global configuration.

Also audit for regressions introduced by the repair bundle, especially API
call compatibility, zero-width increments, fail-closed clearing, unit
conversion, and source/particle mass semantics.  You may run builds and unit or
startup-negative tests, but do not run RAMSES time integration.

Known scientific blockers are intentionally not papered over: no immutable
population/fate source is approved, the provisional 40--120 Msun interval has
winds but no terminal fate, and PISN/PPISN remnant semantics await F-P3.  Judge
whether these are honestly reported and kept fail-closed.  Do not invent a
physical choice.

Report two explicit verdicts:

- `D1--D5 ENGINEERING VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.
- `OVERALL F-P1 SCIENTIFIC VERDICT`: PASS, CONDITIONAL PASS, or BLOCK.

Give file/line evidence for any remaining defect.  Distinguish an actionable
in-scope implementation defect from an acknowledged source/model decision.
Keep RT, dust, AGN, SNIa DTD implementation, PISN implementation, and generic
infrastructure out of scope unless they directly invalidate an F-P1 claim.
