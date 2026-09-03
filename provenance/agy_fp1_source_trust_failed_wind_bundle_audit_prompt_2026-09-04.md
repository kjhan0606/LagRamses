# AGY bundle-end audit prompt — F-P1 source trust and LC18 cross-check

You are conducting an independent, read-only implementation audit of the
completed bundle in `/gpfs/kjhan/LRD_JWST`.

Use model `gemini-3.8-flash-high`. Do not edit files, write generated artifacts,
download data, contact authors, launch RAMSES or other jobs, or treat the
current dirty worktree as committed source. You may inspect files, hashes, git
diffs, and run narrowly scoped read-only tests if the CLI sandbox permits it.

## Final project purpose

lagRamses is being developed as a production-ready and publication-ready
high-level hydrodynamics stack focused on radiative transfer, stellar/AGN
feedback, and dust. The current bundle concerns the stellar-feedback source
admission boundary and evidence quality. A review-only result must not be
presented as physical approval or runtime readiness.

## Audit scope

Audit the implementation delta from the declared base commit `296dd0c`, with
special attention to:

1. `simulation/snrt/tools/validate_fp1_source_identity_rights.py` and
   `simulation/snrt/tools/fp1_gate_validator_registry.py`: code-owned trust
   root, exact candidate/release identity, inventory and hash binding, license
   evidence, type/date/path validation, symlink/non-regular-file rejection,
   malformed-input and exception handling, and resistance to a coherent
   manifest/contract/bytes rewrite.
2. `simulation/snrt/tools/audit_fp1_fate_admission.py`, its sidecar and physical
   package contract: exact repository-relative artifact binding, publication
   invariant, validator registry semantics, and whether an uncommitted or
   mutable artifact can be promoted.
3. `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py` and its
   test/JSON evidence: reuse of the existing parsers, exact 108/108 join,
   52/56 classification, CDS comparison, duplicate collapse, phase monotonicity,
   96/12 table7 handling, no silent correction, rights status, and preservation
   of the four physical hard blockers.
4. Adversarial tests and generated evidence. Check that the tests actually
   prove the claimed invariants and cannot pass merely because the fixture is
   self-consistent.
5. The expected fail-closed state: zero physical nodes and false
   `production_ready`, `publication_ready`, `canonical_conversion_allowed`,
   and `runtime_deposition_allowed` are intentional. Do not call that state a
   defect unless an implementation path can bypass it.

## Required report

Return a self-contained audit with:

- verdict: PASS, CONDITIONAL PASS, or FAIL;
- findings ordered by severity, each with file/line or exact symbol, concrete
  evidence, exploit/reproduction path, and an actionable disposition;
- independent numerical/physical sanity checks, explicitly separating
  release-internal control discrepancies from cross-source discrepancies;
- any claims that are not reproducible from the staged local data;
- residual risks that are correctly deferred to the next bundle;
- a clear statement of whether the bundle is complete enough to let the driver
  design the next implementation bundle.

Do not infer author confirmation, source licensing, age-resolved wind,
per-node energy mapping, canonical momentum/deposition, physical nodes, or
production approval from schemas or synthetic fixtures.
