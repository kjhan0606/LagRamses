# Claude Opus 5 R3 implementation-stage audit

Act as the independent scientific, algorithmic, wiring, and implementation
auditor for completed F-P1R step R3 in
`/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Work read-only. Do not edit
files, run shell commands or jobs, select physical sources, contact authors,
or redistribute data. AGY is retired and must not be called or treated as an
active approval authority.

Audited implementation commit: `2921584` (`Enforce same-run high-mass
evidence freshness`). Read:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/audit_governance_amendment_2026-09-04.md`
- `simulation/snrt/tests/run_fp1_population_fate_contract.sh`
- `simulation/snrt/tests/fp1_high_mass_freshness.py`
- `simulation/snrt/tests/fp1_high_mass_seam.py`
- `simulation/snrt/tools/audit_fp1_high_mass_seam.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`
- `simulation/snrt/tools/audit_fp1_fate_admission.py`
- `simulation/snrt/data/fp1_high_mass_seam_review.json`
- `simulation/snrt/data/fp1_physical_package_admission_audit.json`
- `simulation/snrt/data/fp1_fate_admission_audit.json`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. R3 is a narrow evidence-execution gate: it
must prove that high-mass review evidence is regenerated before both
physical-package and fate admission consume it. It must not select a source,
create physical nodes, enable canonical conversion/runtime deposition, or
broaden into checkpoint, AMR, MPI, or generic hydro work.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` for stale evidence that can make the runner green, a missing
consumer binding, a hash comparison that can be bypassed, or a changed
fail-closed physics/admission state. Use `CONDITIONAL PASS` only for a
material but non-blocking evidence gap. Do not treat byte-identical
regeneration as a defect; it is the intended production state.

## R3 contract

1. `run_fp1_population_fate_contract.sh` must regenerate
   `data/fp1_high_mass_seam_review.json` before either admission consumer runs.
   The high-mass test/audit must occur before the population/fate admission
   sequence, and the old end-of-run-only placement must be gone.
2. In the same shell invocation, capture the post-regeneration SHA-256. After
   both admission JSONs are regenerated, assert that their high-mass evidence
   records bind that exact hash. The physical-package report stores the record
   under `evidence_artifacts.high_mass_review`; the fate-admission report
   stores it under `physical_package_contract.evidence_artifacts.high_mass_review`.
   Verify path, `sha256`, code lock, and contract-declared hash, not just a
   printed label.
3. The freshness check must be read-only, fail closed on missing/malformed or
   stale reports, and must not add a nonce or rewrite evidence. It must also
   preserve the expected terminal state: physical nodes zero,
   production/publication/conversion/deposition false, unresolved fate seams,
   and final `G2_PREFLIGHT_BLOCKED` in the broader preflight.
4. A second standalone high-mass regeneration must be byte-identical to the
   tracked review artifact. The driver reports that this `cmp` check passed,
   the full runner returned exit 0 with `FP1_HIGH_MASS_FRESHNESS_TEST_OK` and
   `FP1_POPULATION_FATE_CONTRACT_OK`, and the G2 preflight remains blocked.
   Assess these claims from the files/provenance and current lock structure;
   do not rerun the artifact-writing commands in this read-only audit.
5. No production export, source selection, runtime activation, or unrelated
   infrastructure changes may be present. AGY must appear only as retired or
   historical provenance; current implementation auditing is Opus 5.

Inspect for subtle ordering/wiring errors: high-mass audit after a hidden
consumer, a freshness script that validates only one report, wrong nesting for
the fate report, expected hash captured before regeneration, a hash field not
bound to actual bytes, path-only validation, stale data accepted because the
runner overwrites the wrong artifact, or a test that can pass under a changed
terminal state. Separate optional portability/style issues from gate failures.

End with the verdict, checks actually established, mandatory fixes (if any),
non-blocking findings, and a direct statement on whether the F-P1R bundle is
complete and whether a new bundle may be planned. AGY is retired and must not
be called.
