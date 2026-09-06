# AGY F-P1 identity/publication closure bundle audit request

You are the independent scientific and software auditor for the completed
F-P1 identity and publication closure bundle in the lagRamses/SNRT project.
Work read-only in `/gpfs/kjhan/LRD_JWST`. Do not edit files, write repository
artifacts, launch RAMSES/HPC jobs, select a physical yield source, contact an
author, or redistribute CDS data. Inspect the current checkout at HEAD
`5aeb6d3`; implementation code is in `25bd05f`.

The project's final purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. This bundle is an integrity/admission
boundary, not physical source activation. The current review-only blocked
state is intentional: zero physical source nodes, unresolved fate seams at
`[0.8,1.0]` and `[40,120] M_sun`, and false production/publication/runtime
deposition flags must remain fail-closed.

Return one verdict exactly: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
`PASS` means no remaining F-P1 blocker. A conditional pass may identify
non-blocking follow-up work. Use `BLOCK` only for an implementation defect or
missing F-P1 acceptance evidence that makes the boundary unsound. Do not
turn intentionally deferred physical-source, license, high-mass-fate, or
runtime-activation work into a defect of this bundle unless the bundle claims
to have solved it.

Audit algorithmic correctness, data-flow/wiring, provenance binding, claim
scope, scientific/numerical legitimacy, fail-closed behavior, and the
reproducibility of the evidence. Independently inspect the implementation and
run only bounded, read-only/fixture tests as appropriate. If a command would
write tracked files, inspect its logic instead of running it or use an
isolated temporary output. Do not rely only on the checked-in report labels.

## Required F-P1 acceptance checks

1. Validator admission: a passed identity validator must carry a valid
   package fingerprint; positive selection must require all nine
   identity-matched passed validator reports, bind the selected package hash
   to executable source-identity evidence, and preserve the blocked real
   selection state.
2. Source-node mapping: the shared canonical serializer must define the exact
   bytes hashed, normalize numeric values deterministically, reject invalid or
   duplicate coordinates, and bind schema, row count, source-node coverage,
   source-node contract hash/approval, package hash, and asset hash. The
   converter must compare the generated mapping to the admitted mapping before
   writing; proposal mode must be non-writing.
3. Publication gate: derived LC18 artifacts must be blocked unless the
   code-locked terms path and bytes, verified production license, source and
   derived redistribution permissions, attribution, explicit approval
   identity, artifact/candidate identity, and `review_use_only=false` all
   pass. A mutable report label must not open the gate. Check that the current
   terms evidence is not silently upgraded.
4. LC18 diagnostics: verify the symmetric endpoint accounting and its meaning:
   48 positive/4 zero successful controls, 53 positive/3 zero failed models,
   101 positive/7 zero total, without changing the unresolved anomaly policy.
5. Evidence integrity: verify deterministic JSON regeneration, fixture
   hash-invariance, compilation/tests, exact contract hash propagation, and
   that no production runtime path has been enabled accidentally.

## Files to inspect first

- `provenance/fp1_bundle_plan_identity_publication_closure_2026-09-04.md`
- `provenance/fp1_identity_publication_closure_bundle_2026-09-04.md`
- `simulation/snrt/tools/fp1_source_node_mapping.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`
- `simulation/snrt/tools/convert_yield_rows_to_canonical.py`
- `simulation/snrt/tools/fp1_publication_rights.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- the corresponding F-P1 tests and `simulation/snrt/data/*.json` reports
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/config/g2_source_selection_matrix_v1.json`

Report findings in severity order. For every non-PASS item give exact
file/line evidence, physical or algorithmic impact, whether it is an F-P1
blocker or a later-gate deferral, and a concrete remedy. Explicitly state
which of the five acceptance checks were independently reproduced. End with
the safe next-bundle constraint: the next bundle must wait for driver planning
and Fable approval under the project governance.
