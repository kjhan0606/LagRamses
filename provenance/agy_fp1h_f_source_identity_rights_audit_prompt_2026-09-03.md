# Independent audit request: F-P1H-F source identity/rights validator bundle

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit files, commit, or launch a
simulation. Use model `gemini-3.8-flash-high` and inspect the implementation
directly rather than trusting the implementation note or generated reports.

Scope only the newly completed executable-validator bundle relative to commit
`296dd0c`:

- `simulation/snrt/tools/validate_fp1_source_identity_rights.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- related F-P1 tests, generated audits, roadmap, and bundle note

The intended result is narrow: only the Boccioli--Roberti 2026 candidate may
pass `source_identity_and_rights`. It must remain production-unqualified with
eight missing gates, four explicit physical blockers, zero physical nodes, and
all canonical-conversion/runtime/production/publication approvals false.

Independently verify:

1. The executable registry cannot be bypassed with a self-declared validator,
   hash-only artifact, unknown ID, wrong gate, malformed result, or altered
   contract validator list.
2. Candidate identity cannot be substituted across the admission candidate,
   acquisition manifest, candidate source contract, staged release root,
   version-specific Zenodo record, and use-terms record.
3. The validator actually re-hashes every manifest-scoped local source file,
   verifies byte counts and contract/manifest SHA256 parity, and computes the
   stated composite fingerprint without trusting a stale audit JSON.
4. Absolute paths, `..`, duplicate paths, symlinks, missing files, malformed
   JSON, altered DOI/version metadata, license downgrade, and redistribution
   downgrade fail closed. Identify any concrete untested bypass.
5. CC-BY-4.0 and version-record evidence justify this rights-only gate, without
   being mistaken for physics approval or resolving the failed-wind anomaly.
6. Sidecar/hash wiring, generated reports, and F-P1/G2 runners are consistent.
7. The implementation remains reproducible and appropriately fail-closed when
   the external staged package is absent from a clean checkout.

Run focused read-only tests or construct temporary adversarial fixtures under
`/tmp` if useful. Report:

- verdict: PASS, CONDITIONAL PASS, or FAIL;
- findings ordered by severity, with exact file/line references and a concrete
  exploit or failure mode;
- whether the five rights requirements are genuinely established;
- whether production/publication remains correctly blocked;
- minimal required remediation for every nontrivial finding.

Do not expand into unrelated RT, RAMSES topology, or previously closed gates.
