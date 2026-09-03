# Independent Grok audit: F-P1H-F source identity/rights bundle

Audit `/gpfs/kjhan/LRD_JWST` read-only with xAI Grok. Do not edit files,
commit, or launch simulations. Independently inspect the uncommitted bundle
relative to commit `296dd0c`; do not rely on another model's findings.

Primary scope:

- `simulation/snrt/tools/validate_fp1_source_identity_rights.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- associated tests, reports, roadmap, and bundle note

The intended result is rights-only. Boccioli--Roberti 2026 may pass only
`source_identity_and_rights`; eight gates, four physical blockers, zero
physical nodes, and all runtime/production/publication approvals must remain
closed.

Try concrete `/tmp` mutations. Audit identity binding, registry bypasses,
live source-byte hashing/composite fingerprinting, DOI/version/license/rights
checks, absolute and `..` paths, duplicate records, internal/external symlinks,
missing and malformed inputs, wrong types, unexpected exceptions, absent
external assets, hash sidecars, generated evidence, and runner behavior.

Return PASS, CONDITIONAL PASS, or FAIL; findings by severity with exact
file/line references and reproducible failure modes; whether the five rights
requirements are established; whether production/publication remain blocked;
and minimal remediation. Do not expand into unrelated RT or old RAMSES work.
