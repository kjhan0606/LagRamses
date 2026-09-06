# Independent re-audit request: F-P1H-F source identity/rights bundle

Use `gpt-5.6-sol` to audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit,
commit, or launch simulations. This is an independent re-audit: inspect the
code and construct your own adversarial cases without relying on another
model's findings.

Scope the uncommitted executable-validator bundle relative to commit
`296dd0c`, especially:

- `simulation/snrt/tools/validate_fp1_source_identity_rights.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- related tests, generated audits, roadmap, and bundle note

The intended result is rights-only: `boccioli_roberti2026_lc18` may pass only
`source_identity_and_rights`, while eight gates, four physical blockers, zero
physical nodes, and every runtime/production/publication approval remain
closed.

Audit candidate identity binding, registry bypass resistance, actual source
byte hashing and composite fingerprinting, DOI/version and machine-readable
license checks, redistribution evidence, path traversal and every symlink
case, malformed/missing/wrongly typed input handling, absent external assets,
sidecar hashes, generated evidence, and F-P1/G2 runner behavior. Verify that
all failures are fail-closed rather than uncaught exceptions. Try concrete
temporary mutations under `/tmp` where useful.

Return PASS, CONDITIONAL PASS, or FAIL; list findings by severity with exact
file/line references and reproducible failure modes; state whether the five
rights requirements are truly established; confirm whether production and
publication remain blocked; and give minimal remediation. Do not expand into
unrelated RT or previously closed RAMSES infrastructure work.
