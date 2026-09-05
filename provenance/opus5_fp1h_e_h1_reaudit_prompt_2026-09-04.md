# Claude Opus 5 H1 remediation re-audit — F-P1H-E

Perform a read-only re-audit in `/gpfs/kjhan/LRD_JWST`. Do not edit files,
run jobs, build RAMSES, or launch a simulation. Inspect only the H1 files
listed in the previous prompt and the new H1 audit record:

- `simulation/snrt/tools/fp1_gate_validator_blocks.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/tests/fp1_physical_package_admission.py`
- `simulation/snrt/data/fp1_physical_package_admission_audit.json`
- `simulation/snrt/data/fp1_fate_admission_audit.json`
- `provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md`
- `provenance/opus5_fp1h_e_h1_audit_2026-09-04.md`

The first Opus audit returned `CONDITIONAL PASS` and requested: (a) an
explicit equality assertion between the imported source-identity validator ID
and `GATE_VALIDATOR_IDS["source_identity_and_rights"]`, and (b) lowercase-only
SHA-256 validation. Both fixes are now present in
`fp1_gate_validator_registry.py`. The driver also reran, on GPFS, Python
compilation, the focused admission test (`FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK`),
and `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh`, which
ended `FP1_POPULATION_FATE_CONTRACT_OK`; the generated state remains
fail-closed with four candidates, zero physical nodes, null selection, and
false production/publication/conversion/runtime-deposition flags.

Verify specifically that the requested fixes are correct, that the generated
code hashes and sidecar are fresh, and that the original H1 conclusions still
hold: exact nine-gate registry/contract identity, eight never-passing
unavailable adapters, preserved LC18 rights result, no declarative-status
bypass, and no physical promotion. H5 will add the remaining adversarial
tests for unregistered, mis-bound, stale-hash, and malformed-report paths.
Do not demand those H5 tests as a condition of this narrow remediation
re-audit unless the current H1 change itself introduces a defect.

Return a concise severity-ranked report with file/line evidence and one clear
verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`. No edits.
