# Claude Opus 5 implementation-stage audit — F-P1H-E H2/H3

Perform a read-only audit in /gpfs/kjhan/LRD_JWST. Do not edit files, run
jobs, build RAMSES, or launch a simulation. Inspect the following exact files
and their direct imports only:

- simulation/snrt/tools/audit_fp1_source_node_contract.py
- simulation/snrt/tools/validate_fp1_source_identity_rights.py
- simulation/snrt/tools/fp1_gate_validator_registry.py
- simulation/snrt/tools/audit_fp1_physical_package_admission.py
- simulation/snrt/config/fp1_source_node_contract_v1.json
- simulation/snrt/config/fp1_physical_package_admission_contract_v1.json
- simulation/snrt/data/fp1_source_node_contract_audit.json
- simulation/snrt/data/fp1_physical_package_admission_audit.json
- simulation/snrt/data/fp1_fate_admission_audit.json
- simulation/snrt/tests/fp1_source_node_contract.py
- simulation/snrt/tests/fp1_physical_package_admission.py
- provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md

This is the combined H2/H3 stage of the fail-closed physical stellar-source
admission gate for the production/publication-ready lagRamses high-level
hydro project (RT, stellar/AGN feedback, dust, and coupled source terms).
The current repository must remain review-only: four review candidates, zero
physical nodes/canonical rows, null selection, false production/publication/
conversion/runtime-deposition flags, and unresolved [0.8,1.0) and
[40,120] M_sun seams.

## H2 claims to verify

1. An approved_physical_nodes source-node contract cannot become eligible
   from node-local research_use_status or redistribution_status strings.
   Every node must map its source_id uniquely through the code-owned
   LOCKED_CANDIDATE_PROFILES to the registered
   fp1.source_identity_and_rights.v1 validator.
2. The executed rights validator must pass for the same candidate, expose
   verified source-file bytes, and return a package fingerprint identical to
   the node's fingerprint. A missing/blocked validator, unknown source, or
   different package must fail before promotion.
3. Review-only/empty node contracts must remain valid without falsely
   claiming an approved rights binding. Existing LC18 rights and license
   boundaries must not be widened into redistribution authorization.

## H3 claims to verify

1. Every current candidate admission record names exactly all nine validator
   IDs; the generated candidate report records each gate as one of
   missing_evidence, validator_blocked, validator_error, or pass.
2. validator_error and missing evidence cannot qualify a candidate, and only
   nine executed pass results plus the existing package/node/upstream guards
   can select a physical package. A blocked adapter is not silently treated
   as a pass.
3. Current LC18 remains pass only for source identity and validator_blocked
   for the other eight; the other candidates remain blocked for all nine.
   Generated evidence is fresh and the real physical-package state is
   unchanged.

The driver reports these GPFS checks as passed:

- python3 simulation/snrt/tests/fp1_source_node_contract.py
  -> FP1_SOURCE_NODE_CONTRACT_TEST_OK
- python3 simulation/snrt/tests/fp1_physical_package_admission.py
  -> FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK
- bash simulation/snrt/tests/run_fp1_population_fate_contract.sh
  -> FP1_POPULATION_FATE_CONTRACT_OK

Review algorithmic/wiring correctness, trust-root boundaries, exception and
fingerprint semantics, status classification, selection guards, test
adequacy, and final-purpose fit. Do not demand unrelated AMR/HDF5 work or
physical source data in this stage. Return severity-ranked findings with
file/line evidence and one verdict: PASS, CONDITIONAL PASS, or BLOCK.
No edits.
