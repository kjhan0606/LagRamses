# Claude Opus 5 targeted H2/H3 audit retry — F-P1H-E

Read-only audit in /gpfs/kjhan/LRD_JWST. Do not edit files, run commands,
run jobs, build, or launch a simulation. Give the verdict quickly after
checking the specified ranges; do not recursively search or read large files
in full.

Read these focused ranges:

- simulation/snrt/tools/audit_fp1_source_node_contract.py:1-215 and 690-825
- simulation/snrt/tools/validate_fp1_source_identity_rights.py:1-40,
  356-370, 599-668
- simulation/snrt/tools/fp1_gate_validator_registry.py:25-185
- simulation/snrt/tools/audit_fp1_physical_package_admission.py:78-115,
  307-375, 440-490, 720-770
- simulation/snrt/config/fp1_source_node_contract_v1.json:1-25 and 195-208
- simulation/snrt/config/fp1_physical_package_admission_contract_v1.json:
  candidate_qualification section and gate_validation section only
- simulation/snrt/data/fp1_source_node_contract_audit.json:1-30
- simulation/snrt/data/fp1_physical_package_admission_audit.json:
  top-level state and candidate_qualification section only
- simulation/snrt/tests/fp1_source_node_contract.py:40-190
- simulation/snrt/tests/fp1_physical_package_admission.py:395-570

H2 must establish that approved nodes cannot rely on local rights strings:
source_id must map uniquely through code-owned LOCKED_CANDIDATE_PROFILES to
the registered source-identity validator; that validator must pass, expose
verified source bytes, and return a package fingerprint equal to the node's.
Unknown source, blocked/absent validator, and package mismatch must fail before
promotion. Empty/review-only nodes must remain non-approved.

H3 must establish that the checked-in four candidates each name exactly all
nine validator IDs. The generated report must include a status for every gate
using only missing_evidence, validator_blocked, validator_error, or pass.
Only all nine executed pass statuses plus the existing node/upstream/selection
guards can promote; current LC18 must be one pass plus eight
validator_blocked, and the other candidates nine validator_blocked.

Project purpose is production/publication-ready lagRamses high-level
hydrodynamics: RT, stellar/AGN feedback, dust, and coupled source terms.
Current state must stay review-only: zero physical nodes/canonical rows, null
selection, all promotion flags false, and unresolved [0.8,1.0) and
[40,120] M_sun seams. Do not demand unrelated AMR/HDF5 work or physical data.

The driver already reports on GPFS:
fp1_source_node_contract.py -> FP1_SOURCE_NODE_CONTRACT_TEST_OK
fp1_physical_package_admission.py -> FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK
run_fp1_population_fate_contract.sh -> FP1_POPULATION_FATE_CONTRACT_OK

Return only a concise severity-ranked audit and one verdict PASS, CONDITIONAL
PASS, or BLOCK. Note any limitation from the read-only scope. No edits.
