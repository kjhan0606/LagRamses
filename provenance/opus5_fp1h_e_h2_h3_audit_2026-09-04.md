# Claude Opus 5 H2/H3 implementation audit — F-P1H-E

Date: 2026-09-04
Project: /gpfs/kjhan/LRD_JWST (kjhan0606/LagRamses)
Auditor: Claude Opus 5 CLI, read-only
Prompt: provenance/opus5_fp1h_e_h2_h3_audit_prompt_2026-09-04.md

## Verdict

PASS. No H2 or H3 defect or unsafe promotion path was identified.

## H2 — executed rights binding

- For every node in an approved_physical_nodes contract, source_id is
  resolved through the code-owned LOCKED_CANDIDATE_PROFILES and must map to
  exactly one candidate.
- The registered fp1.source_identity_and_rights.v1 validator is executed
  for that candidate. Its report must pass registry identity/requirement/
  fingerprint checks, expose non-empty verified source-file records, and
  return a package fingerprint equal to the node fingerprint.
- Unknown sources, absent/blocked/mis-registered validators, and package
  mismatches fail before promotion. Node-local research or redistribution
  strings are not sufficient.
- Empty/review-only node contracts remain non-approved, and the checked-in
  state has physical_node_count: 0 and rights_bindings: [].
- Tests cover an unknown source, a package mismatch, and a blocked executed
  rights validator.

## H3 — complete candidate evidence and statuses

- All four candidates name exactly the nine validator IDs, enforced by the
  admission audit.
- Each generated candidate report carries a status for every required gate
  from the closed vocabulary missing_evidence, validator_blocked,
  validator_error, and pass.
- Qualification requires all nine pass statuses and no hard blockers.
  Selection additionally requires all nine verified reports, all-pass status,
  no validation errors, non-empty physical nodes, and the existing upstream
  guards.
- LC18 remains one pass for source identity plus eight validator_blocked; the
  other three remain nine validator_blocked.
- The generated repository state remains review-only: four candidates, zero
  physical nodes/canonical rows, null selection, false production/publication/
  conversion/runtime-deposition flags, and unresolved [0.8,1.0) and
  [40,120] M_sun seams.

## Low-severity observations

- missing_evidence is intentionally supported by the subset-tolerant helper
  and mutation tests, but the checked-in production contract rejects a
  candidate record that omits one of the nine IDs before that status can
  appear in the normal generated report.
- Review nodes are schema-checked but do not execute rights binding until the
  contract enters the approved physical-node state; they cannot promote.
- The eight unavailable adapters share one code file/hash; they remain
  explicit never-passing blockers and this is not a promotion defect.

## Driver verification

On GPFS, after H2/H3 implementation and artifact regeneration:

- python3 simulation/snrt/tests/fp1_source_node_contract.py ->
  FP1_SOURCE_NODE_CONTRACT_TEST_OK
- python3 simulation/snrt/tests/fp1_physical_package_admission.py ->
  FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK
- bash simulation/snrt/tests/run_fp1_population_fate_contract.sh ->
  FP1_POPULATION_FATE_CONTRACT_OK

Opus did not edit files, launch jobs, build RAMSES, or run simulations.
