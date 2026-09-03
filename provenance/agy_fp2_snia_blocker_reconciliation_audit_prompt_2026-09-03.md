# AGY re-audit prompt — F-P2 blocker reconciliation

Audit `/gpfs/kjhan/LRD_JWST` read-only with model `gemini-3.8-flash-high`.
Do not edit files, select a physical source, activate runtime, commit, or push.
This is an independent re-audit of the latest grouped change.

The previous bundle runner passed before the change. The change now:

- adds the missing activation prerequisites to
  `simulation/snrt/config/fp2_snia_dtd_contract_v1.json`;
- declares a review-only `promotion_requirements` schema in
  `simulation/snrt/config/fp2_snia_event_source_approval_sidecar_v1.json` and
  validates it in `audit_fp2_snia_event_source_admission.py`;
- classifies the HESMA `n300c` gross discrepancy as
  `source_data_anomaly_requires_quarantine`, while retaining `n1600c` as a
  profile warning;
- adds source physical-review status and quarantine flags to the audit and
  selection packet;
- adds the two negative-path Python tests to
  `simulation/snrt/tests/run_fp2_snia_dtd_contract.sh`.

Read the relevant configs, tools, tests, generated reports, DTD/event-ledger
sources and mirrors. The runner
`bash simulation/snrt/tests/run_fp2_snia_dtd_contract.sh` now passes all tests,
including `FP2_SNIa_EVENT_SOURCE_ADMISSION_TEST_OK` and
`FP2_SNIa_DTD_CONTRACT_TEST_OK`. Runtime activation remains false and no model
or physical event value is selected.

Assess algorithm/wiring correctness, fail-closed behavior, schema consistency,
HESMA warning/quarantine semantics, and whether the new declarations are
honest about what is still unimplemented. Pay particular attention to whether
WD remnant-reservoir debit, momentum deposition, portable provenance/commit
binding, and canonical promotion remain correctly blocked rather than merely
declared. Return a concise structured verdict (`PASS`, `CONDITIONAL PASS`, or
`BLOCK`), severity-ranked findings with exact file/field references, and the
next required actions. Do not modify anything.
