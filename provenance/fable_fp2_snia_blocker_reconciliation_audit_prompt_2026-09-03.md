# Fable re-audit prompt — F-P2 blocker reconciliation

Audit `/gpfs/kjhan/LRD_JWST` read-only with the actual `fable` model. Do not
edit files, select a physical source, activate runtime, commit, or push. Do not
consult AGY or any other audit result.

Review the latest grouped F-P2 change: the DTD/event-ledger sources and exact
mirrors; `fp2_snia_dtd_contract_v1.json`; the event-source sidecar and its
admission validator; HESMA source audit, model comparison, profile estimator
comparison, and selection packet; the contract audit; and directly relevant
tests/tools. The change adds explicit activation prerequisites, a review-only
`promotion_requirements` schema, `n300c` source-anomaly quarantine
classification, physical-review status, and negative-path tests to the shell
runner. The full `bash simulation/snrt/tests/run_fp2_snia_dtd_contract.sh`
passes, including both negative-path Python tests. No model or physical event
values are selected; runtime activation and canonical conversion remain false.

Independently determine whether the change is algorithmically and physically
sound for review-only use, whether it improves or weakens fail-closed behavior,
whether `n300c`/`n1600c` treatment is justified, and whether the declared
promotion requirements honestly cover the unresolved WD-remnant debit,
momentum convention/deposition, IMF and realization policy, thermal coupling,
metallicity, portable provenance, and commit binding. Check for false claims,
stale hashes, and test coverage gaps.

Return a structured report with verdict (`PASS`, `CONDITIONAL PASS`, or
`BLOCK`), severity-ranked exact file/field findings, what is resolved versus
only declared, and required next actions. Do not modify anything.
