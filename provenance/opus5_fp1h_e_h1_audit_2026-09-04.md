# Claude Opus 5 H1 implementation audit — F-P1H-E

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5 CLI, read-only
Prompt: `provenance/opus5_fp1h_e_h1_audit_retry_prompt_2026-09-04.md`

## Verdict

**CONDITIONAL PASS.** H1 is accepted as a fail-closed prerequisite and H2
may proceed after the following H1 conditions are addressed:

1. Add the registry rejection-path tests required by H1/H5 for an
   unregistered gate, gate/validator mis-binding, stale validator code hash,
   and malformed validator report. The existing focused test covers only a
   subset of these paths; the controls themselves appear correct.
2. Add an explicit equality assertion that the imported source-identity
   validator ID equals `GATE_VALIDATOR_IDS["source_identity_and_rights"]`.
3. Restrict `_valid_sha256` to lowercase hexadecimal, matching the producer
   contract and avoiding a misleading later equality failure.

The two one-line hardening fixes are treated as H1 remediation before H2.
The broader rejection-path tests remain scheduled for H5, where the complete
fixture matrix will exercise the registry guards.

## Findings confirmed by the audit

- The nine gate IDs, nine validator IDs, and all requirement sets match the
  contract exactly; code hashes in the generated artifacts and sidecar match
  the current working-tree bytes.
- The eight unavailable adapters are visibly controlled blockers: they return
  `status: "blocked"`, `passed: false`, all requirements false, a null package
  fingerprint, and non-empty blockers. They cannot select a package or enable
  runtime deposition, and the source-identity gate is excluded from the
  placeholder adapter.
- All four candidates carry all nine validator identities. LC18 preserves the
  real source-rights result, while the remaining physical gates remain
  blocked. The generated state remains zero physical nodes/canonical rows,
  null selection, false approval flags, and unresolved `[0.8,1.0)` and
  `[40,120] M_sun` seams.
- H1 is appropriately bounded to executable admission wiring and does not
  invent physical values or perform AMR/HDF5/runtime work.

## Audit execution note

Opus inspected only the files named in the retry prompt and did not execute
the focused test or generator under its read-only tool restriction. The
driver independently ran the focused H1 test and the full F-P1 runner before
this audit; their results remain recorded in the bundle work log and generated
artifacts. No files, jobs, builds, or simulations were changed or launched by
the auditor.
