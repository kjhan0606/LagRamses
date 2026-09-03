# F-P1H-F source identity and rights validator bundle

Date: 2026-09-03
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Scope: first executable physical-package admission validator

## Provisional implementation outcome

The package-admission gate now has a code-owned validator registry. Candidate
contracts can reference only validator identifiers that are both present in
the Python registry and listed in the admission contract. A hash-only or
unexecuted evidence artifact cannot pass a gate.

**Post-implementation audit status: superseded by the 2026-09-04 remediation
bundle; bundle-end audit pending.** AGY found a passing candidate-identity
substitution, and an independent Codex `gpt-5.6-sol` re-audit reproduced it and
found additional trust-root bypasses. The generated rights-gate PASS in this
original bundle was provisional test output, not valid approval evidence. The
trust-root and adversarial-test remediation is recorded in
`fp1_source_trust_failed_wind_bundle_2026-09-04.md`; see
`agy_fp1h_f_source_identity_rights_audit_2026-09-03.md`,
`gpt56sol_fp1h_f_source_identity_rights_reaudit_2026-09-03.md`, and
`fp1h_f_source_identity_rights_audit_comparison_2026-09-03.md`.

`fp1.source_identity_and_rights.v1` independently re-reads the staged
Boccioli--Roberti 2026 source package and verifies all five declared gate
requirements:

- article citation and version-specific Zenodo record/DOI identity;
- manifest, source contract, byte counts, per-file SHA256, and composite
  package fingerprint;
- machine-readable `cc-by-4.0` metadata from the staged Zenodo record;
- explicit CC-BY-4.0 redistribution permission and attribution terms;
- a confined local source mirror whose bytes still match the acquisition
  manifest.

The computed package fingerprint is
`3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`.
The uncorrected implementation reports that only
`boccioli_roberti2026_lc18` passes this first gate. That pass is not accepted
because candidate identity and release bytes are not independently pinned. It remains
non-qualified because eight executable gates and its existing physical
blockers remain open. No source node was added, and canonical conversion,
runtime deposition, production, and publication flags remain false.

## Initial fail-closed tests

The automated tests reject:

- an unknown or contract-unapproved validator identifier;
- a validator record with undeclared hash-only evidence fields;
- a contract that disables or diverges from the code registry;
- altered source-package bytes;
- downgraded redistribution terms;
- substituted version-specific data DOI;
- unsupported candidates without a code-registered rights profile.

The audits showed that this list is insufficient: coherent manifest/contract/
byte rewrites, empty inventories, internal symlinks, candidate/root mismatch,
wrong scalar types, malformed dates, mutable rights assertions, external fate
sidecar artifacts, and unexpected exceptions require additional coverage.

## Executed evidence

- `simulation/snrt/tests/fp1_source_identity_rights.py`: PASS
- `simulation/snrt/tests/fp1_physical_package_admission.py`: PASS
- `simulation/snrt/tests/run_fp1_population_fate_contract.sh`: PASS with the
  expected physical-package block
- `simulation/snrt/tests/run_g2_preflight.sh`: all subordinate tests PASS and
  terminates with the expected `G2_PREFLIGHT_BLOCKED`
- `git diff --check`: PASS

## Remaining promotion blockers

The next executable-validator bundles are coordinate hull/population, fate and
remnant structure, lifetime/wind history, terminal mass/species closure,
decay projection, energy/momentum/deposition, pair instability, and runtime
invariance/reproduction. Boccioli--Roberti also remains quarantined pending a
resolution of the failed-model wind-table anomaly and complete age-resolved
wind, per-node injected-energy, momentum, and deposition evidence. Physical
node inventory remains empty and the RAMSES source-node consumer is absent.

Per the audit cadence, AGY `gemini-3.8-flash-high` audited this bundle and the
negative result was independently re-audited with Codex `gpt-5.6-sol`. The
matching findings are mandatory inputs to the next bundle plan. Grok was
excluded because the xAI service was unavailable and the user directed the
project to continue without it.
