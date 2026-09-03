# Claude Opus 5 F-P1 source-cell/admission bundle audit — 2026-09-02

Auditor: Claude Opus 5 CLI. Scope was read-only review of the source-cell,
unresolved-bucket, and admission-sidecar bundle. No RAMSES time integration
was launched.

## Empirical-verification limitation

The auditor's session had Bash disabled for the auditor and its subagents, so
it did not execute tests, builds, or SHA256 commands. The project-side
verification was run independently before and after the audit:

- `run_fp1_population_fate_contract.sh`: PASS;
- `run_g2_population_ledger.sh`: PASS;
- `tests/run_stellar_population_contract_unit.sh`: PASS;
- `P0_DIAGNOSTIC=1 run_g1_native_contract.sh`: PASS, including the JAX 0.11.1
  CPU differential;
- resolver, candidate-source, Python/JSON, and `git diff --check`: PASS.

The recorded production-linked build evidence predates this bundle and is not
used as a current production-linked PASS.

## Auditor findings and disposition

The initial Opus review returned:

`F-P1 SOURCE-CELL/ADMISSION ENGINEERING VERDICT: CONDITIONAL PASS`

`OVERALL F-P1 SCIENTIFIC VERDICT: BLOCK`

The scientific BLOCK is correct and remains active: the fate resolver has zero
physical nodes, the 0.8–1 and 40–120 M☉ intervals are unresolved, and no
source package/license/checksum/approval has been promoted.

The actionable engineering conditions were:

1. **Order-dependent closed edge in the resolver.** Fixed in
   `simulation/snrt/tools/fp1_fate_resolver.py`: the closed upper edge is now
   assigned to the matching cell with the largest upper bound, independent of
   JSON record order. A reversed-node regression test was added.
2. **Resolver approval fields had no reader.** Fixed in
   `simulation/snrt/tools/audit_fp1_fate_admission.py`: resolver approval,
   canonical-conversion, runtime-deposition, and approval-id fields are now
   read and checked against the sidecar/map.
3. **Fortran unresolved intervals were not cross-checked against JSON.** Fixed
   in the same admission audit: both the production and native
   `stellar_enrichment_config.f90` mirrors are parsed and compared to the
   audited F-P1 intervals.
4. **Production/native table-reader error-order divergence.** Fixed by making
   the age conversion occur after physical range validation in both mirrors;
   doubly-invalid rows now have the same error precedence.

The auditor also confirmed that the declared native/production source split
matches `stellar_source_identity_v1.json`; the three RAMSES-only deposition
adapters are intentionally absent from the native mirror and are not a parity
defect.

## Current decision

The bundle is an engineering **CONDITIONAL PASS pending a future Opus audit of
the post-disposition tree**, with all actionable conditions above addressed
and independently tested. It is not a scientific pass. The next in-scope
work is licensed source-package selection and immutable physical source-node
staging; no 40–120 M☉ physical fate, yield, remnant, lifetime, decay, energy,
or momentum value may be fabricated meanwhile.
