# AGY bundle audit prompt — F-P2 SNIa/HESMA source admission

Audit `/gpfs/kjhan/LRD_JWST` read-only using your own scientific and code
judgment. Do not edit files, select a production model, activate runtime code,
commit, or push. This is an independent bundle audit, not a request to rerun
the whole project.

Scope:

1. `simulation/snrt/native/phase0/stellar_snia_dtd.f90`
2. `simulation/snrt/native/phase0/stellar_snia_event_ledger.f90`
3. Matching production mirrors under `patch/lagRamses/`
4. `simulation/snrt/config/fp2_snia_dtd_contract_v1.json`
5. `simulation/snrt/config/fp2_snia_event_source_approval_sidecar_v1.json`
6. `simulation/snrt/data/fp2_snia_hesma_source_audit.json`
7. `simulation/snrt/data/fp2_snia_hesma_n100_review_normalized.json`
8. `simulation/snrt/data/fp2_snia_hesma_model_comparison.json`
9. `simulation/snrt/data/fp2_snia_hesma_profile_estimator_comparison.json`
10. `simulation/snrt/data/fp2_snia_hesma_source_selection_packet.json`
11. `simulation/snrt/data/fp2_snia_event_source_admission_audit.json`
12. `simulation/snrt/data/fp2_snia_dtd_contract_audit.json`
13. The directly relevant builders, adapters, and tests in
    `simulation/snrt/tools/` and `simulation/snrt/tests/`.

The runner `bash simulation/snrt/tests/run_fp2_snia_dtd_contract.sh` has just
passed. The HESMA package contains 15 models. The selection packet identifies
13 profile-consistent review candidates and two physical-warning candidates
(`n300c`, `n1600c`), but intentionally selects no model or mixture and keeps
all physical event fields null. Runtime activation and canonical conversion
must remain disabled.

Audit these questions:

- Is the interval-integrated DTD and event ledger algorithm physically and
  numerically justified for the current review-only stage?
- Is the wiring from config → artifacts → source sidecar → contract audit
  complete, hash-bound, deterministic, and fail-closed?
- Does the HESMA inspection distinguish format/integrity evidence from the
  missing physical closure (decay convention/horizon, isotope aggregation,
  returned mass/remnant, energy, momentum, and population weighting)?
- Are the 5% profile-consistency screen and the two warnings scientifically
  defensible as diagnostics rather than approval criteria? Is the profile
  estimator comparison used safely?
- Identify omissions, false assurances, or overclaims. Separate blockers for
  production activation from medium/long-term improvements.

Return a structured report with: verdict (`PASS`, `CONDITIONAL PASS`, or
`BLOCK`), findings ranked by severity, exact file/field references, required
actions before production approval, and a short assessment of whether the
current review-only fail-closed status is appropriate. Do not make changes.
