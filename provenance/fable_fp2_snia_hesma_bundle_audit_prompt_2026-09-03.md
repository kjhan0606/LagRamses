# Fable bundle audit prompt — F-P2 SNIa/HESMA source admission

Audit `/gpfs/kjhan/LRD_JWST` read-only with the actual `fable` model. Do not
edit files, select a production model, activate runtime code, commit, or push.
This must be an independent audit; do not consult another model's result.

Review only the F-P2 bundle:

- the interval-integrated SNIa DTD and event ledger in
  `simulation/snrt/native/phase0/` and their exact mirrors in
  `patch/lagRamses/`;
- `simulation/snrt/config/fp2_snia_dtd_contract_v1.json` and
  `simulation/snrt/config/fp2_snia_event_source_approval_sidecar_v1.json`;
- the HESMA source audit, n100 normalized review fixture, all-15 model matrix,
  profile-estimator comparison, selection packet, source-admission audit, and
  DTD contract audit under `simulation/snrt/data/`;
- directly relevant builders/adapters/audits/tests under
  `simulation/snrt/tools/` and `simulation/snrt/tests/`.

The complete runner `bash simulation/snrt/tests/run_fp2_snia_dtd_contract.sh`
has passed. The HESMA review package has 15 models: 13 are profile-consistent
under the existing 5% diagnostic screen and `n300c`/`n1600c` carry physical
warnings. No source model, mixture, estimator, or physical event values are
selected. Canonical conversion and runtime activation are deliberately
disabled pending physical approval.

Assess independently:

1. algorithmic correctness and physical justification of the DTD/event ledger;
2. configuration and artifact wiring, provenance/hash binding, determinism,
   and fail-closed behavior;
3. whether source-format evidence is correctly separated from the unresolved
   physical contract (decay convention/horizon, isotope→project-element
   policy, returned mass/remnant, energy, momentum, and population weights);
4. whether the profile screen, estimator comparison, and warning treatment
   are suitable diagnostics without silently becoming selection criteria;
5. missing blockers, false assurances, and medium/long-term improvements.

Return a structured report with verdict (`PASS`, `CONDITIONAL PASS`, or
`BLOCK`), severity-ranked findings with exact file/field references, required
actions before production approval, and a concise assessment of whether the
review-only state is correctly enforced. Do not modify anything.
