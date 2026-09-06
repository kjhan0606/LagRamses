# Claude Opus 5 implementation-stage audit — F-P1H-E H4

Perform a read-only audit in /gpfs/kjhan/LRD_JWST. Do not edit files, run
jobs, build RAMSES, or launch a simulation. Inspect only these files and the
direct shared-module imports:

- simulation/snrt/tools/fp1_limongi_phase_history.py
- simulation/snrt/tools/audit_g2_limongi_phase_mass_history.py
- simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py
- simulation/snrt/tests/g2_limongi_phase_mass_history.py
- simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py
- simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json
- simulation/snrt/data/g2_limongi_phase_mass_history_audit.json
- simulation/snrt/config/g2_limongi_phase_mass_history_contract_v1.json
- simulation/snrt/config/fp1_physical_package_admission_contract_v1.json
- provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md

This is H4 of the fail-closed F-P1H-E physical stellar-source admission
bundle for the production/publication-ready lagRamses high-level hydro
project: RT, stellar/AGN feedback, dust, and coupled source terms. H4 is
review evidence only. It must not approve or reconstruct physical values.
The current state remains zero physical nodes/canonical rows, false promotion
flags, and unresolved [0.8,1.0) and [40,120] M_sun seams.

Verify:

1. The LC18 failed-wind cross-check exits with a controlled non-zero result
   when cumulative age, total-mass, negative-cumulative-wind, or terminal-PSN
   invariants fail. Its normal unresolved anomaly result must remain blocked
   for production/publication/conversion/deposition and must not infer a
   physical zero from parsed zeros.
2. Both the G2 phase-history audit and the F-P1 cross-check use the same
   fp1_limongi_phase_history.py implementation, and their generated reports
   attest its SHA256. Confirm the common aggregation preserves 108 models,
   845 unique phase rows, 19 exact collapsed duplicate rows, and the
   existing 52/56, 48/4, 53/3, and 101/7 accounting.
3. Phase order, positive duration, total-mass monotonicity, missing PSN,
   source precision, cross-source residuals, failed-model anomaly, and
   review-only publication/rights semantics remain explicit. No failed-wind
   correction, terminal energy inference, momentum inference, or cross-source
   reconciliation is introduced.
4. Test coverage is sufficient for the shared invariants, controlled CLI
   failure, differential count checks, deterministic source/config protection,
   and the existing admission blockers. Do not demand unrelated AMR/HDF5
   work or missing physical source data in H4.

The driver reports these GPFS checks as passed:

- python3 simulation/snrt/tests/g2_limongi_phase_mass_history.py
  -> G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK
- python3 simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py
  -> FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK
- bash simulation/snrt/tests/run_fp1_population_fate_contract.sh
  -> FP1_POPULATION_FATE_CONTRACT_OK

Return severity-ranked findings with file/line evidence and one verdict:
PASS, CONDITIONAL PASS, or BLOCK. Note read-only limitations. No edits.
