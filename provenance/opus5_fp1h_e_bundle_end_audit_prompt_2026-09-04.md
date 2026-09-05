# Claude Opus 5 bundle-end audit — F-P1H-E

Perform a read-only final audit in `/gpfs/kjhan/LRD_JWST`. Do not edit files,
run jobs, build RAMSES, launch a simulation, or modify generated artifacts.
The auditor is Claude Opus 5; Grok and AGY are not active auditors for this
bundle.

## Project and decision boundary

This is the bundle-end audit of F-P1H-E, the executable validator and
physical-package admission closure for the production-ready and
publication-ready lagRamses high-level hydro project: radiative transfer,
stellar/AGN feedback, dust, and coupled source terms.

F-P1H-E is a prerequisite and review/evidence bundle, not a physical-source
approval. The correct repository state is fail-closed: zero physical source
nodes and canonical rows, null package selection, false production,
publication, canonical-conversion, and runtime-deposition flags, and
unresolved `[0.8,1.0)` and `[40,120] M_sun` seams. No missing lifetime, wind,
fate, remnant, decay, energy, momentum, or deposition values may be inferred.

## Scope

Inspect the implementation, contracts, tests, generated evidence, and
provenance records for H1 through H5, especially:

- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/tools/fp1_gate_validator_blocks.py`
- `simulation/snrt/tools/audit_fp1_source_node_contract.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`
- `simulation/snrt/tools/fp1_limongi_phase_history.py`
- `simulation/snrt/tools/audit_g2_limongi_phase_mass_history.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tools/audit_stellar_yield_asset.py`
- `simulation/snrt/tools/convert_yield_rows_to_canonical.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_source_node_contract_v1.json`
- `simulation/snrt/data/fp1_physical_package_admission_audit.json`
- `simulation/snrt/data/fp1_source_node_contract_audit.json`
- `simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`
- `simulation/snrt/data/g2_limongi_phase_mass_history_audit.json`
- `simulation/snrt/tests/fp1_physical_package_admission.py`
- `simulation/snrt/tests/fp1_source_node_contract.py`
- `simulation/snrt/tests/g2_limongi_phase_mass_history.py`
- `simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tests/stellar_yield_asset.py`
- `simulation/snrt/tests/yield_converter.py`
- `provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md`

## Required audit questions

1. H1: Do the contract and code registry name exactly the same nine gates in
   both directions, with exact requirement sets, code-file hashes, candidate
   identity, and strict report validation? Do the eight unavailable adapters
   remain explicit never-passing blockers rather than synthetic physical
   validators?
2. H2: Can an approved node pass by self-declaring rights, using an unknown
   source, a blocked/substituted validator, or a mismatched package
   fingerprint? Verify that the code-owned source-rights validator is actually
   executed and that review/empty states cannot promote.
3. H3: Does each candidate name all nine executable gate reports and statuses,
   and does qualification/selection require all nine pass reports, complete
   physical nodes, exact fingerprints, and upstream gates? Verify the current
   four-candidate state remains blocked with no selected package.
4. H4: Is the shared phase-history aggregator fail-closed on source order,
   duplicate handling, missing PSN, positive duration/cumulative age,
   nonnegative cumulative wind, and nonincreasing mass? Are 108 models, 845
   unique phase rows, 19 collapsed extras, and 52/56, 48/4, 53/3, 101/7
   accounting preserved? Are parsed zeros, failed-wind anomalies,
   cross-source residuals, source precision, and no-inference semantics
   explicit?
5. H5: Are direct registry adversarial tests present for unregistered,
   gate/validator mis-bound, malformed-report, and stale-code-hash cases? Do
   G2 and LC18 stored artifacts compare against live reports, and does the G2
   report carry the caveat that intermediate-burning order is a project
   contract assumption rather than source-attested data? Are synthetic test
   seams confined to tests and restored without changing production behavior?
6. Does the whole bundle preserve the high-level hydro scope without pulling
   in unrelated AMR/HDF5/restart work or launching a RAMSES calculation?

## Driver evidence

The driver ran on GPFS and obtained:

- `python3 simulation/snrt/tests/fp1_physical_package_admission.py` →
  `FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK`
- `python3 simulation/snrt/tests/fp1_source_node_contract.py` →
  `FP1_SOURCE_NODE_CONTRACT_TEST_OK`
- `python3 simulation/snrt/tests/g2_limongi_phase_mass_history.py` →
  `G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK`
- `python3 simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py` →
  `FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK`
- `.venv/bin/python simulation/snrt/tests/stellar_yield_asset.py` →
  `STELLAR_YIELD_ASSET_TEST_OK`
- `.venv/bin/python simulation/snrt/tests/yield_converter.py` →
  `YIELD_CONVERTER_TEST_OK`
- `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh` →
  `FP1_POPULATION_FATE_CONTRACT_OK`
- `bash simulation/snrt/tests/run_g2_preflight.sh` → terminal
  `G2_PREFLIGHT_BLOCKED`
- Python `compileall` and `git diff --check` → pass

The expected blocked state is a success of the admission gate, not a failed
test: no authoritative physical package is currently qualified.

Return severity-ranked findings with file/line evidence and exactly one
verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`. A low-severity improvement
must not be promoted to a condition unless it can permit unsafe physical
promotion or invalidate the bundle's stated evidence. Note read-only
limitations. Do not edit files or run jobs.
