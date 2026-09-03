# Independent audit request — F-P1 40–120 M☉ internal implementation

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit files, commit, push,
download or alter source data, launch RAMSES, or enable production/runtime
deposition. Inspect the actual checkout and current dirty diff rather than
trusting this prompt or the provenance summary. Use generous reasoning time.

This is a post-implementation audit of the internal F-P1H-A--E controls. It is
not a claim that the 40–120 M☉ physical fate gap has been solved. A correct
audit may PASS the engineering controls while BLOCKING physical promotion.

Review at minimum:

- `provenance/fp1_high_mass_internal_gate_implementation_2026-09-03.md`
- `provenance/fp1_high_mass_required_data_comparison_2026-09-03.md`
- `provenance/feedback_population_dtd_active_roadmap.md`
- `simulation/snrt/config/fp1_population_fate_map_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/config/fp1_source_node_contract_v1.json`
- `simulation/snrt/config/fp1_terminal_deposition_contract_v1.json`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/g2_candidate_grid_coverage_contract_v1.json`
- their audit tools, negative tests, generated JSON evidence, Fortran config
  mirrors, converter/asset-audit wiring, and shell runners.

Independently verify:

1. F-P1H-A: a namelist-supplied digest/approval id cannot self-authorize a
   production build; the compiled identity is exact and review builds remain
   blank/fail closed.
2. F-P1H-B: the source-node sidecar losslessly retains every resolver axis and
   distinguishes missing from physical zero. Failed/direct-collapse nodes must
   be explicit wind/zero-terminal/remnant records. Check whether the canonical
   conversion and asset audit really preserve/bind the new identity fields.
3. F-P1H-C: channel 3 is an 8–120 M☉ candidate domain filtered by source-node
   fate, not a universal explosion interval. Check mutually exclusive ownership,
   energy-kind semantics, scalar radial momentum storage, deposition geometry,
   and exactly-once rules. Identify any missing runtime consumer or bypass.
4. F-P1H-D: branch-specific nodes are retained for F23 single/binary, LC18,
   WH07, Z9.6, W18/N20, Limongi set-R, and Heger-Woosley Pop III. Verify that a
   flattened union cannot be treated as continuous/interpolable coverage and
   that 40–120 M☉ nodes are not silently clipped. Recompute or inspect the
   rounded-source mass residual, failed-node counts, W18/N20 wind comparison,
   radioactive reference-epoch warnings, and K-40 duplicate handling.
5. F-P1H-E: the nine-gate physical-package contract is complete enough for a
   production/publication admission decision, all evidence is checksum-bound,
   no candidate can self-declare qualification, and zero physical nodes cannot
   lead to canonical conversion or runtime deposition.
6. Test wiring: confirm the tests are actually invoked, mutation tests exercise
   meaningful negative cases, generated evidence is reproducible, and the
   source/native/production mirrors agree. Treat the reported source-parity
   status `blocked=production_linked_build_evidence` as an open item and assess
   whether it blocks only promotion or invalidates an internal-control claim.
7. Scientific validity: identify every required datum still absent before the
   40–120 M☉ gap is physically resolved. In particular assess source rights,
   multi-Z/rotation/population coverage, failed-collapse remnant/wind records,
   lifetimes/age-resolved release, decay epoch/network, terminal mass/species
   closure, injected-energy mapping, radial momentum/deposition, PPISN/PISN,
   LC18 failed-model Wind anomaly, and source-hull boundary behavior.
8. Look for overclaims, schema-only checks that do not prove runtime behavior,
   circular checksum/evidence dependencies, stale generated artifacts, unsafe
   future promotion logic, or current implementation defects.

You may run these bounded checks and other read-only diagnostics:

- `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh`
- `python3 simulation/snrt/tests/g2_candidate_grid_coverage.py`
- `bash tests/run_stellar_feedback_policy_unit.sh`
- `python3 simulation/snrt/tools/validate_stellar_source_parity.py`
- `git diff --check`

Do not run the full RAMSES simulation.

Return:

- top-level verdict: PASS, CONDITIONAL PASS, or BLOCK;
- separate verdicts for engineering controls, physical gap resolution,
  production readiness, and publication readiness;
- a severity-ranked findings table with implemented/verified/partial/missing
  status and exact `file:line` evidence;
- independently recomputed key numerical evidence where feasible;
- precise mandatory fixes before F-P1H-F and before runtime activation;
- explicit assessment of whether the local wording is honest about the gap
  remaining physically unresolved.
