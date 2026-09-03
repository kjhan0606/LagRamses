# Independent re-audit request — F-P1 40–120 M☉ internal controls

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit files, commit, push,
download or alter source data, launch RAMSES, create a physical source package,
or enable production/runtime deposition. Inspect the current checkout and
dirty diff directly. This is a confirmation audit after remediation of the
first AGY/Fable findings; do not trust the remediation claims without tracing
the code and tests.

The claim under review is deliberately narrow: the internal F-P1H-A--E
engineering controls are complete and fail closed, while the 40–120 M☉ physical
gap, F-P1H-F, production readiness, and publication readiness remain blocked.

Read at minimum:

- `provenance/fp1_high_mass_internal_gate_implementation_2026-09-03.md`
- `provenance/fp1_high_mass_internal_gate_dual_audit_comparison_2026-09-03.md`
- `provenance/feedback_population_dtd_active_roadmap.md`
- `simulation/snrt/config/fp1_population_fate_map_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/config/fp1_source_node_contract_v1.json`
- `simulation/snrt/config/fp1_terminal_deposition_contract_v1.json`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/g2_candidate_grid_coverage_contract_v1.json`
- their Python audit tools, mutation tests, generated JSON evidence, converter,
  asset auditor, aggregate runners, Fortran patch/native mirrors, and policy
  tests.

Independently verify the first-audit remediations:

1. High-mass seam conclusions are derived from source records rather than
   trusted literals, and mutations can falsify remnant, isotope-duplicate, and
   closure conclusions.
2. A real source-node record validator enforces all 84 fields, undeclared-field
   rejection, resolver identity, null/zero distinctions, domain checks,
   cumulative-wind constraints, direct-collapse zero-terminal plus baryonic-
   remnant semantics, PISN zero-remnant semantics, and unique node ids.
3. Canonical conversion and asset audit resolve, hash, and inspect actual
   source-node contract/mapping/asset files. Verify coordinate-to-node mapping,
   row counts, contract hashes, asset hashes, and strict energy/momentum
   vocabulary. Look for path substitution or digest-only bypasses.
4. Physical-package qualification is derived from nine candidate-specific,
   checksum-bound evidence and validator pairs. A candidate must not
   self-declare passed gates, and zero qualified packages/zero physical nodes
   must remain impossible to promote.
5. The F-P1 runner is invoked by G2 preflight and meaningful 121 M☉ and
   high-mass-runtime negative tests are included.
6. Fate policy, map SHA256, and approval id travel through the actual namelist
   parser transactionally. Runtime strings cannot self-authorize a build with
   blank or mismatched compile-time identity.
7. Until a real channel-3 source-node fate/deposition consumer exists, the
   driver refuses any channel-3 upper mass above 40 M☉. Check every cumulative
   and interval entry point for bypasses and patch/native mirror agreement.
8. Wording no longer claims source parity passes, an 8–40 M☉ configured window
   is physically resolved, or source-node filtering/deposition exists at
   runtime.

Also independently recompute the W18/N20 numerical evidence and inspect for new
defects, circular evidence, stale generated artifacts, unsafe future-promotion
logic, or test-only behavior. Distinguish a correct fail-closed review system
from an implemented physical feedback model.

You may run these bounded checks and other read-only diagnostics:

- `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh`
- `python3 simulation/snrt/tests/g2_candidate_grid_coverage.py`
- `bash tests/run_stellar_feedback_policy_unit.sh`
- `python3 simulation/snrt/tools/validate_stellar_source_parity.py`
- `git diff --check`

Do not run the full RAMSES simulation.

Return:

- top-level verdict: PASS, CONDITIONAL PASS, or BLOCK;
- separate verdicts for internal engineering controls, physical gap resolution,
  F-P1H-F eligibility, production readiness, and publication readiness;
- for every first-audit finding, `VERIFIED FIXED`, `PARTIAL`, or `OPEN`, with
  exact `file:line` evidence;
- any new findings ranked by severity and independently reproduced where
  feasible;
- exact remaining requirements before physical source admission and before
  runtime activation;
- an explicit answer whether the narrow internal-control claim is confirmed
  without weakening the continuing physical/production/publication BLOCK.
