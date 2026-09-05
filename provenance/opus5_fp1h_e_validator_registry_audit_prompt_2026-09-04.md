# Claude Opus 5 implementation-stage audit prompt — F-P1H-E H1

Audit the H1 implementation in `/gpfs/kjhan/LRD_JWST` as a read-only,
independent reviewer. Do not edit files, run jobs, build RAMSES, or launch a
simulation. Use only repository inspection and, if needed, read-only test
commands.

The project final purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, dust, and coupled source terms. F-P1H-E is the
fail-closed physical stellar-source admission gate. The current scientific
state must remain review-only: zero physical nodes and canonical rows,
unresolved `[0.8,1.0)` and `[40,120] M_sun` seams, and false production,
publication, conversion, and runtime-deposition flags.

The user explicitly authorized implementation despite the unavailable Grok
start-plan audit; no Grok approval is claimed. AGY is retired. This is the
first Claude Opus 5 audit of H1.

## H1 changes to inspect

- `simulation/snrt/tools/fp1_gate_validator_blocks.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/tests/fp1_physical_package_admission.py`
- generated `simulation/snrt/data/fp1_physical_package_admission_audit.json`
  and `simulation/snrt/data/fp1_fate_admission_audit.json`
- parent/plan records, especially
  `provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md`,
  `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`,
  `provenance/feedback_population_dtd_active_roadmap.md`, and
  `provenance/production_publication_readiness_plan.md`.

## Questions

1. Does the registry enumerate exactly the nine contract gates with exact
   requirement sets, correct gate/validator identity, code-hash binding, and
   controlled result validation?
2. Are the eight currently unavailable adapters honestly fail-closed and
   clearly distinguished from real physical validators, or could their
   contract/registry/evidence wiring be mistaken for scientific validation?
   Confirm that they cannot pass, cannot manufacture a package fingerprint,
   and cannot enable selection or runtime deposition.
3. Does the current candidate matrix execute all nine identities and preserve
   the real LC18 source-rights result while reporting the other physical gates
   as blocked? Check that the update does not silently treat declared JSON
   fields as authority.
4. Are imports, path confinement, exception behavior, generated artifact
   determinism, and native/production test coverage sound? Identify any
   absolute-path or circular-evidence issue that is material to this H1.
5. Is H1 a justified next step for the final purpose, and is the scope
   appropriately limited without pretending to resolve missing physical data?

Run the focused H1 test and the full F-P1 population/fate contract runner if
available. Report exact commands/results only if executed. Do not modify any
file. Return a clear verdict of `PASS`, `CONDITIONAL PASS`, or `BLOCK`, with
severity-ranked findings, file/line evidence, and mandatory fixes if any.
