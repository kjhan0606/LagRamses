# Claude Opus 5 H1 audit retry — F-P1H-E validator registry

Perform a read-only implementation audit in `/gpfs/kjhan/LRD_JWST`. Do not
edit files, run jobs, build RAMSES, or launch a simulation. Inspect only the
following exact files; do not recursively search the repository:

- `simulation/snrt/tools/fp1_gate_validator_blocks.py`
- `simulation/snrt/tools/fp1_gate_validator_registry.py`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`
- `simulation/snrt/tests/fp1_physical_package_admission.py`
- `simulation/snrt/data/fp1_physical_package_admission_audit.json`
- `simulation/snrt/data/fp1_fate_admission_audit.json`
- `provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md`

H1 is an explicitly fail-closed registry/admission wiring step for the
production/publication-ready lagRamses high-level hydro project (RT,
stellar/AGN feedback, dust, and coupled source terms). It must not approve a
physical source package: the current state must remain zero physical nodes,
zero canonical rows, unresolved `[0.8,1.0)` and `[40,120] M_sun` seams, and
false production/publication/conversion/runtime-deposition flags.

Assess these points:

1. The registry and contract enumerate exactly the same nine gate and
   validator identities, exact requirement sets, and code-hash binding.
2. The eight unavailable adapters are visibly controlled blockers: they can
   never pass, cannot manufacture a package fingerprint, cannot select a
   package, and cannot enable runtime deposition. They must not be mistaken
   for physical validators.
3. The candidate matrix contains all nine identities, preserves the real LC18
   source-rights result, and does not promote candidates from declarative JSON
   status strings.
4. Review imports, exception/result semantics, path behavior, generated
   artifact consistency, and the focused test's assertions. Flag only
   material H1 defects; do not demand unrelated AMR/HDF5/runtime hardening.
5. Decide whether H1 is a justified, appropriately bounded prerequisite for
   the project's final purpose.

Give a concise severity-ranked report with file/line evidence and one clear
verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`. If a command cannot be run
under the read-only tool restriction, say so rather than inferring its result.
No edits.
