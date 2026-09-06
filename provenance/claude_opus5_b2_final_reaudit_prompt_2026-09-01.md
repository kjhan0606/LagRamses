# Claude Opus 5 final B2 re-audit request

Date: 2026-09-01
Project root: `/gpfs/kjhan/LRD_JWST`
Mode: read-only; do not edit or execute calculations
Scope: only the two blockers in
`provenance/claude_opus5_b2_reaudit_2026-09-01.md`

Verify whether B2-1 and B2-2 are closed in the current files:

- `simulation/snrt/tools/validate_multiphysics_b2.py`;
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`;
- `simulation/snrt/tests/b2_multiphysics_artifact.py`;
- `simulation/snrt/B2_PRODUCTION_SOLVER_VALIDATION.md`.

Changes:

1. Artifact provenance now includes `snrt_core_sha256`, a sorted mapping for
   every `simulation/snrt/snrt_core/*.py` file (32 modules in the current
   artifact). The artifact test reconstructs and requires exact mapping
   equality, so edits, additions, or removals in the core invalidate it. The
   validator itself remains separately hashed. Documentation narrows the claim
   to validator and `snrt_core` edits.
2. The result table now reports worst-of-five H-ledger L1 `5.60382e-5`, owned
   by `secondary_200ev_on`, rather than the baseline value.

Additional non-blocking improvements from your second audit were also made:
the dust-run H ledger is normalized by gas-absorbed photons, each of the five
run dictionaries is asserted to have at least 20 iterations, the one-iteration
inert shadow is disclosed, and zero-H conservative-primordial coverage was
added.

Return one verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED`. State explicitly
whether B2-1 and B2-2 are closed. Introduce a new B2 blocker only for a material
algorithm, wiring, gate, or provenance defect. Separate later-gate improvements.
This verdict is H-only B2, not overall production/publication readiness. Keep
the response under 1000 words.
