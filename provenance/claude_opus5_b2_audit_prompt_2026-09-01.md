# Claude Opus 5 audit request: B2 production-solver validation only

Date: 2026-09-01
Project root: `/gpfs/kjhan/LRD_JWST`
Mode: independent read-only audit; do not edit files or run new expensive
calculations.

Audit only B2, the production multiphysics RT transport/chemistry validation.
Do not expand this into a whole-project, yield-table, thermal-atlas, or live
RAMSES audit. Overall production/publication readiness is known to remain
blocked by later gates and unapproved assets.

The original findings and acceptance contract are in:

- `provenance/claude_opus5_rt_audit_2026-09-01.md`, especially B2 and M6;
- `simulation/snrt/B2_PRODUCTION_SOLVER_VALIDATION.md`.

Inspect at minimum:

- `simulation/snrt/snrt_core/multiphysics.py`;
- `simulation/snrt/snrt_core/implicit.py`;
- `simulation/snrt/snrt_core/conservative_hydrogen.py`;
- `simulation/snrt/snrt_core/conservative_primordial.py`;
- `simulation/snrt/snrt_core/thermochemistry.py`;
- `simulation/snrt/tools/validate_multiphysics_b2.py`;
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`;
- `simulation/snrt/tests/b2_multiphysics_artifact.py`;
- `simulation/snrt/tests/p2_p3_validation.py`;
- `simulation/snrt/tools/p4_run_transport_pilot.py`;
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`.

Recorded local results include Solver-A analytic radius ratio `0.9559179`,
Solver-A/B xHII L1 `1.15449e-6`, cap-active cell-step fraction `0`, maximum
fixed-point residual `2.39611e-5`, H-ledger L1 relative error `1.94513e-6`,
dust mean-xHII delta `-0.00400357`, secondary delta `+0.0122545`, and Solver-A
S8/A192 shadow relative difference `0.00439426`. The source-bound artifact
test passes under JAX 0.11.1 CPU.

Evaluate, with file/line evidence:

1. whether the former atom-inventory attenuation cap and dead/direct branch
   are genuinely removed from the production Solver-A path;
2. whether the absorbed-count-to-rate construction, time-averaged H opacity,
   He backward-Euler opacity state, recombination/collisional terms, dust
   partition, secondary terms, and species ledgers are mathematically and
   physically consistent at the converged fixed point;
3. whether 20 fixed iterations with 0.5 under-relaxation and the recorded
   residual are adequately fail-closed by the validator and P4/P5 runners;
4. whether the Solver-A/Solver-B differential is independent and configured
   fairly enough to satisfy B2 rather than comparing two aliases;
5. whether the dust, secondary, and full-Solver-A shadow controls exercise the
   claimed wiring and whether their acceptance checks could pass vacuously;
6. whether the cap-activity, fixed-point, non-finite, positivity, and ledger
   diagnostics are correctly propagated and normalized;
7. whether the zero-helium underflow repair is safe and introduces no new
   divide-by-zero path;
8. whether any material bug, missing acceptance check, stale claim, or test
   weakness prevents B2 closure.

Give exactly one B2 verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED`.
Distinguish mandatory B2 fixes from later-gate improvements. PASS means no
material B2 algorithm/wiring/validation defect remains. Keep the report under
2500 words. Do not treat B2 PASS as overall production or publication PASS.
