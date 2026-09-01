# Claude Opus 5 helium-recombination re-audit request

Date: 2026-09-02
Project root: `/gpfs/kjhan/LRD_JWST`
Mode: read-only; do not edit files, launch simulations, or regenerate artifacts
Scope: closure of the four conditions in
`provenance/claude_opus5_he_recombination_audit_2026-09-01.md`

Re-audit roadmap stage 1 in the current dirty working tree. Verify these fixes:

1. `simulation/snrt/snrt_core/primordial_cooling.py` now implements
   `beta_HeIII,B(T) = 8 beta_HII,B(T/4)` by calling one shared H II cooling
   function, so the explicit temperature prefactor is evaluated at `T/4` and
   the former 4x cooling error and duplicate formula are both removed.
2. `simulation/snrt/tests/b1_thermal_coupling.py` checks the first-principles
   `beta/alpha < 1.5 k_B T` bound for the radiative H II, He II, and He III
   channels at `1e4, 2e4, 4e4, 1e5 K`. It separately checks the He II
   dielectronic cooling/rate coefficient ratio because the radiative bound does
   not apply to the autoionizing-level energy.
3. `simulation/snrt/tests/helium_recombination.py` no longer compares the
   production He III function to itself or to its own direct dependency. It
   uses a separately written NumPy evaluation of the Hui--Gnedin formulas and
   literature-locked reference arrays whose source is named in the artifact.
4. The original three-recombination-time run is explicitly labeled a
   cross-module consistency gate. A second run uses one literal physical
   interval, `2.0e12 s`, for all temperatures, produces non-degenerate final
   fractions, and remains below the predeclared 2% error limit.

Also inspect:

- `simulation/snrt/data/helium_case_b_recombination_validation.json`
- `simulation/snrt/tests/helium_recombination_artifact.py`
- `simulation/snrt/HELIUM_RECOMBINATION_VALIDATION.md`
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`
- `simulation/snrt/tests/b2_multiphysics_artifact.py`
- `simulation/snrt/B2_PRODUCTION_SOLVER_VALIDATION.md`

Recorded local validation after the final refactor:

- helium one-zone test: PASS; fixed-time maximum relative errors
  `0.00168485` (He II), `0.00932733` (He III);
- B1 thermal coupling: PASS;
- B2 zero-density: PASS;
- P2/P3 implicit/sharding: PASS on two CPU devices;
- full canonical B2 rerun: PASS with unchanged physics payload;
- both helium and B2 fail-closed artifact tests: PASS;
- compileall and `git diff --check`: PASS.

Recompute or independently inspect whatever the read-only tool set permits.
Cite exact `file:line` locations for any remaining blocker. Do not block stage
1 on the coupled H+He front/convergence work intentionally assigned to stage
6. Separate future improvements from stage-1 defects.

Return concise Markdown with the exact model identifier and exactly one final
verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED`. State explicitly whether
all four original conditions are closed. Introduce a new blocker only for a
material algorithm, wiring, test, or provenance defect. This is a stage-1
verdict only, not full RT production/publication readiness.
