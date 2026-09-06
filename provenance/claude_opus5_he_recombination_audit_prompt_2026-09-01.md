# Claude Opus 5 helium-recombination audit request

Date: 2026-09-01
Project root: `/gpfs/kjhan/LRD_JWST`
Mode: read-only; do not edit files, launch simulations, or regenerate artifacts
Scope: roadmap stage 1, helium atomic physics

Independently audit the current dirty working tree for the completion of stage
1 of the RT production/publication-readiness roadmap:

1. He III recombination must use the hydrogenic case-B relation
   `alpha_HeIII(T) = 2 * alpha_HII,B(T / 4)`.
2. He II radiative recombination must consistently use a Hui & Gnedin-style
   case-B coefficient. Its dielectronic contribution must be physically and
   algebraically well-defined, without omission or double counting.
3. Temperature-resolved one-zone reference tests must cover at least
   `1e4, 2e4, 4e4, 1e5 K` and meaningfully test the implementation rather than
   merely reproducing the same function on both sides of an assertion.
4. Every chemistry/thermal/transport caller must be wired to the corrected
   coefficients with no stale mixed case-A/case-B formula left in executable
   source.
5. Recorded machine artifacts, provenance hashes, documentation, and existing
   B2 regression evidence must support the stage claim.

Inspect at minimum:

- `simulation/snrt/snrt_core/primordial.py`
- `simulation/snrt/snrt_core/implicit.py`
- `simulation/snrt/snrt_core/multiphysics.py`
- `simulation/snrt/snrt_core/photon_coupling.py`
- `simulation/snrt/snrt_core/conservative_primordial.py`
- `simulation/snrt/tests/helium_recombination.py`
- `simulation/snrt/data/helium_case_b_recombination_validation.json`
- `simulation/snrt/HELIUM_RECOMBINATION_VALIDATION.md`
- `simulation/snrt/P3_IMPLICIT_TPU.md`
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`

Audit requirements:

- Check the equations and coefficient conventions against the cited Hui &
  Gnedin (1997) source, including units and the meaning of the He II
  dielectronic term.
- Enumerate all relevant call sites and identify any stale formula or semantic
  mismatch.
- Assess whether the hard-coded reference values in the one-zone test are
  sufficiently independent, and whether the 3-recombination-time, 512-step,
  less-than-2-percent criterion is a valid numerical check.
- Check that the artifact hashes correspond to the current files and that the
  B2 rerun did not hide a physics-payload regression.
- Cite exact `file:line` locations for every blocking or conditional finding.
- Do not block stage 1 merely because the coupled H+He opacity/front and
  timestep-convergence gate is intentionally assigned to roadmap stage 6,
  unless that deferred work makes this local implementation invalid.
- Separate future improvements from defects that must be fixed before stage 1
  closes.

Return concise Markdown with the exact model identifier and exactly one final
verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED`. A conditional or blocked
verdict must list the minimum concrete fixes required for re-audit. This is a
stage-1 verdict, not a claim that the full RT stack is production/publication
ready.
