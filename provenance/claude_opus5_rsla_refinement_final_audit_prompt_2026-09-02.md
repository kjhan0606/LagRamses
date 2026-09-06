# Final stage-3 audit: RSLA and refinement

Audit the current dirty worktree in `/gpfs/kjhan/LRD_JWST` read-only with
`claude-opus-5`. Do not edit files. This is the final mandatory audit after the
first re-audit in
`provenance/claude_opus5_rsla_refinement_reaudit_2026-09-02.md` returned
`CONDITIONAL PASS` with new findings N1–N4.

Return exactly one verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`. PASS
requires the original F1–F4 and new N1–N4 to be closed, the canonical v3
artifact to be source-bound and internally correct, and no new material
stage-3 defect. Later dynamic-source, radiation-hydrodynamic, helium, and dust
gates remain out of scope.

## Remediation to verify

- N1: `simulation/snrt/tools/validate_rsla_refinement.py` now extrapolates the
  three high-speed points in two coordinates: inverse reduced-light fraction
  and directly measured photon-storage fraction. In each coordinate it uses
  the two adjacent linear fits and a three-point quadratic fit, pads the maximum
  by the full fit-order spread, and hard-gates the larger coordinate-specific
  upper bound. The storage coordinate is selected: upper radius ratio
  `1.0093227854966493`, production RSLA term `0.017269837730221192`, total
  envelope `0.019764600042460632` (98.8% of the 2% gate). Verify every value
  independently from `data/rsla_refinement_validation.json`, and judge the
  scientific interpretation and narrow margin.
- N2: `RSLA_REFINEMENT_VALIDATION.md` now says old P1 used `4*size` at `0.01c`
  and already had fixed physical duration; the fix raises ĉ to `0.03c` with
  `12*size` while adding an analytic assertion.
- N3: the escape roundoff tolerance `1e-4` of emitted photons is declared in
  the report, validator constant, JSON acceptance thresholds, and test.
- N4: `tests/rsla_refinement_artifact.py` pins the mesh allowance to exactly
  `0.005`, pins escape tolerance to `1e-4`, independently recomputes both
  linear and the Lagrange quadratic intercepts in both coordinates from the
  matrix, checks coordinate selection, and recomputes the production RSLA
  error.

Also recheck original F1–F4 and all hash/document provenance. Primary files:

- `simulation/snrt/tools/validate_rsla_refinement.py`
- `simulation/snrt/tests/rsla_refinement_artifact.py`
- `simulation/snrt/data/rsla_refinement_validation.json`
- `simulation/snrt/RSLA_REFINEMENT_VALIDATION.md`
- `simulation/snrt/tests/p1_validation.py`
- `simulation/snrt/snrt_core/ionization_front.py`
- `simulation/snrt/P1_CONVERGENCE.md`
- B2 and P5 canonical reports/artifacts

Run from `simulation/snrt`:

```bash
.venv/bin/python tests/rsla_refinement_artifact.py
JAX_PLATFORMS=cpu .venv/bin/python tests/p1_validation.py
.venv/bin/python tests/b2_multiphysics_artifact.py
.venv/bin/python tests/p5_secondary_ionization_artifact.py
```

Inspect the diff and do not trust the prompt. Give a compact F1–F4/N1–N4
closure table, list any new findings with severity, and finish with the single
overall verdict.
