# Focused stage-3 re-audit: RSLA and refinement

Audit the current dirty worktree in `/gpfs/kjhan/LRD_JWST` read-only. Do not
edit any file. Use the `claude-opus-5` model's own scientific and code judgment;
do not merely accept this prompt's claims.

This is the mandatory re-audit after the stage-3 audit in
`provenance/claude_opus5_rsla_refinement_audit_2026-09-02.md` returned
`CONDITIONAL PASS`. Return exactly one overall verdict: `PASS`,
`CONDITIONAL PASS`, or `BLOCK`. A PASS requires all four findings below to be
closed without introducing a new stage-3 defect. Later source-variability,
radiation-hydrodynamics, helium, and dust limitations are out of this gate and
must not be confused with defects in this fixed H-only benchmark.

## Findings to verify

1. **F1, false helium attribution:** At the 20 eV P1 fixture He I has zero
   opacity and remains neutral. Verify that
   `simulation/snrt/snrt_core/ionization_front.py` and
   `simulation/snrt/P1_CONVERGENCE.md` now say this correctly and no longer
   attribute the residual radius deficit to helium.
2. **F2, finite-0.03c error reference:** Verify the new implementation in
   `simulation/snrt/tools/validate_rsla_refinement.py`, its artifact, test, and
   report. It fits the three high-speed points in inverse reduced-light
   fraction using two adjacent-pair linear intercepts and one quadratic
   intercept. The three infinite-light radius-ratio estimates are
   `1.0073980`, `1.0070156`, and `1.0069731`. It takes the largest estimate plus
   the complete fit spread as a one-sided upper bound, `1.0078229`. The
   production-to-bound term is `0.0158073`; adding direct same-0.01c mesh
   `0.0021720`, angular `0.0002655`, and Solver A/B `0.0000572` sensitivities
   gives `0.0183020`, hard-gated below 2%. Determine whether this is a
   scientifically defensible conservative, benchmark-specific bound. Check
   code and JSON calculations independently, including the exclusion of the
   visibly non-asymptotic 0.001c point.
3. **F3, hidden mesh-degradation allowance:** Verify the 0.005 allowance is now
   present in the documented acceptance contract, JSON thresholds, validator,
   and artifact assertion, with the correct limited interpretation (guard
   against material degradation, not convergence-order evidence).
4. **F4, stale B2 provenance:** Verify the B2 report now records artifact HEAD
   `ca90a391296e4fbd99d183df3850de10c537cef4` and that B2/P5/RSLA canonical
   artifacts have been regenerated or rebound to the current full-core hash.

## Primary evidence

- `simulation/snrt/tools/validate_rsla_refinement.py`
- `simulation/snrt/tests/rsla_refinement_artifact.py`
- `simulation/snrt/data/rsla_refinement_validation.json`
- `simulation/snrt/RSLA_REFINEMENT_VALIDATION.md`
- `simulation/snrt/tests/p1_validation.py`
- `simulation/snrt/snrt_core/ionization_front.py`
- `simulation/snrt/P1_CONVERGENCE.md`
- `simulation/snrt/B2_PRODUCTION_SOLVER_VALIDATION.md`
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`
- `simulation/snrt/data/p5_secondary_ionization_validation.json`

Run from `simulation/snrt`:

```bash
.venv/bin/python tests/rsla_refinement_artifact.py
JAX_PLATFORMS=cpu .venv/bin/python tests/p1_validation.py
.venv/bin/python tests/b2_multiphysics_artifact.py
.venv/bin/python tests/p5_secondary_ionization_artifact.py
```

Independently recompute the extrapolation and envelope from the canonical JSON,
check every new criterion and hash binding, and inspect the source diff. Report
a concise F1-F4 closure table, any new stage-3 findings with severity, and one
final verdict. Do not grant PASS for numerical correctness alone if a physical
claim, gate implementation, or provenance statement remains materially false.
