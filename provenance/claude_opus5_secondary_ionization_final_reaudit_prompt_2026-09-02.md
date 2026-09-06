# Final focused re-audit: SNRT stage-2 C4 closure

Work read-only in `/gpfs/kjhan/LRD_JWST`. The full stage-2 audit and first
re-audit are preserved at:

- `provenance/claude_opus5_secondary_ionization_audit_2026-09-02.md`
- `provenance/claude_opus5_secondary_ionization_reaudit_2026-09-02.md`

Your first re-audit returned `CONDITIONAL PASS` solely because C4's RK4
fixture used `dt=1e8 s`, had a 33% error relative to the physical step change,
and was inaccurately described as a resolved branch-selection reference.
Inspect the current changes and determine whether C4 is now honestly closed.
Do not edit files.

The remediation is in:

- `simulation/snrt/tests/coupled_photo_collisional_hhe.py`
- `simulation/snrt/SECONDARY_IONIZATION_VALIDATION.md`

The fixture now deliberately separates two claims:

1. The exact electron-free hot neutral multi-root limit verifies that the
   stationary physical branch is preserved exactly.
2. A separate single-root H-only one-step accuracy fixture uses `T=1e5 K`,
   `dt=3e5 s`, photoionization rate `1e-12 s^-1`, and initial xHII `1e-4`.
   A 4096-substep independent host RK4 reference is compared against the local
   implicit update. The measured absolute xHII error is `2.38311e-10`; it is
   `5.74350e-4` of the actual RK4 physical change. Hard limits are `5e-10`
   absolute and `1e-3` relative to the step change.

The docs no longer call this a generic branch-selection test. They explicitly
state it is a resolved single-root accuracy reference, identify the exact
neutral branch test separately, and leave large-step backward-Euler accuracy
to stage 5.

Please run:

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
JAX_PLATFORMS=cpu .venv/bin/python tests/coupled_photo_collisional_hhe.py
.venv/bin/python tests/secondary_furlanetto_stoever_artifact.py
.venv/bin/python tests/b2_multiphysics_artifact.py
.venv/bin/python tests/p5_secondary_ionization_artifact.py
```

Also verify the wording and arithmetic independently. Briefly revisit whether
any advisory from your first re-audit becomes a stage-2 blocker after these
edits. Return exactly one overall verdict: `PASS`, `CONDITIONAL PASS`, or
`BLOCK`. `PASS` only if C4 and all prior stage-2 conditional findings are
closed; explain any remaining boundary precisely.
