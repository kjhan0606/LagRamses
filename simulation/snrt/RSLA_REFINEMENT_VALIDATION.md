# RSLA and refinement validation

Date: 2026-09-02
Stage status: **PASS — Claude Opus 5 closure audit**

## Scope and acceptance contract

This stage quantifies the reduced-speed-of-light approximation for the B2
H-only Strömgren fixture at fixed physical duration `0.5 t_rec`. The source,
density, temperature, domain, one 18 eV group, and S4 transport match B2. The
independent conservative H solver supplies the full RSLA and discretization
matrix. The production multiphysics Solver A is run separately at `0.01c` and
compared cell by cell to the conservative result.

The limits were fixed before the accepted run:

- production `0.01c` analytic-radius error below 2%;
- production and conservative `0.01c` results each within 2% of the `0.03c`
  conservative reference;
- Solver-A/Solver-B mean absolute xHII below `5e-5`;
- 32-to-64 mesh radius-ratio change below 0.03;
- S4-to-S8 radius-ratio change below 0.02;
- a linear-sum production radius-error envelope below 2%, using the larger
  padded upper bound from the predeclared `1/(ĉ/c)` and photon-storage
  extrapolations;
- at `0.003c`, mesh refinement may worsen the absolute analytic-radius error by
  no more than 0.005 (0.5 percentage points); this is a guard against material
  degradation, not a convergence-order claim;
- every run independently passes `1e-4` fixed-point, `1e-3` H-ledger, finite,
  and production electron-root gates; inferred escape may be negative only by
  the declared accumulation-roundoff allowance `1e-4` of emitted photons.

These are validation limits for this fixed H-only front, not a universal RSLA
guarantee for arbitrary source variability, hydrodynamic fronts, helium,
dust, or live feedback.

The extrapolation-coordinate set for this v3 gate is closed and consists
exactly of `1/(ĉ/c)` and photon-storage fraction: the expected leading RSLA
scaling and its directly measured physical mechanism. Adding or substituting a
coordinate is a new gate definition requiring a schema/version change and a
fresh predeclared run; coordinates are not added post hoc to move this narrow
98.8%-utilized decision boundary.

## RSLA matrix

All rows use 32 cubed cells, S4, 32 opacity iterations, and exactly the same
physical duration.

| ĉ/c | Steps | Rion/Ranalytic | Analytic error | Photons still in domain / emitted | Ionized-volume deficit |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `0.001` | 39 | `0.845542` | `15.45%` | `43.63%` | `39.55%` |
| `0.003` | 116 | `0.955902` | `4.41%` | `16.19%` | `12.65%` |
| `0.010` | 385 | `0.991949` | `0.805%` | `5.02%` | `2.40%` |
| `0.030` | 1153 | `1.001993` | `0.199%` | `1.68%` | `-0.60%` |

The radius converges monotonically toward the infinite-light-speed analytic
solution as `ĉ` increases. The old 4.41% B2 discrepancy is quantitatively
consistent with finite photon storage: at `0.003c`, 16.19% of emitted photons
remain in flight and the ionized-volume deficit is 12.65%. The slight
`0.03c` overshoot is below the 2% reference criterion and reflects the
remaining discretization/chemistry difference rather than an asserted exact
finite-c solution.

## Production `0.01c` result

Solver A gives `Rion/Ranalytic=0.991892`, an absolute analytic error of 0.811%.
Against the conservative `0.01c` field, mean absolute xHII is `4.3528e-6` and
the radius-ratio difference is `5.7189e-5`. Relative to the conservative
`0.03c` run, the production radius differs by 1.008%. That finite-ĉ comparison
is retained as a diagnostic, not used as the infinite-light-speed error
estimate.

The infinite-light limit is estimated from the three highest-speed runs in two
coordinates: `1/(ĉ/c)`, the expected leading RSLA scaling, and measured photon-
storage fraction, the diagnosed physical cause of the lag. In each coordinate,
linear fits to the `0.003–0.01c` and `0.01–0.03c` pairs and a quadratic fit
through all three are compared. The largest intercept plus the complete
fit-order spread gives a one-sided bound for that coordinate; the hard gate
uses the larger coordinate-specific bound.

| Coordinate | Linear pair intercepts | Quadratic intercept | Padded upper bound |
| --- | --- | ---: | ---: |
| `1/(ĉ/c)` | `1.0073980`, `1.0070156` | `1.0069731` | `1.0078229` |
| Photon-storage fraction | `1.0081344`, `1.0070696` | `1.0069459` | `1.0093228` |

This makes the earlier audit's `1.00813` storage-coordinate estimate explicit
and prevents fit-coordinate choice from understating the gate. The `0.001c`
point is excluded because it is visibly outside the high-speed asymptotic
regime, with 43.63% of emitted photons still in flight. Including it in higher-
order fits gives compatible limits but not a reliable adjacent linear slope.

The adopted error envelope is a linear sum, not a quadrature uncertainty:

| Contribution | Radius fraction |
| --- | ---: |
| Production `0.01c` vs cross-coordinate infinite-light upper bound | `0.0172698` |
| `0.01c` 32-to-64 mesh sensitivity | `0.0021720` |
| `0.01c` S4-to-S8 angular sensitivity | `0.0002655` |
| Solver A/B radius difference | `0.0000572` |
| **Linear envelope** | **`0.0197646` (1.98%)** |

Thus `0.01c` passes the declared 2% radius-error gate for this benchmark, but
uses 98.8% of the allowance. The margin is narrow and is not a basis for
generalizing `0.01c` to a different source history or gas configuration.

## Mesh and angular refinement

The original B2 `0.003c` sensitivity is retained because it directly explains
the historical 4.41% radius deficit.

| Grid / quadrature | Rion/Ranalytic | Change from 32³ S4 |
| --- | ---: | ---: |
| 32³ S4 | `0.955902` | reference |
| 64³ S4 | `0.951651` | `-0.004251` |
| 32³ S8 | `0.954583` | `-0.001319` |

The factor-two mesh change is 0.425 percentage points and the angular change
is 0.132 percentage points. The 64³ value does not move monotonically toward
the analytic radius; therefore this is a measured sensitivity bound, not a
Richardson-extrapolated convergence order. The predeclared degradation guard
allows at most 0.5 percentage points of increased analytic error; the measured
0.425 percentage points passes but uses 85% of that allowance. At `0.003c`,
RSLA photon storage dominates both discretization changes.

The production error envelope uses a second, same-`ĉ` refinement family at
`0.01c`, avoiding transfer of a discretization estimate between RSLA regimes.

| Grid / quadrature at `0.01c` | Rion/Ranalytic | Change from 32³ S4 |
| --- | ---: | ---: |
| 32³ S4 | `0.991949` | reference |
| 64³ S4 | `0.989777` | `-0.002172` |
| 32³ S8 | `0.992215` | `+0.000266` |

Both refined values remain within 2% of the analytic radius. Their measured
sensitivities, 0.217 and 0.027 percentage points, are the mesh/angular terms in
the 1.98% production envelope above.

## P1 correction

The older P1 benchmark already used a grid-independent physical duration via
`4*size` steps at `0.01c`, but asserted only grid-to-grid consistency while the
front was about 20% below the analytic radius. The corrected test keeps that
physical duration, raises the light speed to `0.03c` with `12*size` steps, and
explicitly asserts every grid is within 10% of the time-dependent analytic
radius. The corrected 32³/48³/64³ analytic ratios are `0.9351`, `0.9296`, and
`0.9261`; any reduced-light choice still requires a measured reference.

## Reproduction

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
JAX_PLATFORMS=cpu .venv/bin/python tools/validate_rsla_refinement.py \
  --output data/rsla_refinement_validation.json
.venv/bin/python tests/rsla_refinement_artifact.py
JAX_PLATFORMS=cpu .venv/bin/python tests/p1_validation.py
```

The canonical JSON binds the validator, every current `snrt_core/*.py`, the B2
artifact, JAX/backend/device identity, and repository HEAD. Large-step source
variability, radiation hydrodynamics, and science-observable convergence remain
later gates.

- Canonical JSON: `data/rsla_refinement_validation.json`, SHA256
  `899ee41325e2c2ec935db4ccf58d38c777aa7e23aecab8121e2cf34ab5fe1431`.
- Validator SHA256:
  `e3972142762216fbf8968ec6a6e84091873ab67ddbef941502c0e112d46ab92c`.

## Independent audit

Claude Opus 5 performed one full audit, two full re-audits, and a focused
closure audit. The successive verdicts were `CONDITIONAL PASS`, `CONDITIONAL
PASS`, `CONDITIONAL PASS`, and final **PASS**. The audit trail is preserved in:

- `provenance/claude_opus5_rsla_refinement_audit_2026-09-02.md`;
- `provenance/claude_opus5_rsla_refinement_reaudit_2026-09-02.md`;
- `provenance/claude_opus5_rsla_refinement_final_audit_2026-09-02.md`;
- `provenance/claude_opus5_rsla_refinement_closure_audit_2026-09-02.md`.

The final audit independently reproduced the v3 extrapolation and envelope,
verified all original and follow-up findings closed, checked full hash closure,
and confirmed that the two-coordinate definition is both documented and
enforced by the artifact test.
