# Stage-3 re-audit — RSLA and refinement

Model: `claude-opus-5`
Mode: read-only
Repository HEAD during audit: `ca90a391296e4fbd99d183df3850de10c537cef4`

All four requested commands passed:

```text
RSLA_REFINEMENT_ARTIFACT_OK ratios=0.845542,0.955902,0.991949,1.00199 production=0.991892 error_envelope=0.018302 mesh=0.00425058 angular=0.00131851
P1_VALIDATION_OK spatial_spread=0.0096 analytic_ratios=0.9351,0.9296,0.9261
B2_MULTIPHYSICS_ARTIFACT_OK radius=0.955685 A_B_L1=1.52855e-05
P5_SECONDARY_IONIZATION_ARTIFACT_OK delta_mean_xhii=2.24653e-08
```

The auditor independently reproduced the inverse-reduced-light-fraction
intercepts `1.0073980090245733`, `1.0070156249072968`, and
`1.0069731377831546`, their spread `4.248712414187672e-4`, upper bound
`1.007822880265992`, RSLA term `0.01580727705654007`, and total envelope
`0.01830203936877951` exactly. It also confirmed the independent B2 Solver-B
radius agrees with the stage-3 `0.003c` fixture to `2.3e-7`.

## F1–F4 closure

| Finding | Result |
| --- | --- |
| F1, false helium attribution | **CLOSED.** The 20 eV fixture is below the 24.59 eV He I threshold; code and P1 report now state that helium remains neutral and does not explain the residual. |
| F2, finite `0.03c` reference | **CLOSED with new caveat N1.** `0.03c` is now diagnostic and the inverse-ĉ envelope is 1.8302%, but coordinate sensitivity was not included. |
| F3, hidden mesh allowance | **CLOSED.** The 0.005 degradation guard is in the contract, JSON, validator, and artifact test, with its limited interpretation disclosed. |
| F4, stale B2 provenance | **CLOSED.** The report has HEAD `ca90a391…`; B2 was regenerated and P5 rebound to the current full-core hash. |

The exclusion of `0.001c` was accepted: normalized photon storage is 14% away
from the leading-order asymptote there, versus within 1% for the retained
points; fits including it do not reveal hidden contrary information.

## New findings

- **N1 — MEDIUM:** the adopted padding spans fit order only at fixed `1/ĉ`.
  Repeating the same two-linear-plus-quadratic rule in photon-storage fraction
  gives intercepts `1.008134`, `1.007070`, `1.006946`, upper bound `1.009323`,
  and total envelope about `0.019765` (98.8% of the 2% gate). This reproduces
  the prior audit's estimate and exceeds the new `1/ĉ` upper bound. The report
  therefore understated extrapolation-coordinate uncertainty.
- **N2 — LOW:** the RSLA report incorrectly says old P1 kept step count fixed.
  It used `4*size` and the same grid-independent physical duration; the actual
  correction is the increase from `0.01c` to `0.03c` at fixed duration.
- **N3 — LOW:** the escape-sign tolerance `-1e-4 * emitted` is absent from the
  contract and JSON thresholds.
- **N4 — LOW:** the artifact test reads the mesh allowance from the payload
  instead of independently pinning `0.005`, and does not independently
  recompute extrapolation intercepts from the matrix.

## Verdict

**CONDITIONAL PASS**

F1–F4 are genuinely closed and the scientific conclusion passes under every
coordinate tested, but N1 is the same class of uncertainty the gate is meant
to expose. The auditor required coordinate sensitivity and the corrected P1
history before final PASS; N3 and N4 were recorded as low-severity hardening
items.
