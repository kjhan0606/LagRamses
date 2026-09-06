# Independent stage-3 audit: RSLA and refinement

Audit the current dirty worktree in `/gpfs/kjhan/LRD_JWST` read-only. Do not
edit files. This is the mandatory Claude Opus 5 gate after stage 3 of the SNRT
production/publication-readiness sequence.

Return exactly one overall verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
`PASS` requires that the RSLA matrix, production `0.01c` error bound,
mesh/angular wiring, P1 analytic correction, and provenance are numerically
and scientifically justified with no stage-3 conditional finding left open.
Separate later hydro/source-variability limitations from defects in this gate.

## Original audit finding being closed

Read M4 and its acceptance test in
`provenance/claude_opus5_rt_audit_2026-09-01.md`. It found that P1's old
`0.01c` run was about 20% below the finite-time analytic radius while checking
only grid-to-grid spread, and requested:

1. correct the false front-speed docstring;
2. add an analytic radius assertion to P1;
3. run `ĉ/c = {1e-3, 3e-3, 1e-2, 3e-2}` and show convergence;
4. add mesh/angular refinement data for the B2 4.41% radius discrepancy;
5. determine a defensible production `0.01c` error range.

## Implementation and canonical evidence

Primary files:

- `simulation/snrt/tools/validate_rsla_refinement.py`
- `simulation/snrt/tests/rsla_refinement_artifact.py`
- `simulation/snrt/data/rsla_refinement_validation.json`
- `simulation/snrt/RSLA_REFINEMENT_VALIDATION.md`
- `simulation/snrt/tests/p1_validation.py`
- `simulation/snrt/snrt_core/ionization_front.py`
- `simulation/snrt/P1_CONVERGENCE.md`

The B2 and P5 canonical reports were regenerated/rebound after the core
docstring correction so their source-closure tests remain current.

### Physical/numerical contract

- Fixed B2 H-only Strömgren problem: 32 cubed, S4, 18 eV, nH=0.01 cm^-3,
  T=1e4 K, Q=1e49 s^-1, duration exactly 0.5 recombination times, domain width
  4 Strömgren radii, Courant 0.4, 32 opacity iterations.
- Only `ĉ` and the resulting explicit step count vary across the RSLA matrix.
- Independent conservative H Solver B supplies all four matrix values.
- Production multiphysics Solver A is separately run at `0.01c`; secondary,
  helium, and dust are off to isolate RSLA and match B2.
- Every run hard-gates fixed-point residual, H ledger, finite values, inferred
  escape sign; production additionally hard-gates the electron-root flag.

### Measured matrix

At `ĉ/c = 0.001, 0.003, 0.01, 0.03`, respectively:

- radius ratios: `0.845542, 0.955902, 0.991949, 1.001993`;
- photon-storage fractions: `0.43627, 0.16186, 0.05015, 0.01684`;
- fixed physical duration, step counts `39, 116, 385, 1153`.

The radius is required to be monotone nondecreasing. The `0.03c` reference
must lie within 2% of the infinite-light-speed analytic radius.

### Production `0.01c`

- Solver A radius ratio `0.991892` (0.811% analytic discrepancy).
- Solver A/B mean absolute xHII `4.3528e-6` and radius-ratio difference
  `5.7189e-5`.
- Solver A `0.01c` vs conservative `0.03c` reference: `1.00814%`.
- Production electron-root failures: zero.

The declared production radius error envelope is a deliberately conservative
linear sum, all using direct same-`ĉ` discretization controls:

- production `0.01c` vs conservative `0.03c`: `0.0100814`;
- conservative `0.01c`, 32-to-64 mesh: `0.0021720`;
- conservative `0.01c`, S4-to-S8 angular: `0.0002655`;
- production-vs-conservative radius difference: `0.0000572`;
- total `0.0125762` (1.26%), hard-gated below 2%.

Check whether this is a defensible conservative benchmark-specific bound,
including whether any terms are mixed or misinterpreted. It is explicitly not
claimed as a universal dynamic-source/hydrodynamic RSLA error.

### Refinement data

At original B2 `0.003c`:

- 32 cubed S4 `0.955902`;
- 64 cubed S4 `0.951651`, absolute change `0.004251`;
- 32 cubed S8 `0.954583`, absolute change `0.001319`.

At production `0.01c`:

- 32 cubed S4 `0.991949`;
- 64 cubed S4 `0.989777`, absolute change `0.002172`;
- 32 cubed S8 `0.992215`, absolute change `0.000266`.

The report explicitly says the non-monotone 0.003c mesh result is sensitivity
data, not a Richardson convergence order.

### P1 correction

`ionization_front.py` no longer claims a universal front-speed margin.
`p1_validation.py` uses fixed 22,265.6 yr duration at `0.03c` by using
12 steps per linear grid count and directly asserts each grid is within 10%
of the time-dependent H-only analytic radius. Measured 32/48/64 analytic
ratios are `0.9351, 0.9296, 0.9261`; spatial spread is 0.96%. The document
marks the old self-consistent-but-wrong `0.01c` values superseded.

## Required reproduction

Run from `simulation/snrt`:

```bash
.venv/bin/python tests/rsla_refinement_artifact.py
JAX_PLATFORMS=cpu .venv/bin/python tests/p1_validation.py
.venv/bin/python tests/b2_multiphysics_artifact.py
.venv/bin/python tests/p5_secondary_ionization_artifact.py
```

The full RSLA calculation is source-bound in the canonical JSON; re-run
focused independent calculations if needed. Verify the analytic formula,
fixed-duration construction, source normalization, radius estimator,
photon-storage interpretation, refinement comparison, all criteria, and hash
closure rather than trusting this prompt. Finish with a concise closure table
for M4's five actions and the single overall verdict.
