# B2 production-solver transport validation

Date: 2026-09-02
Project root: `/gpfs/kjhan/LRD_JWST`
Gate verdict: **PASS — source-bound FS2010/coupled-electron rerun**

## Scope

B2 tests the H-only production multiphysics path in `snrt_core/multiphysics.py`,
called Solver A in the 2026-09-01 Opus RT audit. It does not promote the whole
SNRT project, the physical dust model, source SEDs, the thermal atlas license,
or live RAMSES coupling. Those remain later gates.

The gate implements the audit's fixed benchmark: 32 cubed cells, S4, one
18 eV group, `n_H=0.01 cm^-3`, `T=10^4 K`, `Q=10^49 s^-1`, reduced light speed
`3e-3 c`, and duration `0.5 t_rec`. The baseline has zero dust and disables
secondary ionization. Solver B is the independent H-only conservative solver
in `snrt_core/conservative_hydrogen.py`. The baseline, dust, and Solver-B
ionization-front runs set `n_He=0`. The matched 200 eV secondary ON/OFF control
uses `n_He/n_H=0.079`, because the FS2010 table is a primordial H/He closure.
B2 still makes no front/convergence validation claim for coupled helium; that
path remains for a later H+He gate.

## Algorithm repair

The former Solver-A path capped local gas attenuation by the initial neutral
atom inventory and returned the excess photons to the radiation field. B2
removes that cap. Absorbed primary and secondary counts are converted to rates
at the current opacity iterate. H uses an analytic neutral-fraction relaxation
and its solved time-averaged H I fraction; helium uses its backward-Euler end
state. In Solver A, local photoionization, collisional ionization,
recombination, and the electron density are closed by a nearest-root scalar
bisection; this preserves the exact neutral solution for a hot cell with no
photons while avoiding the unstable electron-density Picard loop once a cell
is ionized. The bracket-success diagnostic is a hard gate. The backward-Euler
root remains timestep dependent, so large-step branch accuracy belongs to the
declared stage-5 timestep-convergence gate. Twenty under-relaxed outer
iterations synchronize chemistry and opacity. The secondary target-
availability mask is fixed from the start-of-step species inventory, so newly
created He II cannot switch its own secondary channel on halfway through that
fixed point. Unabsorbed photons remain in
the field through the converged opacity, not through an atom-count clipping
operation.

The compatibility diagnostic `gas_absorption_scale` is now identically one.
P4 and P5 outputs additionally record the cap-activation fraction, minimum
scale, and maximum fixed-point residual; their production CLIs reject fewer
than 20 opacity iterations. They now also fail on an unbracketed electron root
or a photoelectron-energy ledger violation. If an FS2010 target species is
absent, its energy is returned to heat. The zero-helium path therefore neither
creates He ionizations nor relies on overflow through a density floor.

Solver B intentionally retains its simpler H-only closure: electron density is
evaluated from the mean H II opacity state, collisional ionization is omitted,
and no coupled nearest-root routine is shared with Solver A. The two solvers
share transport, cross sections, the analytic H relaxation, and case-B
recombination, but the chemistry iteration remains independent. Their field
difference is accepted below `5e-5` absolute mean xHII (about 0.25% of this
fixture's mean); the analytic Strömgren radius is the fully independent
physical reference.

## Recorded result

| Check | Result | Threshold | Status |
| --- | ---: | ---: | --- |
| Solver-A analytic radius ratio | `0.9556850` | relative error `< 0.05` | PASS |
| Solver-A vs Solver-B mean absolute xHII | `1.52855e-5` | `< 5e-5` | PASS |
| Worst H ledger L1 relative error | `2.85782e-5` (`secondary_200ev_off`) | `< 1e-3` | PASS |
| Maximum opacity fixed-point residual | `2.39611e-5` | `< 1e-4` | PASS |
| 200 eV secondary-ON H/He II/He III residuals | `5.36e-7 / 1.28e-5 / 7.51e-6` | each `< 1e-4` | PASS |
| Retired-limiter structural guard | all runs: active fraction `0`, scale `1` | exact invariant | PASS |
| Solver-A S8 vs A192 shadow difference | `0.00439426` | `< 0.02` | PASS |

The S8 and A192 shadow transmissions are `0.1333715` and `0.1339602`,
respectively. This test uses the full Solver-A dust transport with inert gas,
a 48 cubed mesh, a 12-cell clump radius, clump absorption 8, and 150 steps.
Because all gas cross sections are zero, one opacity iteration is sufficient;
this is a transport-wrapper check, not a reactive-chemistry shadow test.

## Controlled physics deltas

These are wiring controls, not approved LRD physical prescriptions.

- Adding a deliberately strong `1e-20 cm^2/H` dust cross section at 18 eV
  absorbs `20.18%` of all absorbed photons and changes mean xHII from
  `0.0224783` to `0.0184787` (`delta=-0.00399958`).
- In the primordial-composition 200 eV control, enabling FS2010 changes mean
  xHII from `0.00364859` to `0.0227100` (`delta=+0.0190614`) and produces
  `0.603538` secondary H ionizations per emitted photon.

The gate uses fixture regression bands, not only signs: dust absorbed fraction
must lie in `[0.10,0.30]`, its mean-xHII delta in `[-0.006,-0.002]`, secondary
yield in `[0.50,0.75]` per emitted photon, and its mean-xHII delta in
`[0.015,0.025]`. Baseline, dust, secondary-off, secondary-on, and Solver B
must all independently satisfy their fixed-point and H-ledger thresholds.

## Reproduction and artifacts

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
JAX_PLATFORMS=cpu .venv/bin/python tools/validate_multiphysics_b2.py \
  --output data/b2_multiphysics_transport_validation.json
.venv/bin/python tests/b2_multiphysics_artifact.py
```

- Passing JSON: `data/b2_multiphysics_transport_validation.json`, SHA256
  `dabb8ca58d27dfb0f0ed06ddb485f511f3e3bd67737b5ac2971a1f1d7c700360`.
- Validator SHA256:
  `16f53d809f79f77092e07d06a80073e9e5b2a9ab2fd83f3e42906f7b3cb42fef`.
- JAX `0.11.1`, CPU backend, repository HEAD recorded in the artifact as
  `f792791069810c92c3e23c2339cc16370cef15b3` with source-file hashes because
  the worktree is intentionally dirty.
- The `failed_fp12` and `failed_fp16` JSON files preserve the convergence and
  zero-He underflow failures encountered during implementation. They are not
  passing gate artifacts. The `pass_pre_api_cleanup` file preserves the
  numerically identical pass before removal of two obsolete no-op API knobs;
  the `conditional_pre_gate_tightening` file preserves the first pass reviewed
  by Opus 5. The canonical JSON above is the current source-bound artifact.

The artifact contract test recursively rejects non-finite values and verifies
the validator plus all 32 Python modules currently in `snrt_core`, including
transport, chemistry, cross sections, dust, secondaries, quadrature, shadow,
and source deposition. Any validator or `snrt_core` edit, addition, or removal
therefore requires a fresh B2 run.

## Independent audit

Claude Opus 5 first returned a conditional pass, then verified each repair in
two focused re-audits. The final H-only B2 verdict is PASS with no remaining
B2 blocker:

- `provenance/claude_opus5_b2_audit_2026-09-01.md`;
- `provenance/claude_opus5_b2_reaudit_2026-09-01.md`;
- `provenance/claude_opus5_b2_final_reaudit_2026-09-01.md`.

The final audit explicitly leaves coupled helium, a reactive-gas shadow,
resolution refinement, physical dust/SED approval, and live RAMSES coupling to
later gates. B2 PASS is not overall production or publication readiness.
