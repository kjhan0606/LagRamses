# B2 production-solver transport validation

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Gate verdict: **PASS — independently re-audited by Claude Opus 5**

## Scope

B2 tests the H-only production multiphysics path in `snrt_core/multiphysics.py`,
called Solver A in the 2026-09-01 Opus RT audit. It does not promote the whole
SNRT project, the physical dust model, source SEDs, the thermal atlas license,
or live RAMSES coupling. Those remain later gates.

The gate implements the audit's fixed benchmark: 32 cubed cells, S4, one
18 eV group, `n_H=0.01 cm^-3`, `T=10^4 K`, `Q=10^49 s^-1`, reduced light speed
`3e-3 c`, and duration `0.5 t_rec`. The baseline has zero dust and disables
secondary ionization. Solver B is the independent H-only conservative solver
in `snrt_core/conservative_hydrogen.py`. Every ionization-front run in this
gate sets `n_He=0`; B2 therefore makes no validation claim for coupled helium
transport/chemistry. That path remains for a later H+He gate.

## Algorithm repair

The former Solver-A path capped local gas attenuation by the initial neutral
atom inventory and returned the excess photons to the radiation field. B2
removes that cap. Absorbed primary and secondary counts are converted to rates
at the current opacity iterate. H uses an analytic neutral-fraction relaxation
and its solved time-averaged H I fraction; helium uses its backward-Euler end
state. Twenty under-relaxed fixed-point iterations synchronize chemistry and
opacity. Unabsorbed photons remain in the field through the converged opacity,
not through an atom-count clipping operation.

The compatibility diagnostic `gas_absorption_scale` is now identically one.
P4 and P5 outputs additionally record the cap-activation fraction, minimum
scale, and maximum fixed-point residual; their production CLIs reject fewer
than 20 opacity iterations. The zero-helium path also evaluates count per time
before division by its density floor, preventing an underflow-generated `0/0`
in code-unit shadow tests. The same strict-positive floor and divide ordering
are regression-tested for the sibling conservative solvers.

Solver A and Solver B intentionally share the transport operator, H
relaxation, cross sections, and case-B recombination coefficient. Their
differential is therefore a wrapper/wiring cross-check; the analytic
Strömgren radius is the independent physical reference.

## Recorded result

| Check | Result | Threshold | Status |
| --- | ---: | ---: | --- |
| Solver-A analytic radius ratio | `0.9559179` | relative error `< 0.05` | PASS |
| Solver-A vs Solver-B mean absolute xHII | `1.15422e-6` | `< 1e-5` | PASS |
| Worst H ledger L1 relative error | `5.60382e-5` (`secondary_200ev_on`) | `< 1e-3` | PASS |
| Maximum opacity fixed-point residual | `2.39611e-5` | `< 1e-4` | PASS |
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
  absorbs `20.19%` of all absorbed photons and changes mean xHII from
  `0.0224948` to `0.0184912` (`delta=-0.00400357`).
- At 200 eV, enabling the current high-energy secondary prescription changes
  mean xHII from `0.00420259` to `0.0164571`
  (`delta=+0.0122545`) and produces `0.387678` secondary H ionizations per
  emitted photon.

The gate uses fixture regression bands, not only signs: dust absorbed fraction
must lie in `[0.10,0.30]`, its mean-xHII delta in `[-0.006,-0.002]`, secondary
yield in `[0.20,0.60]` per emitted photon, and its mean-xHII delta in
`[0.008,0.018]`. Baseline, dust, secondary-off, secondary-on, and Solver B
must all independently satisfy their fixed-point and H-ledger thresholds.

## Reproduction and artifacts

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
JAX_PLATFORMS=cpu .venv/bin/python tools/validate_multiphysics_b2.py \
  --output data/b2_multiphysics_transport_validation.json
.venv/bin/python tests/b2_multiphysics_artifact.py
```

- Passing JSON: `data/b2_multiphysics_transport_validation.json`, SHA256
  `dbc21fee12278f1288db585276e403357015d5f737e2aee0f5924af4eb59038e`.
- Validator SHA256:
  `e38cac82753c42ca2e9cad1a669187b765fe0836deea810db924b337574f8dbc`.
- JAX `0.11.1`, CPU backend, repository HEAD recorded in the artifact as
  `45cee683788fdfb0e9aa66a978d1e56c1af2c10c` with source-file hashes because
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
