# Claude Opus 5 plan re-audit — F-P1.5-R — 2026-09-04

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Auditor: Claude Opus 5 (primary; read-only)  
Prompt: `provenance/opus5_fp1_5_agn_effective_efficiency_convention_plan_reaudit_prompt_2026-09-04.md`  
Result: completed; no files modified, jobs run, or runtime activated  
Model cost reported: approximately USD 1.72

## Verdict

**APPROVE WITH CHANGES**

Opus confirmed that the amended plan remains the correct next bounded
engineering task, but found seven additional mandatory clarifications.

## Mandatory clarifications integrated into the plan

1. `dMBH_coarse` is Bondi supply, not the writer's supplied inflow.  The
   driver must account the per-step increment of
   `min(dMBH_coarse,dMEd_coarse)` and use `dMsmbh` only for a one-sided
   retained-mass check.
2. The retained-mass relation is not an equality because `accrete_bondi`
   clips by the gas floor and suppresses/changes supply in saved-energy
   states.  The plan now uses the one-sided bound
   `dMsmbh <= (1-epsilon_eff)*dM_inflow*(1+tol)`, names the admitted slack,
   and fail-closes unexplained excess or singular effective efficiency.
3. `spin_bh` defaults to `.true.`, `eps_sink` starts at zero, and the
   spin-enabled first-step accretion read is unguarded.  The plan now records
   the deliberate `.1` ledger/photon fallback as a visible divergence rather
   than claiming exact accretion parity; `spin_bh=.false.` remains the
   deliberate `.1` branch.
4. The helper must define the `eps_sink >= 1` case and preserve the asymmetric
   raw `(0,1)` versus effective `[0,1)` contract.
5. The AGN audit tool must accept an explicit helper source path and record its
   SHA256 so the one-helper/two-consumer claim is machine-verifiable.
6. Both stale mismatch assertions are named for update: the audit tool's
   limitation and `simulation/snrt/tests/agn_ledger_transaction.py:157–158`.
7. `snrt_agn_source_smoke.f90` is recorded as the third positional consumer
   of the photon-budget API.

## Additional source-confirmed guidance

- The driver source loop's photon-budget mass convention must be supplied
  inflow. `dMBH_coarse` accumulates Bondi supply, `dMEd_coarse` supplies the
  Eddington cap, and `dMsmbh` is the retained BH increment.
- The `0.5 * snrt_group_energy_fraction` spectral split, the latched
  `SNRT_RT_ENABLE` gate, and all-or-nothing group transaction are outside the
  physics change and must remain explicitly unchanged/unapproved.
- The helper belongs in unconditional `MODOBJ`, while direct consumer edges
  and `snrt_ramses_driver.o`→`amr_commons`/`pm_commons` dependencies are needed
  for parallel and default-build correctness.  The hydro module is already
  transitively reached through `pm_commons.o`.
- The deferred-`Esave` mass carry-over case should be tested after switching
  accounting from retained to limited supplied inflow.

The amended plan now incorporates these conditions.  Implementation remains
paused pending explicit user approval, as required by the project governance.
