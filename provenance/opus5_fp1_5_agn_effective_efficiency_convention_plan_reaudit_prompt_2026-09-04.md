# Claude Opus 5 pre-implementation plan re-audit — F-P1.5-R

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Read-only plan re-audit. Do not edit files, run jobs, activate a runtime,
commit, or push. The initial Opus plan audit returned `APPROVE WITH CHANGES`;
its mandatory amendments are now integrated into the plan. No implementation
has started and explicit user approval is still required.

Read:

1. `provenance/fp1_5_agn_effective_efficiency_convention_bundle_plan_2026-09-04.md`
2. `provenance/opus5_fp1_5_agn_effective_efficiency_convention_plan_audit_2026-09-04.md`
3. `patch/lagRamses/sink_particle.kjhan.f90`
4. `patch/lagRamses/snrt_ramses_driver.f90`
5. `patch/lagRamses/snrt_agn_source.f90`
6. `patch/lagRamses/pm_commons.f90`
7. `patch/lagRamses/amr_parameters.jaehyun.f90`
8. `bin/Makefile`
9. `simulation/snrt/tools/audit_agn_coarse_ledger.py`
10. `simulation/snrt/tests/agn_ledger_transaction.py`

Check whether the amended plan is now complete, feasible, and correctly
bounded for the final production/publication objective of RT, stellar/AGN
feedback, and dust. In particular verify:

- the helper has explicit `spin_bh` semantics and reproduces the actual
  `spin_bh=.false.` default rather than reading a zero/uninitialized
  `eps_sink`;
- raw/base efficiency is strict `(0,1)` while effective efficiency is
  `[0,1)` and MAD/zero-Eddington/floor handling is identical;
- `dMBH_coarse` is the supplied inflow increment and `dMsmbh` is the retained
  BH increment, with the photon-budget API and `(1-epsilon_eff)` consistency
  contract stated without ambiguity;
- the shared helper's placement in unconditional `MODOBJ`, direct consumer
  prerequisites, missing driver→`amr_commons`/`pm_commons` edges, and default
  non-SNRT build evidence are correctly specified;
- raw provenance remains separate, the unapproved `0.5` spectral split stays
  outside scope, and no runtime activation or physical SED/dust/yield claim is
  implied.

Return a concise verdict within 1200 words using exactly one of:

- `APPROVE` — plan is complete and implementation may begin only after
  explicit user approval;
- `APPROVE WITH CHANGES` — list mandatory bounded changes;
- `REJECT` — identify the blocker and replacement boundary.

Cite paths and approximate line numbers. End with an ordered acceptance
checklist. Do not modify files or run commands.
