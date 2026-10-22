# Grackle nine-species integration — 2026-10-06

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001`; origin:
`git@github.com:kjhan0606/LagRamses.git`. Base commit:
`1b527042810d3fae09b69f26611a88449883b435`.

## Implemented scope

Independent `GRACKLE_NEQ=1` build and `cooling_method='grackle_neq'`
selector. Modern Grackle local API evolves HI, HII, HM, HeI, HeII,
HeIII, H2I, H2II and electrons, with equilibrium metal cooling.
Dust and transported RT are disabled; UV background is independently
selectable and requires matching cooling data.

Each OpenMP worker owns a chemistry context. Exact H/He composition
grouping prevents the library's default hydrogen abundance from changing
enriched element inventories. Updates are transactional and check
finite states, element conservation and charge neutrality. Nine species
are advected as additional hydro fields, initialized only on fresh runs,
and protected by restart field-index/NVAR checks. Feedback element
inventory changes are reconciled before chemistry integration.

Namelist generator and mkrun selections were updated together.

## Verification obtained

The minimum 8^3 smoke has since been extended. Current detailed evidence,
metrics and caveats are in
[`grackle_neq_validation_2026-10-06.md`](grackle_neq_validation_2026-10-06.md).
The added checks include dynamic refine/coarsen, exact OMP1/OMP2 and MPI1/MPI2
AMR comparisons, thermal-response and HM2012 UV scale-factor tests, a
delayed-cooling on/off control, a synthetic post-feedback inventory restart,
and the initial one-step 128^3 GRAFIC cosmological IC run. The follow-on Slurm
job 414735 on `syn06` completed matched six-step UV-on/UV-off trajectories,
ending at `a=0.01851` (z≈53), plus a one-step no-cooling hydro control. The
cosmological UV comparison is an HDF leaf-array-order comparison; coordinates
were not cross-matched. These checks passed within their stated fixture scope;
the synthetic inventory test does not invoke the actual stellar feedback
depositor, and the cosmological trajectory does not reach reionization.

## Qualification limits

This is still **not production-qualified**. The real stellar/SN mass return
and delayed-reservoir deposition path has not been exercised through an actual
stellar event with production-approved yields. Cosmological evolution reaches
only `a=0.01851`, not the UV/reionization epoch. A stiff-regime timestep-
convergence check remains open. The post-Godunov/pre-cooling negative-energy
diagnostic is reproduced for the initial step with Grackle disabled, but the
later deterioration is not fully isolated from prior-step cooling history;
this blocks qualification of the coupled hydro run. The six-step 128^3 run
spent 82.7% of wall in `cooling_fine` at one MPI rank and eight OpenMP threads;
this is not a production scaling/performance pass. Hydro uses fixed gamma=5/3,
and metal cooling is equilibrium-table based, not a non-equilibrium metal
reaction network. Dust and transported RT coupling are not part of this
route's verified scope. No production-ready or publication-ready claim is
made.

No global-context changes were made. A read-only Astra audit and its
conditional-pass opinion are recorded in
`grackle_neq_astra_audit_2026-10-06.md`.
