# Approved Grackle NEQ accuracy and timestep policy

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001`,
`kjhan0606/LagRamses`.

Status: **APPROVED NUMERICAL POLICY**, adopted by the operator on 2026-10-07
("자네가 말한 수치 기준을 확정하자"). The limits below govern prospective
qualification. Their adoption is not production-readiness approval: the
specific-energy/temperature refinement comparison remains to be performed.
Astra's latest production verdict is CONDITIONAL; the completion email
remains unsent.

## Purpose and supported use

This policy allocates a numerical-error budget for the CPU Grackle NEQ
galaxy-formation route: nine primordial species, equilibrium metal cooling,
fixed gamma=5/3, the tested HM2012 UV/reionization switch and the independently
tested stellar mass-return connection. Dust and transported RT remain off.
The present live accuracy evidence is the fixed level-7, 128^3 cosmological
trajectory ending at a=0.1065. Acceptance would qualify this declared
configuration and numerical treatment, not all densities, epochs,
refinement patterns, or feedback-event timescales.

The adopted limits are operational accuracy budgets for bulk thermal and
chemical evolution. They are not universal literature tolerances, an
observational calibration, or estimates of absolute error relative to an
exact solution. Existing results inform their practicality, so a new,
prospective refinement pair must test them. Applications
requiring accurate individual-cell helium transitions or line emission
cannot rely solely on the bulk/tail allowances below.

## Quantities and adopted limits

Compare candidate and smaller-step trajectories at exactly the same epoch
and on the same leaf mesh. Define thermal energy density
`E_th = E_total - |momentum|^2/(2*rho)` and specific thermal energy
`epsilon = E_th/rho`. Obtain temperature using the Grackle temperature API
with each run's actual species state and the verified cosmological units.
Temperature must not use a fixed molecular weight. Retain E_th and rho
differences as explanatory diagnostics.

Relative differences for positive epsilon and T are
`abs(q_candidate-q_reference)/q_reference`, without a floor that hides cold
cells. Species differences are absolute differences of each run's species
fraction of its own elemental inventory. H2 fractions count hydrogen
nuclei; `e/H` counts electrons per hydrogen nucleus.

| Quantity | Adopted acceptance limit |
|---|---|
| epsilon and T relative differences | Both volume-weighted and mass-weighted p95 <= 2% |
| Global inventory ratios HII/H and HeIII/He | Absolute difference <= 1e-4 for each ratio |
| Each elemental species fraction, and e/H | Both volume-weighted and mass-weighted p95 absolute difference <= 1e-3 |
| Cells with epsilon OR T relative difference > 10% | Union occupies <= 0.1% of total volume AND <= 0.1% of gas mass |
| Cells with ANY elemental species or e/H absolute difference > 0.01 | Union occupies <= 0.01% of total volume AND <= 0.01% of gas mass |
| State and inventories | Finite, rho/epsilon/T/H/He positive, species nonnegative; maximum H/He/charge residual normalized by rho < 1e-6 |

The 2% bulk thermal budget and 1e-4 global ionization budget express the
approved project accuracy target. The latter is 0.01 percentage points;
the 1e-3 local-composition budget is 0.1 percentage points. The tail limits
allow localized discrepancies while requiring their affected volume AND
mass to remain small. Maximum errors and cell counts must still be reported.
The owner adopted these numerical choices as the project's accuracy budget;
small tails alone are not physical validation of their individual cells.

Use symmetric gas-mass weights
`0.5*(rho_candidate+rho_reference)*cell_volume` for paired percentiles and
tail fractions. Volume percentiles use cell volumes; equal-cell quantiles
are equivalent only on a uniform grid. Global HII/H and HeIII/He are ratios
of separately integrated species and elemental inventories in each run.
Tail masks must be unions, rather than checking species independently and
potentially missing the aggregate affected population.

## Adopted candidate timestep treatment

The candidate supported production profile uses `courant_factor=0.125`,
with the tested `nsubcycle=1`. These are profile settings, not a change to
the repository's general hydro defaults. Grackle's existing internal
chemistry/cooling substeps and failure/transaction checks remain active;
no unimplemented chemistry subcycling parameter is implied.

Qualify this candidate against `courant_factor=0.0625`, keeping other physics
and refinement parameters identical. Prior 0.25-versus-0.125 results
characterize the coarser configuration; they do not by themselves establish
the accuracy of 0.125. A global Courant limit also does not guarantee a bound
on the chemistry timestep relative to cooling or reaction times.

If an adopted limit fails, inspect the corresponding states/intervals and
adjust the actual timestep or local chemistry treatment, then repeat only
the affected verification. Do not relax adopted limits merely to obtain a
pass. New physical regimes need a relevant bounded accuracy check; the
existing stiff fixed-state discrepancies remain an explicit limitation until
the chosen treatment is shown to control them in its intended regime.

## Authorized qualification follow-through

Extend the existing exact-epoch analyzer with epsilon, Grackle temperature,
weighted percentiles, and the two union-tail masks above. Test the analyzer
with density-rescaling and composition-dependent-temperature controls before
running the pair. Reuse the corrected executable, IC and UV table; no new
benchmark framework or broad resolution campaign is required.

Run only the candidate/reference pair on H100 or H200 through Slurm, one MPI
rank/eight OpenMP threads, with explicit memory/GRES and sufficient wall time.
Use a step ceiling >= 500 to accommodate the smaller Courant factor, plus
the existing exact epoch/stopping checks. Each case has one scheduled dump
at aout=0.1065, match_aout=true; suppress periodic dumps/backups. Before
submission, verify effective namelists, estimate the two dump sizes, and
check free storage. Preserve snapshots through complete evaluation and any
needed follow-up analysis, then remove only the evaluated raw outputs and
retain inputs, logs, executable/source identities, compact results and a
cleanup manifest.

Evaluate every adopted limit once, then request Astra's final read-only
verdict on the scoped evidence. Only an explicit APPROVE triggers the
authorized completion email to `kjhan0606@gmail.com`, including the audit
opinion. Policy approval authorizes this qualification; an email still
requires Astra's explicit production APPROVE.

## Existing evidence versus this policy

Job 415379 completed its three trajectories and inventory checks. Its
0.25-versus-0.125 energy-density p95 was 1.765%, and global HII/H and
HeIII/He differences were 1.313e-5 and 1.106e-5. However, epsilon/T,
aggregate union tails and a 0.125-versus-0.0625 comparison were not measured.
The existing evidence therefore cannot be relabeled as passing this policy.

See [validation and Astra findings](grackle_neq_validation_2026-10-07.md).
