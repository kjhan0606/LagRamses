# Grackle NEQ validation addendum — 2026-10-07

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001`; origin:
`git@github.com:kjhan0606/LagRamses.git`. This addendum supplements
[`grackle_neq_validation_2026-10-06.md`](grackle_neq_validation_2026-10-06.md).
It records the newly closed UV tests, the still-open time-integration gate,
and the cosmological unit-contract correction now under validation.
**This is not production-readiness approval.**

## UV/reionization wiring

The runtime Grackle NEQ initializer now receives RAMSES `z_reion`; for
cosmological runs, RAMSES applies its existing start-redshift cap before
initialization. The bridge switches the configured HM2012 background off for
`z >= z_reion` and on below the transition. The focused regression used
`z_reion=8.5` and verified:

- At `a=0.1` (`z=9`), configured UV-on matched the explicit no-UV sentinel
  exactly in energy and all nine species.
- Forcing UV on at that pre-reionization epoch changed energy by
  `3.9274e-5` relative and produced a maximum absolute species difference per
  total density of `0.44897` relative to the configured gate.
- At `a=0.12` (`z=7.333`), configured UV-on differed from the no-UV sentinel
  by `1.8946e-5` relative in energy and `0.44897` in maximum absolute species
  difference per total density.

The five-case regression completed in Slurm job `415084` on A40 and passed.
The new bridge/test wiring is in `patch/lagRamses/grackle_neq_bridge.c`,
`patch/lagRamses/grackle_neq_runtime.f90`, and
`patch/lagRamses/init_time.kjhan.f90`.

## Matched 128^3 cosmological UV pair

Fresh `GRACKLE_NEQ=1` build and matched UV-on/off RAMSES runs were executed in
job `415084`, one MPI rank and eight OpenMP threads, on `syn06` (A40 allocation;
the solver is CPU-only). The fixed level-7 mesh contained 2,097,152 leaf cells.
Subsequent unit-contract review found that this build passed RAMSES'
supercomoving `scale_t ∝ a²` into Grackle, whose comoving time unit must be
constant. It also passed RAMSES specific energy without converting to
Grackle's `LU/(a TU)` velocity-squared unit. Therefore the run is retained as a
historical integration diagnostic only; its thermal/chemical trajectory is
not valid production evidence and must be repeated after the correction.
Both runs completed 62 steps to `a=0.12218196` (the requested output target was
`a=0.12`) with the same coordinates and scale factor to `3.6e-12`.

Both final states were finite and positive, with maximum normalized H/He/charge
residuals of `7.601e-11`, `2.400e-11`, and `3.90e-16`; no negative
post-hydrodynamics energy diagnostic appeared under the tested
`pressure_fix=.true., beta_fix=0` configuration. Mass-weighted `HII/H` was
`0.99999536` with UV and `0.81895931` without UV. Wall times were `826.568 s`
and `798.116 s`; Slurm reported peak batch RSS `20,279,644 KiB` under the
20-GiB request. This establishes the tested trajectory and late-UV response,
not feedback-event coupling or a broad performance/scaling qualification.

After the comparison metrics were written to
`.grackle-neq-production-qualification-415084/summary.json`, both evaluated
HDF5 `output_00001` directories were removed. Effective namelists, logs,
summary, build log and hashes remain. The job's terminal marker was emitted
before its scope was clarified; its PASS means only that the matched
cosmological trajectory gates passed, not that Grackle is production-ready.
The script now states this scope explicitly.

## Stiff chemistry interval convergence

The timestep fixture covers 16 primordial-network states over code density
`1e-4..1e2` and energy-array values `1e10..1e16`, with UVB active at `a=0.1`.
For this fixture's `density_unit=m_H/cc`, Grackle's comoving conversion gives
proper hydrogen number density `0.76*rho_code cm^-3` (not `/a^3`). Its
`length_units=time_units=1` and `a=0.1` imply Grackle velocity unit
`LU/(a TU)=10 cm/s`; hence the energy-array range represents
`1e12..1e18 erg/g`. Time intervals
are `1e10`, `1e12`, `1e14`, `7e14`, and `1e15 s`, with uniform external
subdivision through 131,072 calls. The latter intervals bracket the live
trajectory's approximately `7e14 s` physical Grackle call near `a=0.12`.

At `dt=7e14 s`, one bridge call differs from the 131,072-subdivision result by
`13.49%` in maximum relative internal energy, `1.286%` in maximum species mass
per total density, `1.466%` in the normalized H inventory L1 norm, and
`10.73%` in the normalized He inventory L1 norm. At `dt=1e14 s`, corresponding
one-call discrepancies are `21.78%`, `1.048%`, `1.485%`, and `8.746%`. The
`dt=7e14 s` 65,536-to-131,072 adjacent differences are `0.0300%` energy,
`0.0374%` species-per-density, `0.0254%` H L1 and `0.3122%` He L1. Thus the
fine reference is stable enough to establish that the one-call production
interval is materially different; the test does not assert which trajectory
is the physically correct answer.

Grackle 3.4.1's solver chooses internal steps with a nominal 10%-of-energy
change limiter (`dtit <= 0.1*energy/edot`) in
[`solve_rate_cool_g.F`](https://github.com/grackle-project/grackle/blob/grackle-3.4.1/src/clib/solve_rate_cool_g.F#L674-L677).
The observed external-step sensitivity still requires an explicit scientific
decision and an accuracy-qualified production treatment. No arbitrary
acceptance tolerance or mitigation has been introduced after observing the
result.

The extended matrix completed in Slurm job `415093` on A40 before the adapter
unit fix. Its log is
`.grackle-dt-matrix-415093.log`; its `assets.sha256` records the exact bridge,
UV data and test-source hashes used. The corrected physical interpretation
supersedes the earlier addendum text; the numerical matrix must be rerun with
the fixed bridge before it is used to set a production timestep policy.

Run identities: job `415084` RAMSES binary SHA-256
`64278763c528fd77a0d44ef98ee737787604b12701ce602ec2df87ddacb8b3c0`; bridge
object `4fdc0853179c1c20a2a80c55c11b9a2a7d392f4cdf6df98f82834ab3e7bb8979`;
HM2012 table `8715f1b39e90a7296ec2adcd442fa13a3d45d2ad021c6fa2fae9e4ab7a4700b2`.
The matrix used test-source hash
`d856b9392f189a25dc9634ebb5ab167de723addd02e0ee45533b0048b74579e7`; the
label/documentation-only current source hash is
`cad5dbeb760c309e9e93d3ef3b9951b576b2dcfc1d3049d70d1bf8a2bcdaf599`.

## Interim production-blocker assessment (historical; superseded below)

The following blocker list and verdict described the state before jobs `415107`,
`415121`, and `415186`. Read it as an intermediate audit record only; current
evidence and status are recorded in the final section below.

- **Cosmological unit contract:** the RAMSES adapter now uses a fixed Grackle
  time unit `scale_t_RAMSES/a²`, the required comoving velocity unit
  `length_units/(a*time_units)`, converted RAMSES specific energy, and
  `dt_Grackle=dt_RAMSES*a²`. The focused unit regression passed; a fresh
  cosmological run is still required before this correction is qualified.
- **Stiff-interval accuracy:** the tested RAMSES chemistry interval is in a
  regime where the single-call map differs by up to 13.5% (energy) and 10.7%
  (He inventory norm) from the refined integration. The acceptable physical
  error and any subcycling/timestep treatment have not been established.
- **Physical source event:** the cosmological pair has star formation and
  stellar feedback disabled. Existing synthetic inventory-restart checks are
  not a real stellar/SN `thermal_feedback` event through the Grackle route.
- **Supported-model limits:** this route evolves nine primordial species,
  uses equilibrium metal cooling, fixed `gamma=5/3`, and has dust and
  transported RT disabled. These limits must remain explicit in any supported
  production configuration.
- **Performance:** the 128^3 paired trajectory ran in about 13 minutes per
  case with 8 OpenMP threads, but no controlled scaling comparison or GPU
  implementation was tested; the measured batch RSS was close to the 20-GiB
  request.

The cosmological UV trajectory gate is **PASS**. Overall Grackle NEQ
production qualification remains **NOT APPROVED** pending the stiff-interval
accuracy decision, a physical source-event integration test, and review of
the resource/performance margin. Astra has been asked to independently review
this scoped evidence. No completion/approval email is sent unless Astra
explicitly approves production readiness for the declared scope.

## Corrected-unit regression evidence

The adapter initializes Grackle with the constant cosmological time unit
`scale_t_RAMSES/a²`; each call uses that fixed time unit, the required
`velocity_units=length_units/(a*time_units)`, converts RAMSES specific energy
to Grackle units (`T2_unit_Grackle=a²*T2_unit_RAMSES`), and converts the
RAMSES timestep (`dt_Grackle=dt_RAMSES*a²`). Non-cosmological calls retain
their units at `a=1`. The bridge rejects a changed Grackle time unit after its
rate context has been initialized.

Slurm job `415105` passed the focused regression at `a=0.1` and `a=0.2`:
proper-density, velocity-unit, and temperature-unit relative differences were
zero; adapter and direct Grackle Local API outputs matched exactly at each
epoch; inventories remained closed. Cross-epoch chemistry states differed by
`2.114e-5`; this is a diagnostic, not an acceptance failure, because
cosmological processes may depend explicitly on redshift. An earlier test
design that required cross-epoch state equality was discarded. The
cooling-enabled variant also included redshift-dependent Compton cooling and
was not used as a unit gate.

This pending statement predates the fresh build/trajectory, source-event
integration, and three-level timestep comparison. Jobs `415084` and `415093`
remain historical evidence only; the subsequent evidence and release status
are below.

## Current production-scope evidence — 2026-10-07

### Corrected units and cosmological UV trajectory

- Focused unit regression `415105`: at `a=0.1` and `0.2`, proper-density,
  comoving velocity-unit, and temperature-unit conversions matched. Adapter
  and direct Grackle Local API outputs matched exactly at each epoch, and
  inventories closed.
- Fresh Fortran build plus qualification `415107`: fixed level-7 mesh,
  2,097,152 leaf cells, one MPI rank/eight OpenMP threads, CPU Grackle solver.
  Matched UV-on/off runs completed to `a≈0.121235`. Both states were finite,
  positive and budget-closed (maximum normalized H/He/charge residuals about
  `7.60e-11`, `2.40e-11`, and `3.64e-16`). No negative post-hydro energy
  diagnostic was recorded. Mass-weighted `HII/H` was `0.9995303` with HM2012 UV
  and `0.0012877` without UV. Runtime was 526 s UV-on and 486 s UV-off.
- Slurm batch peak RSS for `415107` was `21,003,500 KiB` during the fresh
  build; the simulation steps used about `1.73–1.74 GiB`. Later runtime jobs
  peaked around `1.74–1.77 GiB`; the runtime-only convergence job requested
  8 GiB. Build and simulation memory envelopes are distinct.

### Stellar mass return through the Grackle route

Job `415121` ran the user-selected Kroupa + KL16/LC18 + effective-SSP SNIa
configuration with Grackle NEQ, then checked only RAMSES `PTYPE_STAR` records.
Across four snapshots, stellar-particle count grew from 512 to 2,048 and 3,072
persistent stellar-mass-decrease records were measured. All four gas states
were finite, positive, and H/He/charge budget-closed (largest normalized
residual about `2.15e-14`). This verifies the tested mass-return/chemistry
coupling; it is not an SNIa event count, yield calibration, or universal
approval of the selected source tables. The short run does not claim that a
delayed SNIa explosion occurred.

### Three-level hydro/chemistry timestep comparison

The final comparison combines the base and half-step states retained by
`.grackle-neq-cosmo-timestep-convergence-415175` with the quarter-step state
from job `415186`. All used the same level-7 leaf mesh and targeted the UV-on
epoch near `a=0.1063`. Actual output scale factors were `0.1063093`, `0.1062713`,
and `0.1064630` for Courant factors `0.5`, `0.25`, and `0.125`; this is a
near-matched, not mathematically identical-time, comparison. All states were
finite/positive and had closed H/He/charge budgets.

| Pair | p95 internal-energy relative difference | p95 species difference / rho | mass-weighted `HII/H` difference | Finer-run wall time |
|---|---:|---:|---:|---:|
| 0.5 vs 0.25 | 2.475% | `2.827e-4` | `1.240e-4` | 811 s |
| 0.25 vs 0.125 | 1.270% | `3.418e-4` | `6.201e-6` | 1,605 s |

The p95 energy difference decreases by a factor `0.513` between adjacent
refinement pairs; the p95 species metric increases by a factor `1.209`, while
the global mass-weighted HII difference decreases by about 20x. Maximum local
relative differences remain large in a small tail (up to 94% energy for
0.5-vs-0.25 and 246% for 0.25-vs-0.125); they are reported, not hidden behind
an invented threshold. This is a mixed but improving global/thermal trend, not
an independently declared scientific tolerance. The corrected fixed-state
matrix remains in
`.grackle-neq-production-qualification-415107/timestep_convergence.log` and
shows material single-call sensitivity near the live interval.

The first three-level harness run (`415175`) stopped before the quarter-step
target because `nstepmax=160`; its base and half-step outputs were valid and
retained. A quarter-only wrapper attempt (`415184`) failed before RAMSES due to
a missing-file branch. Corrected job `415186` completed the quarter path,
combined all three states, wrote the summary, and then removed the three
evaluated raw HDF5 output directories. These were harness/configuration
failures, not solver failures.

### Supported scope and release status

The candidate scope is CPU `GRACKLE_NEQ`: nine primordial species,
equilibrium metal cooling, fixed `gamma=5/3`, the tested HM2012
UV/reionization gate, and the tested stellar mass-return coupling. Dust,
transported RT, Grackle GPU execution, multi-rank production performance, and
broad physical yield approval are excluded.

The evidence bundle is complete for the declared tests. Production-ready is
**not declared until Astra's final read-only audit issues an explicit verdict**.
No completion email is sent unless Astra explicitly approves this stated
scope.

### Exact-epoch timestep follow-up — 2026-10-07

Job `415210` completed all three cosmological trajectories with finite,
positive states and closed chemistry inventories. Its first analyzer pass
correctly refused to compare the base result at `a=0.1063093333` to requested
`a=0.1063`, because the nominal target fell just below the first snapshot.
The six HDF5 snapshots were retained, not discarded. Analysis-only job
`415236` then evaluated their linearly interpolated conserved fields at the
common `a=0.1065`; it passed the comparison and removed exactly those six
evaluated snapshot directories. The trajectory wall times were 530.110,
856.779 and 1647.185 s for Courant factors 0.5, 0.25 and 0.125; step peak RSS
was 1.66, 1.70 and 1.70 GiB. The interpolated states were finite/positive and
had maximum normalized H/He/charge residuals below `7.61e-11`, `2.41e-11` and
`6.61e-16`.

These are **interpolation-only diagnostics**, not the final accuracy gate.
At `a=0.1065`, p95 internal-energy differences were 2.5865% (0.5 vs 0.25),
1.7321% (0.25 vs 0.125), and 1.0848% (0.5 vs 0.125). The corresponding p95
species-abundance differences normalized by total density were
`2.779e-4`, `4.878e-4`, and `6.065e-4`; mass-weighted `HII/H` differences
were `1.129e-4`, `3.904e-4`, and `2.776e-4`. Species-resolved p95 abundance
differences (`|Δ species|/rho`) were:

| Pair | HI | HII | HeI | HeII | HeIII | e |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 vs 0.25 | `1.137e-4` | `4.687e-4` | `3.950e-5` | `1.495e-4` | `3.249e-8` | `5.061e-4` |
| 0.25 vs 0.125 | `2.988e-4` | `7.307e-4` | `1.526e-4` | `2.869e-4` | `1.188e-7` | `8.021e-4` |
| 0.5 vs 0.125 | `2.377e-4` | `1.033e-3` | `1.301e-4` | `3.674e-4` | `1.512e-7` | `1.123e-3` |

The adjacent-pair convergence is not monotonic for abundance or global
ionization metrics; local maxima are also large in sparse tails. Since this
diagnostic linearly interpolates between outputs, it cannot distinguish that
behavior from interpolation error. A direct equal-epoch run is queued as job
`415237`: the existing RAMSES `match_aout=.true.` path will land all three
Courant trajectories on exactly `aout=0.1065`, with one output per case. The
one-snapshot exact-time comparison will supersede this diagnostic if it
completes. No scientific tolerance has been invented, and production
readiness remains pending the explicit timestep-accuracy decision and Astra's
final verdict.

### Direct exact-epoch timestep comparison — 2026-10-07 (supersedes interpolation)

Job `415237` reran the same three Courant trajectories with RAMSES
`match_aout=.true.` and one scheduled output at `aout=0.1065`. The final
cosmological step was clipped to land on the requested epoch; HDF5 headers
were `0.10650000000000019`, `0.10650000000000023`, and
`0.10650000000000004`. This supersedes the linearly interpolated `415210`
diagnostic and removes output-epoch interpolation from the timestep comparison.
All three runs completed with finite/positive states and closed chemistry
budgets (maximum normalized H/He/charge residuals `7.61e-11`, `2.41e-11`,
`3.64e-16`). Their wall times were 477.137, 861.377 and 1696.423 s at
Courant factors 0.5, 0.25 and 0.125; per-step peak RSS was 1.67, 1.66 and
1.66 GiB. The three exact-epoch HDF5 outputs were deleted only after the
comparison succeeded; compact summaries, effective inputs, logs and hashes
remain under `.grackle-neq-cosmo-timestep-convergence-415237/`.

| Pair | p95 internal-energy difference | p95 species difference / rho | mass-weighted `HII/H` difference | mass-weighted `HeIII/He` difference |
|---|---:|---:|---:|---:|
| 0.5 vs 0.25 | 2.543% | `2.544e-4` | `1.776e-5` | `1.636e-5` |
| 0.25 vs 0.125 | 1.765% | `3.080e-4` | `1.312e-5` | `1.108e-5` |
| 0.5 vs 0.125 | 1.009% | `5.524e-4` | `4.643e-6` | `5.282e-6` |

Species-resolved p95 `|Δ species|/rho` for `(HI,HII,HeI,HeII,HeIII,e)` were:

| Pair | HI | HII | HeI | HeII | HeIII | e |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 vs 0.25 | `1.686e-5` | `4.452e-4` | `7.836e-6` | `1.406e-4` | `2.732e-8` | `4.803e-4` |
| 0.25 vs 0.125 | `1.210e-5` | `5.467e-4` | `1.082e-5` | `1.731e-4` | `1.169e-7` | `5.899e-4` |
| 0.5 vs 0.125 | `8.691e-6` | `9.745e-4` | `5.606e-6` | `3.074e-4` | `1.441e-7` | `1.051e-3` |

The global mass-weighted ionization differences shrink with finer timesteps,
and p95 thermal differences improve. However, the p95 aggregate species
metric and p95 HII/e abundance differences increase across adjacent
refinements; local maxima remain large in sparse tails (energy maxima 94%,
244%, 114%; maximum species difference per total density 0.412, 0.240, 0.480
for the three pairs). Thus the exact-time result still does not demonstrate
uniform monotonic state convergence. It quantifies the tradeoff but does not
choose an acceptable application error or production timestep policy. No
post-hoc threshold has been applied. Astra's final verdict is **CONDITIONAL**;
production-ready status remains pending the owner's accuracy-policy decision.

### Astra final read-only audit — 2026-10-07

Astra answered **Q-GOAL: yes** and **Q-LEAN: proportionate**. It issued a
**CONDITIONAL** verdict, not production approval. The reviewer found that the
exact-epoch comparison resolves the output-time mismatch, all states remain
finite/positive with closed budgets, and the bounded CPU Grackle route
materially advances the project. Dust, transported RT, GPU, multi-rank
performance, delayed SNIa occurrence and broad yield calibration need not be
added to this scoped approval gate.

The remaining blocker is the owner-defined numerical-accuracy contract:
which quantities govern acceptance, their scientifically justified errors,
and the supported production timestep/subcycling policy. The p95 species
comparison currently uses species density difference divided by the finer
run's total density, so it mixes chemistry and gas-density differences; large
local maxima also need interpretation before attributing them to chemistry.
These mixed convergence metrics are not proof of solver failure, but
positivity and inventory closure alone do not establish temporal accuracy.
No further large benchmark campaign is recommended before this policy is
specified. The requested email to `kjhan0606@gmail.com` was **not sent**:
Astra did not explicitly approve production readiness.

### Corrected elemental-abundance and local-tail follow-up — 2026-10-07

At the operator's direction, the exact-epoch comparison was repeated using
per-cell species fractions normalized by the corresponding elemental mass
density, rather than `|Δ species density|/rho`. H2 species are reported as
fractions of hydrogen nuclei carried by that species. Local-error tails are
reported both by volume and by gas mass, with pair mass weighted symmetrically
as `0.5*(rho_coarse+rho_fine)*cell_volume`. Exceedance bins are descriptive
only, not acceptance thresholds.

Slurm job `415379` completed on H100 at the same exact `a=0.1065` and level-7
leaf mesh. All trajectories were finite and positive; maximum normalized
H/He/charge residuals remained `7.61e-11`, `2.41e-11`, and `3.64e-16`, while
per-cell elemental abundance closure was within `5.6e-16`. Wall times were
451.101, 809.241 and 1589.493 s (Courant 0.5, 0.25, 0.125); peak RSS was
1.69, 1.66 and 1.66 GiB.

| Courant | mass-weighted HII/H | mass-weighted HeIII/He |
|---:|---:|---:|
| 0.5 | 0.99859949 | 0.0009647493 |
| 0.25 | 0.99861726 | 0.0009811049 |
| 0.125 | 0.99860414 | 0.0009700424 |

| Pair | p95 thermal energy-density difference (volume / mass weighted) | p95 `|Δ(HII/H)|` | p95 `|Δ(HeIII/He)|` | p95 `|Δ(e/H)|` | global HII/H, HeIII/He differences |
|---|---:|---:|---:|---:|---:|
| 0.5 vs 0.25 | 2.543% / 2.567% | `2.194e-5` | `1.140e-7` | `2.449e-5` | `1.777e-5`, `1.636e-5` |
| 0.25 vs 0.125 | 1.765% / 1.752% | `1.559e-5` | `4.873e-7` | `1.904e-5` | `1.313e-5`, `1.106e-5` |
| 0.5 vs 0.125 | 1.009% / 1.048% | `1.037e-5` | `6.010e-7` | `9.565e-6` | `4.642e-6`, `5.293e-6` |

For the three pairs, the fractions of gas mass above a 5% thermal difference
were 0.0999%, 0.0625%, and 0.1249% (volume fractions 0.0235%, 0.0080%, and
0.0291%). Thus the maximum thermal differences (94%, 245%, 115%) occur in a
small tail, while the p95 thermal difference decreases with timestep
refinement. HII/H abundance p95 also decreases across these pairings. The
HeIII/He p95 is very small but is not monotonic, reflecting the sparse,
weakly populated HeIII component.

Large local composition maxima are now quantified rather than treated as
unweighted failures. The HeIII/He absolute-difference maxima are 0.944, 0.994,
and 0.996; cells exceeding an absolute difference of 0.01 number 49, 27, and
49, respectively, or about 0.0013–0.0023% of volume and 0.0032–0.0054% of
pair-weighted gas mass. These tails are small but remain visible and should
not be concealed by percentile summaries. The exact per-species percentiles
and tail distributions are in
`.grackle-neq-cosmo-timestep-convergence-415379/timestep_convergence_summary.json`.

The corrected result supports the same scoped interpretation: the adjacent
p95 thermal-density and HII/H differences decrease with refinement, while
HeI/HeII/HeIII adjacent-pair percentiles are mixed/nonmonotonic. The fraction
of pair-weighted gas mass with a thermal-energy-density difference above 1%
is 99.925%, 99.976%, and 7.834% for the three pairings. Thus the >5% thermal
discrepancy is a sparse tail, but the >1% discrepancy is not; the selected
energy-density metric must not be summarized as a tail-only effect. This
quantifies timestep sensitivity without choosing an owner acceptance
tolerance or declaring universal temporal convergence.
The job removed only the three evaluated HDF5 snapshots after successful
analysis; the exact submitted batch script, effective namelists, logs, hashes,
summary and cleanup manifest are retained under
`.grackle-neq-cosmo-timestep-convergence-415379/`. The corrected-metric
Astra re-audit and verdict are recorded below; no email has been sent.

### Astra re-audit of corrected metrics — 2026-10-07

Astra answered **Q-GOAL: yes** and **Q-LEAN: proportionate**, and returned
**CONDITIONAL**, not production approval. It confirmed the former
species-density/total-density confound is resolved by per-run elemental
abundance fractions, the symmetric mass weighting is appropriate, and
unweighted volume percentiles are valid for this fixed uniform level-7 mesh.
It also confirmed the submitted script and recorded executable/source hashes
match.

The audit emphasized two limitations that must remain explicit. First,
HeI/HeII percentiles are not monotonic either, and the nearly unit HeII/HeIII
local maxima, although limited to 27–49 cells, are retained as visible tails.
Second, the current thermal comparison is **internal-energy density**, not
specific internal energy or temperature; it shows >1% differences across
about 99.9% of gas mass for adjacent Courant comparisons, so only the >5%
discrepancy is a sparse tail. Astra did not interpret these observations as
solver failure, but said they cannot establish application-level acceptance
without an owner-defined accuracy/timestep contract specifying quantities,
bulk and tail tolerances, and the supported Courant/subcycling choice. The
fixed-state interval sensitivity remains relevant.

The requested notice to `kjhan0606@gmail.com` was **not sent** because Astra's
verdict remains conditional. Production readiness is pending the owner's
numerical-accuracy policy; broader tests are not justified until that policy
is set.

### Concrete accuracy/timestep proposal prepared — 2026-10-07

The operator authorized preparation of the proposed accuracy/timestep
contract. Its concrete quantities, limits, candidate Courant profile and
minimal prospective comparison are recorded in
[`grackle_neq_accuracy_timestep_policy_2026-10-07.md`](grackle_neq_accuracy_timestep_policy_2026-10-07.md).
The owner subsequently adopted its numerical limits and candidate timestep
treatment on 2026-10-07 ("자네가 말한 수치 기준을 확정하자"). Existing
results have not been relabeled as a pass; the prospective specific-energy/
temperature comparison is authorized and still pending. Astra's production
verdict remains CONDITIONAL until that qualification is evaluated.
