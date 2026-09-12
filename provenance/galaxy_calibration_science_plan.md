# Galaxy calibration and convergence — small-box pilot

Status2026-09-12: **operator authorized preparation/execution of the small-box
pilot, including parallel-performance checks; bounded actual-IC coupled
startup job540015 and performance comparison job541239 completed, not a
calibrated/source-active science run**. See the
[performance record](snrt_performance_2026-09-12.md): RAMSES wall
1000.53 -> 947.75 seconds, actual CPU 28.56 -> 27.43 core hours, with
node/load caveats. CUDA stream-count tuning is stopped per operator direction.
A new128^3 small-box IC is now [generated](galaxy_new_ic_2026-09-11.md);
RAMSES hydro/element reader checks and the bounded cosmological CMB/dust/gas
coupling test have passed; see the active coupling record.
This is not full cosmological RT/feedback/dust admission. The operator has deferred pc-scale
zoom to long-term work; halo selection/refined zoom is not an active prerequisite.
The former volume-first execution matrix below
is superseded, not authorized for launch.
Baseline `66d58a81dd069041a912389ae43f066851675ea1`, pushed to
`kjhan0606/LagRamses/main`. Separate science campaign, not a reopening of
completed implementation gates. See [repair evidence](chimes_long_interval_repair_2026-09-11.md).

## Active operator decision: pc zoom deferred

Latest operator instruction: move pc-scale zoom to long-term work and proceed
with the current small-box campaign. Do not generate a halo-selection campaign
or demand few-pc attainment as a condition for this pilot. Retain the existing
12.5cMpc128^3 gas/DM IC; its coarse mass resolution is explicitly documented.

For future matched resolution comparisons, keep the same volume,
large-scale phases, physical prescriptions and diagnostic definitions. Report
four distinct quantities: actual dense-gas cell widths, gas cell masses,
stellar particle birth masses and high-resolution DM particle masses. Maximum
AMR level alone is not evidence that star-forming gas reaches that resolution.
Use physical pc at each epoch; for a comoving box,
dx_phys=L_com/[2^level*(1+z)]. AMR gas refinement does not increase the
initial DM particle sampling or add missing initial small-scale modes.

For this campaign, choose a computationally feasible resolution
from the existing IC and measured workload.
Do not call this intermediate milestone attainment of the ultimate target,
and do not invent a fixed few-pc runtime/cost before selecting the IC.
At finer scales, check how the existing SF/feedback closures apply without
automatically adding new physical models or an audit gate per refinement.

Local comparisons prioritize star formation, gas phases, feedback response
and dust evolution/conservation. A selected zoom is **not** a volume-complete
sample: the GAMA mass function and cosmic SFR density remain later population
validation targets, not likelihoods computed from a handful of chosen haloes.
Conditional galaxy relations can be compared with declared selection; no
claim of population calibration from an unrepresentative zoom sample.

Approval covers this direction and scoped preparation. Exact IC identity,
halo/sample choice, achieved resolution and revised resource/output plan
must be stated before science execution. Do not launch the former25cMpc
23--29-run ensemble, or retain its startup-only pilot matrix as the new plan.
No IC ownership transfer is implied. Existing physics applicability and retention rules
below remain in force; old resource numbers are proposals, not reservations.

## Active pilot preparation and parallel performance

Prelaunch inspection found a material applicability blocker:
`patch/lagRamses/read_hydro_params.f90` explicitly sets dust_ok=false when
`dust_mass_enabled` and `cosmo` are both true (line404 at inspection).
The subsequent diagnostic requires noncosmological periodic metal hydro.
This is a runtime admission restriction, not an IC-format failure. The existing
integrated evidence was noncosmological; cosmological RT/feedback/dust tuning
was not qualified by it. At that inspection no guard had been removed; no dust-off calculation
has been represented as full-model calibration. Subsequent bounded dust-off
readers are recorded below; no full-model cosmological science run was submitted.
Supporting the requested cosmological combination requires a scoped coupling
assessment (expansion/unit conversions, material/thermal conservation,
primordial initialization and source-table domain), then justified changes
and a bounded expanding-box execution. The guard alone does not establish
which of those mechanisms actually need code changes. The operator has now
approved that scoped repair and preapproved routine implementation/tests.
See [active coupling progress](cosmological_coupling_progress_2026-09-11.md).
Update2026-09-12: the scoped CMB/expansion repair and a uniform4^3 MPI2/OMP2
live test passed (job539357). The guard now admits only coadvected C/silicate
DL01/D03/CHIMES kind7/NENER=0 with an explicitly ledgered optically thin
`2.727/a` CMB bath; Fe/PAH, CR/SGS, MHD, drift, sublimation, SN shocks and
cosmological sinks/AGN remain excluded. No 10K cosmological dust floor or
instantaneous grain-energy projection. Fixed-group IR dilutes as a^-3 without
frequency redshift. This closes the bounded CMB/material coupling item, not
the full cosmological source/AMR campaign or galaxy calibration. The test used
manufactured Z=.02 dusty gas, not the pre-enriched galaxy IC below.
The operator subsequently approved [trace pre-enrichment at Z=1e-10](galaxy_pilot_preenrichment_2026-09-11.md).
The [actual-IC coupled startup pilot](galaxy_coupled_pilot_2026-09-12.md) now
uses that same128^3 input on four normal nodes, retaining the full spectral
and angular layout, with a two-step/20-minute bound and no scheduled dumps.
Submission is not evidence of completion or a performance result.
Use its new explicit passive IC files for the active pilot, preserving the
original zero-metal IC view. This is initial pre-enrichment, not runtime clamping.

Operator instruction (2026-09-11): start work and track CPU time during
parameter fine tuning. First execution venue is the grammar cluster's `debug`
partition, node `grammar-debug`; submit through `grammar`, not syntax Slurm.
Read-only recheck: 64 physical cores, 257647 MiB configured memory, no GPU
GRES; 16 cores and 64000 MiB allocated at inspection. GPFS is shared, with
approximately 95 TiB free (not a project quota or storage reservation).
Initial candidate layout: MPI8 x OMP2, 16 physical cores, at most128 GiB.
Recheck availability and size the selected IC before submission.

Native pre-enriched initial-step comparisons are complete:8x2 versus4x4 at
16 cores, same input/binary. Actual compute CPU432.315s versus415.706s,
RAMSES elapsed23.9596s versus25.3040s; no statistically established optimum.
Both are hydro-only costs, not full RT/CHIMES scaling measurements.

Full dense IR memory cannot inherit the reader's allocation:136 frequencies
x80 directions x128^3 cells x8 bytes is170GiB for ONE field alone. Current
native live/operator code has persistent, trial, transported and candidate
fields, at least680GiB globally before halos/capacity growth/other physics.
Consequently grammar-debug's252GiB is appropriate for small coupling tests,
not this dense128^3 full-IR run. Normal grammar nodes report515697MiB each;
the eventual128^3 resource request must size aggregate/per-rank peaks and
communication for multiple nodes. Do not weaken physical resolution or claim
full-IR performance from a dust-off reader. This is resource sizing, not a new
physics implementation gate or authorization for an unconstrained ensemble.

The registered `p4_pilot_zoom_agn_candidates.json` and
`p4_high_density_manifest.json` describe a 32^3 RT extraction from an external
Run0 snapshot, NOT self-consistent cosmological zoom ICs. They are not admitted
as calibration ICs. The external legacy restart and its binary/chemistry layout
are not interchangeable with the current coupled model. No compatible matched
three-resolution IC family has been identified in the inspected project records.
The operator subsequently authorized new IC generation. See the linked new-IC
record for the completed small-box input; the external extraction is not used.

Use existing phase/direct-cooling/SNRT timers and Slurm accounting; do not add
a profiling framework or modify physics solely for timing. Record separately:

- wall seconds for initialization, evolution and output;
- allocated core-hours = AllocCPUS * ElapsedRaw /3600;
- measured process CPU hours = Slurm TotalCPU /3600 (user + system over ranks
  and threads, using the compute step without double-counting parent records);
- CPU utilization = measured CPU seconds / allocated core-seconds, peak RSS,
  MPI/OpenMP layout, binding, co-tenant/load information and exit state;
- coarse steps, simulated time advanced and AMR leaf-cell updates, alongside
  chemistry/RT/feedback/dust phase costs where existing timers expose them.

`CPUTimeRAW` is allocated CPU time, not measured process CPU time. Timing the
launcher alone does not measure all MPI ranks. Missing accounting remains
missing, not zero. MPI spin-wait can consume process CPU: high utilization
does not by itself demonstrate useful work or good parallel efficiency.

After the selected IC's initial-step check, reuse a bounded source-active
window for MPI16xOMP1, MPI8xOMP2 and MPI4xOMP4 at16 cores. Fix binary, physical
parameters, initial state, stopping interval, output policy and binding policy.
This tests rank/thread layout, not strong scaling. Compare the best layout
with an8-core counterpart for strong scaling; speedup is T8/T16 and parallel
efficiency is (T8/T16)/2 for genuinely matched work. Repeat the selected layout
once if timing differences are comparable to run-to-run variability. Avoid a
full performance matrix for every physical parameter point. Recheck a late
source-active/refined state because startup timing does not price that regime.
For all tuning runs retain wall time AND allocated/measured CPU hours; compare
physical runtime changes with actual work counts rather than confusing changes
in SF/AMR/chemistry workload with changes in parallel efficiency.

Per-dump size, total storage, effective absolute namelist, all output clocks,
runtime/step cap and model applicability will be reported with the chosen IC
before launch. No new job has been submitted during this preparation record.

## Historical volume-first proposal — superseded execution design

Approve this design and **preparation/pilot only**: register compatible ICs
and observational data; at most3 bounded runs, total5,000 CPU-core-hours,
512GiB aggregate memory per job,1TiB peak storage. Each run is limited to
64 allocated physical cores,24 wall-hours and100 coarse steps. No GPU budget.
These are expenditure ceilings, not measured forecasts or an allocation on
the interactive host. Use allocated compute nodes.

After pilot measurements, submit **one** priced full-ensemble proposal for
approval; no approval/audit per run. Missing required capability is a scope
decision, not permission to invent another implementation bundle.

Q-GOAL: constrain observable predictions of the implemented RT/feedback/dust
effective model, with independent predictions separated from fitted data.
Q-LEAN: four existing controls, finite run ceiling, existing timers/output
tools; no generic gate framework or microscopic physics programme.

## Scientific method and frozen model

Target nearby ordinary galaxies first, nominal Mstar=10^9.5--10^11 Msun,
narrowed by actual resolution/completeness. COLIBRE uses mass/size observables
and emulators at multiple resolutions. Adopt its design/fit/direct-validation
strategy, **not its fitted parameters or particle-code resolution rules**.
[Chaikin et al.](https://arxiv.org/abs/2509.04067),
[COLIBRE](https://colibre.strw.leidenuniv.nl/).

Freeze evaluated C/silicate two-size coadvected grains, CHIMES157, kind7 RT,
v5 stellar source + empirical DTD/N100, CPU/OpenMP/NENER0; fixed reaction and
spectral banks, reduced-c, source energy/yields, DTD, grain optics and size
exchange/destruction laws. No joint CR/MHD/Fe/PAH/drift/MAD parameter search.

Preparation must explicitly resolve applicability, without bypassing guards:

- Combined evidence is **noncosmological, existing-sink** evidence, not
  cosmological BH seeding qualification. Do not insert arbitrary z=99 BHs
  or label an existing-sink test a cosmological AGN population calibration.
  A compatible physically justified seed path/IC is needed before AGN
  population fitting; otherwise report the limitation before the ensemble.
- Retain the explicitly normalized **Kroupa** effective-Ia reference for
  this proposal; repository Chabrier defaults remain unchanged. Match
  observational SPS/IMF inference with uncertainty, not a guessed offset.
- Verify mass/Z/age coverage and primordial gas treatment. Do not extrapolate
  stellar tables to Z=0 or absent ages merely to finish a run.
- No science IC is selected yet. Register actual paths/hashes, cosmology,
  source assets and effective namelist before launch. Do not take over
  unrelated simulations. Startup-only timing cannot establish late-time cost.

## Observables: fit versus independent prediction

| Role | Dataset/statistic | Required comparison |
| --- | --- | --- |
| Primary fit | GAMA DR4 low-z total stellar mass function | Published completeness/volume selection, IMF/SPS and aperture matching; .25dex bins |
| Primary fit | GAMA Lange mass--half-light-radius relation | Same band, projected light, morphology and surface-brightness cuts; not raw 3D half-mass radius |
| Dust fit | DustPedia dust/gas versus gas-phase O/H, training subset | One metallicity calibration; consistent dust emissivity and HI/H2/He conventions |
| Independent holdout | Remaining DustPedia objects; xGASS HI fractions vs mass and sSFR | Galaxy-identity split; upper limits/censoring and survey weights |
| Predictions only | SFR, quenched fraction, gas/stellar metallicity, z=1/2 evolution, supported AGN trends | Publish discrepancies; do not silently add these to the fit |

Sources: [GAMA DR4](https://arxiv.org/abs/2203.08539),
[GAMA releases](https://www.gama-survey.org/),
[Lange et al.](https://arxiv.org/abs/1411.6355),
[DustPedia De Vis et al.](https://arxiv.org/abs/1901.09040),
[CDS catalogue](https://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J%2FA%2BA%2F623%2FA5),
[xGASS](https://xgass.icrar.org/).
DustPedia's main site did not load in this check; use paper/CDS and verify
the required dust/gas fields rather than assuming one table has everything.
No data download/checksum is claimed completed in this planning turn.

Preassign DustPedia70/30 train/holdout, stratified by mass/metallicity, seed
20260911. Account for shared GAMA objects/covariance. Apply selection to
synthetic observables. If matched photometry is unavailable, report that
before fitting sizes; an intrinsic radius is not an equivalent substitute.

## Four existing parameters

These are proposed search bounds, **not literature-calibrated priors**.

| Namelist control | Centre; range | Meaning |
| --- | --- | --- |
| `eps_star` | .02; .005--.05 | SF efficiency, fixed SF prescription/threshold |
| `eAGN_T` | .15; .05--.30 | Thermal AGN coupling; `eAGN_K=1` fixed |
| `dust_sticking` | .30; .10--1 | Growth sticking probability |
| `dust_condensation(1:3)` | common q=.10; .01--.30 | Tied wind/AGB/SNII condensation, one freedom |

Sample log-uniformly. Controls exist in `read_hydro_params.f90` and frontend
definitions. SN/SNIa energy, IMF, DTD and yield normalization stay fixed;
this is not a claim to tune every SN feedback parameter. If the sample has
no thermal-mode AGN response, freeze `eAGN_T` and label it unidentifiable,
not successfully fitted. BH seeding is not an extra nuisance parameter.

## Historical resolution and ensemble matrix — not the active run plan

Proposed25 comoving Mpc periodic box, z_init99 to0, matched large-scale
phases. Proposed fixed cosmology Omega_m=.315, Omega_b=.049, h=.674,
sigma8=.811, n_s=.965; actual IC metadata must agree. Numbers are derived
from this proposal, not selected ICs or copied COLIBRE resolution values.

| Resolution | Base grid / DM count | Initial gas cell mass | DM mass | Proposed max AMR level / z=0 cell width |
| --- | --- | --- | --- | --- |
| Coarse | 128^3 / 128^3 | 4.60e7 Msun | 2.50e8 Msun | 14 /1.53kpc |
| Intermediate | 256^3 / 256^3 | 5.75e6 Msun | 3.12e7 Msun | 15 /.763kpc |
| Fine | 512^3 / 512^3 | 7.19e5 Msun | 3.90e6 Msun | 16 /.381kpc |

Masses use rho_crit=2.775e11*h^2 Msun/Mpc^3, component mass/N^3. Keep
dimensionlessly matched refinement rules; record actual leaf/star masses.
Require >=100 star particles for mass, >=1000 for size, Re>=4 cell widths,
and >=30 eligible objects/bin. Common resolved range may be narrower than
the nominal range. Insufficient bins are insufficient evidence, not permission
to loosen cuts. A25Mpc box cannot qualify rare massive AGN or JWST/LRD counts.

After approval of measured full-ensemble cost:

- Coarse: centre +8 axial perturbations +8 space-filling combinations
  +3 reserved direct candidate checks = **20 full runs**.
- Same chosen parameters at intermediate/fine: **2 runs**, strong convergence.
- Second coarse IC phase: **1 run** for variance sensitivity; two phases
  do not estimate a full cosmic-variance covariance.
- Optional at most**6 intermediate recalibrations** if needed, fine held out.
- Total23 primary,29 maximum; no20-point design at every resolution and no
  automatic enlargement of box/sample or repeated search until success.

Three resolutions measure a trend, not a validated parameter-interpolation
law. Interpolation would require fitted endpoints plus an independently
held-out intermediate test; that additional study is not in this allocation.
Strong convergence (fixed parameters) and weak convergence (refit) stay
separate throughout.

## Historical budget estimates and retained launch-safety rules

Pilot covers startup and at least one cold/source-active state if physically
compatible ICs/checkpoints are available. Do not manufacture cross-resolution
restarts by rewriting fields. If suitable states are absent, report before
spending the allocation. Pilot approval does not approve an unsupported
physical combination or an arbitrary new source model.

| Pilot slot | Proposed input and purpose | Limit |
| --- | --- | --- |
| P1 | 12.5cMpc/64^3, compatible cosmological startup; verify units/source coverage and initial memory | 100 steps,64 cores,24h |
| P2 | 25cMpc/128^3, same mass scale; measure volume/communication and memory scaling | 100 steps,64 cores,24h |
| P3 | A registered source-active checkpoint in the same admitted model, no more than128^3 base grid | 100 steps,64 cores,24h |

P3 is conditional on a legitimate available checkpoint. If absent, do not
substitute another empty early-time run and call late-time cost measured.
All three share the stated total resource ceiling; none is a new fit sample.

The64-cell source-active regression took920.690s/4steps, cooling99.8%.
This warns against naive box extrapolation; no cosmological wall time or GPU
speedup has been measured. Use existing timers for cell-steps/s, leaf growth,
chemistry load and peak memory; price the full ensemble with a forecast range.

Planning allowance16KiB/leaf for hydro/chemistry/RT alone implies
32/256/2048GiB at128^3/256^3/512^3, **before** AMR, solver staging, particles
or MPI duplication. Real memory can be several times larger. Fine requires
a distributed allocation; the observed64-core host is not that allocation.

Pilot: <=128GiB/full dump, <=2 dumps/run,6 total,1TiB peak including restart
copy and products. Stop/re-scope if measurements exceed this ceiling.
Full-ensemble storage proposal32TiB, **not authorized**; CPU budget and even
the feasible volume/resolution await pilot measurements. Current GPFS free
space~95TiB is shared free space, not our project reservation.

Before every launch report absolute effective namelist, run purpose,
noutput/aout/tout/foutput/fbackup, stopping condition, dump/total size and
free space. Full-run proposed schedule: z=2,1,0 only, noutput3,
aout=1/3,1/2,1; suppress periodic dumps and allow at most one rolling restart
if needed. Pilot schedule must match its actual starting state, not blindly
reuse that schedule. Do not use the69Myr stress-test timestep as a production
accuracy prescription; use the physical/numerical timestep controls.

Delete evaluated **pilot/test** raw after needed restart evaluation, retaining
inputs, hashes, metrics and cleanup manifest. Full scientific snapshots are
not test trash: retain final/restart products within budget until validated;
obtain explicit approval for destructive scientific-archive cleanup.

## Analysis and completion criteria

Use observational, sampling and measurement uncertainty, covariance where
available, censored gas likelihoods and documented systematics. Missing
covariance means an explicitly approximate likelihood/sensitivity analysis,
not a precisely constrained posterior inferred from invented independence.

Inspect direct responses first. Use a small GP only if withheld simulation
errors are <half observational error and roughly .05dex in populated fitted
bins. Otherwise keep direct comparisons; do not build a larger framework.

Proposed science targets, **not new engineering acceptance gates**:
mass-function/size residual <=.15dex and dust/gas <=.20dex in >=80% eligible
bins, reporting uncertainties/outliers. Strong-convergence differences
<=.10dex for mass/size, <=.15dex for dust/gas in the common resolved range.
Publish holdout failures; do not weaken conservation, move holdouts into
training or restart completed implementation gates to obtain a favorable fit.

Deliver one report with parameter degeneracies, fit/holdout figures,
strong/weak convergence, actual cost and applicability limits. Stop at the
approved run ceiling; insufficient counts or failed fits are reported results.
No publication-ready galaxy prediction claim without matched selection,
resolved scales and independent observational support.

## Historical approval request — replaced by the active decision above

Requested now: design and preparation/pilot only,5,000 core-hours/1TiB.
No simulations launched. The priced23--29-run ensemble and any material
physics/IC ownership change require separate approval. Microscopic physical
calculations stay in the long-term backlog.
