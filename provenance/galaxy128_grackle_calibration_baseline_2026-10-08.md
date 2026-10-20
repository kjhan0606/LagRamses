# Small-volume Grackle calibration baseline — 2026-10-08

Operator authorized return to parameter tuning after closing Grackle
qualification on the owner-verified numerical policy. No final Astra audit
or completion email is required. This does not certify galaxy predictions.

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001`,
origin `git@github.com:kjhan0606/LagRamses.git`.
Read-only inspection confirmed `/gpfs/kjhan/LRD_JWST` is another worktree
of the same repository. Its former tuning run and pinned source assets are
used as references/inputs, not modified or assumed to be current code.

## Submission and effective inputs

**Job 415961**: one H100, MPI4 x OpenMP2, 8 CPUs, 64 GiB, 24-hour limit;
at most 192 allocated core-hours and 24 GPU-hours for this segment.
CUDA-enabled Grackle build: `GRACKLE_NEQ=1 HDF5=1 USE_CUDA=1`,
`SNRT=CHIMES=DUST_LIVE=USE_FFTW=0`, `NVECTOR=32`, `NVAR=27`.
Poisson GPU enabled; hydro and Grackle chemistry remain CPU calculations.
The preceding Grackle qualification used CPU Poisson; startup verification
of this new execution profile is therefore performed before science evolution.

Deployed input bundle:
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-inputs.wsfIJi`.
New run:
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-baseline-415961`.
Effective inputs are `preflight/effective.nml` and `science/effective.nml`
under that run. Build, actual executable hashes, source/input manifest and
logs are retained there. The job validates the 22-file IC manifest and checks
frozen source/input hashes before and after compilation. A failed startup
stops the script; it does not launch the baseline anyway.

## Model and scientific intent

Existing 12.5 cMpc / 128^3 gas+DM IC, trace Z=1e-10 pre-enrichment,
levels 7--13 and original mass/Jeans refinement (8/4), auto capacity growth.
Cosmology is read from the original IC with omega_b consistency enforced.
Density-based SF: `eps_star=.02`, `t_star=0`, `n_star=1`, `m_star=.1`.
Courant=0.125 and nsubcycle=1 follow the owner-adopted numerical profile;
its fixed-level early-epoch qualification is not universal AMR/late-time
accuracy proof.

Nine Grackle primordial species + equilibrium metal cooling, HM2012 UV
active with z_reion=8.5, fixed gamma=5/3. RT, dust, CR, sinks and AGN are off.
The earlier CHIMES/dust baseline was UV-off; this explicit background choice
is now fixed across prospective parameter comparisons, not hidden as a
mere performance change. No claim of active ISRF, shielding or radiation
pressure is made.

Retain the previous effective user-selected channel-resolved source:
Kroupa IMF=1, 0.08--600 Msun, version-5 mixed-lowmass/terminal proxy history,
effective SSP plus empirical Maoz/N100 SNIa, wind/AGB/SNII/PISN enabled.
No replacement by the narrower KL16/LC18 diagnostic population and no
implicit conversion of its DTD normalization to Chabrier. Known source
approximations remain declared; source admission guards are not bypassed.

This is ONE reference run, not a sweep or fitted result. First obtain
resolved galaxies, stellar mass/SFR/metallicity and HI statistics, sample
counts, actual cell/star masses and wall/CPU costs. Then price paired eps_star
perturbations and compare observables with matched IMF/SPS/aperture/selection.
The existing calibration plan supplies finite-volume/resolution limits.
Inactive dust and AGN parameters cannot be tuned here; colors require a
separate declared measurement-compatible SPS/attenuation treatment.

## Output and stop policy

Startup: nstepmax=1, future aout=1.1, tout=1e100, periodic clocks=1e6;
must create no raw snapshot. No full code-test campaign is added.
Science: noutput=3, aout=1/3,1/2,1 (z=2,1,0), tout=3*1e100,
match_aout=true, foutput=fbackup=1e6, nstepmax=100000.
Automatic walltime dump is disabled. Slurm sends USR1 one hour before
the 24h limit, requesting one clean checkpoint-stop through jobcontrol.
A partial checkpoint/step ceiling is not z=0 completion. Preserve scientific
outputs and the necessary restart; no automatic ensemble or restart jobs.

Base measured Grackle dump: ~0.72 GiB. Adaptive/star-active output is unknown;
planning envelope is <=128 GiB/dump, <=512 GiB for three snapshots plus one
checkpoint. Latest /gpfs free space ~260 TiB, not a private reservation.
CPU cost must use Slurm compute-step accounting and existing timers, not
launcher CPU alone. Late-time completion time is not inferred from startup.

The prior temporary /gpfs write restriction was lifted and directory creation
verified before deployment. Submission is not evidence that initialization,
source-active evolution or observation matching has completed.

## Startup outcome and corrected continuation

Job 415961 completed compilation and its compute preflight (step 415961.0
COMPLETED, exit 0, 3m32s). The source/IC checks passed. Grackle and
HIGH_MASS_HISTORY_READY were present; GPU MG converged on levels 7 and 8,
and logged NaN checks were zero. The run advanced from a=0.01000 to
approximately 0.01010 and reported `Run completed`.

The batch nevertheless exited 1 because its `output_*` directory check
mistook `output_00000/stars_00000.out00001` through `00004` for a full
snapshot. The cuRamses star-formation implementation opens these birth
record files whenever sf_birth_properties is enabled, even without births.
These were empty files, not AMR/hydro/particle dump data. This was a driver
classification error, not a chemistry or physics failure.

Submitted continuation job 415971 with dependency afterany:415961 using
`.galaxy128-grackle-inputs.wsfIJi/resume-baseline.sbatch`. It requires the
completed preflight and all existing physics/error checks, permits only
regular stars_*.out* diagnostic files in preflight output directories,
rechecks frozen source/input and executable/asset hashes, and refuses to
overwrite a science run.log. It reuses the exact binary and effective
science input; no rebuild, repeated preflight, changed physics or bypass
of source-admission checks. H100 x1 / MPI4 / OMP2 / 64GiB / 24h and the
three-snapshot plus clean checkpoint-stop policy are unchanged.

Initial preflight cost is not a late-time forecast: its routine timer total
was 164.012s, refinement 41.051s average, cooling 1.678s average. GPU MG
base and AMR timers overlap the enclosing poisson timer and must not be
summed as independent costs. Adaptive grid capacity grew from 150016 to
643904 per-rank slots; the printed `mem=****%` is formatting overflow, not
by itself OOM evidence. Actual step MaxRSS was 3575276K; batch/build peak
was 12151744K. No science completion or calibration success is claimed.

Continuation log confirms all seven asset hashes matched and
`BASELINE_START job=415971 time=2026-10-08T01:53:09+09:00 MPI=4 OMP=2`.
Science run.log was created and began reading the pinned GRAFIC IC fields.
This proves baseline launch, not yet successful first science evolution.

At elapsed 2m06s, job 415971 was RUNNING on syn08. The science log reports
GPU MG level-7 convergence (residual 1.164e-5), zero logged NaN counts,
step-0 coarse interval 39.43s, direct cooling 1.479s and 1.1GB reported
per-rank memory before capacity growth. These are startup timings, not
late-time cost or a CPU/GPU speedup measurement. The next step is live
after growth to 643904 grid slots; no completed step-1 yet at this check.

Completed accuracy-test raw HDF5 snapshots were removed (1.44GiB total)
only after owner recheck and exact-file SHA256 recording. Compact results,
metadata and logs remain; see
`../.grackle-neq-accuracy-415811/raw_cleanup_2026-10-08.md`.
No calibration scientific snapshot or active restart was removed.

## AMR thermal-state failure: baseline stopped

At elapsed 4m14s job 415971 was live and had reached coarse step 1,
a approximately 0.01010. The newly populated level 8 has 2,097,152 grids
(16,777,216 cells). Its entry diagnostics count 2,313,593 negative-internal-
energy cells across four ranks. More importantly, the subsequent
godunov+sync diagnostics still count 2,272,775 negatives, with minimum
internal-energy density approximately -4.656e-10. These are finite values
and therefore were not detected by the NaN-only startup gate. The completed
nstepmax=1 preflight had already printed negative level-entry diagnostics;
it did not exercise the following fine-level hydro/cooling update.

Cancelled own baseline job 415971 at 01:58:21 KST after confirming the
post-hydro failure. Do not accept this as a valid thermal trajectory or use
its state to calibrate parameters. Preserve its logs/inputs/build; no valid
science snapshot or restart is asserted. The fixed-level qualification is
not contradicted within its scope, but does not cover this AMR failure.

Submitted bounded isolation job 415972, H100 x1 / MPI4 / OMP2 / 64GiB /
30min, reusing the exact existing binary with no rebuild. Two original-IC,
two-coarse-step controls in `.galaxy128-amr-isolation-20261008`:

- `interp_internal`: cooling unchanged; existing `interpol_var=1` internal-
  energy AMR interpolation instead of default total-energy interpolation.
- `no_cooling`: original interpolation unchanged; cooling disabled to test
  whether prior-step cooling is necessary for the negative entry state.

Both retain original refinement/Jeans thresholds, source population,
Courant factor, gravity and all inactive-module choices. Each has an 800s
timeout and future-only output clocks; only small stellar birth records are
expected. No full snapshot is requested, no tolerance/floor is changed.
Existing internal-energy interpolation is a diagnostic control, not an
approved production fix: its conservation implications must be evaluated
before adopting it. Source, prolongation and gravity source-term ordering
remain candidate causes pending measured results. No Grackle-only root
cause is inferred from a negative pre-cooling state.

First isolation submission 415972 did not evolve: interpol_var was
incorrectly placed in HYDRO_PARAMS rather than REFINE_PARAMS, and the
no-cooling control retained cooling_method='grackle_neq', whose admission
correctly requires cooling enabled. Both driver mistakes were corrected
without changing source code or admission guards. Rejected logs are retained
as rejected-415972.log per case; input hashes from that attempt remain.
The no-cooling control now explicitly selects cooling_method='original'
with cooling=false (not a different active cooling solver). This control
does not exercise Grackle species initialization and is used only to
isolate the hydro energy producer, not as a chemical accuracy comparison.
Corrected isolation submission is job 415975.

At elapsed 4m36s, 415975 was RUNNING on syn08 in the interp_internal
control. Level-8 MG converged and the case reached coarse step 1 at
a approximately 0.01010 without a logged negative level-entry diagnostic,
unlike the original baseline at the same boundary. No completed fine-level
hydro/cooling outcome is claimed yet. Step accounting reported MaxRSS
3572472K and nonzero accumulated CPU. The no-cooling control has not started.

Source inspection: CPU synchro_hydro_fine explicitly removes old kinetic
energy, updates momentum, and adds new kinetic energy, preserving internal
energy algebraically for that isolated kick. This does not prove the whole
gravity/AMR ordering correct. In contrast, interpol_var=1 interpolates
internal energy then adds each child's kinetic energy, without a final
constraint matching the parent total-energy average. Therefore a positive
internal-energy control does not by itself justify production adoption.
Any actual repair must retain declared conservation and must not simply
replace negative energy by a temperature floor.

### Isolation evidence and conservative repair candidate

415975 interp_internal completed two steps at a approximately 0.01020,
exit 0, with no negative entry/post-hydro/post-cooling diagnostic. Its routine
timer total was 326.142s, Godunov average 121.310s (37.2%), cooling 16.212s
(5.0%). This is a startup/adaptive cost, not a calibrated late-time forecast.
The no_cooling control subsequently reproduced negative level-8 entry
states on all four ranks (minimum approximately -4.688e-10). Hence active
Grackle cooling is not necessary for this failure. The control is still
finishing; these observations isolate dependence, not every producer detail.

Added a scoped conservative prolongation limiter in
patch/cuRamses/interpol_hydro.kjhan.f90 for grackle_neq + interpol_var=0.
It leaves admissible stencils unchanged. For an inadmissible child of an
admissible parent, a single common theta contracts every child conserved
field about its parent using the existing symmetric offsets. Child averages
therefore retain parent density, momentum, total energy and passive-field
inventories to roundoff. The thermal-admissible set is convex; bisection
selects a valid common factor. This adds no energy floor and does not turn
internal-energy interpolation into the default. An inadmissible parent is
rejected instead of repaired silently. Scope excludes MHD/NENER Grackle
models; this is not a claimed general MHD positivity method.

Submitted job 415981 afterany:415975 with the same H100 x1 / MPI4 / OMP2 /
64GiB layout, bounded 30min. A separate build directory preserves the
original binary; inherited object mtimes are preserved so Make recompiles
the changed source. The actual routine is tested against a cold stencil
that fails unlimited interpolation, conserved child averages, positive
thermal energy, and unchanged warm/uniform stencils. Only after these
checks does the two-step original-IC run execute. Input/source/executable
hashes and logs are retained in `.galaxy128-amr-repair-20261008`.
No production baseline is relaunched before measured repair verification.
Cancelled pending 415981 before execution and resubmitted as 415983 after
correcting object copying to preserve mtimes. Slurm captures the submitted
script, so editing its on-disk copy alone would not correct a queued job.
415983 is the authoritative repair job; no duplicate calculation is running.

Stopped 415975 after the no-cooling control established the same negative
entry, post-Godunov and unchanged post-cooling states. It was beginning a
further capacity expansion 643904 -> 4970304 slots, with observed MaxRSS
22215580K for one task; continuing an already-invalid thermal trajectory
would add resource cost without changing this isolation result. The
no-cooling control is therefore a reproduced failure, not a completed
two-step PASS. The internal-energy case did complete two steps. Pending
415983 uses afterany and will execute after this cancellation; no test result
is fabricated from cancellation. Future source mutations are on hold while
the repair build records its input hash.

415983 compiled/relinked the changed prolongation source successfully, but
the standalone actual-routine checker failed to link because its explicit
library list omitted /home/kjhan/local/lib64 (the main Makefile already
uses that location). No check or IC evolution ran in this attempt. Corrected
the checker link path and submitted 415995 reusing the separate existing
repair build; it refuses to overwrite run.log and retains per-job build and
input-hash logs. This is a driver/link error, not a verdict on conservation
or physics. The original repair build.log and repair-inputs.sha256 remain
historical evidence; no repeated full compilation is needed.
415995 actual-routine checker passed all three stencils: cold unlimited
failure repaired with positive thermal energy, conserved per-field child
averages (measured maximum error zero), and warm/uniform results unchanged.
Original-IC srun subsequently started. Integrated AMR repair verification
is still in progress; these local checks alone are not production evidence.

415995 reached level-8 entry, but reported residual negatives of order
5e-24--8e-24 in 69,917 cells in total. The original approximately 4.7e-10
failure was removed, but an exact-positivity claim would still be false.
The original limiter bisected to thermal>0 with essentially zero arithmetic
headroom. Componentwise kinetic subtraction/gravity synchronization can
change this cancellation at a few ULPs. Stopped 415995 before treating that
state as a pass or allowing further invalid evolution.

Adjusted the admissibility predicate to require thermal energy greater than
128 ULPs of max(abs(total energy),abs(kinetic energy)). This selects a
slightly smaller common slope factor: no total energy is added, no physical
temperature floor is introduced, and child averages remain conservative.
An unresolvable parent still fails closed rather than receiving artificial
energy. New job 415996, same resources and two-step original-IC bound,
uses `.galaxy128-amr-repair-roundoff-20261008`; earlier artifacts remain
unchanged. The same actual-routine checker now also includes a cold parent
with thermal energy 2e-13 relative to kinetic energy 0.5 to exercise
near-cancellation. The roundoff explanation is a numerical inference to
be tested by integrated results, not a proven complete root-cause diagnosis.
415996 checker passed all four cases with measured maximum conserved-field
average error zero. The updated binary compiled and the original-IC
integration started; at elapsed 6s the job was RUNNING on syn08. Full
integration remains pending.

### Verified bounded repair and baseline restart

415996 COMPLETED, exit 0, batch elapsed 9m27s, compute step 9m22s,
MaxRSS 3622540K. Four actual-routine stencil checks passed with zero measured
parent-average error. Original IC reached a approximately 0.01020 at step 2
and `Run completed`; no negative entry/post-Godunov/post-cooling diagnostic,
nonzero NaN count, ERROR or MG nonconvergence. The runner emitted
BOUNDED_REPAIR_PASS. This closes the reproduced early AMR thermal-state
failure within this measured scope, not arbitrary AMR/late-time accuracy.

Existing timer total 519.520s: refine average 156.339s (30.1%), Godunov
115.841s (22.3%), cooling 20.837s (4.0%). This conservative repair costs
more than the nonconservative internal-energy diagnostic (326.142s total),
but they do not evolve the same numerical energy representation. Not a
matched CPU/GPU speedup claim. Runtime/load and finer-level costs remain
uncertain; no z=0 ETA is inferred. Refinement/limiter efficiency is a visible
cost, not grounds to silently remove conservation or numerical headroom.

Updated mkrun.py's Grackle choice to explicitly select interpol_var=0,
interpol_type=1 before the optional editor; other modes retain their old
defaults. The namelist generator describes energy representations correctly
and warns, rather than silently rewriting an explicit nonconservative
choice. Python syntax and targeted diff whitespace checks passed. README
records the measured scope. No new Fortran namelist key was added.

Submitted baseline 416003 using the exact repair-qualified binary, new run
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-conservative-416003`,
effective namelist `effective.nml` there. H100 x1 / MPI4 / OMP2 / 128GiB /
24h. Host-memory headroom increased from 64GiB in view of observed adaptive
capacity growth; physical refinement criteria, Courant, source population
and module choices are unchanged. Existing three full dumps at z=2,1,0
plus at most one clean walltime checkpoint retained (128GiB/dump planning
envelope, <=512GiB; /gpfs approximately 260TiB available at launch).
This is the same one-baseline calibration plan, not a new parameter sweep.

Runner verifies repair PASS, four stencil PASS records, executable/assets,
and the original frozen source/input manifest except exactly the measured
prolongation source change. All checks passed. Log confirms BASELINE_START
job=416003 at 2026-10-08T02:36:35+09:00; science log initialized Grackle and
began AMR construction. No repeated build/preflight and no invalid old
checkpoint is used. This establishes launch, not low-redshift completion or
an observational fit. Dust/RT-on readiness is not inferred from this run.

Baseline live check: 416003 RUNNING on syn08 at elapsed 10m26s; reached
coarse step 2, printed a=0.01020 (z approximately 97). No negative thermal
or nonzero NaN diagnostic, ERROR, or MG failure is recorded. Last completed
coarse interval 343.52s and reported memory 3.4GB per rank. Step-2 mcons
is 4.28e-15; econs=-2.09e-5 is only a native evolution diagnostic in a
cooling/expanding/gravitating model, not proof of a closed energy budget.
The new baseline has traversed the formerly failing early AMR state; it is
not at the first z=2 science snapshot and no galaxy catalogue is available.

Stable-cost sample: 416003 RUNNING at elapsed 32m49s, step 8,
a approximately 0.01084. Last five completed coarse intervals are
239.36, 238.53, 238.82, 239.48, 240.82 seconds, followed by 240.94s;
reported per-rank memory remains 3.4GB with no additional capacity growth.
No negative thermal/NaN/MG error is recorded. Slurm step sstat sampled
AveCPU=00:13:50, AveRSS=3608060K, MaxRSS=3611740K. Four-task average CPU
time implies approximately 55m20s accumulated task CPU (0.922 core-hours)
at that sample; this rounded live accounting is not the final charged
allocation or a backend speedup measurement. Retain final sacct TotalCPU
and allocation elapsed separately when the segment terminates.

The existing four-case executed checker is now retained in the repository
as `simulation/snrt/tests/test_grackle_conservative_prolongation.f90`.
cmp confirmed byte identity with the job-415996 input; SHA256 is
1757c02462eb382d4979799ccc6dc1d7f403f51fd39b070ecca0e7373a601217.
README records its actual-routine linking requirements. This preserves the
already-executed regression, not a new benchmark, job or acceptance gate.
The running science executable is unchanged.

Extended live sample: job 416003 remains RUNNING on syn08 at elapsed
2h28m15s, coarse step 36, printed a=0.01435. Recent completed coarse
intervals remain approximately 245--247 seconds (latest 246.33s).
Slurm compute-step sstat reports four tasks, AveRSS=3612388K and
MaxRSS=3657636K (approximately 3.49GiB maximum per task). The log's
`mem=****%` formatting is not an OOM diagnosis; measured RSS remains
close to the earlier sample. A full current-log scan found no ERROR,
negative-state DIAG, fine-MG failure, or nonzero NaN_CHK uold/f/d0 count.
No restart, source/binary change, additional benchmark or parameter sweep
was performed during this wait. This is early-trajectory evidence only,
not low-redshift completion or RT/dust-on qualification.

Five-hour live sample: job 416003 is RUNNING on syn08 at elapsed
5h00m54s, coarse step 73, printed a=0.02079. Latest completed coarse
interval is 248.63 seconds. Four-task sstat reports AveRSS=3593009K,
MaxRSS=3666528K (approximately 3.50GiB maximum per task). The latest
12000 log lines contain no ERROR/DIAG/negative-state/fine-MG failure
match or nonzero NaN_CHK uold/f/d0 count. This bounded log inspection
does not replace final whole-run evaluation. The reference remains
early in its cosmological trajectory; neither low-redshift calibration
nor full RT/dust readiness is established. Inputs and executable remain
unchanged; no additional job or restart was submitted.

AMR capacity transition: at coarse step 164 (a=0.04949), job 416003
expanded grid capacity from 643904 to 4970304 per rank. At elapsed
11h25m41s the job remained RUNNING; sstat's preceding sample reported
MaxRSS=22280160K (approximately 21.25GiB per task). The GRID_GROW message
is emitted only by rank 1, after that rank's array replacements; missing
messages from the other ranks do not establish a stalled allocation.
The collective target is the maximum required grid count across ranks,
with at least 25 percent growth, not an unconditional 7.7-fold multiplier.
Actual post-transition mesh counts are not yet logged.

Read-only eu-stack samples on syn08 locate PID 3109821 in
hilbert3d -> cmp_ordering -> cmp_cpumap -> make_grid_fine -> refine_fine,
and PID 3109824 in interpol_hydro -> make_grid_fine -> refine_fine.
Their OpenMP worker threads are waiting. These samples demonstrate
serial refinement work on two ranks rather than MPI waiting at those
instants; they do not prove completion, rule out every subsequent stall,
or quantify the dominant cost over the entire transition. No running
binary, input, job allocation or checkpoint was changed.

At elapsed 11h44m53s, the log advanced beyond GRID_GROW: post_maketree,
post_kinfb, post_agnfb, density and first synchronization NaN counts were
all zero. Level-7 GPU MG used the expanded capacity (ncell=39762433,
ngm=4970304) and converged in five iterations to 1.134e-5; post_force
NaN counts were zero. This establishes that the observed refinement
phase was not permanently stuck. Step 165 and the post-transition
physical mesh counts were still pending; ncell here is allocated extent,
not a measured leaf-cell count. MaxRSS at that sample was 22495048K
(approximately 21.45GiB per task). No source or runtime contract changed.

At elapsed 11h49m41s the synchronized diagnostic reached step 165,
a=0.04987, dt=0.02929. The actual mesh now contains 16777216 level-9
grids (4194304 per rank), corresponding to a fully refined 512-cubed
cell level, despite the unchanged 128-cubed initial conditions.
Level-9 GPU MG converged in four iterations to 2.075e-5. Mass diagnostic
mcons=3.02e-14; native econs changed from +4.04e-5 to -2.28e-2.
The source in patch/cuRamses/update_time.f90 defines econs from kinetic
and gravitational potential energy with an expansion integral; it does
not include thermal energy or a cooling/feedback exchange ledger.
The jump coincides with the new gravitational mesh and must not be
described as a measured 2.28-percent thermal interpolation loss or
dismissed as harmless without further evidence. Post-transition timing
and subsequent diagnostic trend remain pending. This is not a science
calibration PASS or a full energy-accounting qualification.

The transition interval completed at elapsed 12h12m39s: 3197.36 seconds
(53.29 minutes), versus 253.45 seconds for the preceding interval.
This 12.62-fold increase includes refinement creation and is NOT the
steady-state level-9 cost. Direct cooling was 148.671 seconds versus
18.499 seconds previously (8.04-fold increase, approximately matching
the eightfold cell growth), only 4.65 percent of this transition interval.
SNRT timers remained zero; particle sub-timers totaled 8.107 seconds.
The remaining time cannot all be assigned to hydro from these partial
timers. Live stacks separately observed serial refinement and OpenMP
Godunov work, without an integrated time split. The next post_maketree
NaN counts were zero. Source, inputs and executable remain unchanged;
the subsequent steady-state cost and econs trend are still pending.

First unchanged-mesh interval: at elapsed 12h43m36s the step-166 timer
reported 1834.48 seconds (30.57 minutes), 7.24 times the pre-transition
253.45 seconds for eight times as many fine cells. This is one interval,
not a converged steady-state cost distribution or an OMP scaling test.
Direct cooling was 148.205 seconds (8.08 percent of the interval), with
SNRT timers zero. Printed timestep approximately halved after refinement,
so wall-clock cost per unit evolved time also rises; the original early
mesh cost cannot forecast z=0 completion. Native econs at step 166 was
-2.26e-2 versus -2.28e-2 at 165, and mcons=-4.32e-14. One subsequent
sample does not resolve the energy jump. The next post_maketree and
feedback NaN checks were zero. No physics thresholds, input contract,
allocation or binary were changed, and no continuation was submitted.

The second unchanged-mesh interval (step 167) reported 1834.49 seconds,
with direct cooling 148.564 seconds (8.10 percent), SNRT timers zero,
and particle sub-timers totaling 8.098 seconds. The two intervals agree
within 0.01 second, supporting a repeatable approximately 30.57-minute
cost for this mesh and runtime configuration, not a resolution-independent
cost or a parallel scaling claim. Native econs=-2.23e-2 and
mcons=-3.57e-16; the energy jump remains unresolved despite the modest
subsequent recovery. Post-interval tree/feedback/density NaN checks were
zero. The active run and its physics contract remain unchanged.

At elapsed 23h24m16s on 2026-10-09, step 187 completed in 1830.78
seconds, with direct cooling 149.695 seconds (8.18 percent), particle
sub-timers 8.038 seconds and SNRT timers zero. Native econs=-1.81e-2
and mcons=1.14e-14; the earlier gravitational-mesh energy jump remains
unresolved. The launcher delivered its one-hour-before-limit stop
request, and the code acknowledged `stop requested at step 188`.

By elapsed 23h27m42s the log reported `HDF5 output done.` and a valid
Poisson warm-start checkpoint. The retained checkpoint is
`.galaxy128-grackle-conservative-416003/output_00001/data_00001.h5`
(40007155072 bytes, approximately 37.26 GiB). Its COMPLETE marker
contains 00001, and info_00001.txt records nstep_coarse=188,
aexp=0.0540802071628654 (z approximately 17.49), ncpu=4 and
ngridmax=4970304. POISSON_PHI_VALID records the same step and epoch.
These are output-completion metadata, not a successful restart test.
At that observation the Slurm job was still RUNNING and `Run completed`
had not been seen. Subsequent tree/feedback diagnostics were zero.
Normal exit and restart qualification remain pending. No continuation
was submitted, no scientific output was deleted, and no source or
physical refinement threshold was changed. This epoch is not a low-z
galaxy-observable calibration result.

Terminal confirmation: Slurm job 416003 and its batch/compute steps
COMPLETED with ExitCode 0:0 after 23:33:03. The native log contains
`Run completed` at step 188 and the launcher records BASELINE_SEGMENT_EXIT
at 2026-10-09T02:09:37+09:00. Final mcons=-3.85e-14 and econs=-1.79e-2;
the earlier energy jump is still unresolved. This proves a clean segment
exit with a retained checkpoint, not low-z completion or restart fidelity.

The terminal whole-segment timer reports rank-average Godunov time
45022.386 seconds (53.1 percent of the average TOTAL), cooling 6325.743
seconds (7.5 percent), and reverse hydro ghost-zone time 4911.284 seconds
(5.8 percent). Reverse ghost-zone rank times range from 467.236 to
11793.878 seconds (std/average 0.852); Godunov ranges from 38320.226
to 49316.264 seconds (std/average 0.090). These are cumulative timer
categories, not a new strong-scaling measurement, and communication
waiting may reflect upstream compute imbalance. They identify hydro
transport and rank imbalance as the next cost investigation rather than
justifying another chemistry optimization. The 128-cubed IC evolved to
a fully refined level-9 mesh, so it cannot be costed as a fixed 128-cubed
run. Checkpoint restart qualification and a bounded performance review
are required before another long calibration segment; preserve the
checkpoint and do not weaken physical refinement to reduce its cost.

Restart reader attempt 417258 failed after six seconds: the checkpoint
header and level counts were read, but the original cold-start
ngridtot=600000 allocation exhausted free grids at level 7. The native
MPI step exited 0 through clean_stop despite the ERROR, while the runner
correctly rejected the log and the batch exited 1. This is a restart
capacity admission issue, not an OOM or demonstrated checkpoint damage.
The reader runner now derives ngridmax=4970304 per rank from the saved
info file instead of reusing the cold-start global grid capacity. All
physics controls, rank count and frozen binary remain unchanged. Retry
417259 was submitted with 192 GiB, one H100, four ranks and two threads
per rank for at most two hours, requesting no additional evolution/dump.
Automatic capacity admission in the restart reader remains a code-level
limitation until independently corrected; this input workaround is not
claimed as its implementation fix.

Attempt 417259 also exhausted capacity at level 9 because the legacy
same-ncpu reader replicated every grid on every rank. The source now
routes multi-rank Hilbert restarts through the existing distributed
reader; K-section handling is unchanged. Attempt 417264 exposed an
incorrect global Hilbert extent: cmp_minmaxorder was passed half the
coarse-cell width. Correcting that alone in 417267 was insufficient:
the call also supplied x_tmp(1,3), while the explicit-shape callee uses
x(NVECTOR,NDIM). With NVECTOR=32 this gives invalid coordinate strides.
Both width and array extents are corrected in restore_hdf5.f90. The
failed 417264 and 417267 jobs were canceled after native ERROR records
to release ranks still waiting in collectives; checkpoints were untouched.

Attempt 417272 rebuilt only the affected reader/provenance objects in
an isolated build. At elapsed 1:15 its rank-1 active counts at levels
2 through 8 are exactly one quarter of the global counts; level 8 is
524288/2097152 rather than the 1835009/2097152 seen in 417264 or the
2097152/2097152 seen in 417267. This supports corrected distribution,
but level-9, hydro/chemistry/particle restoration and normal exit are
still pending. It is not yet a restart PASS or continuity comparison.

417272 subsequently reached a new failure in virtual-grid creation:
`refinement pre-count mismatch 16 12`, followed by MPI_Abort. Unlike
the previous failures, distribution through level 8 was balanced.
The reader directly removed grid slots from the free list without
updating used_mem until the end of the entire AMR restore. The
refine_fine capacity reservation therefore used stale occupancy.
Both direct-allocation sites now update used_mem=ngridmax-numbf,
matching make_grid_fine's existing invariant; capacity/count assertions
remain intact. Attempt 417274 rebuilt and started with the same
checkpoint and unchanged physics. At elapsed 1:22 it was restoring
level 9 with no new ERROR record. This remains a pending verification,
not a claimed solution or restart qualification.

At elapsed 4:10, 417274 completed AMR restore: level-9 active grids on
rank 1 were 4194304/16777216, exact Hilbert per-domain counts ranged
4793490--4793491 (imbalance 1.000000), and Morton self-lookup mismatch
was zero. Capacity grew collectively from 4970304 to 6212928 to cover
virtual grids; local occupied grids were 4970340. Grackle initialization
and HDF5 hydro restoration then began. Normal exit remained pending.

Read-only stacks subsequently placed all four ranks in the streaming
hydro reader's MPI_Bcast, not chemistry integration. Two /proc I/O
samples separated by 60 seconds showed rank-0 rchar increasing from
5673999605 to 6680640061 bytes (1006640456 additional bytes); other
ranks increased only about 7.5 KiB. This confirms continued root-read
progress during the apparently quiet log interval, not a deadlock.
The implementation hdf5_read_dataset_chunk_dp performs independent
rank-0 reads followed by replication to all ranks. These counters are
logical read bytes and not physical disk throughput (read_bytes was
zero on rank 0), nor a whole-file performance forecast. The job was
not canceled or restarted based on the quiet log.

At elapsed 55:29, hydro, Poisson and particle restore had completed;
initialization took 3265.932 seconds. Subsequent forced load balance
completed with post_varcpu_load_balance NaN counts all zero, and
level-9 GPU MG converged in four iterations to 2.134e-5. The script's
original zero-evolution description was incorrect: restore sets
nstep_coarse_old=nstep_coarse, and update_time checks stopping criteria
only after that counter changes. Consequently nstepmax=188 does not
stop before the first resumed coarse step. This job must be evaluated
as a bounded restart-and-evolution smoke test, not reader-only; its
native first Fine step=188 log retained the saved epoch before advancing.
No running input or binary was modified. Future runner wording was
corrected; any original RESTART_READER_PASS wording in the frozen
417274 runner must not be interpreted as proof of zero evolution or
restart continuity. Final step, epoch, diagnostics and exit remain pending.

At elapsed 1:13:07, Slurm confirmed 417274 RUNNING on syn08; the first
resumed coarse step was still calculating. The Fine-step diagnostic's
negative `Ekin_tot=-4.39e-6` is not negative kinetic energy: the existing
update_time.f90 format 779 labels arguments mass_tot, epot_tot, ekin_tot,
eint_tot as Mtot, Ekin_tot, Epop_tot, Eint_tot respectively. Interpret the
second value as gravitational potential energy and the third (2.85e-6)
as kinetic energy. This is a diagnostic-label defect, not evidence of a
negative physical state or restart failure; the running binary is unchanged.
The source format 779 now labels the unchanged argument sequence as
Mtot, Epot_tot, Ekin_tot, Eint_tot. Static inspection confirmed the sole
format-779 call and git diff --check passed. This output-only correction
has not yet been rebuilt or observed in a new runtime; 417274 retains
the original labels. No numerical expression, physics setting or saved
checkpoint was changed by the correction.

At elapsed 1:18:04, 417274 remained RUNNING. Two process samples about
75 seconds apart showed CPU-time increases of 42, 41, 35 and 32 seconds
on the four MPI ranks, with RSS stable near 23 million KiB per rank.
A bounded rank-0 stack snapshot placed both OpenMP computation threads
inside Godunov (ctoprim/trace3d -> unsplit -> godfine1), not a chemistry
solver, restart-read collective or stream wait. This is evidence of
ongoing hydrodynamic work, not completion or a measured scaling result.

Subsequent affinity inspection found all four main computation threads
bound to CPU 8 and all four OpenMP workers bound to CPU 9, despite
Slurm allocating eight CPUs. The job step's effective cgroup CPU set
was 8-11,32-33,58-59. Runtime environment confirmed OMP_NUM_THREADS=2,
OMP_PLACES=cores, OMP_PROC_BIND=close and I_MPI_PIN=0. Thus the observed
execution was not a properly separated four-rank/two-thread performance
sample; the low CPU-time increments cannot establish intrinsic hydro cost.
Without changing physical state, binary or inputs, the existing main/worker
threads were rebound within that allocated set: rank 0 to 8/9, rank 1
to 10/11, rank 2 to 32/33 and rank 3 to 58/59. taskset confirmed the
changes. Timing spans both placements and is not a clean benchmark.
The future restart launcher now passes --cpus-per-task explicitly to
srun and enables verbose core binding. This launcher change still requires
an actual affinity check on its next execution; the exact binding-collapse
origin is not yet established, and the original baseline's completed
thread placement cannot be inferred from this restart observation.

### Restart-and-evolution smoke result: 417274

Slurm COMPLETED, ExitCode 0:0, elapsed 01:23:28 (native step 01:23:25),
MaxRSS 23072900 KiB. Native Run completed at Main/Fine step 189,
t=-11.2423 and a=0.05426 after restoring step 188. Final mcons=1.75e-14,
econs=-1.77e-2, epot=-4.43e-6, ekin=2.88e-6, eint=1.25e-9.
Level-9 GPU MG converged in four iterations to 2.137e-5; all reported
NaN checks were zero. Log scan found no ERROR, abort, Fortran runtime
failure, MG nonconvergence or explicit negative/invalid-state diagnostic.
No output_00002 exists; the original checkpoint remains 40007155072
bytes. Binary checksum verification passed. This is a PASS for this
bounded distributed restart-and-one-step smoke only, not a comparison
against uninterrupted evolution or all-physics production readiness.

Direct resumed-step cooling was 37.902 seconds, compared with 149.695
seconds in the earlier baseline step 187. The native coarse interval was
1655.08 seconds, but includes mixed CPU placement and restart overhead;
do not infer a clean speedup from these measurements. Rank-average
terminal timers still put Godunov first (936.171 seconds, 54.2%) and
reverse ghost-zone imbalance at std/av=0.948. Grackle cooling was
38.360 seconds (2.2%). The native Total elapsed time=1160.120 seconds
does not include full startup: update_time initializes its private timer
on first entry. Use Slurm elapsed for total restart wall cost.

The next safe action before another long calibration segment is a short
launch-affinity check with explicit srun task width and the actual Intel
MPI/OpenMP runtime. Validate disjoint two-core rank placements, rather
than assuming the launcher edit proves the fix. Preserve the science
checkpoint for continuation; there is no new raw simulation dump to remove.

### Launch placement qualification

Short Intel MPI/OpenMP probes 417333, 417335 and 417336 completed with
ExitCode 0 in 3, 5 and 6 seconds respectively. The original launch and
explicit srun --cpus-per-task=2 both reproduced all four ranks sharing
the same two CPUs. I_MPI_PIN=1 with I_MPI_PIN_DOMAIN=omp also reproduced
the overlap. scontrol show config reported TaskPlugin=task/cgroup,
without task/affinity; thus core-binding flags alone were not an adequate
verified placement mechanism here.

The rank wrapper reads its allowed CPU set, assigns each SLURM_LOCALID
a disjoint SLURM_CPUS_PER_TASK slice, restricts the process with taskset,
and sets explicit OMP_PLACES before Intel MPI/OpenMP initialize. It fails
when there are insufficient allowed CPUs rather than sharing cores.
In 417336 this produced rank pairs 16/17, 18/19, 40/41 and 42/43;
both OMP_DISPLAY_AFFINITY and sched_getcpu output confirmed all eight
computation-thread placements. This passes launch placement for the
tested four-rank/two-thread runtime configuration, not a speed benchmark
or general multi-node placement guarantee. The restart runner now uses
this verified wrapper, copying it into the run directory for retention.
No simulation was run by these probes and no raw outputs were created.

### Approved baseline continuation submission: 417338

Submitted .galaxy128-grackle-resume-20261009.sbatch to H100: one GPU,
four MPI ranks, two OpenMP threads/rank, 192 GiB and 24 hours. The verified
rank wrapper is copied into the new run directory and OMP_DISPLAY_AFFINITY
is enabled for actual-runtime placement inspection. The original science
checkpoint remains a read-only input by convention, linked as output_00001;
no original run file or binary is replaced. An isolated incremental build
uses the verified 417274 objects with the output-label correction.

Unchanged frozen sources/assets were checked before submission; only the
previously measured prolongation repair, verified restart reader changes
and output-label correction are admitted relative to the original source
contract. The effective input changes only nrestart=1 and startup capacity
ngridmax=4970304. Physics and the three scheduled aout values (z=2,1,0)
are retained, periodic dumps suppressed with foutput=fbackup=1000000.
A Slurm signal one hour before expiry requests a clean checkpoint-stop.
At the existing 37.26 GiB dump size, three scheduled dumps plus one stop
checkpoint estimate 149 GiB of new output, excluding the retained source
checkpoint; AMR growth may increase this. The shared filesystem had 260 TiB
free at submission; the existing 512 GiB project storage budget is unchanged.

First check confirmed 417338 PENDING (Resources). No run directory or
effective runtime copy exists yet, so actual input/build/placement checks
remain pending. Submission is not execution or progress to low redshift;
inspect the new run's effective.nml once it starts. Do not change physics
or claim a speedup merely because this continuation has been queued.

417338 subsequently started on syn08 at 2026-10-09 04:13 KST. Actual
effective input /gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-resume-417338/effective.nml
was compared against the original baseline: exactly nrestart 0->1 and
ngridtot=600000->ngridmax=4970304 changed; physics/output controls match.
The isolated build recompiled only update_time plus build/git metadata;
binary/source/asset checksum checks passed. Actual simulation runtime
OMP_DISPLAY_AFFINITY confirmed disjoint pairs 8/9, 10/11, 32/33, 58/59
on MPI PIDs 684329--684332, with workers 684498/684497/684500/684499.
AMR restore and Morton/MPI checks completed, exact domain imbalance was
1.000000, and capacity grew to 6212928 as in the passing smoke. Grackle
reported nine species with equilibrium metals, two private owners and
no dust/RT; hydro HDF5 restore had started. First evolved-step correctness
and uncontaminated wall-clock timing remain pending; do not equate this
startup confirmation with successful low-redshift calibration.

At elapsed 2:55, 417338 had completed hydro/Poisson/particle restore;
native initialization time was 81.684 seconds versus 3265.932 seconds
in 417274. Forced load balance took 33.527 seconds. Reported NaN checks
remained zero, level-9 GPU MG converged to 2.134e-5 in four iterations,
and all four MPI processes were live. Rank-0 logical reads totaled
38799499981 bytes with read_bytes=0. This approximately 40-fold observed
startup-time difference follows corrected disjoint CPU placement on the
same node/checkpoint, but is not an isolated attribution experiment;
do not call it disk throughput or whole-simulation speedup. First-step
wall timing and evolved conservation checks remain pending.

At elapsed 4:42, actual runtime also verified the format-779 correction:
Fine step 188 printed Epot_tot=-4.39e-6 and Ekin_tot=2.85e-6, matching
the original unchanged argument values. This closes the output-label
runtime check without treating negative gravitational energy as a
negative-state failure. The first resumed step was still computing.

The first resumed coarse interval in 417338 completed in 433.78 seconds
with direct cooling 37.877 seconds. This includes restart/load-balance
overhead and must not substitute for a steady-state timing sample.
For context, the prior baseline's step 187 was 1830.78 seconds with
cooling 149.695 seconds, but these are different steps and the completed
baseline's thread placement was not measured. At elapsed 9:56, Main step
189 reported mcons=1.75e-14, econs=-1.77e-2, epot=-4.43e-6,
ekin=2.88e-6 and eint=1.25e-9, matching the rounded diagnostics of the
bounded 417274 smoke at the same step. This scalar agreement does not
prove field-level or uninterrupted-restart continuity. No new ERROR or
nonzero NaN diagnostic was found. Step 189 is now evolving; retain the
same inputs and collect subsequent normal coarse intervals before
claiming a production cost improvement.

The subsequent normal coarse interval completed in 398.21 seconds
(6.64 minutes), with direct cooling 37.840 seconds and native memory
21.9--22.0 GiB/rank. This is about 4.60 times shorter than the prior
baseline's 1830.78-second step, without changing resolution or physics.
It remains a cross-step observation rather than an isolated same-state
benchmark; the prior baseline's affinity was not directly observed.
At elapsed 16:38, Main step 190 had mcons=-1.43e-14,
econs=-1.75e-2, epot=-4.48e-6, ekin=2.91e-6, eint=1.27e-9 and
a=0.05445. No new error/nonzero NaN diagnostic was found. The inherited
cosmological energy diagnostic offset persists and is not erased or
declared resolved by the speedup. Step 190 continues on the same contract;
follow the approved low-redshift baseline rather than opening another
chemistry optimization or performance gate.

A second normal coarse interval was 398.20 seconds, direct cooling
37.574 seconds, confirming two consecutive approximately 6.64-minute
intervals at unchanged level-9 mesh size. Main step 191 reported
mcons=4.18e-14, econs=-1.73e-2 and a=0.05463, with no new nonzero NaN
or error record. The local cost improvement is sustained across these
two intervals, not proof of cost at later AMR refinement or low redshift.
The calculation remains active toward the approved observable-calibration
baseline; do not create additional performance gates merely to restate
these already recorded samples.

## 2026-10-10: terminal baseline and bounded failure repair

The preceding statements describe observations while 417338 was live, not
its final state. Slurm now confirms 417338 TIMEOUT, elapsed 24:00:20;
417338.0 FAILED (1:0). The final log reports step 391, a=0.08510
(z approximately 10.75), mcons=-1.29e-6, econs=-3.79e-3, followed by
`Grackle NEQ solve/budget failed, level/status: 9 2`. Recent completed
coarse intervals were about 418 seconds, cooling about 40 seconds.
These are not low-redshift observable-calibration results.

Status 2 is adapter INPUT rejection, before local_solve_chemistry, not a
Grackle nonlinear convergence verdict. It combines invalid scalars,
negative/nonfinite species, electron mismatch and nucleus/metal closure;
the old log cannot identify which predicate failed. No tolerance is relaxed.

The launcher wrote `jobcontrol.txt` containing `0 -1` at 03:12:51 KST,
but the runtime has no corresponding job-control acknowledgement and no
new output_00002. Preserve the original step-188 output_00001 checkpoint.
Do not describe 391 as a restartable saved state.

VPATH selects clean_stop from patch/cuRamses/update_time.f90, which calls
finalize_timer and MPI_FINALIZE, NOT amr/end.f90's MPI_ABORT implementation.
The old cooling tile called clean_stop inside an OMP loop. A unilateral
failure can therefore enter collective finalization while other ranks
continue computing/communicating, and it calls MPI from an OMP owner.
This explains a credible hang mechanism but is not a captured stack trace
of the already terminal job; the exact failure timestamp is not available.

Bounded repair: propagate a Grackle tile failure through an OMP max
reduction, exchange status on the main thread at the end of cooling_fine,
and abort the communicator on error. Normal checkpoint/clean_stop behavior
is unchanged. Failure-only adapter diagnostics now report owner, tile cell,
all nine species, density, energy, metals, electron and mass residuals.
The caller reports rank/coarse step/level and native cell IDs. Input and
post-solve admission tolerances and successful chemistry are unchanged.

Regression script .galaxy128-grackle-failurecheck-20261010.sbatch retains
the same RT/dust-off physics and 4 MPI x 2 OMP affinity, limits evolution
to step 190, and runs the existing mixed-owner transactional/direct-reference
and cosmological-unit tests. No new scheduled dump: the three original
z=2,1,0 aout values are beyond this bounded run; periodic clocks stay large.
The original checkpoint is a read-only symlink, approximately 37.26 GiB;
new full dump storage expected zero, /gpfs free space measured 257 TiB.
Only tiny test/build/log outputs are expected.

418460 compiled the production binary successfully but FAILED before
simulation launch: the C test link omitted GNU Fortran runtime dependencies
of static libgrackle. Add the same HDF5/GFortran/quadmath libraries used by
the production Makefile; replacement job 418471 submitted. This is a test
launcher error, not a physics failure. Do not claim either regression or
the late step-391 problem passed until measured results exist.

418471 exposed a separate reference-test defect: its direct-reference
contexts initialized/read HDF5 tables concurrently inside the OMP loop.
HDF5 attribute errors and its shutdown error accompanied FAIL (3).
Production initializes all Grackle owners serially before the OMP solve;
serialize only the reference test's local_initialize_chemistry_data calls
to match this supported setup. Keep the reference solves parallel.
Replacement 418473 reports private-owner/mixed delayed-cooling/transactional
PASS (0), and cosmological-unit contract PASS, adapter-vs-direct max
difference 0.0. Its full RAMSES bounded restart has begun reading the
checkpoint; integrated evolution is not yet declared passed. Preserve all
failed-attempt logs instead of suppressing these failures.

At 18:53 KST, 418473 is RUNNING on syn08, has restored step 188 and
completed the initial level-9 GPU MG solve. Wrapper CPU sets are disjoint:
0/1, 2/3, 32/33, 58/59. Do not infer a deadlock from a quiet restart log.

Reproduction job 418474 is PENDING (Dependency), afterok:418473. It checks
the regression's PASS marker and binary/source hashes before starting;
never restart merely because an observation expires. It uses the same
4 MPI x 2 OMP, one H100 and 192 GiB, with a 48-hour ceiling to avoid the
previous 24-hour cutoff. Bounded evolution is 188 -> 392, not a new
low-redshift production run. Existing timings imply roughly one day of
elapsed work; resource/AMR changes can alter this estimate.

One explicit `380 1` job-control request retains the pre-failure state
without changing the original three z=2,1,0 scheduled aout values or large
periodic output clocks. During this high-z bounded reproduction those
scheduled outputs are not expected. One approximately 37.26-GiB new dump
is expected; a signal-driven clean-stop could add another (approximately
74.52 GiB total new dump allowance at the present mesh size). Growth can
increase dump size. Free /gpfs space is 257 TiB as measured at 18:53.
The original step-188 checkpoint is preserved and linked, not overwritten.
The existing 512-GiB project dump budget remains the limit.
Effective inputs will be
`.galaxy128-grackle-reproduce-418474/effective.nml`; the submitted script
constructs them from 417338 changing only nstepmax to 392. Physics, output
redshifts and refinement parameters remain unchanged. An error must retain
failure diagnostics and return nonzero, not silently bypass a cell.

### Input-closure hypothesis (not yet attribution of step 391)

Source inspection finds independent minmod slopes for every conserved
variable in interpol_hydro; the existing Grackle limiter checks thermal
admissibility, not element closure. Independent nonlinear slope limiting
does not in general preserve a linear nucleus identity in each child,
although the mean of the eight children still conserves each parent field.
An independent arithmetic example of the inspected minmod formula has
left/center/right rho=(0.90,1.00,1.10), H=(0.70,0.76,0.70),
He=(0.20,0.24,0.40), metals=0. Every input has H+He=rho.
The implemented half-difference minmod gives slopes rho=0.05, H=0,
He=0.02. For child offsets +/-0.5 the closure residual is -/+0.015,
well above the adapter's unchanged 1e-8*rho input allowance.

This arithmetic reproduction is NOT a captured failing simulation stencil
or a native full-prolongation execution. Constant-composition stencils
would not reproduce this example. A similar consistency issue can also
arise from independently limited passive fluxes; do not assign it to AMR
without the actual rejection values and saved pre-failure state. The
reproduction binary remains frozen; do not replace it with a speculative
normalization, looser threshold or untested reconstruction change.

Output counter safety was checked directly with an attribute-only h5dump:
the original checkpoint has `/header/ifout=2`. restore_hdf5 restores that
counter; dump_all forms the directory from ifout before incrementing it.
Thus the planned step-380 output targets new output_00002, not the preserved
output_00001 symlink. This is a source/header contract check, not proof that
the pending reproduction has already written its new checkpoint.

## Integrated regression completed; reproduction started

418473 COMPLETED (0:0), 17:42 Slurm elapsed, at 19:07:37 KST.
The launcher records FAILURECHECK_PASS. Native Main/Fine step 190 has
a=0.05445, mcons=-1.43e-14, econs=-1.75e-2, epot=-4.48e-6,
ekin=2.91e-6, eint=1.27e-9. These printed scalars match the old baseline
at step 190; no field-level equality is inferred. First/second intervals
were 458.23/424.47 seconds. MaxRSS for the Slurm step was 23595608 KiB.
No ERROR, forrtl, nonzero NaN or MG nonconvergence record was found, and
no output_00002 was created. Existing Grackle direct-reference/transactional
and cosmological-unit tests also passed.

This verifies the successful chemistry/restart path after the failure-handling
edit. It does NOT yet exercise an actual distributed fatal input rejection,
prove its abort latency, or resolve the later step-391 error. These are
explicitly left for the bounded reproduction, not marked complete here.

418474 started automatically at 19:07:37 and is confirmed RUNNING.
Actual copied effective.nml differs from 417338 only in nstepmax
100000 -> 392. The runtime jobcontrol.txt is the intended `380 1` request.
All three input hashes passed; the copied executable and successful regression
executable have identical SHA256
`c16b06efca712df30a44fa8c312bab707dc77669db6c29824edb5783e8a58ebe`.
The previous FAILED attempts are retained and not relabeled as passes.

At 00:56 KST on 2026-10-11, Slurm still reports 418474 RUNNING on syn08,
with 5:48:46 elapsed and a 48-hour limit. The native log has advanced through
Main step 237 (scale factor 0.06288); the effective target remains step 392.
The live Slurm step reports 10:42:30 accumulated CPU time and about 23.8 GiB
maximum RSS, so this is active computation, not an idle pending job. From the
observed 189-to-237 progression, a linear estimate leaves about 18.5 hours,
placing step 392 near 19:30 KST on Oct 11. This is an ETA only: timestep and
AMR cost vary, and the reproduced step-391 failure could terminate the job
earlier. The log has no terminal error marker as of this check.

Follow-up at 01:00 KST: 418474 remains RUNNING and has completed Main/Fine
step 238 (a=0.06305). Slurm CPU time increased to 10:50:09 and the run log
advanced during the preceding minute, confirming forward progress through
another coarse step. The simple rate-based estimate remains approximately
18 hours to step 392 (around 19:00--20:00 KST Oct 11); retain the earlier
warning that step 391 may fail first.

### 2026-10-11 bounded step-391 reproduction result

Job 418474 terminated at Main/Fine step 391 with the intended Grackle adapter
input-rejection path, not a scheduler timeout, OOM, or Grackle nonlinear-solve
failure. The C bridge returned status 2; the distributed caller issued
`MPI_Abort`, so `srun` ended with status 137 and Slurm records the job as
FAILED. The reproduced coarse-step scalars match the earlier trajectory at
printed precision (`a=0.08510`, `mcons=-1.29e-6`, `econs=-3.79e-3`).

Failure-only diagnostics identify owner 1, zero-based tile cell 27, which maps
to the 28th `ind_leaf` entry, native cell ID 27911025, at level 9, with
`rho=104.01971842107167`, `metal=0.025668816368954364`,
and nucleus/mass residual `H+He+metal-rho=-1.5168492808470546e-6`. Its
relative magnitude is `1.45823244e-8`, while the unchanged admission bound is
`1e-8*rho=1.04019718e-6` (the residual exceeds it by `4.76652097e-7`, or a
factor 1.458). `invalid_species=0`; the logged electron abundance equals the
charge-derived expectation, and H, He, density, energy, metal, length, and
cooling inputs are otherwise admissible. Thus the observed rejection is
specifically nucleus/metal closure, not a bad species value or electron
charge mismatch.

Correction after direct directory/header inspection: `output_00002` WAS
written at step 380. Its `COMPLETE` marker is `00002`, its HDF5 header has
`nstep_coarse=380`, `ifout=3`, `nvar=27`, and `aexp=0.0838181620328843`.
The HDF5 file is 39,987,305,480 bytes (about 37.24 GiB). The log acknowledges
the step-380 job-control request and HDF5 output. The earlier no-checkpoint
report was incorrect; no restart from step 188 is needed for the next diagnosis.
The `output_00001` input link, new checkpoint and all failure evidence remain
preserved. This confirms the late-step input closure failure and narrows
its symptom; it does **not** yet prove that AMR
minmod prolongation caused it because the failing parent stencil was not
captured. Do not relax the `1e-8` admission tolerance or normalize the species
by guesswork. The next useful step is a targeted reproduction/capture of the
actual failing prolongation state, followed by a conservation-preserving
reconstruction change only if that mechanism is confirmed.

### Targeted diagnosis from the verified step-380 checkpoint

The bounded diagnostic launcher `.galaxy128-grackle-closure-20261011.sbatch`
restarts from output 2 and retains the same step-392 ceiling, four MPI ranks,
two OMP threads/rank, one H100, 192 GiB and unchanged RT/dust-off physics.
The effective namelist changes only `nrestart=1 -> 2` and the job-control
path, whose request is beyond the bounded run. All original scheduled outputs
remain z=2,1,0, with `foutput=fbackup=1000000`; no new dump is expected.
The step-380 checkpoint is linked read-only. Current /gpfs free space is
255 TiB; preserving the checkpoint adds no copied 37.24-GiB file.

`GRACKLE_CLOSURE_TRACE=1` is a compile-only diagnostic, disabled by default.
It records at most the first locally observed admitted-parent/bad-child AMR
stencil and the first affected leaf's old/new conservative state at the
hydro commit on each rank. It changes no field, timestep, tolerance or error
handling. It can distinguish inherited closure from a new hydro defect and
show a real prolongation violation; proximity alone does not assign the
specific rejected cell to either cause. The diagnostic is a scoped response
to the observed failure, not another production qualification campaign.

Submitted job `419596` on 2026-10-11 at 19:59:12 KST; Slurm confirms
RUNNING on H100 node `syn09` from 19:59:13, initially compiling the isolated
diagnostic binary. Effective inputs are
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-closure-419596/effective.nml`.
The estimate is 1--2 hours after execution starts, based on the previous
~418-second coarse intervals and roughly twelve steps, with a three-hour
resource ceiling. This estimates the bounded diagnosis, not low-redshift
calibration, the repair itself, or whole-project completion.

At 20:02 KST the isolated build and source-hash checks passed; both diagnostic
markers are present in the actual executable. Binary SHA256 is
`ab22317d8e0ad62c1bcf785a4d9e91b66d12aa131faf0cf3f37a0d182e01ac2e`;
effective namelist SHA256 is
`077a66f49acf5c3f0a6a73f35c0e4998474e2fce4e4f540f5408d7ee721d032b`.
The actual namelist diff contains only the two declared restart/control-path
changes. The simulation has read the checkpoint, reported time integration
startup and entered load balancing; the four ranks have initialized the
assigned H100. No new coarse-step or closure result is claimed yet.

A read-only Astra consultation (agent `01a12aa3-0d4e-7340-b711-a9a6007a4926`)
confirms the structural risk in independent nonlinear composition slopes,
but does not attribute the specific rejected cell to it. It recommends
closed, positive reconstructed left/right states before the Riemann solve
(rather than dividing signed carrier fluxes at zero net mass flux), and
parent-composition scaling by the final child density for AMR. These preserve
independent conservative element integrals through shared fluxes and child
means. Face admissibility still needs a full native-update check. The advice
also identifies an important limit: inconsistent saved averages cannot be
made exactly closed while preserving those same averages. No polluted
restart state is silently normalized. This is advisory under the standing
Astra consultation permission, not an additional audit gate.

The exact diagnosed Makefile and hydro/Grackle sources are retained under
`.galaxy128-grackle-closure-419596/diagnosed_sources/` before subsequent
corrections; the running binary remains frozen. This preserves the actual
baseline behind the trace, even if source corrections subsequently land.

At 20:12 KST Slurm confirms `419596` RUNNING, elapsed 13m05s. The restarted
simulation has reached the Main/Fine step 381 marker (`a=0.08394`,
`mcons=-9.32e-7`, `econs=-4.06e-3`), demonstrating evolution beyond checkpoint
restore/load balancing. Neither closure-trace site nor the Grackle input
rejection has emitted a record so far. This is forward-progress evidence;
the original late-step failure has not yet been reached or resolved.

### Confirmed native AMR closure defect and bounded correction

While the frozen diagnostic run continues, the existing native prolongation
driver was extended with an optional closed, variable-composition fixture.
It calls the actual routine linked to the diagnosed build, not a reimplemented
formula. At 20:16 KST the retained uncorrected objects reproduce child
mass residuals of +0.05/-0.05 from closed `(rho,H,He,Z)` parents
`(1,.7,.2,.1)`, `(2,1.2,.6,.2)`, `(4,2.8,.8,.4)`. Each independent parent
field mean remains conserved. This confirms a real structural AMR defect;
it does not identify the producer of the particular cosmological rejected
cell. Evidence: `.galaxy128-grackle-closure-419596/native_closure/mass.log`;
native executable SHA256
`a28c3b8af5b4741c96edbd4c0253fb60001622df79e76e27bfc70e04ac76f7e5`.

The scoped correction, compiled only with `GRACKLE_NEQ`, uses the final
child density and parent abundances for `imetal`, `ichem:ichem+10` and
`igrackle:igrackle+8`. It leaves unrelated passives and hydro density,
momentum and energy reconstruction intact. Each child mean preserves its
own parent inventory, with first-order subcell composition. The accepted
common thermal contraction still runs first. An inconsistent parent is
not renormalized: its relative defect is inherited.

The isolated corrected actual-routine build passed all four existing thermal
cases, with exactly zero reported parent-mean residuals. The same composition
fixture now has max absolute mass residual `2.220446049250313e-16`, positive
composition, parent fractions and species budgets preserved. A deliberately
inconsistent parent retains its inventory and relative residual `-5e-7`,
demonstrating that the correction does not hide a bad checkpoint.
Evidence: `native_repair/thermal.log` and `native_repair/mass.log` under the
same diagnostic run; native executable SHA256
`c379bec7111b440ca8fa1bbf9be5490a56d0281ee5b9b4941834113819d72348`.
These tiny native checks ran as additional short CPU-only Slurm steps in
the already allocated H100 job, with isolated build objects; they created no
raw simulation snapshots and did not replace the running executable.

This closes the demonstrated AMR reconstruction defect at native scope.
The frozen cosmological diagnosis is still running; MUSCL face consistency,
the particular step-391 failure and corrected integrated evolution remain
unresolved. No full-project or calibration completion is claimed.

### Confirmed native face defect and corrected runtime candidate

The same native driver now also calls the actual `unsplit` routine over a
periodic 4x2x2 box assembled from two tiles. It exercises variable density,
H/He/Z and ion/molecular fractions, forward/reverse HLLC flow, and LLF
counterflow with 16 zero-net-mass-flux faces. The uncorrected frozen objects
produce maximum face closure errors of `2.8216218879002775e-4`,
`2.5467261894727283e-4`, and `4.8016666666666694e-4`, respectively;
post-update cell errors reach `4.0623019801921956e-4`,
`4.4278074480752494e-4`, and `6.1216666666652486e-4`.
Independent periodic global inventories still conserve to roundoff, which
shows why that check alone would not detect the local composition defect.
Evidence: `native_face_before/mass.log` in diagnostic run 419596.

The correction closes only temporary traced left/right compositions before
the Riemann solve. It contracts inadmissible predicted components toward the
originating cell, applies H+He+Z normalization to the temporary fractions,
aligns predicted primordial species with elemental budgets, and derives
electron charge. The hydro prediction is unchanged and each carrier keeps
its own shared conservative face flux. No division by a signed flux is used,
so the zero-net-mass-flux LLF case is covered. Invalid cell averages are not
bulk repaired, and the Grackle admission tolerance remains `1e-8`.

The corrected actual-routine face test gives maximum errors of
`3.4694469519536142e-18`, `6.9388939039072284e-18`, and
`3.5236570605778894e-18`; maximum updated cell errors are
`4.4408920985006262e-16`, `2.2204460492503131e-16`, and
`1.1102230246251565e-16`. Independent global inventories, positive species,
updated H/He species budgets and electron charge all pass. The existing four
thermal AMR cases and inherited-invalid-parent check also pass.

At 20:36 KST, CPU-only Slurm step 419596.6 completed successfully in 1m31s,
building an isolated **real runtime** candidate with diagnostic tracing
disabled, then linking the same regression against its actual objects.
Candidate SHA256:
`2706069613798f0a095c5fa15c11873a83a2118349a3f1b7838c4c1e2e803a45`.
Its retained input/source hashes, build log, AMR, face and thermal test logs
are under `.galaxy128-grackle-closure-419596/transport_candidate/`.
The frozen diagnostic executable was not replaced.

This correction is in CPU hydro. Grackle-NEQ with `gpu_hydro=true` is now
explicitly rejected rather than silently using uncorrected GPU reconstruction;
GPU Poisson remains supported. The namelist generator enforces the same
constraint, and `mkrun.py` selects CPU hydro for this profile. The actual
effective cosmological namelist already uses CPU hydro/GPU Poisson, so this
does not change its active solver assignment. Equivalent GPU hydro closure
is unverified and is not claimed supported by these tests.

At 20:37 KST, diagnostic job 419596 was RUNNING on H100 node syn09, elapsed
38m19s, and had reached coarse step 384 from checkpoint 380 without a
reported closure rejection yet. A corrected replay launcher is prepared as
`.galaxy128-grackle-transport-20261011.sbatch`; it uses the identical effective
namelist, checkpoint, rank binding, physics, tolerances and output clocks,
and the corrected pinned executable, with four MPI ranks, two OMP threads per
rank, one H100 and 192 GiB. It refuses to overwrite run directories and
requires the existing native passes. It is not submitted before reviewing
the frozen diagnostic result. The bounded target is still step 392; no new
dump is expected (scheduled aout values remain 1/3, 1/2, 1, and periodic
dump clocks are 1000000). The 37.24-GiB step-380 checkpoint is linked, not
copied; /gpfs currently has about 255 TiB free.

The two structural defects are fixed at native-test scope and the real
candidate has built. Attribution of the actual step-391 rejection and the
corrected cosmological replay remain pending. A contaminated saved average
cannot be repaired by these conservative reconstructions; if it triggers
admission, diagnose its origin and restart from an appropriate valid state,
without silently normalizing the checkpoint or relaxing the bound.

Step-reporting clarification from source inspection: the selected VPATH uses
`patch/cuRamses/update_time.f90`. `Main step=N` is printed in `update_time`
before the in-flight interval's hydro/cooling work, not after that whole
interval has completed. `adaptive_loop` advances the coarse counter after
the base-level call returns and then prints the elapsed-time/SFR summary.
Accordingly a live marker means the run has reached that interval; it alone
does not prove its chemistry step passed. The retained target marker 392
together with `Run completed`, successful exit and absence of failure
diagnostics tests evolution through the old failing interval 391. No test
ceiling, tolerance, physics or stopping rule is changed by this reporting
correction.

### Captured late-step hydro producer and corrected replay submission

Diagnostic job 419596 terminated at 21:29:02 KST with `FAILED` after 1h29m49s
(Slurm recorded ExitCode `9:0`; the runner records the application-step
status 137 after collective `MPI_Abort(2)`). It again fails in interval 391,
not from OOM, walltime expiry or Grackle nonlinear iteration failure.

The first hydro trace is native rank 3, level 9, local cell **27376990**,
which is also the 28th tile-cell entry in the subsequent owner-1/cell-27
Grackle rejection. Local storage IDs are not assumed invariant across the
step-380 restart/load balance; the earlier run reported 27911025.
The traced pre-update state has
`rho=112.40205663358053`, `H=85.385950335808303`,
`He=26.988291444979861`, `Z=0.027814864020920340` and
`H+He+Z-rho=+1.1228564744669711e-8`, relative `9.989643500273051e-11`.
It is inside the unchanged `1e-8` admission bound. The post-update state has
`rho=104.01971977591144`, `H=79.018430627334382`,
`He=24.975618814953670`, `Z=0.025668816775090385` and
residual `-1.5168482860872246e-6`, relative `-1.4582314674130582e-8`.
The absolute residual change is `-1.5280768508318943e-6`; the relative defect
grows by about 146 times.

The chemistry adapter logs those same post-update values, finite/nonnegative
species, and exactly charge-derived electron abundance. This ties the real
rejection to a new mismatch across the complete conservative hydro-update
window, before chemistry integration, rather than an already inadmissible
old cell average. No prolongation trace was emitted. It is not a proof that
AMR can never contribute, or a decomposition of every term in the hydro
window, but the independent native MUSCL reproduction already demonstrates
the non-closing face reconstruction and its correction. The failing physical
state closely reproduces 418474 despite the restart-dependent local IDs;
bitwise equality is not claimed.

After reviewing this terminal trace and rechecking the candidate binary and
source hashes, the corrected bounded replay was submitted as **419694** at
about 21:32 KST using `.galaxy128-grackle-transport-20261011.sbatch`.
Its intended effective namelist is
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-transport-419694/effective.nml`.
It copies the diagnostic effective namelist unchanged, links the preserved
step-380 checkpoint, and uses the pinned trace-free runtime candidate
`2706069613798f0a095c5fa15c11873a83a2118349a3f1b7838c4c1e2e803a45`.
No additional source model, tolerance, stopping criterion or output request
was introduced. It requests H100x1, MPI4 x OMP2, 192 GiB and 3 hours.
The retained aout schedule is z=2,1,0 and both periodic dump clocks are
1000000; no new dump is expected in this replay. /gpfs free space remains
about 255 TiB. Submission is not startup or a repair pass; actual execution
and progression through the failing interval still require observation.

At 21:34:53 KST, Slurm confirms 419694 RUNNING on the same H100 node syn09,
started at 21:31:58. The actual effective file is present and `cmp` against
419596's effective namelist succeeds; SHA256 is
`077a66f49acf5c3f0a6a73f35c0e4998474e2fce4e4f540f5408d7ee721d032b`.
The runner checked the pinned binary and began the restart at 21:31:59.
The simulation has entered GPU MG; its reported level-8 residual reaches
`4.421e-5`, below the retained `1e-4` tolerance. Startup is verified, but no
repaired late-step outcome is claimed yet.

### Used-face optimization and first corrected interval timing

The trace-free correction initially checked all 384 predicted face states
per 3-D grid stencil. Only 72 of those states are consumed by `cmpflxm` for
the eight active cells. The source now restricts the same closure helper to
those interface slabs: the left state uses `qm` at the coordinate one lower
along the face normal, and the right state uses `qp` at the interface
coordinate. The closure formulas, conserved cell averages, physics switches
and chemistry admission tolerance are unchanged. This reduces helper calls,
not the physical stencil or the number of evolved interfaces.

An isolated build in
`.galaxy128-grackle-closure-419596/transport_face_pruned/` completed as a
CPU-only overlapping Slurm step of allocation 419694, without rebuilding or
replacing the running executable. Its production binary SHA256 is
`986c751ff3ac9ee95e014f0e018c659f7b3ae151952bf20af3f91bf520759efe`.
The expanded existing native fixture covers forward/reverse HLLC and
zero-net-mass LLF counterflow along each of x, y and z (nine cases). All pass:
maximum face closure error `6.9388939039072284e-18`, updated-cell error
`4.4408920985006262e-16`, and global absolute inventory discrepancy
`1.1368683772161603e-13`. Each counterflow case exercises 16 zero-net-mass
faces. The native AMR test and all four thermal prolongation tests also
pass, including preservation of the deliberately invalid inherited parent
residual. Logs and build-time hashes are retained in that isolated build;
the four relevant source snapshots are in its `corrected_sources/`.
These results do not assert bitwise whole-run equivalence or measured
cosmological walltime improvement.

Job 419694 continues to run the frozen **unpruned** candidate SHA256
`2706069613798f0a095c5fa15c11873a83a2118349a3f1b7838c4c1e2e803a45`.
At 21:59 KST it is RUNNING on syn09 and has completed interval 380 and
entered 381. The completed interval reports `1067.69 s`, with direct
cooling `54.595 s`, SNRT time zero, no NaN flags, and approximately 26.6 GiB
reported memory. The GPU MG level-9 residual reaches `9.723e-5` below the
retained `1e-4` bound. This proves advancement, not passage of the original
391 failure. If subsequent intervals retain this cost, the original
three-hour limit ending around 00:32 KST will be insufficient; the first
interval alone is not used to declare a stall or cancel the live replay.

At about 22:03 KST, an attempt to extend this same allocation to 4h30m
(`scontrol update JobId=419694 TimeLimit=04:30:00`) was rejected with
`Access/permission denied for job 419694`. No scheduler state was changed;
the three-hour limit and original 392 stopping criterion remain in force.
The verified live replay is not restarted merely because this request was
denied. Subsequent completed-interval timings must determine whether the
remaining window fits; the independently built used-face candidate is
available if measured cost justifies an optimized replay. Its native pass
does not substitute for the actual late-step cosmological verification.

For comparison, the frozen diagnostic's same first interval took `519.78 s`
with direct cooling `50.534 s`, versus the corrected unpruned replay's
`1067.69 s` and `54.595 s`. The full interval is about 2.05 times slower,
while the separately reported cooling cost changes by about 4 seconds.
This is not a fine-grained hydro profile, but it does not support blaming
the extra first-interval cost on the chemistry integration itself. The
used-face change targets the unnecessary reconstruction-helper calls
without reducing chemistry accuracy or relaxing physical admission.

The existing `.galaxy128-grackle-transport-20261011.sbatch` launcher is now
prepared to use the retained `transport_face_pruned` artifact by default;
`GRACKLE_TRANSPORT_CANDIDATE=transport_candidate` explicitly selects the
original unpruned artifact. Other candidate names are rejected. The launcher
retains the same checkpoint, effective namelist, resource layout and target
392, checks the selected binary hash and native test results, and records
the selected name at startup. `bash -n` passes. This does not submit a new
job or change 419694: its original submitted script and binary identity
remain frozen in its run directory. At 22:05 KST 419694 is still RUNNING
in interval 381, with no reported input rejection or error.

### Measured steady cost and optimized replay — 22:07 KST

The unpruned run's next completed interval, 381, reports `906.84 s` and
total runtime `2053.610 s`. The matching diagnostic interval took
`428.47 s`; the corrected unpruned interval is about 2.12 times slower.
At that observed cost the remaining ten intervals would give approximately
3h05m total runtime, beyond the unextendable 3h allocation. This measurement,
not log silence or an observation timeout, motivated switching to the
already built used-face implementation. Its expected performance benefit
remains an estimate until measured in the real replay.

After both completed intervals and the absence of rejection were captured,
the driver deliberately cancelled its own job **419694** and submitted
**419714** using the prepared launcher with `--time=04:00:00`. This is an
operator cancellation, not a chemistry failure or a passed interval-391
verification. The old logs, effective inputs, original submitted script,
binary and native checks remain available; the saved step-380 checkpoint is
not modified or deleted.

The new run's intended effective namelist is
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-transport-419714/effective.nml`.
It is the same 380-to-392 test, with H100x1, MPI4 x OMP2 and 192 GiB,
using used-face binary SHA256
`986c751ff3ac9ee95e014f0e018c659f7b3ae151952bf20af3f91bf520759efe`.
No physics, source, tolerance, resolution or output schedule changes.
As before, `noutput=3`, aout z=2,1,0, periodic clocks 1000000; no new full
dump is expected before 392. The existing checkpoint is 37.24 GiB and is
linked read-only for restart; current /gpfs free space is 255 TiB.
The requested four-hour cap supplies headroom, not a longer integration.
Submission alone does not prove startup or passage of interval 391.

At 22:08 KST, Slurm confirms 419714 RUNNING on H100 node syn08 with
the four-hour limit; startup is recorded at 22:07:13. The runner explicitly
reports `candidate=transport_face_pruned` and verifies its pinned binary
hash. `cmp` confirms the actual effective namelist is unchanged from 419596;
SHA256 remains `077a66f49acf5c3f0a6a73f35c0e4998474e2fce4e4f540f5408d7ee721d032b`.
The native log has entered HDF5 restart and is reading the level-9 grids.
419694 is now terminal CANCELLED after 35m13s; its failed application-step
status reflects that deliberate cancellation. No `output_00003` was
created by it. The optimized replay is on another node of the same GPU
type, so subsequent timing comparisons are operational walltime evidence,
not a controlled same-node microbenchmark.

### First optimized completed interval — 22:19 KST

419714 has completed interval 380. Its elapsed-time summary reports
`572.77 s` and total runtime `652.8400 s`, with direct cooling `40.858 s`.
The matching unpruned completed interval cost `1067.69 s`: operational
walltime is reduced by approximately 46.4% (1.86x speedup), not inferred
from the helper-call count. The node changed from syn09 to syn08, so this
is not an isolated same-node performance attribution. No chemistry input
rejection is reported, all printed NaN counters are zero, and no new full
dump exists. Initial GPU MG reached residual `8.392e-5`, within `1e-4`.
Slurm still reports RUNNING with average CPU time increasing and MaxRSS
`28471760K` per recorded task at the 22:19 observation.

A provisional operational finish window is 23:45--00:15 KST, based on the
first completed interval and the original restart's cold-versus-steady cost.
This is not a promise of terminal success: the next completed interval
must refine the estimate, and passing the old failure still requires
evolution through 391, marker 392, `Run completed`, normal exit and no
failure diagnostics. The four-hour deadline is around 02:07 KST, so the
new observed pace has substantial headroom without extending the physics.

### Steady-interval cost and generated backend choice — 22:31 KST

The optimized replay completed interval 381 in `460.79 s`, with total
runtime `1113.620 s` and direct cooling `40.952 s`. The same interval in
419694 took `906.84 s`; operational walltime decreases by approximately
49.2% (1.97x). The original diagnostic cost was `428.47 s`. At the 22:29
observation, 419714 is RUNNING in interval 382, without a reported input
rejection. Extrapolating ten more intervals at the measured steady cost
gives approximately 95 minutes total from the 22:07 start, or around
23:43 KST; the provisional reporting window is 23:40--00:00 KST.
This remains an estimate conditional on continued numerical acceptance.

While that fixed-input run proceeds, existing frontend regressions exposed
an actual generated-input omission: the setup database defines
`gpu_hydro=False`, and the formatter previously discarded a supplied False
as a default. The selected CUDA source instead initializes
`gpu_hydro=.true.` (`patch/lagRamses/amr_parameters.jaehyun.f90`), with the
non-CUDA disable in `read_params.jaehyun.f90` guarded by `!HYDRO_CUDA`.
Thus a generated Grackle namelist could pass setup validation but inherit
GPU hydro and fail native admission. The current replay is unaffected:
its effective namelist already explicitly writes `gpu_hydro=.false.`.

The shared formatter now retains every supplied GPU hydro choice, including
False, and explicitly writes False for `grackle_neq` even if that key was
absent. Both direct generation and `mkrun.py`/GUI use this formatter; GPU
Poisson selection is preserved. Two cases were added to the existing GUI
test module for direct round-trip/admission and the actual shared wizard.
The existing grouped-form fixture's count was also corrected from seven to
eight fields, reflecting the previously added cooling-selector prompt;
the production GUI layout was not changed by this fixture correction.

The initial full frontend run **completed** in 24.898 s (it did not timeout):
the two new assertions exposed the omission and one existing grouped-form
expectation was stale. After correction, the headless existing suite reports
60 tests in 1.390 s, 59 PASS and one display-dependent skip; all three
existing radiation-selector tests also PASS. The focused three-case rerun
passed before the full suite. No simulation was launched by these tests,
and no native binary or running effective input was changed. README records
the CUDA/default distinction; `git diff --check` and the running binary's
pinned checksum both pass. These are backend/input consistency checks, not
additional physical qualification gates.

At 22:34 KST, 419714 remains RUNNING and has also completed interval 382:
`460.04 s`, total runtime `1573.660 s`, direct cooling `41.214 s`. The two
post-startup completed intervals are therefore stable at about 7m40s,
supporting the 23:40--00:00 provisional finish window. The original
interval-391 failure is still ahead; no terminal repair pass is recorded.
The Makefile explicitly names `amr_parameters.jaehyun.o` and
`read_params.jaehyun.o`; conventional unsuffixed `amr_parameters.f90`
paths are absent in this checkout and are not used as source evidence for
the generated-backend fix.

Read-only preparation for the subsequent restart confirms that the existing
HDF5 reader distinguishes the saved file counter `ifout` from schedule index
`iout`. `patch/lagRamses/restore_hdf5.f90` restores the counters but then
recomputes `iout` using the existing `restart_output_index` in
`amr_parameters.jaehyun.f90`; the namelist's aout/tout arrays remain
authoritative. For checkpoint `aexp=0.0838181620328843` and the retained
schedule `aout=(1/3,1/2,1)`, that helper selects the first future target,
not the third target merely because `ifout=3`. No output-clock code or
effective input was changed and no new clock test was launched. This source
inspection is not a claim that the no-output bounded replay qualifies the
later z=2/1/0 dumps.

### Prepared continuation, not submitted — 23:14 KST

The existing October-9 resume script refers to an earlier reader build and
checkpoint, so it is not the continuation entry point for the composition
repair. `.galaxy128-grackle-calibration-resume-20261011.sbatch` now prepares
the previously approved reference model using the pinned face-pruned binary
and preserved step-380 checkpoint. It refuses launch unless job 419714 is
Slurm `COMPLETED 0:0` and has the normal step-392 acceptance marker. Bash
syntax passes; an in-memory effective-input diff contains only
`nstepmax=392 -> 100000` and the stop-control path changing to the new run's
local `jobcontrol.txt`. That path is essential: a signal-written local stop
request would not affect the bounded replay's external control file.

The prepared resource/output policy is unchanged from the approved science
segment: one H100, MPI4 x OMP2, 192 GiB, 24 hours, z=2/1/0 dumps and at most
one clean checkpoint-stop; no automatic restarts or sweep. It checks the
four actually selected physical input hashes rather than old unrelated
executables. It has **not** been submitted, and neither its syntax check nor
the continuing repair test certifies a z=0 trajectory or high-level RT/dust
production. At 23:14, 419714 is RUNNING in interval 388 after completing
380--387, most recently in 458.68 s, without reported rejection.

### Corrected replay closure — normal exit at 23:43:51 KST

Job 419714 and its actual MPI compute step are Slurm `COMPLETED 0:0`.
The retained log reaches `Main step=392` and `Run completed`; the runner
emits `CLOSED_TRANSPORT_REACHED_392`. All twelve intervals 380--391 finish,
including the independently diagnosed failure interval. A complete log
search finds no Grackle input rejection, ERROR/abort/Fortran runtime error,
MG nonconvergence, or nonzero printed NaN counter. Effective input and
binary checksums still match their frozen records. No new raw snapshot
exists: the run has only the linked `output_00002` checkpoint, which remains
needed for the approved calibration restart and is preserved.

The compute step lasts 1h36m38s, consumes 11h56m31s total CPU, and reports
MaxRSS `28471760K` per recorded task. Eight allocated CPUs correspond to
12.884 core-hours (approximately 92.7% CPU utilization), not the launcher's
CPU usage. Godunov remains the dominant recorded timer (3597.921 s average,
63.0% of the enclosing timer total); cooling is 508.780 s (8.9%), Poisson
478.990 s (8.4%). Nested MG timers must not be added to Poisson again.

The driver closes the conservative AMR/face-composition repair for this
unchanged cosmological Grackle profile. No chemistry tolerance was relaxed,
no cell-average composition normalization was added, and GPU Poisson stayed
enabled while hydro remained the verified CPU implementation. This is not
z=0 calibration, MHD/RT/dust admission, or universal production approval.
The independent coupled M5 result and outstanding physical-source/release
identity concerns retain their separately documented scopes.

### Approved baseline resumed — job 419850

After the terminal repair evaluation, the prepared continuation was submitted
as 419850. Slurm confirms RUNNING on H100 node syn08 at 23:48 KST; the
launcher records start at 23:47:51 and all binary, selected-source and four
physical-input hash checks PASS. The actual effective input is
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-calibration-419850/effective.nml`,
SHA256 `5d646ebb2ae93908a04941aa3b39c4890c474beb882811677a16529fcbb80667`.
Its diff from the evaluated replay contains exactly the two prepared run
controls, not a physics/source/refinement/backend change. It reuses the
step-380 checkpoint because the bounded no-output replay did not save a new
step-392 state. The repeated short interval is necessary restart evolution,
not a second qualification campaign. GPU Poisson remains enabled; CPU hydro
and nine-species Grackle are intentionally retained.

Resource ceiling: one H100, eight CPUs (MPI4 x OMP2), 192 GiB, 24 hours:
at most 192 allocated core-hours and 24 GPU-hours for this segment. Output
policy remains three science dumps at z=2/1/0, plus at most one checkpoint
on the one-hour-before-limit stop request; the planning envelope is 512 GiB,
not a measured late-time dump size. `/gpfs` currently has 255 TiB available.
No future jobs, ensemble or automatic repeated restart has been submitted.
Successful launch is not evidence of later redshift attainment, observable
calibration, or high-level RT/dust physical release.

### Later closure rejection — 2026-10-12, job 419850

The continuation is FAILED, not still running: the application aborts at
06:26:58 KST after 6h39m09s. At coarse step 431, level 9, rank 3, Grackle
rejects input status 2. This is not an OOM or GPU Poisson failure. For the
rejected cell, rho=0.46749856869861578, H=0.35529891128816726,
He=0.11219966206259736, Z=4.674985795925399e-11, and the printed
H+He+Z-rho residual is 4.6988987145724082e-09, slightly above the existing
1e-8*rho admission bound. Charge agrees with its species-derived value;
all reported species are nonnegative. The scale factor is
0.089483972134108747 (z approximately 10.175), not a calibration endpoint.

The native caller aligns Grackle H/He species to the separate elemental
carriers before solving and commits only chemistry and thermal state.
Thus the elemental carrier/density producer needs investigation; the
closeness to the numerical bound does not establish round-off as the cause.
No tolerance relaxation, cell-average normalization, or speculative physics
change has been applied. The bounded 380--391 replay PASS remains valid for
its actual interval and must not be promoted to full-trajectory approval.

The prepared `.galaxy128-grackle-late-closure-20261012.sbatch` reuses that
evaluated face-pruned Grackle build's Makefile, objects and unchanged native
sources, rebuilding only to enable the existing closure traces. It restarts
the preserved step-380 checkpoint and stops at 432, retaining all physics
and output controls except the finite step limit. Resources: one H100,
MPI4 x OMP2, 192 GiB, at most 10 hours (80 allocated core-hours, 10 GPU-hours).
No new raw dump is requested; the existing 37.24 GiB checkpoint is preserved.
This is a first-origin diagnostic replay, not a new qualification campaign
or a repaired-run claim. Submission and results are recorded separately.

The bounded diagnostic is submitted as **420107**. At 07:18:45 KST Slurm
reports RUNNING on H100 node syn08; the MPI launcher starts at 07:18:10.
Its build log contains exactly two Fortran compilations (interpolation and
Godunov) with `GRACKLE_CLOSURE_TRACE`; the before/after checksum comparison
confirms every other evaluated object is unchanged. Native source, input
and Makefile checks pass before launch. Only `nstepmax=100000 -> 432`
changes in the effective namelist. This is successful diagnostic launch,
not a reproduced failure verdict, repaired trajectory, or z=0 completion.

### Residual-transport diagnosis and operator sequencing — 2026-10-12

The read-only diagnostic 420125 completes 0:0 in seven seconds. Its native
contrast links the evaluated transport objects, not a rewritten Python
transport model. An exactly closed state stays closed to 1.17e-16. An
admitted relative residual of 9.95e-9 stays constant in absolute mass to
2.22e-16 while density falls from 1 to 0.95; its new relative residual is
1.0473684375929344e-8 and would fail the unchanged Grackle admission bound.
This proves the existing face-only projection can strand inherited residual
mass in a cell rather than advect it consistently with density.

The same job streams only the four required hydro fields and leaf flags from
the preserved checkpoint in bounded chunks. All 134,141,080 level-9 leaves
and 9,581 level-8 leaves are within the existing admission bound at step 380.
Their maximum relative residuals are 1.0001104381525136e-10 and
1.0000593018532103e-10, respectively; none exceed 5e-9. The actual passive
IC stores H=0.7599999904632568, He=0.24000000953674316 and
Z=1.000000013351432e-10. H+He already sums to one before adding Z, as the
original pre-enrichment provenance explicitly records. No checkpoint, IC,
input tolerance, or production source is modified by this diagnostic.

Evidence is retained under `.grackle-closure-diagnosis.mEOoNv/`; the earlier
420124 probe failed to link omitted native support objects and is not a
physics verdict. Its corrected 420125 probe supplies the results above.

The operator now requires final confirmation in the actual failed trajectory
BEFORE the consistent-transport correction. At the current observation,
420107 has reached coarse step 392 without a closure trace or input rejection;
the actual first offending cell is still pending. Do not label the native
contrast alone as final proof for the step-431 cell, apply a speculative
production patch, normalize cell averages, or weaken chemistry admission.
After confirmation, retain the existing independent global inventories,
advect the admitted residual consistently with density, preserve the exactly
closed control, and compare the original/corrected native interval including
the previous step-431 rejection. This is a bounded numerical repair, not a
new physical model or qualification program.

### Actual first failure confirmed; consistent transport repair — 2026-10-12

420107 aborts at 13:54:28 KST after 6h36m26s, reproducing the same native
Grackle rejection at coarse 431, level 9, rank 3. Before the chemistry call,
the hydro trace identifies cell 22284867. Its density decreases from
0.47752029133641127 to 0.46749856869861250 (2.0987%), while the elemental
H+He+Z-rho residual is EXACTLY the same printed
4.6989026003529943e-9 before and after the conservative hydro update.
The relative residual consequently crosses from 9.84021555859421e-9 to
1.0051159329607035e-8. The subsequent chemistry input rejects that same
state. This confirms the stranded-residual mechanism in the actual failed
trajectory, not just the manufactured native contrast. It is not OOM, a
new chemistry integration error, or an AMR-prolongation creation event;
the trace is emitted after hydro exchange and before the uold commit.

Following the operator's confirmation-before-repair instruction,
`grackle_close_face` now preserves the originating admitted parent's
H+He+Z composition sum instead of unconditionally setting it to one.
All 21 carrier fields retain the common reconstruction/positivity factor;
the existing H/He nucleus, electron and tracked-metal checks remain. No
cell-average normalization, admission-tolerance relaxation, new field,
physics parameter, or namelist change is introduced. Nonadmitted parent
states still bypass repair and remain subject to the existing rejection.

The existing native fixture is extended, not replaced by a new framework:
constant admitted residuals +9.95e-9, -9.95e-9 and +1e-10 exercise the same
three axes, positive/negative HLLC advection and zero-mass-flux LLF cases.
It checks residual/density consistency, independent global inventories,
species/charge closure and positivity, retaining the exactly closed tests
and AMR/thermal controls. The new test must also reject the old object.

Build/regression 420772 is submitted; its isolated build freezes every
evaluated production object except `umuscl.kjhan.o`. At the current
observation it is PENDING (`QOSGrpGRES`), so no PASS or corrected-trajectory
claim is made yet. The prepared replay uses the exact effective input from
420107, preserved step-380 checkpoint and the gated pinned new binary,
with one H100, MPI4 x OMP2, 192 GiB and a ten-hour ceiling. It stops at 432.
The unchanged scheduled aout values are 1/3, 1/2 and 1, beyond this bounded
interval near a=0.09; foutput/fbackup are 1000000, with no initial/requested
dump. Expected new raw dumps: zero. The 37.24 GiB restart is preserved;
current /gpfs free space is 281 TiB. This is not a full-z=0 restart or release.

The corrected replay is submitted as **420773**, with Slurm
`afterok:420772` dependency and `GRACKLE_RESIDUAL_BUILD_JOB=420772`.
It cannot launch if the native regression/build fails and independently
checks the native acceptance markers and pinned source/binary before MPI.
Current states: 420772 PENDING (`QOSGrpGRES`), 420773 PENDING (`Dependency`).
The bounded replay is therefore submitted but not running; no test PASS or
corrected cosmological result has yet been obtained. The retained production
source SHA256 is
`9b4a309153634fa14cfdc09c83f6e1f0e0ff4b34cf3dc056c8e8a0c0172f5365`;
the extended existing fixture SHA256 is
`d30685eb8d40d7d009cc7e6dbb552b40f9f9012818ddb62b3ff6717a24df92c9`.

### Residual repair accepted — 2026-10-13, 22:09:26 KST

Native build/regression 420772 completes 0:0 in 18 seconds at 09:57:27.
The unchanged thermal/prolongation tests, closed-face tests and new inherited
residual advection tests pass. Linking the same new fixture against the old
transport object fails explicitly with
`transport amplified admitted relative residual past gate` (status 128),
providing the intended negative control rather than a linker/runtime failure.

Corrected cosmological replay 420773 and its MPI step complete 0:0 at
22:09:26. The retained log reaches `Main step=432` and `Run completed`,
and the launcher emits `GRACKLE_RESIDUAL_REPLAY_REACHED_432`. Thus interval
431, which aborted both 419850 and the traced reproduction 420107, completes
with the consistent residual transport. A full retained-log search finds
no Grackle input rejection, ERROR/abort/Fortran runtime error, MG
nonconvergence or nonzero printed NaN counter. The effective namelist is
byte-identical to the failed traced replay; the sole production object
change remains the confirmed transport correction.

The actual MPI compute step lasts 7h05m49s and records 53h21m18s total CPU,
approximately 94% of eight allocated CPUs. Its reported task MaxRSS is
28627684K; this is not an aggregate four-rank peak. Godunov remains 66.7%
and cooling 8.5% of the native timing total. The accepted binary SHA256 is
`feffbd862cd23bb659a961c55bcc54ef27c0d9f0342a303e24c00885c639f7fa`.
Post-run binary/source/input checksum checks and live-source comparison pass.

The run directory contains only the linked existing `output_00002` restart,
whose COMPLETE marker remains `00002`; no new raw simulation dump exists
to clean up. Preserve that checkpoint for the calibration continuation.
The driver closes this numerical repair and its 380--431 integration test.
The z=2/1/0 science trajectory and observable parameter tuning remain pending.

### Parameter calibration resumed — job 421829, 2026-10-13

The operator explicitly directs parameter tuning after accepting the residual
repair. The reference trajectory is submitted as 421829 and starts on H100
syn08 at 23:01:01 KST; MPI launch is 23:01:02. Binary/source and all four
selected physical-input checks PASS before restart restoration. At the first
observation the job is RUNNING and reading the level-9 checkpoint, not yet
qualified as a completed new evolution interval or an observational fit.

Effective input:
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-calibration-421829/effective.nml`,
SHA256 `5d646ebb2ae93908a04941aa3b39c4890c474beb882811677a16529fcbb80667`.
The only input difference from accepted replay 420773 is
`nstepmax=432 -> 100000`. It retains eps_star=.02, the 12.5 cMpc/128^3 IC,
Grackle nine-species chemistry/HM2012, original stellar sources, RT/dust/AGN
off, CPU hydro/chemistry and GPU Poisson. The pinned executable is the
accepted 420772 build with SHA256
`feffbd862cd23bb659a961c55bcc54ef27c0d9f0342a303e24c00885c639f7fa`.
Since the bounded replay intentionally saved no step-432 dump, restart is
from the preserved step-380 checkpoint. This requires repeating the verified
interval; its state cannot be recovered from a completed process's log.

Allocation: one H100, MPI4 x OMP2, 192 GiB, at most 24h/192 allocated
core-hours/24 GPU-hours. Output is three science snapshots at z=2/1/0,
aout=1/3,1/2,1 and tout=3*1e100, periodic clocks=1000000. Slurm USR1 at
23 hours requests one clean checkpoint-stop via the existing native
jobcontrol path. Budget envelope: <=128 GiB per late-time dump, <=512 GiB
for three science dumps and the checkpoint, against observed /gpfs free
space of 280 TiB. The retained step-380 dump is approximately 37.24 GiB;
late-time size remains a forecast. Preserve scientific and restart products.

The active first calibration control is eps_star, centred on .02 within the
existing proposed .005--.05 range. Obtain the reference resolved galaxy
counts, stellar mass/SFR, metallicity and HI statistics at low redshift,
then price paired eps_star responses and select subsequent values using
measurement-compatible observations and finite-volume uncertainties. Inactive
dust/AGN controls are not fit targets in this run. No parameter sweep,
resolution enlargement or automatic restart chain has been launched.

### Reference trajectory observation — 2026-10-14 10:02 KST

Job 421829 remains RUNNING on H100 syn08 after 11h01m. The latest
completed main-step line observed is 460 (a=.09229); the subsequent native
SFR diagnostic reports z=9.82, Nstar=61, Mstar=3.951e7 Msun,
SFR=.7000 Msun/yr and SFRD=3.604e-4 Msun/yr/Mpc3. These are native
simulation-wide diagnostics, not measured galaxy-bin observables or evidence
of a fitted parameter value. Latest printed NaN counters are zero; level-9
MG converged in eight iterations to 7.110e-5. The reference evolution still
precedes the first requested science output at z=2.

The current runner is specifically anchored to the old output_00002
(step 380). After its planned clean stop, select and verify the newly written
COMPLETE checkpoint and its native epoch before constructing a continuation;
resubmitting this runner unchanged would repeat step 380. Neither a normal
process exit alone nor the scheduled USR1 time proves checkpoint completion.

Diagnostic interpretation checked against the copied build source
`build/patch/cuRamses/adaptive_loop.jaehyun.f90`: the printed Mstar uses
`mstar_tot`, the accumulated formed-mass counter, rather than a sum of current
stellar particle masses after mass return. The printed SFR differences that
counter between diagnostic calls, dividing by the elapsed code time converted
with the current units; with ncontrol=1 this is a coarse-step diagnostic, not
an observational 10/100 Myr tracer average. SFRD divides it by the code's
comoving box-volume convention. Use actual particle masses, birth masses and
ages with a declared aperture and averaging interval for galaxy observations;
do not fit a stellar mass function directly from this global log counter.

Live Slurm accounting at 2026-10-14 10:23:33 KST: 421829.0 has four tasks,
elapsed 11:22:32, AveCPU=21:24:00 per task, AveRSS=27,949,043 KiB,
MaxRSS=28,724,056 KiB. This corresponds to approximately 85.6 CPU core-hours
over 91.0 allocated core-hours (94% utilization), a sampled task-RSS sum
of 106.6 GiB and maximum single-task RSS of 27.4 GiB. These are live Slurm
samples, not a final simultaneous job-wide peak or a strong-scaling result.
The 192 GiB allocation has headroom at this epoch; future AMR growth remains
unmeasured. The printed RAMSES mem percentages describe array capacity and
must not be substituted for these physical-memory measurements.

### New failure at coarse step 530 — 2026-10-14

On the operator's one-shot status request after monitoring was stopped,
Slurm reports 421829 FAILED at 19:46:21 KST after 20:45:20; the simulation
step exits 46:0 in accounting and the batch runner records status 174.
The first explicit error is `conservative Grackle AMR parent is thermally
inadmissible`, followed by SIGSEGV in MPI collectives. Last printed coarse
step is 530, a=.09922 (z approximately 9.08); level 10 has 31,282,121 grids.
The last native SF diagnostic reports 137 stars and accumulated formed mass
1.013e8 Msun. Slurm's maximum single-task RSS is 64,699,000 KiB; this is not
a job-wide aggregate or an OOM verdict. No new output exists in the run:
output_00002 remains the symlink to the step-380 checkpoint. The failure
preceded the scheduled 22:01 clean-stop request.

The operator directs the next repair. Address-to-symbol resolution of the
actual pinned executable gives `finalize_timer -> clean_stop ->
interpol_hydro -> godfine1 -> godunov_fine OpenMP region`. The selected
cuRamses clean_stop invokes MPI timer collectives and MPI_FINALIZE from a
worker-thread-local error path. This explains the secondary collective crash;
it does not establish the producer of the invalid parent state. The numerical
predicate also rejects positive internal energy below 128 ULP of total or
kinetic energy, so the old message alone cannot identify negative energy,
insufficient arithmetic margin, nonfinite fields or low density.

The scoped patch replaces only this fatal clean_stop call with serialized
stderr diagnostics (rank, step, parent density, E/kinetic/thermal/margin and
the complete local input stencil), flushes them and executes ERROR STOP 17
without MPI. The launcher kills remaining tasks on failure. Numerical
interpolation, admission tolerances, chemistry and physical parameters are
unchanged. Build/regression job 422621 checks the existing thermal/composition/
face/advection regressions and deliberate negative-energy and positive-but-
below-margin rejection from OpenMP worker 1 under a two-task Slurm launcher.
Do not claim the thermal-state producer is repaired from these error-path tests.

The prepared bounded replay is
`.galaxy128-grackle-thermal-replay-20261014.sbatch`: restart 380, stop 531,
save step 520 via foutput=520, retain original physical inputs. One H100,
MPI4 x OMP2, 384 GiB (headroom for the observed level-10 growth), 30-hour
ceiling/240 allocated core-hours, optional clean stop at 29 hours. Original
noutput=3 at z=2/1/0 remains beyond this short interval; expected new output
is the step-520 checkpoint plus at most one walltime-stop checkpoint,
estimated <=128 GiB each and <=256 GiB total. Measured free space at submission
is 275 TiB. Admission requires a successful native build; the replay is queued
with an afterok dependency and rechecks all native PASS markers and the pinned
binary/source hashes in its prologue. Preserve the new nearby checkpoint
for diagnosis rather than repeatedly replaying from step 380. Automatic
continuous monitoring remains stopped per the operator's instruction.

### 64-ULP thermal-admissibility follow-up — 2026-10-20

The bounded replay 422623 reproduced the step-530 rejection with a finite,
positive parent thermal energy of `1.0757318430129614e-19`; the former
128-ULP cutoff was `1.0842021724855044e-19`. Thus this was a representable
positive state narrowly below an implementation margin, not negative thermal
energy. The slope-contraction criterion is now 64 ULP (still purely numerical,
not a physical temperature floor); total energy and physical state are not
floored or renormalized. The cold-state native fixture now starts at 128 ULP
and verifies conservative prolongation keeps every child above the new margin.

Build/regression job 429929 completes 0:0 on syn103 without requesting a GPU.
All four native thermal/conservation cases, composition mass closure, face
closure, inherited-residual advection, and both intended OpenMP-worker fatal
diagnostic cases pass. The latter deliberately exit 17; the Slurm build wrapper
checks their diagnostics and itself completes successfully. This qualifies the
numerical interpolation gate only, not the cosmological trajectory.

The targeted integration replay 429930 is configured to restart from the
verified step-520 `output_00003` written by failed replay 422623 and stop at
step 531, with no scheduled raw output. It requires one H100, MPI4 x OMP2 and
384 GiB. At this record update it is pending with `QOSGrpGRES`; no post-change
simulation result is available and the thermal producer is not yet considered
resolved. The step-520 checkpoint is retained.

Submitted bounded replay **422623**, `afterok:422621`, using
`GRACKLE_THERMAL_BUILD_JOB=422621`. Intended effective input is
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.galaxy128-grackle-thermal-replay-422623/effective.nml`;
it is generated, checked and hashed on the allocated node before launch.
At submission the native build is Priority-pending, so neither its regression
verdict nor execution of the replay is yet claimed. No recurring monitoring
process or automatic further restart chain was installed.

### Eight-ULP arithmetic-bound follow-up — 2026-10-20

The 16-ULP replay 430369 again reached step 530 and rejected three finite,
positive stored states. Independent extended-precision evaluation of the
printed binary64 state gave residuals of 14.2042, 14.8052 and 15.1651 ULP.
The rejection diagnostic had retained a stale 64-ULP display value while the
predicate was 16 ULP; both now report the selected 8-ULP reserve. For
`K=0.5*sum(p_i**2)/rho`, the six rounded operations give a first-order
`gamma_6` bound of about 3 ULP near `E ~= K`; subtraction `E-K` is exact in
this cancellation regime by Sterbenz's lemma. Eight ULP is a numerical
arithmetic reserve, not a physical thermal floor. Native job 430585 passed
six conservative-prolongation cases, exact 8-ULP boundary equality, the
14-ULP unchanged-parent regression, conservation/face/advection checks, and
both expected OpenMP-worker rejection cases. The bounded step-520-to-531 A100
replay 430588 terminated at step 530 with status 17. Four finite positive
parent states were rejected because their computed double residuals were
6.1-7.6 ULP, below the imposed 8-ULP reserve. Independent decimal evaluation
of the round-tripped binary64 values gives positive exact stored-state
residuals of 6.0829, 6.7414, 6.8214, and 7.5793 ULP. This falsifies the hard
8-ULP admission threshold as overly restrictive for the observed trajectory;
it is not evidence for negative thermal energy, an OOM, or Grackle cooling
instability. The four source stencils and Slurm accounting are retained in
`.galaxy128-grackle-thermal-margin-replay-430588/run.log` and the matching
Slurm log. The repair keeps the 8-ULP double fast path and evaluates the
stored-state sign in extended precision only inside the near-cancellation
region; it also requires the downstream double-precision residual itself to
be strictly positive. Native job 430767 passed the prior six prolongation
cases, all four captured runtime states unchanged, the exact 8-ULP boundary,
a positive 4-ULP parent, a negative 1-ULP parent rejected from an OpenMP
worker, and mass/face/advection conservation gates. The integrated A100 replay
430769 restarted from step 520 and is running; production trajectory
qualification remains pending its step-530/531 result.

### Step-530 total-energy cancellation blocker — 2026-10-21

Integrated replay 430769 terminated at step 530 with status 17 after 1:19:48;
it did not reach the requested step 531 and produced no new RAMSES dump. Four
finite parents (rank/vector 1/11, 4/4, 3/11, 2/14) had `E-K == 0` in the
binary64 arithmetic used by RAMSES/Grackle. Re-evaluation of the exact stored
binary64 inputs in extended precision gives positive residuals of only
0.4291–0.5628 ULP. The current guard deliberately rejects these states because
the downstream production path also computes `E-K` in binary64; this is not a
rejection by the 8-ULP reserve alone. Relaxing that reserve cannot resolve a
zero binary64 thermal state.

An offline reproduction from the four logged seven-cell stencils, applying the
current minmod slopes and common conservative contraction, found no tested
contraction that yields eight strictly positive binary64 child thermal
energies. At ordinary contraction factors, only four children are positive;
the remaining four are negative, while zero contraction leaves all children
at zero. Thus preserving the conserved parent total energy and momentum does
not provide admissible positive-pressure children for these stored states.
The failure is a high-Mach total-energy cancellation/thermal-resolution limit,
not evidence of an OOM or a Grackle rate-table failure. The native tests that
passed on job 430767 exercised nearby but representably positive binary64
states and did not cover this zero-residual case.

The effective input has `pressure_fix=.true.` and `beta_fix=0`. Existing
`enew` is initialized from the same binary64 `E-K` residual at each hydro step;
it is not a persistent independent thermal/entropy state and therefore does
not recover the lost thermal information here. Do not weaken the rejection,
inject an arbitrary energy floor, or call the simulation qualified. The next
physics/numerics bundle must design and test a genuinely independent thermal
state (dual-energy/internal-energy or entropy evolution), including shock
heating switches, Grackle source updates, AMR prolongation/restriction, restart
and feedback coupling, before another long replay. This is the established
high-Mach use case for dual-energy schemes (e.g. [Bryan et al. 1995](https://doi.org/10.1016/0010-4655(94)00191-4)); the project still needs a RAMSES-specific conservation and coupling design.

The separate A40 job 430823 passed the bounded current-source M5+dust+CHIMES
wiring smoke; it does not qualify this cosmological Grackle/thermal gate.
