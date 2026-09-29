# Cosmological coupling — active work, not admitted yet

Operator approved the coupling repair and preapproved normal implementation,
regression, bounded execution and performance checks. New scientific model
choices must be identified explicitly; no repeated approval for routine work.
Few-pc zoom remains long-term. Bounded cosmological IC/initial-step jobs have
run below; the full RT/feedback/dust cosmological model is not yet admitted.

## Completed bounded correction: CHIMES CMB epoch

The native receiver had a fixed2.727K Compton background. Added a validated
serial-boundary expansion setter; `chimes_prepare_level` supplies aexp for
cosmology and1 for noncosmo before cell callbacks. Each existing cell solve
copies that configuration privately. Do not call the setter concurrently with
OpenMP cell solves. Invalid factors leave the previous temperature unchanged.
The epoch is derived from the existing clock, not another checkpoint field or
namelist selector; frontends require no new option. Dust's IR background is a
separate unresolved contract, NOT fixed by this gas-Compton change.

Compiled bridge, ISO-C interface and RAMSES runtime adapter in the existing
NVAR199/NENER0/SNRT/CHIMES/DUST_LIVE/DUST_DYNAMICS CPU build. Extended the existing
long-interval fixture rather than creating a new test framework.

Allocated grammar-debug tests in `.cosmo-coupling.g8Heia`:

-538991: a=.01 ->272.7K, a=0 rejection, a=1 exact restoration; existing
  long interval passed.13s job wall,2 cores,26 allocated core-seconds;
  compute-step TotalCPU12.119s, sampled MaxRSS107236KiB.
-538995: additionally exercises real native Compton cooling/heating of the
  SAME dilute ionized H at100K for1e12s. At a=1 it reaches99.9999984456566K;
  at a=.01 it reaches228.278739890876K, with accepted nucleus/charge budgets.
  Restoring a=1 leaves the existing68.65Myr photon-chemistry regression passing:
  element1.73875e-10, charge6.88275e-14, photon5.55599e-12, energy1.32895e-11.

These are native receiver tests, not full expanding-box or galaxy tests.
No simulation raw outputs were created. Preserve logs and binaries; no cleanup
of ICs or unrelated runs. Slurm CPUTimeRAW is allocation, not measured CPU.

## Remaining in the approved coupling scope

- Pressureless grain specific enthalpy must stay constant absent exchanges:
  with scale_v proportional to a^-1, its code energy/mass must scale as a^2.
  Check the AMR in-flight/reflux clocks before choosing the runtime hook.
- Photon code number density already converts via an a^-3 density unit.
  Spectral redshift/group migration and the physical dust IR/CMB bath need an
  explicit consistent treatment; do not insert an extra a^-3 dilution.
- Fresh GRAFIC initialization currently has a hardcoded near-zero metal floor
  and H/He defaults, not a specified primordial electron/molecule population.
  Actual source table starts above Z=0; this applicability question remains.
- Check integrated proper-time source intervals, initialize admitted state,
  run bounded expansion/reader checks, then compare identical-work MPI/OMP
  layouts. Keep the cosmological dust admission guard until this is supported.

## Implemented conversions awaiting coupled execution

Added `dust_expansion_level` and a DUST_LIVE-only `update_time` hook, matching
the existing MHD in-flight-state epoch: scale only unew(idust_energy) by
(a_new/a_old)^2 on active and reception cells of levels1..ilevel. No dust mass,
gas total energy or photon dilution changes. Added explicit object dependency;
VPATH ordering is unchanged. This compiles; live expansion/reflux verification
is still required before admission, not established by dimensional analysis.

Native initial H/He indices now come from elem_h/elem_he in the enrichment
configuration, independently of an optional legacy yield table. Legacy mode
keeps its original path. No abundance floor or initial-metal mixture was
changed; that is a scientific initial-condition choice. Compiled successfully.
No new namelist option was introduced, hence no frontend option change.

## IC-only baseline execution (not coupled-model validation)

`.cosmo-coupling.g8Heia/reader-run/run.nml`: cosmological hydro+DM,128^3,
levelmin=levelmax=7; cooling, star formation, sinks, RT and dust all inactive.
Existing pinned `ramses_chem_fixed3d` used, not the newly patched modules.
Legacy enrichment mode is explicit solely to satisfy this binary's required
namelist while all stellar activity is off. No legacy stellar yields run.
GRAFIC omega_b=.049 matches; all NaN_CHK counters zero; main step1 reports
mcons=0 and econs=0, a=.01058. nstepmax=0 nevertheless performs the documented
initial coarse evolution. No output directory was created: only input/logs.

Job539003 exited0 after a namelist error (missing required enrichment group):
NOT a pass. Job539005 failed in unsplit's OMP worker with missing stack setup.
Restored KMP_STACKSIZE=512M and unlimited main stack. Job539006 completed:
29s job wall,16 CPUs,464 allocated core-seconds; compute step28s and403.509
actual CPU seconds. Sampled maximum per-task RSS2205228KiB, NOT aggregate RSS.
Hydro Godunov phase21.128s rank mean,83% of the reported25.456s phase total.
This does not measure late-time cooling/RT cost or establish parallel scaling.
Reader wrapper rejects errors even if RAMSES returns0 and reads NaN counters,
not the literal diagnostic label. Logs of failed attempts retained.

Second CMB job538995:13s wall,26 allocated core-seconds, compute CPU12.075s,
sampled RSS107816KiB. No raw simulation outputs need deletion.

## Fable review and driver disposition

Fable read-only plan review completed via Claude model fable, no --bare,
Read/Grep/Glob only. Requested Q-GOAL then Q-LEAN and essential corrections;
output `fable_cosmological_coupling_plan_2026-09-11.json` is the review record.
One review covers this coupling design; do not create per-fix audit gates.
Accepted the required dust-energy conversion and native H/He mapping, with
the admission guard retained. Do NOT adopt its claim that an evolving dust
CMB bath only matters outside the pilot: this IC begins atz99. Also, a
perturbed cooling run is not a homogeneous adiabatic expansion test; do not
expect exact a^-2 gas temperature while Compton/chemical losses are active.
No-redshift local RT, fixed IR bath and absent UV background are omissions,
not completed cosmological radiation physics. Its proposed Z~1e-10 seed
metallicity is pre-enrichment, not a numerical conversion; placing all metals
in one element is also a chemical-pattern assumption. Request the operator's
scientific choice rather than silently adopting these to satisfy the tables.

Reference context: [Rosdahl et al.2013](https://academic.oup.com/mnras/article/436/3/2188/1247446)
describes RAMSES-RT's supercomoving formulation. It is context, not proof that
this project's S_N/CHIMES/dust implementation already has those corrections.

## Approved pre-enrichment: native input and CPU comparison

The operator subsequently approved explicit Z=1e-10 pre-enrichment. See
`galaxy_pilot_preenrichment_2026-09-11.md`; no all-metals-in-one-element shortcut.
Both runs below use the exact same effective namelist SHA256
`645ffe01e7f3c93578587dacd2f627d22596b723be30b940744074ae1bf3eb73`
and pinned binary `.cosmo-coupling.g8Heia/ramses_cosmo_preseed3d`, SHA256
`94660e06243d00ae61e750eec5a0d6e196ddb1660ccd272ea44e954d307bddd9`.
This binary includes CMB epoch, native H/He mapping and grain expansion hook,
but dust/cooling/RT/star formation are deliberately inactive in these readers.
Native channel-resolved element layout is active; all 12 explicit pvar files
were read. This establishes input consumption, not chemical evolution.

| Job / directory under `.cosmo-coupling.g8Heia` | MPI x OMP | Job wall | Compute TotalCPU | Per-task MaxRSS |
| --- | --- | --- | --- | --- |
|539024 / `preseed-reader`|8 x 2|31s|432.315s|2258668KiB|
|539031 / `preseed-reader-4x4`|4 x 4|36s|415.706s|3585880KiB|

The compute steps both completed. The batch wrappers failed: grammar-debug
does not have rg, and list-directed Fortran wraps long input paths onto two
lines. Replaced checks with standard awk and reran postprocessing on both
preserved logs, successfully. Do not relabel the Slurm FAILED wrapper states
as COMPLETED. No expensive calculation was repeated merely to repair a check.
Reader validation checks every pvar, zero NaN counters, clean completion and
absence of output directories. Both main-step summaries have mcons=econs=0,
epot=-3.23e-8, ekin=2.09e-8, eint=3.76e-11, final a=.01058 (printed precision).
These are not cellwise equality measurements. Initial evolution is one coarse
step even though nstepmax=0. No raw outputs/checkpoints were created or deleted.

Allocated job CPU seconds:496 and576; actual compute CPU above excludes batch
overhead. Compute elapsed31 and35s; RAMSES elapsed23.9596 and25.3040s. Godunov
rank mean21.335 and21.419s, respectively. Poisson rank mean1.436 vs2.863s.
8x2 is faster in this one initial-step sample, whereas4x4 uses slightly less
actual CPU. This is not strong scaling, a repeated benchmark, or the cost of
late-time CHIMES/RT. Keep the baseline allocation until coupled costs exist.

## Native grain expansion and off-knot bath work

`dust_expansion_native` links the actual production objects. Its compact
state test exercises two levels, active and reception cells, two successive
expansion ratios, a disabled hook and unit ratio. Only unew grain energy gains
(a_new/a_old)^2; all gas/mass variables, uold and inactive cells are unchanged.
Native test passed. This tests the hook, NOT an AMR reflux evolution or live
cosmological dust admission. Reject a nonpositive ratio before squaring it.

The native IR initializer now accepts a background between temperature knots,
using the material solver's same log(T)-linear emissivity interpolation; exact
knots retain their original powers. An off-knot equilibrium check at sqrt(20*50)
K passes, with no out-of-range extrapolation. Existing native transient, halo,
rollback and coupling fixture checks also pass (fixture balance7.38116e-10).
This is an interpolation capability only: do not turn on cosmological dust
merely because the initializer now accepts a time-dependent temperature.

The persistent IR field is in physical erg/cm3, unlike the primary photon's
code number-density convention. Its cosmological conversion must be addressed
separately. Moreover, the regular material operator represents excess over a
fixed bath and requires grains at/above that bath. Updating only the bath would
reject cold injected grains or silently change the radiation zero point.
Adding CMB heating to primary power without changing the IR bookkeeping would
double-count bath-derived re-emission/reabsorption. None of these shortcuts was
implemented; the cosmological dust admission guard remains intact.

Final off-knot build and existing native OpenMP backend regression passed,
including material interpolation, gas/dust exchange, stiff/cold absolute
material, scattering, Fe, sublimation and rollback. No GPU validation claimed.
Pinned final binary `.cosmo-coupling.g8Heia/ramses_cosmo_offknot3d` SHA256
`f6d62c41283d728072ff9b2b132d1844a39b959ef195642987e95c4d140352dc`.
Native hook and IR logs are `dust-expansion.log` and `ir-offknot.log` there.

Follow-up design review requested only for the unresolved CMB/RSLA coupling,
not another general audit: `fable_cmb_rsla_design_2026-09-11.json` (completed).
The existing absolute C/silicate/Fe receiver admits zero Fe masses and cold
Planck emission, so another material solver may be unnecessary. However,
physical-c Planck emission with reduced-c absorption means simply loading
physical blackbody energy is not a Tcmb equilibrium. The review is asked for
equations/units and the smallest defensible reuse, comparing absolute radiation,
signed excess and an explicitly limited local bath. It must answer Q-GOAL then
Q-LEAN and distinguish an operator scientific choice from routine repair.
Current contract bath is10K, not2.727K; CHIMES's Compton bath is separately2.727/a.

The same compact grain test also calls production `update_time(2)` twice with
an explicitly manufactured linear clock table, a=.01 -> .011 -> .012. Both
in-flight levels receive the factors through the actual link-selected clock
hook; all other fields stay unchanged. Passed; this is call-site evidence,
not a cosmological Friedmann solution or full AMR/reflux proof. Kept in
`.cosmo-coupling.g8Heia/dust-expansion-clock.log` with
`dust_expansion_clock_native`. The repaired awk postprocessor was additionally
run successfully on grammar-debug itself without rerunning RAMSES.

## Review disposition and further implementation (2026-09-12)

Fable recommends an analytic optically-thin CMB bath with nonnegative excess
IR, rather than storing the full CMB or adding signed radiation kernels.
Accept reuse as the lean design direction, but its recommended upper-table
temperature clipping and unledgered thermalization are NOT adopted. A table
limit is an applicability boundary, not permission to truncate Tcmb. Nor is
the claimed negligible energy fraction established for every AGB/wind/SN
source. Its suggested max(10K,Tcmb) retains an additional reference background
below z=2.667; it is not pure CMB at low redshift (10K corresponds to1+z=3.667).

The one-sided cold-grain issue remains mathematically real. The current
material backend lower-bound check is target >= U(bath)-Qgas (variable-speed
exchange mode). At zero primary, dust initially at bath and colder gas give
Qgas<0, so target=U(bath) violates it even though the true temperature deficit
may be small. Raising injected temperature does not remove this later case.
A one-sided thermostat must account for its external CMB energy, including
energy delivered to colder gas; it cannot simply loosen the energy tolerance.
Treat any density/depression estimates in the review as advisory, not measured
pilot limits. Also, cooling to a changing bath is a finite-time process; do not
require an arbitrary timestep to produce instantaneous exact equilibration.

Independent literature context: [da Cunha et al.2013](https://arxiv.org/html/1302.0844),
equations3--7, separates stellar heating and CMB absorption and obtains the
net modified-blackbody emission relative to the CMB. It motivates an analytic
background but does not prove our transient one-sided material solver valid
for colder gas, nor authorize energy-free projection onto the bath.

Implemented the unambiguous physical-unit IR dilution: production update_time
calls snrt_dust_live_expand ONCE outside the level loop. All slot-indexed
persistent radiation scales by (a_old/a_new)^3; material still scales with a^2.
No extra saved epoch/checkpoint field is introduced, avoiding restoration-time
ambiguity. Values restored at the current epoch are untouched until the next
real clock advance. This is explicitly the fixed-group/no-spectral-redshift
approximation, not full a^-4 radiation evolution. The same native clock test
uses an actual restored live IR slot and verifies two factors compose, invalid
ratios leave it unchanged, and a restored current-epoch field is unchanged by
a unit update. Passed. No cosmological dust admission guard was removed.

Pinned dilution build `.cosmo-coupling.g8Heia/ramses_cosmo_dilution3d` SHA256
`e5e137b3b4c11affd7614fe94a7b584d1cc59d9d18e589716eb8ead131079f58`;
test `dust_expansion_ir_native` SHA256
`beb4835ca55c63154cd7d29a7eec327bb89aec13283e68711ca297ebe35e0084`,
output `dust-expansion-ir.log` in the same directory. Native test and full
RAMSES link succeeded; no new cosmological RAMSES evolution was claimed.

The remaining operator-facing scientific choice is the analytic bath domain
and low-redshift floor, not permission for another routine test. Recommended:
declare an optically-thin analytic CMB bath, explicitly account for any energy
supplied to cold grains/gas, and use physical2.727/a for the cosmological bath
instead of silently retaining the noncosmological10K reference floor. Keep
noncosmological reference runs unchanged. This requires a declared effective
closure and appropriate cold-material support, not a table-temperature cap or
unrecorded injection-energy projection. Do not activate the cosmological dust
guard before this scientific choice and its actual implementation are settled.

## Approved CMB implementation — 2026-09-12

The operator approved the preceding choice. Cosmological C/silicate material
now reuses the native implicit enthalpy/gas-exchange solver with signed net
Planck power relative to physical `2.727/a`. Negative net emission is an
explicit external CMB receipt; only positive excess enters IR transport.
There is no instantaneous injection-temperature projection, hidden 10 K
floor, table-maximum CMB clamp, or c/chat-inflated stored background field.
This is an optically thin analytic bath, not CMB attenuation/distortion or
frequency-redshifting radiation transport. Noncosmological branches retain
their existing semantics. Cold grain decoding and mass-process enthalpy
updates use the same analytic mixture energy as the cosmological receiver.

The native material cases passed: 20 K grains heated by 272.7 K CMB;
equilibrium at 2.727 K below the old table boundary; CMB-supported transfer
to colder gas with grains slightly below the bath; warm-grain positive
excess radiation. A native IR/material/gas transaction using actual 136-band
quadrature closed its external-bath energy ledger to `2.6335e-13` relative;
rejected trials left gas, grain, radiation, photon and receipt state intact.
These manufactured opacity-basis cases test the operator, not the full
live D03 cosmological integration. The existing CPU/OpenMP backend suite
also passed, including absolute cold material, Fe, gas coupling and rollback.
Logs: `.cosmo-coupling.g8Heia/cmb-native.log` and
`cmb-backend-regression.log`. No raw simulation dumps were produced by them.

The live callback and rank/level-local committed `background_erg` receipt
are wired. Startup is narrowed to kind7, NENER=0, coadvected C/silicate
DL01/D03/CHIMES, v4 exchange-enabled IR; sinks, CR/SGS, MHD, Fe/PAH, drift,
sublimation and SN shocks remain outside this cosmological admission.
`mkrun.py` and the full namelist generator document the same restrictions.
Generator checks accept the manufactured cosmological input and reject
drift, PAH, CR and sink variations. This is not full production qualification.

Bounded integrated input: `.cmb-live.iE6Ld9/run.nml`, uniform 4^3 gas/DM,
Z=.02 and dust=.002, deliberately unrelated to the Z=1e-10 galaxy IC.
No stars, AGN or gravity; cosmological expansion, hydro, CHIMES and live dust
are the test. MPI2/OMP2, nstepmax=2, aexp_step_limit=1e-6; output schedule
noutput=1/aout=1.1/tout=1e100/foutput=fbackup=1000000, so no dumps expected.
ICs are preserved inputs, not cleanup targets.

Initial submission 539339 failed before RAMSES started because the driver
duplicated Makefile's `3d` suffix in the executable name: 1 s allocation,
4 allocated core-seconds. Submission 539344 stopped at namelist preflight:
an additional old noncosmological IR guard plus inherited stellar/AGN source
selectors. Its RAMSES process returned zero despite rejecting the namelist;
the wrapper correctly failed because no CMB commit existed. 6 s allocation,
24 allocated core-seconds; compute TotalCPU=9.568 s. Both failures are retained,
not counted as physical tests. Source selectors are unset for the source-free
test and the additional guard now follows the narrow material admission.

Submission 539357 uses binary
`.cosmo-coupling.g8Heia/ramses_cosmo_cmb_v3_3d` SHA256
`cfb8b5148e75dd48e4dbd17c2d20dcf47a7d38d1a2f23c3da28a97add64bfcf5`,
NML SHA256 `c3c145905ce8e3fd767051eada74763211e8ae17c4d12d5e2f95119ae2fa4de0`.
Submission 539357 completed successfully. Two actual evolution steps each
committed on both MPI ranks. Rank-local CMB receipts were
`2.240657158649e52` erg per rank in the first step and `1.397624473316e48`
erg per rank in the second. Both ranks own disjoint equal volumes here;
these are not duplicated global counters. Reported IR balance errors were
`6.0272e-14` and `3.0358e-13`, with zero IR escape. Main-loop econs is not
expected to be zero with cooling/expansion and an external thermal reservoir;
do not substitute it for the explicit material/IR conservation ledger.

Slurm job: COMPLETED 0:0, wall 6 s, 24 allocated core-seconds. Compute step:
wall 5 s, TotalCPU=11.009 s, per-task MaxRSS=401724 KiB; RAMSES evolution
elapsed 0.75308 s. The existing broad cooling timer (includes coupled work)
was mean0.735 s, Godunov0.014 s. This tiny, initialization-dominated case is
not a scaling benchmark. Log `.cmb-live.iE6Ld9/live-539357.log`; no output_*
directory was generated, so raw-output cleanup removes nothing. ICs retained.

Final native fixture pinned as `.cosmo-coupling.g8Heia/dust_cmb_native`, SHA256
`55e41bebaeca5bb1af1a44498e02259d63428428c4c0a66e735998ff673e438d`;
rerun `cmb-native-final.log` passed. Source diff whitespace check and both
namelist frontend syntax checks passed. No GPU run, cosmological stellar/AGN
source run, full128^3 cost measurement or galaxy calibration is claimed.

## 2026-09-25 live D03 CMB recheck (current uncommitted work)

An initial source inspection incorrectly treated the fixed 10 K reference
table background as the live cosmological bath. In git HEAD the live callback
actually computed `2.727/a` itself and passed it to the signed-net
`iron_radiative_cell` C/silicate receiver with Fe mass zero. The receiver
accepts baths below the first 5 K table knot via analytic DL01 cold enthalpy
and per-band Planck ratios. Native tests already covered a 2.727 K bath,
cold gas drawing energy from grains slightly below the bath, and the explicit
CMB receipt. The older integrated z=99 test did exercise that receiver, but
does not prove all later coupling changes or low-z live integration.
The two-size dust mass runtime explicitly checks that the sums of the first
two and last two size bins equal the C and silicate species carriers within
128 machine epsilons; thus the receiver's four-bin mass weights and the
driver's two-species enthalpy curve share the admitted composition.

A proposed replacement with the C++ material backend and an instantaneous
`U(Tcmb)` projection was a physical regression: that backend bounds the
temperature below by the bath and cannot account continuously for CMB
support of heat lost to colder gas. It was removed after a read-only Fable
physics review (Q-GOAL: restore signed-net receiver; Q-LEAN: no extra knot
or background rebasing). The current callback keeps the expanded interface
needed by the newer IR operator but derives physical bath independently of
the 10 K reference table, then uses the original signed-net C/silicate solve.
The native transport fixture now checks that the table bath remains 10 K
while the material receiver uses 272.7 K. With `SNRT=1 DUST_LIVE=1`, changed
modules compiled and the dust backend smoke passed. This is not yet a new
integrated pass.

A full-AMR smoke link in the reused `bin/` directory failed on mixed-option
HDF5/CHIMES/dust-dynamics objects, so no result from that link counts.
Isolated full-option `SNRT=1 DUST_LIVE=1 DUST_DYNAMICS=1 CHIMES=1 HDF5=1`
build plus the source-free cosmological 4^3 two-step fixture is Slurm job
403334. Its effective namelist is `.cmb-epoch-validate-20260925/run.nml`
(SHA256 `c3c145905ce8e3fd767051eada74763211e8ae17c4d12d5e2f95119ae2fa4de0`),
with no full dump expected. The job was held while correcting the regression
and released after the changed modules compiled. Inspect build log, source
hashes, CMB receipts, exit status, and absent `output_*` before changing the
integrated verdict. The originally dependent low-T preparation-only job
403338 was canceled before execution because it did not exercise the live
material solver. Replacement dependency 403359 runs the existing native
CMB material, IR-transport and expansion smoke after a successful 403334.
The in-flight
128^3 H200 job has a copied older binary and cannot validate this patch.
For the two-rank/two-step fixture, the script's minimum two CMB markers is
only a launch sanity check: final evaluation requires four rank-local
`SNRT_CMB_COMMIT` records, both `Main step=1` and `Main step=2`, finite
nonnegative receipts, explicit IR balance/rollback checks, and no full dump.

The first isolated build job 403334 became runnable after broadening its
eligible GPU partitions, but failed before source compilation: Lmod required
`intel/tbb`, `intel/compiler-rt`, and `intel/umf` before the Intel compiler.
Its dependent 403359 never ran and was canceled. Corrected scripts load all
three prerequisites, then Intel compiler/MPI; the sequence was checked on
the login node without running compute there. Replacement build/integration
job 403389 and dependent native smoke 403390 were submitted across
`a10,a40,a100,h100`. Their status and outputs supersede the pending wording
above; the module failure is not a physics or code-test result.

Job 403389 then reached compilation but failed in one second because a clean
`make -j4` raced Fortran module production: `pm_parameters` and
`hydro_parameters` read `amr_parameters.mod` before it existed. This is a
Make dependency/order error, not a source or physics verdict. Its build log
is retained as `build-403389-failed.log`; valid completed objects were left
in the isolated build directory. Dependent 403390 never ran and was canceled.
The corrected sequential build/integration is job 403394, with dependent
native smoke 403395. No previous full simulation outputs were deleted.

## 2026-09-25 resolved integrated CMB check and bounded MPI-cost correction

Jobs 403394 and 403395 completed `0:0` with the isolated full-option build.
The actual two-rank/two-step run emitted four positive CMB receipts, two
`SNRT_RT_CLOSURE_PASS` and `SNRT_DUST_IR_COMMIT_PASS` records, both main steps,
and `Run completed`. IR relative balances were `6.0272e-14` and `3.0358e-13`.
The CMB material, transport/rollback, dilution and expansion native fixtures
also passed. Source hashes were identical before/after the run; no `output_*`
tree was created. This validates the restored signed-net live receiver for
this manufactured source-free z=99 start, not an approved stellar/AGN SED or
all-redshift science production run.

This tiny run exposed an MPI efficiency failure: 80 directions times 136 IR
groups were exchanged in 16-component halo tiles, causing about 53 s of halo
wall per step. `snrt_dust_live.f90` now enlarges tiles only where the field is
small, with at most 8 MiB additional scratch over the previous 16-column
tile. The same effective NML (SHA256
`c3c145905ce8e3fd767051eada74763211e8ae17c4d12d5e2f95119ae2fa4de0`)
in `.cmb-halo-validate-20260925/` gave 1.105/1.115 s halo times. Its
integrated simulation finished both steps and reproduced the old printed
CMB receipts, RT closure, IR balances and main-step values exactly. Batch
403477 ended `FAILED 1:0` only because its post-run native fixture was
invoked without the MPI launcher. A first native-only retry lacked the dust
contract environment (`STOP 9`); corrected job 403519 passed all native
checks. Neither failure is a simulation/physics rejection.

The remaining ~329 s/step came from `snrt_regrid_upload`: it reversed each
component of the combined photon/IR parent payload in a separate MPI call.
It now packs the same per-component parent densities into bounded tiles and
uses the existing reverse-add AMR halo exchange. The parent multiplicity,
collective validation and commit rules remain unchanged. Its Makefile module
dependency was made explicit. Fresh full-option job 403527 completed `0:0`;
both steps and native fixtures passed. The two-step RAMSES wall fell from
783.58 s (403394) to 657.66 s (halo change only, 403477) to 48.18 s
(both changes, 403527). This is an identical 4^3 fixture comparison on
different GPU nodes, not a large-run speedup or an accelerator claim. Source
hashes were identical before/after job 403527. The old and final S_N logs
match exactly at printed precision for CMB receipts, IR balances, RT closure,
and main-step conservation values. Effective NML:
`/gpfs/kjhan/LRD_JWST/.cmb-regrid-validate-20260925/run.nml`.

The same final binary (SHA256
`6c8e8180517d7464345955127d891d9c38328780a82d39329e15e68d14d21543`)
also passed a two-rank/two-step M5 run, job 403540 `0:0`, with only
`snrt_transport_model='tensor_mn', snrt_moment_order=5` changed in its NML.
Both M5 commits have zero printed photon/energy residual, main-loop steps
complete with the same printed `mcons`/`econs` values, and RAMSES wall is
39.52 s. Effective M5 NML:
`/gpfs/kjhan/LRD_JWST/.cmb-m5-validate-20260925/run.nml` (SHA256
`6e4fea14dbca7ef282d04a5cf81c4887f3152faac6c7d2d29966069c0c4056c8`).
All three evaluated comparison directories had zero raw `output_*` dumps;
no deletion was necessary. These are wiring/conservation and tiny-cost
checks only. Reference-control source/group/dust inputs and larger evolving
AMR/illuminated runs remain separate production-readiness limitations.

The pre-existing `snrt_regrid_native` fixture initially stopped before any
upload check (`STOP 2`): its old +/-mean/32 signed test shift violates the
current group's energy-edge admission for some bands. The fixture alone was
reduced to +/-mean/1000, retaining both signs and all existing assertions.
Job 403560 then passed `SNRT_NATIVE_LEVEL_UPLOAD_PASS`,
`SNRT_NATIVE_REGRID_PASS`, the coarse/fine IR conservation/rollback test,
and the defragmentation checks. No production acceptance range or transport
arithmetic was relaxed for this fixture repair. The preceding failed native
job 403553 and its log remain in `.cmb-regrid-validate-20260925/`.
