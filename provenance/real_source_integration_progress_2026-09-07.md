# Real-source integration: implementation progress

First production-completion bundle, still **in progress**. No bundle-end audit
or comprehensive validation is claimed. Scope and updated four-bundle sequence
are in `production_completion_bundles_2026-09-07.md`.

## Implemented

The active `patch/lagRamses/init_sink.f90` previously located `ic_sink` and
`ic_sink_restart` but never read either file. It now reads on rank 1 and
broadcasts each row, checks positive finite mass, finite state, centered
in-box positions and capacity, initializes sink reservoirs to zero, and
appends a canonical PTYPE_SINK particle on the position's owner rank before
init_tree. The normal sink maintenance then builds the accretion cloud.
This adds ordinary IC loading, not a driver-only radiation seed.

The layout follows the repository's 12-column sink IC convention: code mass,
centered position (box center is zero), velocity, gas angular momentum, and
two upstream-only SMBH-mass/drag fields. The last two must be zero because
this patch uses a single BH mass and has no corresponding drag field. Such
values are rejected rather than discarded. Blank and full comment lines are
accepted. Missing/nonfinite fields and capacity overflow are rejected.
No namelist key or generator schema changed in this increment.

## Short native execution

- Project `/gpfs/kjhan/LRD_JWST`, base `78add24`, active lagRamses VPATH.
- Build: `make -C bin -j1 SNRT=1 USE_CUDA=1 USE_FFTW=0 ramses`.
- Binary SHA256: `79de2f3d77125422a15bde94f81a4c93c55df94edc6825378640a76bed1545ff`.
- Run: `.native-source.zxM39a/effective.nml`, `ic_sink`, and `ramses.log`.
- Reusable inputs: `simulation/snrt/config/snrt_agn_accretion_smoke.nml` and
  `snrt_agn_accretion_smoke.ic_sink` (copy the latter as `ic_sink`).
- 512 uniform cells; one initial BH; CPU hydro/Poisson, CUDA primary RT;
  OMP_NUM_THREADS=1; nstepmax=3; no synthetic driver seed or failure injection.
- SNRT_RT_ENABLE=1, SNRT_AGN_MODEL=partition_reference_v1,
  SNRT_REDUCED_C=0.01, SNRT_RT_LEVEL=3, SNRT_ALLOW_REFERENCE_CONTROL=1;
  group and secondary contracts are the existing reference-control files.
- Output policy: noutput=1, aout=2, tout=1e30,
  foutput=fbackup=1000000. Expected and observed dumps: zero.
  Free space at launch: 171 TiB.

Observed `Sink IC loaded: added=1`, actual grow_bondi/AGN maintenance,
`active sources: 1` on the second RT call, transaction commit and closure,
and `Run completed` at coarse step 3, exit 0. Radiation was produced from
accepted accretion, not SNRT_DRIVER_TEST_SEED_SOURCE. Reported final mcons
was 4.34e-16. The generic hydro econs diagnostic becomes 1.00 after feedback;
it is **not** a source-corrected total energy accounting check and is not
claimed as an energy conservation pass here.

## Still in this first bundle

Stellar-source integration, combined nonzero live dust, and restart storage of
unconsumed AGN energy remain unfinished. This run uses a reference spectrum,
feedback_mode=legacy with no stars, and dust ZERO_SCAFFOLD. It does not qualify
stellar/dust physics, MPI source ownership, sink formation or restart.
The existing SNRT restart/MPI guards remain intact until their prerequisites
are actually implemented. Physical input completion belongs to bundle 2;
automatic GPU/OpenMP allocation is now bundle 3 by explicit user instruction.

## Follow-up: AGN HDF5 persistence implemented

The earlier restart limitation above is superseded **for serial HDF5** by
this increment. In `/sinks`, `agn_state_schema=1` binds the saved AGN model
and RT ownership switch. New per-sink datasets persist radiation, thermal,
jet and deferred energies in erg and retained loading mass in code units.
The reader requires all datasets and validates finite nonnegative values
before publishing the reservoirs. Missing versioned state is accepted only
for a legacy/no-live-RT input; live/reference operation cannot initialize a
missing ledger to zero. Schema/model/RT mismatch is a terminal error.

Native admission now permits serial HDF5 restart with a successfully restored
AGN ledger and the existing radiation checkpoint. MPI source ownership is
still rejected. `snrt_agn_reference_config_ok` admits nonnegative restart
indices, while the actual read_params code enforces HDF5 format/build and
the driver checks that the ledger really was restored. No legacy binary
sink layout was changed, and no claim of SNRT binary-format restart is made.

The shared namelist database now exposes the already existing `informat` key
in OUTPUT_PARAMS; mkrun reports the chosen restart format and live-AGN HDF5
requirements. Unrelated local database removals remain untouched.

### Build and short execution

Isolated complete build: `.agn-restart.AD5PAd`, copied current Makefile with
the same relative VPATH, `make -j1 SNRT=1 HDF5=1 USE_CUDA=1 USE_FFTW=0 ramses`.
No DUST_LIVE macro, NVAR=18. The normal bin executable was not overwritten.
Binary SHA256:
`9e1e2bd0e83af8aa466df6ada2c4050829e40f7a102158bbc6c4e7e99e4e7e5f`.

Runs under that directory retain effective.nml, ramses.log and snapshots:

- `continuous`: actual initial BH accretion, 3 steps, checkpoint each step.
- `resumed`: resume checkpoint 1 and finish step 3; nonzero AGN source occurs
  after restart just as in the continuous run.
- `resumed-pending`: resume checkpoint 2, which contains thermal pending
  energy **7.319992397631857e53 erg**, and finish step 3.
- `invalid-ledger`: disposable copy of checkpoint 2 with radiation pending
  energy set to -1; rejects before time integration, exit 1.
- `invalid-schema`: disposable copy with schema set to 99; rejects before
  time integration, exit 1. Original snapshots were not modified.

Both successful resumptions contain `AGN checkpoint ledger restored: sinks=1`
and `Run completed`, exit 0. All five final reservoirs equal the continuous
run exactly, including deferred energy **1.2435327972056231e53 erg**. Other
four final reservoirs are zero. The measured positive checkpoint reservoir
is thermal; this run does not separately exercise positive radiation, jet or
loading-mass checkpoint values (radiation is consumed before each dump).

At checkpoint 3, resume-from-1 vs continuous: max hydro absolute difference
2.168404344971009e-19, max dataset-scaled difference 9.828483696493675e-14;
max SNRT absolute difference 2.220446049250313e-16. Resume-from-2: hydro
absolute 1.3658406274475593e-20, dataset-scaled 1.9810537450744305e-14;
SNRT state exactly equal. These are narrow restart regressions, not a full
physical energy-closure or production qualification.

Output audit: noutput=1, aout=2, tout=1e30; foutput=1, fbackup=1000000,
nstepmax=3. There are 3 continuous dumps, 2 new resumed dumps and 1 new
resumed-pending dump; test copies retain their input snapshots. Each HDF5
file is 3,783,664 bytes; all five run directories occupy about 37 MiB.
Available space at launch was 171 TiB. Environment is the same as the
actual-accretion run above, with driver seed/failure injection/dust contract
unset. No other simulations were launched.

Existing native efficiency/source/deposition regression
`simulation/snrt/tests/run_fp15_agn_efficiency.sh` passed. Shared generator/GUI
tests passed (21 tests, one display skip), including HDF5 informat round-trip.
No new audit gate was introduced. First-bundle remaining work is still actual
stellar-source integration and simultaneous live dust; MPI ownership is not
enabled by this HDF5-only change.

## Continuation through the four implementation bundles

The user preapproved continued work through bundle 4. No new audit gates or
external model audits were launched. Work remains in the same `/gpfs` repository.
New isolated build/run root: `/gpfs/kjhan/LRD_JWST/.production-native.YCek9O`.

### Native changes

- Added a native stellar photon-table consumer, separate from mechanical SN
  energy. It uses actual STAR particles' initial mass, proper birth epoch and
  metallicity, with piecewise analytic age integration and linear Z interpolation.
  IMF type, mass limits, population type/fraction and common transport closure
  must agree. No age/Z extrapolation or physical BPASS promotion is performed.
- Captured each recursive AMR level's starting proper time for the source
  interval. HDF5 binds the enabled switch and exact photon table; MPI startup
  rejects rank-dependent tables. Existing particle HDF5 fields retain feedback
  progress independently of the new radiation source.
- Actual newly formed stars exposed an existing locator defect: old SoA cell
  arithmetic did not match the active blocked layout, and INT on a negative
  fractional grid coordinate falsely assigned a source to an adjacent grid/rank.
  The locator now uses ICELL_OF and half-open coordinates. Stellar grid hints
  avoid a whole-mesh scan for every particle.
- Implemented OpenMP primary transport and automatic CUDA/OpenMP placement.
  CUDA and CPU share the same species/dust cap implementation. GPU selection
  accounts for node-local rank/UUID ownership, free memory and workload; host
  team sizing respects explicit OMP settings or local allocation. No fallback
  is attempted after GPU execution fails. This remains a CUDA-linked build
  even when no physical GPU is present.
- Added MPI photon owner-count publication so accepted sink receipts are
  consumed on every replica only after the global coupled commit; duplicate
  spatial owners fail closed. Existing mechanical receipt publication is reused.
- Bounded 16-component primary and same-level IR halo messages reuse the
  existing grid maps with private SNRT buffers. The old scalar path performed
  720 primary exchanges per substep. Generic AMR/tree/CPU-box code was not edited.
  Coarse/fine physics and the existing reverse-flux arithmetic are unchanged.

### Evidence obtained so far

Full NVAR=30 build: `SNRT=1 DUST_LIVE=1 HDF5=1 USE_CUDA=1 USE_FFTW=0`.
Current executable SHA256 (after primary/IR tiled exchange):
`8535510b92ae55ffa9398f7f175e495235e3d5b337a3d2525ff5489c6fe04910`.
Normal `bin/ramses_final3d` was not replaced.

- Native CPU/CUDA two-opacity-case comparison passes; max difference normalized
  to the initial photon/atom inventory is 1.37311e-7. Both ledgers close and bad
  input leaves CPU photon/atom state unchanged. The existing CUDA regression,
  including zero-dust bitwise behavior and negative cases, passes unchanged.
- Native stellar integral/split-step/domain-bound tests pass. The same small
  native test checks blocked cell indices, a negative-coordinate rejection,
  shared-face ownership and the stellar grid hint.
- Pre-locator-fix baseline directories `cuda`, `openmp`, `auto-hidden` all
  completed three steps. CPU vs CUDA hydro maximum relative dataset difference
  was 4.4507e-11 and SNRT maximum absolute difference 2.2352e-8. These are backend
  evidence, NOT validation of the old spatial deposition, which was subsequently
  corrected.
- After the locator correction, `stellar-native-fixed` and `stellar-final`
  completed with active sources 0, 513, 961 and live IR commits. RAMSES itself
  formed the stars; no synthetic particle insertion was used. Stellar feedback
  channels were disabled in this wiring control; it does NOT certify newly
  developed stellar mechanical feedback. The synthetic SSP masses/rates are
  numerical controls, not a physical galaxy population.
- `stellar-resume-openmp` resumed a copy of `stellar-final/output_00002` with
  GPUs hidden and completed the final step. Both runs retain 961 stars and
  byte-identical numerical SED identity. Jet liability 1.174787943869677e50 erg
  and loading liability 2.4120782632610393e-11 code mass match exactly; the other
  three pending ledgers are zero. Hydro energy relative difference 3.40323e-8,
  dust-energy relative difference 8.31116e-8, SNRT maximum absolute difference
  5.83716e-6. This is successful execution across a backend switch, not a claim
  of bitwise invariance or an approved physical error envelope.
- Existing AGN native helper/source/deposition regressions pass. Shared
  namelist/GUI tests: 21 total, one display skip, no failures.

MPI debugging runs and their logs were preserved. A long scalar-exchange run
was interrupted while both ranks were inside halo communication; this did not
establish a deadlock. A shorter pre-locator-fix run reached the next step and
correctly rejected duplicate AGN ownership. Do not label those runs successful.
The `mpi-tiled` two-rank run subsequently failed at the third live-IR step;
its first two coupled commits are not a complete-run pass. The precise
reproduction and bounded numerical correction are recorded below.

Every new run used its own absolute `.../<case>/run.nml`. Evolution controls:
noutput=1, aout=2, tout=1e30, foutput=1, fbackup=1000000, nstepmax=3; the resumed
case starts at nrestart=2. Typical dump size is 4--5 MiB. Including retained input
copies and planned final checks, the output budget is 150 MiB; free space was
171 TiB. No pre-existing user snapshots were removed.

### Physical input boundary, not a new audit gate

The attempted channel-resolved profile requires an admitted terminal-fate
map/yield package. The currently selected profile is `review_only_unresolved`;
its rejection must not be bypassed by inventing an approval string. The staged
BPASS binary imf135_300 candidate is not the project's Chabrier 0.08--120 Msun
population, and the native stellar control is not a replacement physical SED.
Common stellar/AGN spectral closure and physical dust input admission remain
unresolved. The synthetic g1 yield fixture was never promoted to production.
Previously parked source-data questions are not being replaced with another
test framework. Bundle 2 and the full physical coupling in bundles 1/4 cannot
be marked production-complete on this evidence.

## Latest closeout increment: source format, AMR and near-bath MPI state

All paths below are relative to `.production-native.YCek9O` unless noted.
Final native build `ramses_ready_final3d`, NVAR=30, the same options above:
SHA256 `5ad0bc6407fdf57f95cdc1b851c270b59d6115a5f8a3e7f5aaacdec2b8863e97`.
Build log: `build-ready-final.log`. This supersedes earlier binary identities,
without relabeling their results as executions of this later binary.

Stellar-enabled checkpoint versions are now 3 (primary) and 4 (primary+IR).
Versions 1/2 remain no-stellar formats. Experimental earlier stellar snapshots
in `stellar-final` used version 2 and MUST NOT be resumed with the final
executable. Fresh `ready-star` uses version 4 and completed three steps with
961 stars and live IR. `ready-star-openmp-restart` resumed its checkpoint 2
with GPUs hidden, auto selected OpenMP, and completed the final step. SED
identity and all five AGN pending reservoirs match exactly. Hydro energy
dataset-scaled difference is 3.4032297086747095e-8, dust energy 8.311157656483656e-8,
SNRT maximum absolute difference 5.83716335498391e-6. These runs preceded the
near-bath correction; their binaries are `ramses_ready3d` and
`ramses_diagnose3d`, respectively, not `ramses_ready_final3d`.

`amr-ready` used the existing level-3/4 IR fixture, reducing its Courant factor
to .005, with `ramses_diagnose3d`. It completed two steps with 296 level-3 and
1728 level-4 leaves and six coupled IR commits; maximum reported IR balance
was 1.4660e-12. No stars/AGN were enabled in this specific AMR regression.
One 17 MiB snapshot was written. This exercises the tiled exchange's serial
coarse/fine path, not a joint MPI+stellar+AMR physical run.

`mpi-tiled` failed at step 3 with `dust_err_state=2`. Read-only inspection of
checkpoint 2 found positive material fields and a temperature minimum only
7.1e-15 K below 10 K. After the next hydro advance, the reproduced IR input in
`mpi-resume-precision` was **9.9999999979853200 K** on rank 2: approximately
2.01468e-10 relative below the bath, outside the old 64-ULP admission band.
The guard failure was therefore real, not a halo deadlock; no failed run is
counted as `Run completed` merely because `clean_stop` returned exit 0.

The native transient solver now admits a below-bath relative material-energy
deficit only within its existing solve tolerance (1e-9 in live IR, with a
64-ULP numerical minimum). It does NOT overwrite the incoming material energy:
the entire resulting floor adjustment remains in the material+radiation
closure residual and must satisfy the same tolerance. This is a bounded
numerical correction, not a new physical background heat source or a claim
of exact energy conservation. Larger deficits still reject atomically.
The native regression explicitly checks a charged 2e-12 correction at 1e-10
tolerance and rejection of a 2e-9 deficit. gfortran and ifx both pass the
existing differential/transient/halo suite in `dust-native-regression-final.log`.
The initial attempt used a nonexistent root `.venv` and did not run; the
successful invocation uses the existing `simulation/snrt/.venv/bin/python`.

`mpi-resume-final` uses the final binary and a fresh copy of `mpi-tiled`'s
checkpoint 2, with nrestart=2, nstepmax=3. The previously failing coupled step
now commits (IR relative balance 2.7213e-15); the run's final completion marker
is recorded separately below. AGN injection in this last step is zero; unique
nonzero source ownership was exercised in `mpi-tiled`'s earlier step, not
invented as an additional positive-source claim for this restart.

The actual-star reference input is now reusable at
`simulation/snrt/config/snrt_stellar_agn_dust_reference_smoke.nml` with the
environment/output policy in `simulation/snrt/NATIVE_RUNTIME.md`. Stellar
mechanical channels remain explicitly OFF. The reference SSP/opacity table
does not become physical by being stored in the repository. No new Python
test framework, external audit, production run, commit or push was performed
in this increment. Existing unrelated generator deletions were preserved.

Final two-rank result: `mpi-resume-final/run.log` contains `Run completed`
after writing `output_00003/data_00003.h5`; mpiexec exit 0. Reported total is
255.482 s, including 178.436 s in the broad legacy `cooling` timer (which
contains coupled RT work); this is not a tuned performance claim. Checkpoint
header reports two CPUs and coarse step 3. Density, dust mass/energy, SNRT
state and all five AGN pending reservoirs were inspected for finite,
nonnegative values. Previous failed reproduction outputs remain intact.

## Continued input wiring; two scientific data issues deferred again

The user explicitly deferred the 40--120 Msun fate/yield discrepancy and
KL16/CK22 AGB source-data issues. No data correction, author contact, scientific
approval or renewed source audit was attempted in this increment.

`stellar_snia_population_contract` now checks a valid handoff against the
configured IMF, population and binary fraction. `stellar_ramses_runtime` calls
it before completing SNIa initialization (mismatch error 41). The approved
FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ Kroupa/binary baseline remains approved;
it cannot silently supply default Chabrier/single-star particles. The existing
native fixture adds five population binding checks. The ifx bounds-checked
fixture passes all 13 checks with `PHASE0_SNIA_RUNTIME_CONTRACT` set to the
existing `fp2_snia_runtime_contract_v1.nml`. A later attempted rerun supplied a
positional argument instead of this environment variable and exited at the
missing-input check; it was rerun with the correct environment. This is not
a nonzero physical SNIa evolution run; its existing WD-reservoir requirement
is unchanged.

The full NVAR=30 SNRT/DUST_LIVE/HDF5/CUDA build succeeded in
`.production-native.YCek9O/build-input-binding.log`. New executable:
`.production-native.YCek9O/ramses_input_binding3d`, SHA256
`8cba5cc187cfecf652a1cb96a996b9306b9459772424bb83cf48c70269b444ed`.

The existing Draine thermal builder now exports native v3 input. It uses raw
table `external/draine_wd01_rv31/kext_albedo_WD_MW_3.1_60_D03.all`, SHA256
`b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c`,
with mass per H 1.398e-26 g and the existing pilot AGN ledger (SHA256
`d2326e4ec0320e9e0b34a58460b08bddafa4b4433322691a9dd346a49762909d`).
Primary opacity is monochromatic at each ledger representative energy; the
full raw-domain IR quadrature has 136 nodes, with 65 temperatures including
the 10 K bath. Heat capacity 1e-24 erg/H/K is explicitly a TEST constant;
all native approval fields remain reference-only. Generated input is
`.physical-inputs.A0XtMt/dust_native.nml`, SHA256
`eb26f979f722e8e92a676d00ba607211cf87c1921cb45ba2fcadb696e91fb531`.
Regeneration matches exactly and nonpositive/nonfinite heat capacity rejects.

The initial no-PIC run `.physical-inputs.A0XtMt/draine-live` committed two dust
updates but failed with SIGSEGV in `output_header`: default
`sf_birth_properties=true` accesses `elem_list`, allocated only in particle
initialization. No generic output-system rewrite was made. The no-particle
dust fixture now explicitly disables that particle-only output. Its original
failed directory is preserved.

Fresh run `.physical-inputs.A0XtMt/draine-live-final.IEstJu` uses that setting,
CPU hydro, automatic primary SNRT placement (selected CUDA), fixed level 3,
512 cells and Courant factor .005. Initial dust temperature is 20 K for the
test capacity. Full effective input is `run.nml` in that directory;
`nstepmax=2,noutput=1,aout=2,tout=1e30,foutput=2,fbackup=1000000`.
Before launch, 170 TiB were free and one dump below 60 MB was budgeted.
The run completed in 2.54 s, exit 0, with two IR commits (maximum reported
balance 7.0304e-13). One 54,430,120-byte HDF5 checkpoint plus metadata and
COMPLETE marker was written. Leaf dust mass and energy and SNRT state are
finite/nonnegative; the final derived dust temperature is
10.000000000000004 K. This checks live cooling and output with real optical
data, NOT nonzero AGN absorption, physical heat capacity, or production approval.

The previously omitted `informat` generator field is selectively staged
without including the user's unrelated 88 deleted lines. Existing GUI tests
reported 20 passes and one display skip in
`.physical-inputs.A0XtMt/generator-tests.log`. No external audit, new test
framework or production-scale calculation was launched.

## High-mass source history -> live mechanical feedback -> restart (2026-09-07)

User instruction: implement the connection behind the previously added model
choices. Work remained in `/gpfs/kjhan/LRD_JWST`, LagRamses main at
37c05b5a8be1976b6f3a9aa9122b08dddd55c1f2 plus uncommitted changes. No new
external audit, commit or push. Existing unrelated generator deletions and
its selectively staged `informat` addition were preserved.

Native implementation now binds explicit source nodes, IMF/population, Z and
terminal lifetimes to all three endpoint presets. The user-selected route is
distinct from scientific approval. It requires single-star SSP/binary=0,
wind+SNII, SNIa/PISN off, and HDF5 particles. Canonical raw values are retained;
the resolved remnant, source labels/settings and actual source contents are
saved/compared on restart and checked collectively across MPI ranks. Wind is
piecewise cumulative-linear; terminal release is a lifetime step. Mass-cell
budgets scale by M/M_node with shared analytic IMF mass weights. The exact-Z
restriction is deliberate, apart from relative 32-epsilon roundoff matching.
See `simulation/snrt/NATIVE_RUNTIME.md` for input schema and limitations.

Full NVAR=30 SNRT/DUST_LIVE/HDF5/CUDA, native enrichment build:
`.production-native.YCek9O/ramses_high_mass_live_final3d`, SHA256
`1d0549a4133ab96687465a722c4495034b9bb01ab9a531359d98f8ab12d3e732`.
Build log `.production-native.YCek9O/build-high-mass-live-final4.log`.
Makefile VPATH precedence was unchanged; native runtime object dependencies
were added to particle initialization and HDF5 backup/restore.

The 24-row synthetic source has M=40,60,120 Msun, Z=.01, lifetime=1000 yr,
wind fraction .1, terminal fraction .2 and raw remnant fraction .7. This is
engineering input, NOT physical stellar yields. Source SHA256
`7c2b7c7f245eb33555782a3885646f8f151ca617470e24924f11007c56003ed2`;
history SHA256 `85dca14b211304050c48c7de8e6c99afd2ff445aa08f48f7c24e0c55da23c2d3`.
The `high_mass_history_test.f90` fixture passed with bounds-checked gfortran
and actual optimized mpiifx production objects. It covers all three presets,
before/at/after lifetime, split increments, no repeat release, mass/element/
energy closure, mixed correction and bad input rejection.

All run directories below are under `.high-mass-live.aFGTS4/`; effective
namelists are their `run.nml`. CPU hydro/Poisson, fixed level 3 (512 cells),
natural star formation, no AGN, no active RT/dust evolution. Output schedule
was reported before launch: nstepmax=4, noutput=1, aout=2, tout=1e30,
foutput=1, fbackup=1000000; 170 TB free. Fresh 4 dumps were budgeted <=30 MB,
restart/MPI2 additional outputs <=60 MB. Actual fresh dumps total 1,622,604
bytes including metadata.

- `fresh` and `reproduce.bhweB3` preserve an initial failed run and diagnostic
  restart. Newly formed stars had Z=.010000000000000002--.010000000000000004;
  strict floating equality rejected source queries. Matching now admits only
  relative roundoff, not Z interpolation; out-of-domain tests still reject.
- `final.mhyCBb`: wind-only preset, four steps, exit 0 / Run completed, 2048
  stars. Oldest 512 stars returned 0.00716509 of birth mass once; their mass
  remained fixed in later steps. Final summed stellar loss is
  6.389242644428476e-9 code mass. All saved hydro values are finite, leaf gas
  density/internal energy positive, gas+stellar mass = initial .001 code mass
  to reported FP64 precision. The ordinary screen `mcons` is gas-only here;
  the noncosmological SFRD diagnostic prints NaN, but saved physical fields
  and the explicit gas+star budget are finite. This diagnostic was not patched.
- `restart.mcsFEY`: copied checkpoint 2 into a new directory, resumed to step
  4; exit 0 and STELLAR_SOURCE_RESTART_IDENTITY_PASS. Sorted particle identity,
  mass, mp0, Z and source vector match exactly. Feedback age/birth-time maximum
  difference is 2.7755575615628914e-17 code time. Leaf density, total energy,
  metal, H, He and Fe match exactly. No duplicated mass return.
- `reject-preset.I3SFh9`: valid source-consistent preset replacing wind-only
  on restart rejects at the source identity check, exit 1, no new output.
- `reject-data.mMNSbA`: doubled non-bulk energy in an otherwise valid source
  table rejects at the same check, exit 1, no new output. Checks bind consumed
  numbers, not only user labels.
- `mpi2.JaLAxG`: identical fresh profile on 2 MPI ranks x 2 OMP threads;
  exit 0 / Run completed, 2048 stars, same total fractional mass loss
  .005374941907105354 and conserved gas+star total. Approximately 115 seconds,
  mostly parallel HDF5 I/O; this was a functionality check, not a speed claim.

The reusable input is `simulation/snrt/config/high_mass_feedback_reference_smoke.nml`.
Shared generator/CLI/GUI now expose history path and channel switches, reject
unsupported combinations, and preserve quoted absolute paths on import (the
old regex terminated a group at a slash inside a string). The GUI/setup suite
passes 22 tests with one display skip, log `generator-tests.log` in the test
root. These changes do not certify actual source authenticity or combine the
new single-star model with the separate approved binary/SNIa model. Physical
source choices and the arbitrary-Z prescription remain explicit follow-up
work; no missing data or approval was fabricated.

## SNIa + metallicity interpolation + real event input (2026-09-07)

User requested the new route's SNIa combination, general metallicity
interpolation and physical-data connection. Implemented within the existing
physical-input bundle; no new audit or broad infrastructure work.

- `user_selected_model_v1` now accepts an explicitly matched effective binary
  SSP with SNIa and an AGB WD supplier, in addition to its original single-star
  route. The actual history and SNIa contract must agree on IMF/population/
  binary fraction. The approved Kroupa/binary=.5 N100/Maoz contract is reused;
  default Chabrier/single-star is unchanged. This is not a binary-evolution
  solver and does not relabel single-star yield data as physical binary yields.
- The generic canonical-table grid/coverage audit no longer requires fake
  per-star SNIa rows. The existing DTD and physical event budget own SNIa.
  Only AGB remnants (upper source mass <=8 Msun) fund its WD debit; the native
  test proves NS/BH remnants are unchanged and a zero WD reservoir rejects.
  Runtime accepts only the implemented constant-unity SNIa Z-factor identity.
- History namelist `metallicity_policy` defaults to `exact_nodes`; optional
  `linear_Z_cumulative_mixture` evaluates bracketing source histories at the
  same age and combines resolved cumulative budgets with fixed birth-Z
  weights. No interpolation of individual fates/lifetimes, extrapolation,
  invented primordial yields or cross-engine mixing is introduced.
- Low-mass payloads in the selected source-cell mode scale by M/M_node,
  consistent with high-mass per-initial-mass fractions and IMF integration.
  A split-step test exposed a plateau roundoff artifact in generic AGB energy
  interpolation; identical age payloads now stay exactly constant. Legacy
  tables outside the selected path retain their previous interpolation.
- HDF5 identity adds the Z policy, changed low-mass algorithm semantics and
  all consumed SNIa numerical/label fields (including DTD normalization, event
  energy, yields and thermal coupling). Old exact-Z/high-mass-only/SNIa-off
  identity is retained; a checkpoint using changed low-mass semantics cannot
  silently migrate. Generator validation and the existing CLI/GUI editor
  support the combined population without changing the Chabrier default.

Inputs: `high_mass_snia_z_yields.dat` and `high_mass_snia_z_history.nml` in
`simulation/snrt/tests/fixtures/phase0/` are deliberately SYNTHETIC wind/SNII/
AGB at Z=0,.02; their early WD formation is not physical. The SNIa input is
the actual previously approved `config/fp2_snia_runtime_contract_v1.nml`,
unchanged physically: 1.4004633930489443 Msun and 1.5063100005966762e51 erg per
event, the selected HESMA N100 ejecta, and the Maoz DTD. Its two ampersand
line continuations were removed: Intel accepted them, GNU namelist input did
not. No numbers, approval identifiers or raw HESMA source bytes changed.
New runtime-namelist SHA256:
`e56b1644bf4e08995279a807aad29bde9151b985fe854e443870105453b5b2f0`.

`high_mass_snia_z_test.f90` passes against optimized mpiifx production objects
and GNU bounds-checked modules. It verifies Z mixtures at different lifetimes
and outcomes, out-of-hull rejection, split-step source/DTD equivalence, real
SNIa mass/energy, AGB-only debit and mass closure. Existing three-preset
`high_mass_history_test.f90` still passes. GUI/setup: 22 pass, one display
skip; log `.high-mass-snia-z.dRgxmT/gui-tests.log`. Diff check passes.

Full native binary `.production-native.YCek9O/ramses_snia_z_final3d`, SHA256
`52ba7dbbdf95e155d4eeca80aa3376d2019af54e62f99a997d89ad6f1c57280f`;
build log `build-snia-z-final.log` in the same directory. Flags remain
NVAR=30/SNRT/DUST_LIVE/HDF5/CUDA/native enrichment; VPATH precedence unchanged.

Live evidence is under `.high-mass-snia-z.dRgxmT/`, each input `run.nml`:

- `fresh`: initial successful 4-step run, preserved with its pre-final identity
  extension build. It is not the checkpoint used by the final restart test.
- `final`: four coarse steps, 512 fixed level-3 cells, 2048 naturally formed
  stars, CPU hydro/Poisson, no AGN or active RT/dust evolution. Timescale is
  approximately 1 Gyr/code unit so stellar ages cross the real 40 Myr DTD
  delay. All saved hydro values finite. Birth Z increases from .01 to
  .010000209115435404 inside the supplied Z bracket. Final gas+star mass is
  .001 code units, equal to initial mass to displayed precision. SNIa alone
  returns 2.308207082109318e-10 code mass. Independent analytic Kroupa mass
  fractions plus the power-law DTD reproduce per-particle loss with maximum
  residual 1.3414421921102262e-25 code mass. Screen SFRD is NaN/Infinity in
  this noncosmological control; physical arrays are finite (no unrelated fix).
- `mpi2`: same physics, two MPI ranks x two OpenMP threads, four steps, one
  final output; same total mass and SNIa return (last-digit roundoff only).
  Run completed, exit 0, about 44--46 s including parallel I/O.
- `restart`: copy final's checkpoint 2 into a new directory; resume to step 4.
  Run completed and STELLAR_SOURCE_RESTART_IDENTITY_PASS. Sorted particle
  IDs/mass/mp0/Z/feedback age/birth time and source vector match exactly.
  Leaf gas density, energy, metal, H, He and Fe also match exactly.
- `reject-dtd`: double the otherwise valid event-rate normalization in a
  separate input copy. Identity rejection before evolution, exit 1, no new dump.
- `reject-z`: change only the history policy to exact nodes. Same rejection,
  exit 1, no new dump. Neither test edits approved source files.

Launch reports checked the effective inputs and 170 TB free space. All runs
use noutput=1,aout=2,tout=1e30,fbackup=1000000; final/restart use foutput=1,
MPI2 uses foutput=4. Budgets were 30 MB initial, 40 MB final+MPI2, 20 MB for
restart. Actual final dumps total 1,710,852 bytes; MPI2 dump is 513,786 bytes.
Reusable profile: `simulation/snrt/config/high_mass_snia_z_reference_smoke.nml`.

**Partial physical-input completion, not an all-channel production release:**
real SNIa is connected and exercised. Fully physical LC18 wind/SNII and AGB
inputs are NOT connected by this control. The existing LC18 cross-check still
records missing age-resolved wind composition, injected-energy mapping and
momentum/deposition specification; AGB source issues remain as previously
parked. No omitted physical field was replaced by zero or synthetic data under
a physical label. Completing that part requires a supplied complete source
package or an explicitly selected release/energy/composition approximation.
No external audit, source correspondence, commit or push was performed.

## LC18 actual wind -> native feedback, explicit approximation (2026-09-07)

The operator approved the next LC18 connection. Implemented the actual
40--120 Msun wind input under `user_selected_model_v1` / `wind_only_collapse`,
without reopening unrelated audits or adding a new gate framework.

`tools/build_lc18_native_wind.py` reuses the existing checksum-verified reader
and phase-lifetime aggregation. It does not alter raw data or promote the old
review-only source contract. LC18 Set R table8 supplies integrated wind isotope
ejecta above 25 Msun; table7 supplies the summed evolutionary lifetime. The
native branch has 16 nodes (40,60,80,120 Msun x four source Z values), 18 common
age nodes and 576 canonical wind/SNII rows. The SNII rows contain only the
residual remnant; terminal explosion ejecta/energy are suppressed explicitly.
Total wind mass comes from the isotope sum, not a forced match to rounded
phase masses. Source Z values are 3.236e-5,3.236e-4,.003236,.01345; their
linear cumulative mixture uses the existing native interpolation. No Z
extrapolation or BR26 source merging is introduced.

Required choices: rotation, positive nonrelativistic (<0.1c) effective wind
speed, uniform release until the actual terminal lifetime, fixed average
composition, and isotropic thermalization of E=.5*M*v^2. No source-resolved
wind speed/history is claimed. Isotopes retain their tabulated parent element
without decay. The comparison run chooses rotation=0, vwind=1000 km/s,
Chabrier/single-star, AGB/SNIa off. This speed is not literature-calibrated and
has no default. The existing native bridge adds stellar bulk kinetic energy
separately; directed wind momentum is zero in the selected isotropic model.

Native history gains `net_yield_policy`, default `supplied`. This input selects
`unavailable_diagnostic_zero`: the loader requires zero placeholder net
columns, marks them unavailable and logs the distinction. Actual gross ejecta
remain the physical source. HDF5 source identity binds the new semantics;
all actual input rows and chosen source/model labels remain bound as before.
These are external history-input fields, not new RAMSES RUN namelist fields.
The existing generator's history-file selector therefore needs no new schema.

Verified source SHA256:

- table7: `165201557b9093cab56b978e75b4d5989bc6e12491da54cc2bb26d20d99e6b15`
- table8: `bdfc317ca12a377f545964424dba4a666eab964292d0ebfa8b8c9641f073f218`

Local input `.lc18-native.X4wP2A/input/`:

- yields.dat: `ad2c8376b300fc243c3e337946e510e3b0f4b6bea57dc3846c345a6732117610`
- history.nml: `cddfc79be309df02bf403e6aa428f1fb24969c52e12bc4ffd9230bf4b4a0cc82`

The source/header selection also binds all other verified acquisition files.
Derived numerical files remain local artifacts, not redistributed source data.
Rotation 0 lifetimes range from 3.103519 to 5.393610 Myr. Exporter checks for
rotation 150 and 300 also passed, but those branches were not live-run here.
Doubling selected speed leaves mass/remnant unchanged and multiplies energy
by four; zero/negative/nonfinite or >=0.1c speed is rejected.

Full native build `.production-native.YCek9O/ramses_lc18_wind3d`, SHA256
`eb8ce41174348fa4bdb246be442c8abe854fe87a6d812c59652eb9b91ba543d0`.
Build log: `build-lc18-wind.log` in that directory. Existing lagRamses-first
VPATH and NVAR=30/SNRT/DUST_LIVE/HDF5/CUDA build are retained. The native
`lc18_native_wind_test.f90` passes both optimized Intel objects and GNU
bounds-checked modules: all 16 raw endpoints/11 elements, lifetime timing,
wind+remnant closure, out-of-range Z rejection, SSP split-step equivalence.
For initial SSP mass 10000 Msun at Z=.01, full wind return is
432.087111244174 Msun with 4.295961290478510e51 erg under the selected speed.

Live evidence under `.lc18-native.X4wP2A/`, each effective input `run.nml`:

- `fresh`: single rank x two OpenMP threads, CPU hydro/Poisson, four fixed
  level-3 steps, 2048 naturally formed stars. Run completed, exit 0, final
  time about 18.09 Myr. RT/dust evolution, AGN, AGB and SNIa are off in this
  isolated wind control. Returned stellar mass grows from
  1.2751990779689756e-8 code mass at step 2 to 3.5216440223926617e-8 at step 4.
  Gas+star mass remains .001 code units to displayed precision. Hydro values
  are finite and internal energy is positive. Final maximum Z is
  .010000426641676831, within the supplied metallicity range.
- `restart`: copy checkpoint 2 into a NEW directory and resume to step 4.
  Run completed, exit 0, `STELLAR_SOURCE_RESTART_IDENTITY_PASS`. Sorted particle
  mass/mp0/Z/progress/birth time/position/IDs and source vector match exactly.
  Leaf hydro fields also match exactly except two momentum components:
  maximum differences 6.46235e-27 and 3.23117e-27 code units, each <4e-17 of
  its field's peak. Particle velocity differences are <=6.61745e-24 code
  units (<4e-17 of component peak). Initial bitwise comparison of velocities
  failed at this roundoff level; a separately reported 64-epsilon peak-norm
  comparison passed. No claim of bitwise velocity reproducibility is made.
- `reject-speed`: separately generated vwind=2000 km/s input, same copied
  checkpoint. Native source identity rejects before evolution, exit 1, no
  new dump. Original input/checkpoint remains unchanged.

Output policy is noutput=1,aout=2,tout=1e30,foutput=1,fbackup=1000000.
Fresh four dumps total 2,189,868 bytes. Launch budgets were 30 MB for fresh
and 20 MB additional for restart tests; free space checked at 170 TB.
Reusable profile: `config/lc18_wind_local_smoke.nml` (history path placeholder
must be replaced). No new main-namelist defaults or legacy source files change.

**Completed:** real integrated LC18 wind/element/lifetime data reach native
feedback with explicit release/energy/composition assumptions and protected
restart. **Not completed:** phase-resolved physical wind closure, terminal
explosion physics, full AGB/WD inputs, all-channel science qualification or
source redistribution permission. Next implementation is actual AGB input
and WD supply for the already connected SNIa path. No external audit, commit
or push was performed in this increment.

## Operator-selected KL16 gross normalization (2026-09-07)

Following the explicit instruction to normalize, implemented
`all_listed_elements_to_selected_expelled_mass_v1` in the existing KL16 reader
and recorded that policy in the source-selection matrix. For every active
model, normalize raw gross masses by the sum of all 78 elements, then multiply
by the previously selected total expelled mass. Do not normalize just H..Fe,
modify the remnant, repair source bytes or activate the two commented models.

Each row now contains `selected_ejecta`: original sum, scalar factor, sum
correction, selected total/remnant, all-element gross masses and fractions,
tracked ejecta vector, untracked ejecta and total metals. Old raw fields and
their net diagnostic retain their old meaning; normalized net yields remain
null because initial-composition normalization/full-model matching has not
been selected. The helper always derives from raw values, so repeated calls
do not apply a second correction. Nonfinite/negative/zero-sum input, excluded
nodes, incomplete element payloads, non-closing total/remnant and unsupported
normalization policies reject rather than being silently repaired.

All 62 active nodes pass mass/fraction closure and raw-ratio preservation in
the existing `g2_source_selection_gate.py --include-parked-agb` regression.
Factors range from 0.9941975316583814 to 0.9950586970403992; maximum absolute
gross-sum closure residual is 4.440892098500626e-16 Msun. The 4 Msun, Z=.03
node scales its raw 3.24313662090249 Msun sum by 0.9947160348435394 to the
selected 3.226 Msun return, retaining the .774 Msun remnant and
.0020049972892055502 Msun in untracked elements.
The default source-selection regression also passes. All nine pinned KL16
source-file SHA256 values remain unchanged. No additional audit or simulation
was launched for this input-only change; no commit/push was performed.

The discrepancy's cause is not declared to be rounding: decimal-token analysis
found raw relative excess 0.0049658407--0.0058363335, 44.49--559.04 times the
worst nearest-rounding bound of the printed ejecta columns plus expelled-mass
header. Allowing one full last-place unit per token for one-sided rounding
still leaves a factor 22.24--279.52. Coarser intermediate processing is not
excluded without the original generation code. The chosen normalization is a
mass-conserving processing convention, independent of that unresolved cause.

This completes the requested gross normalization, not the remaining actual
AGB release-time/energy input or all-channel production qualification.

## Actual KL16 lifetime/energy and causal SNIa combination (2026-09-07)

Operator requested all three follow-ups: (1) actual AGB lifetime/envelope/WD
connection, (2) explicit wind energy, (3) combined LC18+AGB+SNIa native check.
Items 1 and 2 are implemented and execute in RAMSES. Item 3's native coupling
and causal rejection execute, but the unchanged physical combination fails;
do not report all three as scientifically complete.

### Source and implementation

- Added `simulation/snrt/data/kl16_stellar_lifetimes.csv`, 71 source rows,
  SHA256 `62e7fd40ce1d63d9024c1755f8226bb501ce92db95831c945ef18eaa2f80ee13`.
  Sources: [Karakas 2014 Table 1](https://arxiv.org/html/1408.5936v1#S3.T1)
  (canonical-helium branches) and
  [Karakas & Lugaro 2016 Table 1](https://arxiv.org/html/1604.02178v1#S2.T1)
  (additional/overshoot models). All 71 numerical durations independently
  matched the downloaded primary HTML. Source HTML SHA256 values are
  `97b4d24ae5b8e00c644e70b447c266e10e5c2a125bed55da7355a4427c85c556`
  and `68aebaa8f1b8a3f9fb53d9ef63f0290e45e9d0c62122f4de98379fb93499d7c9`,
  respectively; local copies remain in `.agb-physical.4LAOTJ/`, not source
  redistribution artifacts.
- The existing reader optionally attaches exact M/Z/Y/overshoot matches to
  all 62 active yields. In particular 1.75 Msun/Z=.014/overshoot=2 uses
  KL16's 1755 Myr, not the no-overshoot K14 value 1756 Myr. KL16 section III.1
  specifies that omitted N_ov denotes no overshoot. M_mix does not change the
  parent stellar track. Durations refer to the source's total evolution to
  its AGB endpoint (K14 includes pre-main sequence), not a resolved wind
  history or WD cooling duration.
- The new offline builder selects 58 ordinary CO-core nodes on common
  1--6 Msun support. Excluded hybrid/ONe cores and commented 8 Msun yield
  blocks are not relabelled CO WDs. Original normalized gross, untracked
  metals and selected remnants are retained; all raw files are unchanged.
- Native `terminal_step` releases envelope and WD together. Every node has
  zero origin, actual terminal age and a constant endpoint through 20 Gyr.
  The terminal-row validator replaces the rectangular-grid requirement only
  for this validated sparse AGB representation. Other grids remain checked.
  Nearest-mass fractions and fixed linear-Z mixtures use native SSP/IMF
  integration; AGB cell boundaries participate in the common quadrature.
- AGB velocity is required, finite and positive, below 0.1c. This comparison
  chooses 15 km/s, with E=0.5*M*v^2 isotropically thermalized and zero net
  directed wind momentum. LC18 retains the explicitly chosen uniform release,
  mean composition and 1000 km/s. These speeds are not source-calibrated.
  Changing AGB 15 to 30 km/s gives exactly four times its energy, unchanged
  mass/chemistry and LC18 rows, and a different source identity.
- Restart binds actual table values, source identity and terminal policy.
  Native SNIa rejects any cumulative WD shortfall at the interval endpoint
  and immediately before every WD formation event. A large dt cannot borrow
  WDs from a later time. No DTD delay, normalization, approved model identity,
  main namelist default or source admission flag was changed.

### Bounded verification

Full NVAR=30/SNRT/DUST_LIVE/HDF5/CUDA native build, lagRamses-first VPATH:
`.production-native.YCek9O/ramses_kl16_lc183d`, SHA256
`c20e416030f7d74279380841117e864bc561f66a79ad9d227d3534a5a60c5721`.
Build log `build-kl16-lc18.log` in the same directory.

`tests/fixtures/phase0/kl16_lc18_native_test.f90` passes optimized Intel and
GNU `-O0 -g -fcheck=all -fbacktrace`: 58 exact endpoint/11-element tests,
pre-event zero WD/envelope, mass and wind-energy closure, split-step SSP
equivalence with 8/64/127 bins, and rejection of unchanged 40 Myr DTD at
50 Myr and even with one large 200 Myr interval. For a 10000 Msun Kroupa SSP
at Z=.01 and age .2 Gyr, total return is 893.188604815424 Msun, including
485.906645104170 Msun AGB ejecta and 109.707363110781 Msun AGB remnants.
Existing `g2_source_selection_gate.py --include-parked-agb` passes including
the 62 exact lifetime attachments. Invalid wind speeds and incompatible
SNIa/Chabrier selection reject during input preparation.

Live evidence root: `/gpfs/kjhan/LRD_JWST/.agb-physical.4LAOTJ/`.
Each run's effective namelist is `run.nml` in its own directory.

- `single`: actual LC18+KL16, Kroupa single-star, SNIa off, single rank/two
  OpenMP threads. Four level-3 CPU hydro/Poisson steps, exit 0, `Run completed`.
  Naturally formed 2048 stars, final time about 181.59 Myr. RT/dust evolution,
  cooling, AGN and sinks intentionally off to isolate stellar return.
  Gas+star mass is .0010000000000000002 at step 2 and
  .0009999999999999998 code units at step 4 (initial .001).
  Returned stellar mass grows from 1.2115979048545474e-8 to
  4.412673338822857e-8 code units. All 30 leaf hydro fields are finite,
  density/internal energy positive; final minimum internal energy
  3.235569233811674e-8 and maximum Z .010000496918701366.
- `restart`: copy checkpoint 2 into a NEW directory, same inputs, resume
  to step 4. Exit 0, `STELLAR_SOURCE_RESTART_IDENTITY_PASS`, `Run completed`.
  Sorted particle mass/mp0/Z/IDs/progress/birth times/positions and the full
  source vector match exactly. Particle velocity differences <=8.73503e-22
  code units (<=5.22e-15 of component peak). Hydro momentum differences
  <=8.78880e-25 (<=5.26e-15 of peak); uold_15 differs by <=1.61559e-27
  (1.40e-16 of peak); all other leaf fields exact. The explicit 64-epsilon
  peak-norm comparison passes, not bitwise velocity reproducibility.
- `snia-reject`: actual source table plus the unchanged approved effective
  binary/Kroupa/N100/Maoz contract. At about 51 Myr, native
  `AGB/SNIa rejected: WD supply cannot fund the DTD causally`, exit 82.
  This is a demonstrated physical incompatibility, NOT successful SNIa-on
  evolution or a technical table/build failure.

Generated single table/history SHA256:
`d01919841578313a35ffa4be13258470948523d85dd244c469d59d31760f8ea1` /
`e1a4be93466972a078ec6bba150d6b1dc86bffa4dffb932e3e2acff99a6177d2`.
SNIa table/history SHA256:
`de15e4e62b757acaac581b4bf536c733cc866d9db49f88a2d4393991cb7c87d2` /
`894339fd9e1b1238b1a658cb5829b20c25cdb94b992c308e405012792091690d`.
Derived numerical packages stay local. Output policy was
noutput=1,aout=2,tout=1e30,foutput=1,fbackup=1000000. Initial two-run budget
20 MB for at most eight dumps; single's four dumps total 2,352,304 bytes.
Restart additional budget 10 MB/two new dumps; free space checked at 169 TB.
No nonconvergent-MG warnings were found in successful runs.

### Remaining physical decision

The approved DTD starts at 40 Myr, but the selected ordinary CO-AGB supplier
first creates a WD at 68.89 Myr (6 Msun, Z=.007; other 6 Msun branches are
69.85/72.55 Myr). Merely extending to the available 7 Msun ordinary CO nodes
would not fix a 40 Myr onset: their lifetimes are 48.44/50.14 Myr, while the
Z=.007 7 Msun model is hybrid. The commented 8 Msun models are not admitted
and must not be fabricated into a prompt CO-WD supply.

Choose a scientifically consistent early binary/WD supplier or a separately
identified, explicitly approved comparison DTD/population before claiming a
physical combined run. Neither is inferred from permission to implement the
three connections. Current approved DTD and strict mass/causal checks remain
unchanged. No new audit bundle, external audit, commit or push in this increment.

## Prompt WD supplier review under unchanged DTD (2026-09-07)

Operator authorized investigating an early WD supplier compatible with the
existing 40 Myr DTD. Literature and native mass accounting were inspected;
no new source was admitted, no DTD/event parameter changed, and no new
simulation, audit framework or binary population-synthesis job was launched.

### What the incompatibility does and does not establish

68.89 Myr is the first WD formation time of the SELECTED common-support
KL16 CO-AGB model, not a universal lower limit on WD formation. Conversely,
the existence of an early binary CO WD is not sufficient to justify an N100
explosion: accreted donor mass, retention, delay to explosion, event frequency
and explosion type also matter.

[Maoz, Mannucci & Brandt 2012, abstract and section 2.1](https://academic.oup.com/mnras/article/426/4/3282/1017965)
recover a prompt bin below .42 Gyr, not resolved measurements of a 40 Myr
turn-on. Their result supports the adopted slope/normalization and Kroupa
mass basis; it does not directly fix a unique binary channel or resolve the
40--69 Myr interval. The 40 Myr cutoff remains an approved project model
choice, not an independently measured onset. It was not changed here.

### Candidate disposition

| Candidate | Verified result | Disposition for this connection |
|---|---|---|
| CO WD + He-star donor, Wang & Han 2010 | Approximately 45--220 Myr SN delays; donor transfer grows the WD toward Chandrasekhar mass; population uses Miller--Scalo primary IMF | Relevant prompt-channel physics, but not a ready Kroupa-normalized 40 Myr supplier or a proof of sufficient event capacity |
| BPASS in Briel et al. 2022 | Their modeled SNIa require at least about 100 Myr, with metallicity/evolution dependence | Not a demonstrated solution for this 40 Myr onset; does not exclude other BPASS versions or channels |
| StarTrack double-detonation study, Ruiter et al. 2011 | Prompt He-star and delayed WD-donor sub-Chandrasekhar branches | Not interchangeable with the approved near-Chandrasekhar N100 mass/energy/yields |
| COSMIC/BSE event histories | Public outputs identify CO versus He/ONe WDs and track stellar masses and binary evolutionary transitions | Practical candidate for a bounded OFFLINE supply calculation; neither a calibrated SNIa population nor an admitted yield table by itself |

Primary references:

- [Wang & Han 2010](https://arxiv.org/html/1003.4050v1), sections 3--4:
  reported onset is an SN delay, not the first WD formation time. Their
  population rate cannot simply be copied into the existing Kroupa contract.
- [Briel et al. 2022](https://academic.oup.com/mnras/article/514/1/1315/6576337),
  section 5.1. The restriction above is specific to the modeled population.
- [Ruiter et al. 2011](https://academic.oup.com/mnras/article/417/1/408/979905):
  distinguish sub-Chandrasekhar double detonations from N100; a short delay
  in a different explosion family does not validate the current event source.
- [COSMIC output documentation](https://cosmic-popsynth.github.io/docs/stable/pages/output_info.html):
  kstar=11 denotes CO WD, 10 He WD, 12 ONe WD. `bpp` records key events;
  `tphys`, component masses/types and binary IDs provide source histories.
  A final remnant count alone does not provide cumulative eligible fuel.
- Recent scope check: [Rajamuthukumar et al. 2025 preprint](https://arxiv.org/abs/2511.11998)
  models a particular hot-subdwarf/WD pathway and reports an integrated
  contribution about 1.69e-5 events per formed Msun, not the full observed
  SNIa rate. This specific pathway is not evidence that every prompt event
  in the current effective DTD can be assigned to He-star donors.

This is a bounded candidate review, not a proof that no physical 40 Myr
channel exists. No mass/age/Z-resolved, consistently normalized early supplier
ready to connect to this KL16/N100 combination was obtained in this review.

### Required mass, independently calculated from the current namelist

For initial SSP mass 1e6 Msun, integrate the unchanged power law with
alpha=-1.07, t_min=.04 Gyr, t_max=13.7 Gyr, N_total/M0=.0013 and
debit/event=1.4004633930489443 Msun. No extra factor .5 is applied: the
contract's binary fraction is already baked into the observed event rate.
For q=alpha+1, cumulative events are
`M0*.0013*expm1(q*log(t/.04))/expm1(q*log(13.7/.04))` within the delay range.

| Age (Myr) | Cumulative expected events | Cumulative debit (Msun) |
|---|---:|---:|
| 45 | 31.827210 | 44.572842 |
| 50 | 60.076115 | 84.134399 |
| 68.89 (just before first selected WD event) | 144.734626 | 202.695546 |

At 68.89 Myr this is 0.0202696% of initial SSP mass but 11.1334% of the
integrated SNIa event count. Small total mass is not permission to manufacture
it. Values are continuous event expectations, not sampled integer events;
the pre-WD-event values are the left limits of the continuous DTD integral.

### Minimal next implementation, conditional on a selected physical model

Native `stellar_ramses_runtime.f90` currently obtains all available SNIa fuel
from `channel_remnant_mass(channel_agb)`, subtracts prior SNIa return, and
debits 1.4004633930489443 Msun per event. This aggregate accounting has no
binary donor-to-WD growth term. It is insufficient for directly inserting
the He-star donor mechanism without changing the population ledger.

1. Before designing a new runtime module, obtain a small, version-pinned
   OFFLINE binary-history sample at the intended Z and declared population
   weights (initial singles + both binary components). COSMIC/BSE is a
   candidate tool, not a selected default physics model. Pin mass-ratio,
   period, common-envelope and retention choices; identify actual CO-WD
   formation, retained donor transfer and eligible explosion histories.
   Test both timing and cumulative available fuel, not merely first WD age.
2. Only with viable histories, add a separately identified binary component
   whose initial mass is removed from the ordinary KL16 SSP component.
   Track living donor -> WD transfer, envelope loss, surviving companion and
   explosion debit conservatively. Do not count the same initial stars again
   through full KL16 AGB return. Separate formation/transfer from gas return;
   include source composition and energy where mass leaves the stars.
   An empirical event normalization must not exceed the eligible binary
   capacity; reweighting needs an explicit population-model decision.
3. Connect the compact cumulative source histories to existing native
   deposition/restart and check one combined run. No binary evolution solver
   belongs in the hydrodynamic timestep, no new gate hierarchy is needed.

Do not implement `WD += missing_DTD_mass`, debit arbitrary living stars at
explosion without a donor model, shorten KL16 lifetimes, or reuse N100 for
an unselected sub-Chandrasekhar pathway. If preserving the exact 40 Myr onset
requires a materially new binary population model, request that scientific
choice rather than representing the present engineering approval as approval
of unspecified binary physics. A separately named delayed-DTD comparison is
an alternative only with explicit approval; baseline N100/Maoz remains intact.

## Executed bounded offline binary feasibility (2026-09-07)

Operator approved the preceding recommendation. Completed the calculation,
not merely another literature review. **Early CO-WD formation and donor
transfer exist in the chosen comparison; a compatible N100 supplier is NOT
established.** The conditional native-source splice was therefore not activated.

### Reproducible comparison and scope

- Workspace `/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`.
  Isolated environment `.prompt-wd.uCcaeJ/venv`, Python 3.13.11,
  `cosmic-popsynth==4.2.0`. Existing SNRT/JAX environments are unchanged.
  Installed dependency versions are in `.prompt-wd.uCcaeJ/requirements.lock.txt`;
  `pip check` passes. This is an offline BSE/SSE calculation, not Python
  replacing native feedback or entering a RAMSES timestep.
- Input physics is the unchanged
  [COSMIC v4.2.0 example Params.ini](https://raw.githubusercontent.com/COSMIC-PopSynth/COSMIC/v4.2.0/examples/Params.ini),
  SHA256 `68b2c0a7a90cd179e4f98c27428935472713e161ea0a19e5c8bf7eef8d3d04e9`.
  In particular SSE, alpha_CE=[1,1], lambda prescription with lambdaf=0,
  ceflag=1, qcflag=5, acc_lim=[-1,-1], epsnov=.001, ifflag=1, eddfac=10.
  These are comparison prescriptions, not newly approved physical defaults.
  The entire BSE/SSE dictionaries are recorded in each `selection.json`.
- `simulation/snrt/tools/run_prompt_wd_feasibility.py` executes 1620
  initially main-sequence binaries: Z=.007,.01,.014,.03; primary midpoint
  masses 4.5--12.5 Msun in nine bins spanning 4--13; five secondary-mass
  quantiles and nine log-period midpoints spanning 1--31623 days, e=0.
  Each Z has 405 systems, evolved to 100 Myr with 1 Myr regular output plus
  exact solver transition records. Two CPU workers, one thread each, timeout
  budget 600 s; all four batches completed normally well within the budget.
- Explicit illustrative weights use a Kroupa-2001 PRIMARY IMF on .08--120,
  binary fraction .5 by SYSTEM count, secondary mass uniform on [.08,m1],
  and uniform log10(P/day) on [0,4.5]. If I_N and I_M are the unnormalised
  IMF number/mass integrals, formed-mass denominator is
  `I_M + .5*(I_M+.08*I_N)/2`, including singles and BOTH binary components.
  No renormalization to the targeted 4--13 subset or observed SNIa rate.
  The sampled binaries represent about .084206 of initial population mass
  under this coarse quadrature. Other primary masses and singles are not
  evolved. A Kroupa primary IMF plus this companion distribution is NOT
  identically the individual-star IMF underlying the approved empirical DTD.
  Weights give a scale comparison, not a calibrated population/event rate.

Run outputs: `.prompt-wd.uCcaeJ/grid-v2/`, log `grid-v2.log`.
Files per Z: original grid/weights, full `initC` (including actual Fortran
random seeds), `bpp` event histories, `bcm` samples and kick histories.
All 405 systems per Z have every requested integer-age sample; 163620 regular
rows are finite and nonnegative in stellar mass, with no sampled total
component mass exceeding its initial binary mass. This check is not a
time-resolved chemical-ejecta conservation proof.
The initial attempt `grid/` failed before evolving a population because
COSMIC appends `bin_num` itself rather than accepting it as a requested
Fortran column. The runner was corrected and a NEW output directory used;
failed evidence was not overwritten.

### Measured results

| Z | First CO WD (Myr) | Components forming CO WDs by 40 Myr | CO mass present at 40 Myr per 1e6 initial Msun |
|---|---:|---:|---:|
| .007 | 27.505616 | 21 | 132.157764 |
| .010 | 27.142089 | 32 | 184.088211 |
| .014 | 26.711868 | 33 | 168.742286 |
| .030 | 21.262627 | 70 | 388.077178 |

Counts are targeted grid components, not observed or random-sample event
counts. CO mass includes wide/ineligible binaries. It is neither an eligible
N100 reservoir nor cumulative fuel delivered to explosions. At Z=.01 the
sampled inventory at 45/50/69 Myr is 283.488/747.460/1593.262 Msun per 1e6
initial Msun. Large aggregate inventory alone cannot establish the required
event capacity or justify pooling unrelated WDs into a single N100 event.

Four selected runs in `.prompt-wd.uCcaeJ/followup/` replay the same physical
flags and actual per-system seeds with .005 Myr regular output; one changes
only the period to 1e9 days as a non-interacting control, not a population
sample. All finish (COSMIC's final output overshoots nominal 100 by .005 Myr).

- `early_CO`: Z=.01, initial 10.5+3.206 Msun, P=177.827941 days. Common-envelope
  stripping at about 23.142 Myr produces a naked He star; CO WD appears at
  27.142114 Myr, mass .64984246 Msun. It remains a WD to the final time.
  `wide_control` with the same masses/Z produces a neutron-star primary,
  not an early WD. Thus binary evolution changes the fate in this comparison.
- `early_disappearance`: Z=.01, 10.5+5.29 Msun, P=5.623413 days. CO WD forms
  at 27.294855 Myr. The companion later supplies He-rich material. Between
  38.890 and 38.905 Myr, WD gains .04677965 Msun, donor loses .05379474,
  and total stellar mass falls .00701510 Msun. These are resolved net
  changes including losses, not a universal retention efficiency.
  The WD disappears at 38.908704 Myr. This is NOT a validated N100 event.
- `donor_growth`: Z=.01, 7.5+5.274 Msun, P=17.782794 days. WD forms at
  53.096783 Myr, .62791744 Msun; donor transfer grows it to .71410023 Msun,
  but the donor undergoes collapse to a neutron star. Growth is not
  synonymous with a SNIa explosion.

Increasing output cadence changes COSMIC integration boundaries: coarse/fine
early-WD ages differ by about 25 yr; the selected disappearance shifts from
38.922140 to 38.908704 Myr. These two cadences confirm the qualitative channel,
not bitwise reproducibility or population/timestep convergence.

### Why the conditional N100 splice remains inappropriate

Inspected [version-pinned evolv2.f](https://raw.githubusercontent.com/COSMIC-PopSynth/COSMIC/v4.2.0/src/cosmic/src/evolv2.f),
lines around 3211--3248. Its CO-WD He-accretion branch destroys the WD when
`mt2 - mass0 >= .15 Msun`, before the separate Chandrasekhar-limit branch.
For the selected disappearance, mass0 is .63119778 Msun: the threshold is
about .7812 Msun, not the native N100 debit 1.40046339 Msun. All selected-grid
CO->massless evol_type=9 transitions before 68.89 Myr have a preceding
type-8 He-star donor (counts 9/10/13/13 by Z); common-envelope mergers are
separately present and must not be counted as SNIa. Neither this classification
nor the generic BSE SN flag supplies an N100 yield/energy/source approval.

**Disposition:** early WD formation is now demonstrated in an explicit
comparison, superseding the absence of executed evidence in the previous
review. Compatible prompt near-Chandrasekhar N100 event supply is still not
demonstrated. Keep KL16, the 40 Myr DTD, N100 and native causal checks unchanged.
Do not turn this inventory into a production source or add it on top of the
full KL16 population. The next scientific choice is a near-Chandrasekhar
accretion/retention population consistent with N100, or a separately selected
sub-Chandrasekhar comparison with its own event yields/mass/energy and rate.
No additional native module, gate hierarchy or large BPS campaign is justified
until that choice. No external audit, commit or push was performed.

Storage: grid data 8.5 MB, selected follow-ups 4.1 MB (compressed CSV plus
small JSON); free GPFS space was 169 TB. Existing RAMSES outputs are untouched.
Package wheel SHA256
`b6a7698acb170e08420d95c50f57a1fc80c691407ed9b55eec31eda7257493f8`;
grid summary SHA256
`781e9547d1064bf0f771ca6798a37abfa2ed0fbc21b6a58c0724ddeb0c593fd9`;
dependency lock SHA256
`e0d65df86f94b10c65726d7e73bb3af76d4a128aa5070cbeb984803c7cec8658`.
The runner pins the upstream settings checksum, refuses existing output
directories and records its own hash. Example invocation uses the isolated
Python, `--params .prompt-wd.uCcaeJ/Params-v4.2.0.ini --output-dir NEW_DIR`;
append `--followup-grid .prompt-wd.uCcaeJ/grid-v2` for the four selected histories.

## Near-Chandrasekhar retention reference selected and screened (2026-09-07)

Operator requested continuing the N100-preserving recommendation. Implemented
`W17_stable_KH04_flash_conditional_v1` as an OFFLINE reference in
`tools/helium_retention_reference.py` and applied it to the existing selected
histories through `run_prompt_wd_feasibility.py --retention-screen`. This is
a selected comparison ingredient, not a calibrated binary population, coupled
binary re-evolution, or native N100 source admission. No native feedback,
DTD, event mass, initial source table or installed COSMIC code was changed.

### Physical prescription and explicit limits

- [Wang, Podsiadlowski & Han 2017](https://arxiv.org/html/1708.07067v1),
  equations 1--2, supplies the lower stable-burning and upper expansion
  rates as functions of CO-WD mass. The implementation uses the stated
  initial-mass domain .6--1.35 Msun and rate grid 1e-8--1e-5 Msun/yr.
  Stable burning has unit retention in this reference, not an automatic
  explosion. Rates above the expansion boundary require a coupled wind/RLO
  model; they are not silently clipped to make the WD grow.
- [Kato & Hachisu 2004](https://arxiv.org/html/astro-ph/0407632v1), equations
  1--6 and the .7 Msun discussion, supplies conditional flash retention.
  Node fits are evaluated at .7,.8,.9,1,1.1,1.2,1.3,1.35 Msun; interpolation
  in WD mass is linear and restricted to the common rate support of both
  adjacent nodes. No extension to .65 Msun or below a fit's stated lower
  rate is made. A missing efficiency stays null, not zero or one. Slight
  overshoots of unity in printed polynomial joins are explicitly capped at
  one. This is not a treatment of erosion in other retention models.
- KH04's flash retention is conditional on the expanded envelope not losing
  additional mass through its Roche lobe. The current compact-WD radius is
  not evidence that its flash-expanded envelope fits. The screen therefore
  does not admit the flash branch solely from the fitted efficiency.
  [Piersanti, Tornambe & Yungelson 2014](https://arxiv.org/html/1409.3589v1)
  explicitly treats this additional binary/thermal sensitivity; KH04 is a
  comparison choice, not asserted to be uniquely correct or most accurate.
- The approximately 2.05e-6 Msun/yr W17 off-centre ignition threshold is a
  risk flag under their hot, nonrotating, solar-composition WD assumptions,
  not an instantaneous universal collapse rule. A later
  [off-centre-burning study](https://arxiv.org/html/1904.05130v1) also finds
  possible OSi-core outcomes. Neither all such objects nor every growing
  WD is labelled a normal N100 event here.
- W17's composition is Z=.02 and accreted He fraction .98. Applying its
  mass/rate reference to the existing Z=.01 histories is explicitly a
  diagnostic comparison, not calibrated metallicity-dependent retention.
  No implicit extension from the fitted WD domain to the exact N100 ejecta
  mass is performed, and no mass is manufactured to reach that mass.

### Execution and results

Evidence root `/gpfs/kjhan/LRD_JWST/.wd-retention.aLENdq/`.
`screen/summary.json` and four small per-case JSON files retain input and
kernel hashes. Only existing fine COSMIC histories were read for this screen.
COSMIC regular output reports donor total mass loss and WD net mass growth,
not the incident helium transfer required by the new prescription. They were
evaluated as TWO DISTINCT RATE PROXIES, not certified upper/lower bounds on
incident transfer or a replacement for a coupled transfer calculation.

The early-disappearance case has 7 sampled CO-WD/He-donor RLO states (14 proxy
evaluations): 3 fall in the reference stable interval, 5 in the unsupported
low-mass flash region and 6 above the expansion boundary. The donor-growth
case has 3 such states (6 evaluations): 2 stable, 1 unsupported flash,
1 expansion and 2 outside the W17 rate grid. The other two cases have no
sampled CO-WD/He-donor RLO states. These are sampled local regimes, NOT
fractions of physical SNIa events or newly integrated retained masses.

Direct checks reproduce eight published KH04 node values, reject nonfinite
or nonpositive rates, preserve null for unsupported interpolation, and test
the W17 stable/expansion boundaries. A bounded 8181-point mass/rate sweep
keeps supported efficiencies in [0,1] and never admits an N100 event. No
new independent gate or external audit was created for these small checks.

An additional paired COSMIC calculation (`ifmr-control/`) checked the
representative 10.5+5.29 Msun, P=5.623413 day, Z=.01 binary with identical
initial conditions, saved seed and .005 Myr output, varying ONLY `ifflag`:

- `ifflag=1` (the original example's Han initial-final-mass prescription)
  reproduces the early CO-WD path and subsequent sub-Chandrasekhar loss.
- `ifflag=0` (standard BSE without that IFMR modification) instead forms a
  neutron-star primary, not an early CO WD for this system.

Thus the earlier 21--28 Myr result remains a result of its explicitly chosen
comparison, not evidence robust to the remnant prescription. This is not a
claim that either option is a COSMIC bug or that all other systems lack early
WDs. The paired raw bpp/bcm/initC histories, full flags and seed are saved;
reproducer `.wd-retention.aLENdq/ifmr_probe.py`. Neither option was promoted
to a project default.

### Disposition / scope control

Completed the reference prescription and its bounded application-domain
screen. NOT completed: a new self-consistent near-Chandrasekhar binary
population. Replacing the old .15-Msun He trigger in saved trajectories and
integrating beyond the original explosion would invalidate donor/orbit/loss
history. A genuine next source calculation must couple incident transfer,
retention/outflow angular momentum and donor/orbit response, with a declared
remnant model and ignition/fate treatment; the current proxy rates cannot
serve as that input. Do not grow a family of validators around this gap or
claim that a formula alone resolves it. Keep the N100 baseline and current
causal WD checks unchanged. No native source splice, commit or push in this
increment.

## Operator-approved scope decision: effective SSP SNIa (2026-09-07)

The operator approved separating microscopic binary evolution from the
immediate RT/feedback/dust objective ("추천한 대로 정리합니다."). The current
implementation priority is now a separately identified effective-SSP SNIa
accounting mode, retaining the empirical DTD and N100 per-event source values
but charging the remaining stellar-particle mass rather than requiring the
selected KL16 CO-WD inventory. The strict causal-WD mode remains a separate
comparison. This explicitly supersedes earlier recommendations that coupled
binary re-evolution is the immediate prerequisite for further native work.

Approval concerns the model direction, NOT a completed native implementation
or an all-channel physical pass. This documentation-only change leaves the
runtime and the measured 40--68.89 Myr strict-mode rejection unchanged. The
new accounting requires its own identity and restart binding, combined-channel
mass/element/energy bookkeeping without double debit or overdraw, and bounded
native/restart verification. An effective SSP does not demonstrate microscopic
consistency between single-star AGB progenitors and the SNIa population.

Empirical SSP-level SNIa delay distributions have a galaxy-simulation precedent
in [IllustrisTNG, Naiman et al. 2018, section 2.1](https://academic.oup.com/mnras/article/477/1/1206/4925009).
That precedent supports the modelling level, not the correctness or validation
of this project's not-yet-implemented accounting choice.

COSMIC/IFMR/retention/transfer/orbit/ignition and metallicity-dependent binary
population work is now medium-term and non-blocking. Existing scripts, raw
histories and negative results are retained, not relabelled as N100 evidence.
No new calculation or audit was launched. Implementation scope and sequence
are in the final operator-decision section of
[the existing completion plan](production_completion_bundles_2026-09-07.md).

## Effective SSP SNIa native accounting implemented (2026-09-07)

Operator "진행 승인." authorized implementation after the scope decision.
Canonical code is in `patch/lagRamses`, with the lagRamses-first Makefile VPATH
unchanged. No new Python runtime, validation framework or external audit bundle.

- External `snia_population_realization` now accepts `mass_accounting`, default
  `strict_wd`, and `accounting_approval_id`, default empty. Effective mode requires
  `effective_ssp` / `SNIA-EFFECTIVE-SSP-2026-09-07`; other combinations reject.
  Separate example: `simulation/snrt/config/fp2_snia_effective_ssp_runtime_v1.nml`.
  All original DTD/N100/thermal values and their three-group source approval
  bindings are unchanged. Kroupa/binary applicability and baked-in rate remain.
- The effective branch computes available mass from actual particle mass minus
  this interval's generic returns. It reuses N100 event-budget arithmetic but
  never sets a fictitious WD reservoir or modifies the generic remnant partition.
  The existing zero-terminal-remnant/equal-debit invariant remains mandatory.
  `close_effective_snia_return` checks aggregate closure, overdraw and consistency
  of reconstructed prior return with the declared cumulative DTD before the
  existing gas/particle/progress transaction commits. No delayed/clipped rate.
- The consumed mode/approval extend HDF5 source identity only for effective mode;
  historical strict layout is preserved. Startup logs the approximation. No main
  namelist field changed, so generator/mkrun changes are not needed here.

### Bounded evidence

Evidence root `/gpfs/kjhan/LRD_JWST/.snia-effective.lL76O2/`.
Incremental builds copy existing objects into this NEW directory, preserving
old binaries and outputs, and rebuild affected modules with the production
Makefile: `HDF5=1 USE_CUDA=1 SNRT=1 DUST_LIVE=1 USE_FFTW=0`, NVAR=30.
Final binary `ramses_snia_effective_final3d`, SHA256
`cd7cafa10c3bcb768e2c48593a3135e5cb469c3f38409e733c21fa19afa9d201`;
build logs `build.log`, `build-final.log`. New example contract SHA256
`3a167ebaf393bd2ac11217cb0a8b278e2e8da793a91c4b36e8a57e4464e4dff5`.

Extended existing Fortran fixtures, compiled against canonical sources:
`fp2_snia_runtime_contract_test`, `fp2_snia_runtime_accounting_test`, and
`kl16_lc18_native_test` pass optimized Intel and GNU
`-O0 -g -fcheck=all -fbacktrace`. Checks cover separate approval/unknown-mode
rejection, unchanged IMF binding, overdraw/nonfinite/missing-or-duplicate-prior
debit rejection, split-interval mass closure, 58 real KL16 endpoints and the
unchanged strict causal deficit. The actual combined source is exercised to
13.7 Gyr at Z=.007,.01,.012,.01345; 10000 initial Msun gives 13 expected SNIa
and 18.2060241096363 Msun total N100 ejecta, without a WD supply claim.
The first trial correctly rejected Z=.014 beyond LC18's .01345 upper bound;
the comparison was restricted to common support, NOT fixed by extrapolation.

Live inputs are each case's absolute `run.nml` under the evidence root:

- `live-final`: one rank / two OpenMP threads, four level-3 CPU hydro/Poisson
  steps using unchanged `.agb-physical.4LAOTJ/snia-input/{yields.dat,history.nml}`
  and the new effective contract. Exit 0, 2048 stars, final time 176.816195 Myr.
  Gas+star mass .001 code units, total stellar return 4.3959007744804755e-8,
  all 30 leaf hydro fields finite, minimum internal energy 4.258844652646133e-8.
  AGN/sinks/cooling/live RT/dust intentionally inactive to isolate stellar return.
- `restart-final`: copied checkpoint 2 in a new directory, same inputs to step 4.
  Exit 0, `STELLAR_SOURCE_RESTART_IDENTITY_PASS`. Sorted particle mass/mp0/Z,
  birth times, positions, IDs, progress and full source identity match exactly;
  velocities differ by at most 8.271806125530277e-22 (3.23609e-15 peak-relative).
  Hydro momentum differs by at most 8.239494382852424e-25 (3.22668e-15 relative);
  other leaf hydro fields are exact. Existing 64-epsilon peak-norm criterion passes.
- `strict`: original contract with the same real source input rejects at the
  early WD deficit, exit 82. This is the retained negative reference, not a pass
  for physically supplied prompt N100. It produced only checkpoint 1.
- `mode-reject`: effective checkpoint 2 + strict contract rejects source identity,
  exit 1. `reverse-mode-reject`: strict checkpoint 1 + effective contract also
  rejects identity, exit 1. Neither switches accounting in an existing history.
- Preliminary `live`/`restart` runs agree with final checks; all are preserved.

Output policy reviewed before launch: noutput=1, aout=2, tout=1e30 (outside
these noncosmological runs), foutput=1, fbackup=1000000, nstepmax=4. Measured
total output storage, including copied checkpoints and preliminary runs, is
10,598,608 bytes, within the 20 MB budget; GPFS free space was 169 TB.
No failed-MG warning in successful final runs. Numerical HDF5 checks read
existing outputs; no additional simulation data framework was introduced.

**Disposition:** native effective-SSP coupling and bounded conservation/restart
verification complete. This is NOT a self-consistent microscopic binary
population, a calibrated wind prescription, a broader-Z source approval, or
simultaneous RT/AGN/dust production qualification. Resume high-level integration
within these limits. BPS remains medium-term. No commit/push in this increment.

## Commit/push and effective SNIa + AGN/RT/dust integration (2026-09-07)

Operator requested committing/pushing and proceeding. Own accumulated source,
input adapters, fixtures and decisions were committed as `3fa4d48` (53 files).
The user's unrelated 88 generator-option deletions and all scratch data/builds
were excluded. Remote commits `e5900c5` (f(R)) and `f93fd53` (DMO restart) were
merged without conflict as `5671108` and pushed to `kjhan0606/LagRamses:main`;
remote tip equality was checked. GUI tests: 23 run, one display skip. The
upstream particle-type regression also passes after the merge. Neither upstream
change was undone or reopened as a high-level physics research task.

The next existing integration increment uses the actual KL16/LC18 source and
effective SNIa alongside accepted BH accretion -> reference AGN primary RT ->
live Draine-optics dust. It does NOT assign the incompatible Chabrier/single-star
SED to Kroupa/binary particles: `SNRT_STELLAR_SED` is absent. AGN spectrum and
constant dust heat capacity retain reference status; no new physical approval.
Reusable input/seed: `config/kl16_lc18_snia_agn_dust_smoke.{nml,ic_sink}`.

Evidence root `/gpfs/kjhan/LRD_JWST/.feedback-rt-dust.dqHCvB/`. Production
Makefile build with lagRamses-first VPATH, NVAR=30, SNRT/DUST_LIVE/HDF5/CUDA,
USE_FFTW=0; merged-base binary `ramses_feedback_rt_dust3d` SHA256
`0bd35a9f8eb6a67f82b3281407b4bfd948d099e04caddf5223644d7ff879146b`.
The runtime uses one rank/two OpenMP threads and forced primary OpenMP.

### Measured execution and narrow restart repair

- `live/run.nml`: four CPU hydro/Poisson steps, exit 0, final time
  102.390489738561 Myr. Real `.agb-physical.4LAOTJ/snia-input` source/history,
  effective SNIa contract, existing reference AGN group/secondary contracts,
  `.physical-inputs.A0XtMt/dust_native.nml`. The BH seed is code mass 1e-5,
  explicitly a numerical stress/control input (2.4527e11 Msun at these units),
  not a galaxy calibration. Cooling/stellar photons/new sink formation off;
  existing sink maintenance, mechanical feedback, primary RT and dust IR on.
  One nonzero accepted AGN source interval, four RT commits and four dust IR
  commits, no chemistry failures. Maximum IR balance residual 3.7746e-13;
  primary dust ledger relative error for the active interval 2.2402e-7.
- `restart/run.nml`: HDF5 stellar, AGN and RT state read successfully, then
  SIGSEGV exit 174 in `kjhan_sync_sink_particle_coordinates` during sink
  maintenance. Address 0x7784b5 maps to its next-particle traversal. At this
  point `create_sink` has gathered particles to level 1 but has not rebuilt
  fine-level lists; the routine unnecessarily walked every level.
- Narrow repair in `patch/lagRamses/sink_particle.kjhan.f90`: update canonical
  BH positions/velocities using the local map freshly filled by cloud creation.
  Zero entries are nonlocal; invalid index/type/ID rejects. No AMR/tree rewrite,
  cloud physics, accretion prescription or namelist change.
- `restart-fixed/run.nml`: copy checkpoint 1 (coarse step 2), nrestart=1, resume
  to step 4 with repaired binary `ramses_feedback_rt_dust_fixed3d`, SHA256
  `df124cd4be8bb3cb4ccc2a44ec74f695abbf001d66df1530ae0b6f23dd056a5c`.
  Exit 0, stellar identity passes. Sorted STAR mass/mp0/Z/IDs/progress/birth
  times, primary+IR SNRT datasets, dust mass/energy, BH mass and all five AGN
  pending reservoirs match the uninterrupted run exactly. Leaf hydro differences
  are <=3.88e-16 peak-relative. BH gas angular momentum, near cancellation,
  differs by <=1.64646e-24 code units (up to 4.27e-10 relative to its tiny value);
  do NOT claim all sink diagnostics bitwise/64-epsilon relative reproducibility.
  This control has spin disabled; spin-sensitive qualification is not established.

Final live/restart STAR count 1733, total stellar return
1.3072717565418019e-8 code mass. All leaf hydro variables finite, internal
energy positive (minimum 2.3046861425049755e-8), dust mass/energy nonnegative.
Gas+STAR+BH mass = .001009905788355202; initial gas+BH = .00101. With epsilon=.1,
BH growth implies radiative mass loss 9.421164479767314e-8 code mass. Adding this
back closes to 4.34e-19 absolute. Cloud tracers are NOT counted as additional BH
mass. The generic mcons/econs lines are not substitutes for this source accounting.

Oldest stellar processed age is 51.363844889651 Myr, before the 68.89 Myr AGB
terminal event. Thus simultaneous nonzero AGB is NOT demonstrated here, although
its input is loaded and standalone AGB already passed. No full physical stellar
SED, physical dust heat capacity, spin-sensitive or combined MPI/AMR promotion.

A smaller-seed exploratory `live-fixed` (1e-7 code mass) was retained but NOT
admitted: its timestep still shrank and star+cloud particles exhausted npartmax
4096 at step 4 (exit 1). It never demonstrated the intended later AGB release.
The reusable seed retains the completed 1e-5 control; do not describe this failed
case as a successful extension. Defer a deliberately age-covering integrated
profile rather than starting a parameter-sweep/gate ladder here.

Output policy: noutput=1, aout=2, tout=1e30 (unreached), foutput=2,
fbackup=1000000, nstepmax=4; checkpoint 1=step 2, checkpoint 2=step 4.
Output storage including copied/failed-case checkpoints: 330,593,437 bytes,
within the revised 450 MB budget; GPFS free space 169 TB. Existing results
are preserved. Next physical input work remains a population-compatible stellar
SED and dust thermal material input, without reviving the parked BPS prerequisite.

## Dust U(T) implementation and stellar-SED compatibility check

User authorized the next physical-input work. Implemented the missing native
temperature-dependent material-energy path, **not a production material-model
approval**. No main RAMSES namelist parameter, IMF, SNIa normalization, stellar
population or existing combined-control selection was changed. Unrelated
88-line deletions in `patch/cuRamses/aux/ramses_nml_generator.py` are preserved.

Native contract v4 carries positive, strictly increasing internal energy per
reference H on the thermal temperature knots and a material-source checksum.
The old constant capacity must be explicitly zero in v4. Runtime uses
`E_volume=nH*relative_dust*U(T)`, linear U(log T), for initial inversion,
implicit emission, and the conserved material state. Primary deposition defers
temperature to the coupled solve rather than applying constant-capacity heating.
HDF5 binds the actual U array; v2/v3 layouts are unchanged. The existing exporter
accepts source-identified specific-energy JSON, preserves its knots, converts
erg/g to erg/reference-H, rejects extrapolation, and grants no approval.

Evidence (all under `/gpfs/kjhan/LRD_JWST`):

- Existing native contract runner passed with ifx and gfortran, now including
  v4 nonlinear heating/cooling, known U(log T) inverse, energy closure,
  range-failure rollback, invalid knot rejection and reference opt-in/reset.
- Existing `dust_thermal.py` passed including v3/v4 export and four invalid
  material cases. Existing `dust_ir_transport.py` passed. Its opt-in
  `native_checks()` passed with `JAX_ENABLE_X64=1`: GNU/Intel relative energy
  differential errors <=1.12e-15 and temperature errors <=4.45e-16. The first
  direct native_checks invocation omitted x64 and was rejected before testing;
  the corrected invocation completed. No new gate framework was introduced.
- Full root-level Makefile build with HDF5=1 USE_CUDA=1 SNRT=1 DUST_LIVE=1
  USE_FFTW=0, NVAR=30, lagRamses VPATH first: executable
  `.dust-material.1056JI/ramses_dust_material3d`, SHA256
  `205c2c203747dd598b2b8c4d3018b0da4ec0dc061d6fd257138b15b99144a1be`.
  Existing same-option object cache was reused; changed modules/dependents were
  rebuilt. Later edits to this module were comments only; added test logic is
  not part of the RAMSES executable.
- `.dust-material.1056JI/material.json` is an explicitly **synthetic** cubic
  knot table. Its exported `dust.nml` uses actual WD01 optical data, 136 IR
  nodes, 67 temperatures, background 10 K, and reference-only status. It is
  not DL01 material data and does not close the physical-input gap.
- First `live` attempt was rejected during namelist admission because its
  inherited `SNRT_AGN_MODEL=partition_reference_v1` requires enabled AGN.
  The process returned 0 but no evolution/dump occurred: not a pass. Preserved.
- Corrected `live-dust-only/run.nml` used `SNRT_AGN_MODEL=legacy` with AGN off,
  SNRT enabled/OpenMP, two threads, no stars/particles/gravity/gas cooling;
  8^3 periodic hydro, four steps. It completed (exit 0, ~4.96 s), with four
  dust IR commits and peak balance 6.1368e-13. Final dust T=10.00000000000001 K,
  positive material energy and unchanged dust mass; no MG failure.
- `restart/run.nml`: copied step-2 checkpoint, nrestart=1, resumed to step 4
  (exit 0, ~2.58 s). All **93 hydro/RT datasets exactly equal** the uninterrupted
  result; the 1640-entry dust identity attribute also matches exactly.
- `restart-changed-U/run.nml`: one U knot changed by 1%, original checkpoint
  linked read-only in intent. Actual HDF5 identity rejection/MPI abort exit 10,
  before any new output. No existing checkpoint was modified.
- Effective schedules: noutput=1, aout=2, tout=1e30 (unreached), foutput=2,
  fbackup=1000000, nstepmax=4. Checkpoint 1=step 2, checkpoint 2=step 4.
  Four real/copy dump directories total **217,852,901 bytes**, below the
  declared 220 MB budget; the negative test uses a symlink, not another copy.
  GPFS free space was 169 TB. No production run was launched.

Stellar SED: directly inspected the local BPASS v2.2.1 HDF5 root/dataset
attributes. It is a binary population with IMF slopes -1.30 over 0.1--0.5 Msun
and -2.35 over 0.5--300 Msun, 51 ages from 1 Myr to 100 Gyr, 13 metallicities,
and Lsun/Hz spectra. It does not match the current Kroupa 0.08--120 Msun plus
binary-fraction convention. Changing an IMF header or multiplying the total
spectrum cannot recover missing track-level population weights. No BPASS SED
was activated, no star younger than 1 Myr was assigned an invented spectrum,
and no synthetic Chabrier control was relabeled as a physical binary source.

Remaining physical decisions, not extra validation gates:

1. Preserve the current feedback population while obtaining compatible SSP
   spectra (or track-level spectra/weights permitting actual re-integration).
   Changing the population to a supplied BPASS preset would be a separate
   model choice, not an adapter fix. The parked microscopic BPS prerequisite
   remains parked.
2. Supply source-identified grain internal energy with explicit composition
   and size approximation for the WD01 optical mixture. [Draine & Li (2001)](https://arxiv.org/abs/astro-ph/0011318) provides vibrational-mode/thermal
   models and explains small-grain stochastic heating. Inferring one bulk
   U(T) from that work still requires a declared mixture approximation;
   the raw [WD01/D03 opacity table](https://www.astro.princeton.edu/~draine/dust/extcurvs/kext_albedo_WD_MW_3.1_60_D03.all)
   contains PAH-like small and graphitic large carbonaceous grains plus silicates,
   not a measured aggregate heat-capacity table. Single-temperature bulk
   closure is not the stochastic PAH solution. The new native path is ready
   for such material input, but **physical dust thermal input remains open**.

## Follow-up: literature-based bulk dust material connected

User requested the remaining connections. The dust material gap above is now
resolved **for an explicitly parameterized bulk comparison closure**, not for
the full stochastic WD01/PAH mixture. New builder
`simulation/snrt/tools/build_dl01_dust_material.py` implements the DL01 Debye
mode prescription via the normalized oscillator density of states and converts
per-atom energy to erg/g. No zero-point energy and no C(T)*T substitution.
Graphite/silicate mass fraction is required at the CLI; the supplied comparison
uses 0.30/0.70 and representative silicate formula MgFeSiO4. It is not an
opacity-derived or calibrated abundance. Source and approximations are stored
in the material JSON. Main namelist/generator code was not modified.

Tracked source/export:

- `simulation/snrt/data/dust_dl01_bulk_030_material_v1.json`, SHA256
  `f9b2109c737d6be003fcb04cf5ae5ec56cdacbf9057759aae7f61fdaff418025`.
  162 material knots (5--300 K, explicitly including 10/20 K).
  U(20 K)=862843.9569500614 erg/g.
- `simulation/snrt/config/dust_dl01_bulk_030_reference_v4.nml`, SHA256
  `8072b09a592a2f13e9c3607f4d3d0a2c5d6e328072b3bfb8429eee55cced077e`.
  Actual existing Draine optics plus the material model, 222 union temperature
  knots, 136 IR quadrature nodes, 10 K bath; reference-only, no approval ID.
- `simulation/snrt/config/kl16_lc18_snia_agn_dl01_dust_smoke.nml` selects the
  existing coupled physics and correctly initializes primitive dust energy
  to 5.43633430456151513e-13 code energy per gas mass at 20 K.

Existing `dust_thermal.py` passes with added checks of material monotonicity,
64/128-point quadrature agreement (<1e-10), and the independent high-T limit
3*kB/atom (<1e-8 relative) for pure graphite and pure silicate. Native U(T)
tests/restart checks from the preceding section were reused, not repeated as
a new gate ladder. No new simulation binary was necessary: the measured native
v4 binary `.dust-material.1056JI/ramses_dust_material3d` directly consumes this
new source, SHA256 unchanged from the preceding section.

Live evidence: `.physical-sed-dust.mkYybi/live-dust-fixed-ic/run.nml` and log.
Four fixed level-3 steps completed, exit 0, elapsed ~12.72 s, with four IR
commits, nonzero AGN absorption and physical wind/SNIa inputs retained. Peak
dust/IR balance=6.6218e-13. Leaf dust T spans 10.000000000000007--71.493316048125 K;
all material energies and masses positive, RT datasets finite; no MG failure.
The HDF5 records v4, zero constant capacity, and the actual U(T) array.
Stellar SED enabled flag remains 0. The run retains the previous high-mass
numerical BH seed and is not an astrophysical galaxy calibration.

The first `live-dust` launch was rejected before evolution because the tracked
shared `.ic_sink` example erroneously had 13 numeric fields, whereas the native
reader requires exactly 12. The previously successful `live/ic_sink` had 12.
Removed only the extra trailing zero in the shared example. Failed input/log
are preserved; the successful retry used a new directory. No reader relaxation
or AMR/sink-algorithm rewrite was made.

Output schedule: noutput=1, aout=2, tout=1e30 (unreached), foutput=2,
fbackup=1000000, nstepmax=4. Two dumps total **110,228,418 bytes**: the approximate
110 MB forecast was low by 0.23 MB; an ad-hoc exact-110-MB assertion consequently
failed, not the simulation or conservation check. Free space 169 TB. There is
no large follow-on run, checkpoint deletion or production approval.

Stellar SED decision remains open. User was asked whether to allow the supplied
BPASS population as an explicitly independent radiation comparison or require
the exact feedback population and obtain/re-integrate appropriate source tracks.
No answer had arrived during this work; the hard population match is unchanged.
Additional source check: the provider's
[BPASS v2.2.1 converter](https://github.com/galacticusorg/galacticus/blob/e8d9c46113eb2639515ecddb8fb29dea98a3989b/scripts/ssps/convertBPASSv2.2.1SSPsToGalacticus.py)
divides raw 1e6-Msun burst spectra by 1e6, converts L_lambda to L_nu, and uses
Zsun=0.02 for its log-solar metallicity coordinate. Thus per-initial-Msun
normalization has provider-code support; this does not solve the IMF/binary
population mismatch or the absent age<1 Myr spectrum. No scalar IMF correction,
young-age spectrum, or common-population claim was invented.

### BPASS independent radiation comparison: user approval and native connection

The user's subsequent **"그렇게 해요"** approves the proposed independent
radiation-population comparison. This supersedes the pending decision above,
not the earlier scientific limitations. Existing Kroupa/binary feedback,
KL16/LC18 channels and approved effective-SSP SNIa remain unchanged.

Implemented in `snrt_stellar_source.f90`: explicit v2 independent-population
reference mode; v1 exact population matching remains the default. v2 retains
the actual BPASS v2.2.1 bin-imf135_300 population (0.1--300 Msun, slopes
-1.30/-2.35, break 0.5 Msun, built-in binary recipe). The sentinel IDs and
binary fraction -1 cannot be mistaken for the feedback's Kroupa/.5 selection.
Reference mode cannot be promoted by setting approved_production. Source SHA,
population, escape fraction and age/tail policies are appended to the numeric
MPI/HDF5 identity for v2 only; v1 identity length/content is preserved.
Startup records the distinct population and approximations. No main RAMSES
namelist parameter was added, so no generator change is required. The unrelated
generator deletions remain untouched. Makefile VPATH is still lagRamses-first.

`tools/build_bpass_native_sed.py` reads the checksum-pinned local BPASS HDF5
(531204704 bytes, SHA256
`b53d7bf4e8c50ae0a02458eae9d5b6dff5d5782f9f2b23dd40514ad85716c9b3`).
The already-normalized Lsun/Hz per initial Msun is integrated as
Q = integral Lnu*3.827e33/(h*lambda_Angstrom) d(lambda_Angstrom).
No second 1e6 division or scalar IMF reweighting. Group intervals are clipped
to actual source wavelengths before integration; this avoids the unmeasured
low-energy triangular tail of the older candidate ledger integration method.
The old ledger was not overwritten or silently reinterpreted.

Export `config/snrt_stellar_sed_bpass_independent_v2.nml`: 247448 bytes,
SHA256 `ebd5e097e0e9eda50ca029747cb12d89f0a1d89ac02dda665f1381fcccbab840`.
52 ages (explicit young extension plus 51 original nodes), 13 Z, 9 groups.
Sidecar `data/snrt_stellar_sed_bpass_independent_v2.json` records normalization,
provider converter link, actual metallicity axis and limitations. Supplied
escape fraction=1, applied once to rates. Birth--1 Myr holds the first 1-Myr
spectrum as an explicit comparison approximation; no other extrapolation.
Wavelength tails outside 1--100000 Angstrom are zero. Transport still uses
the common reference AGN/stellar grey energies/cross sections. This is neither
a self-consistent common-population model nor production science approval.

Full native build: `.bpass-native.v0ZwR6/ramses_bpass_native3d`, SHA256
`5555c1ba1fb428aa643f88332ccafe8eac277c0e37eae3e0ce9c027c2a598dba`.
HDF5=1 USE_CUDA=1 SNRT=1 DUST_LIVE=1 USE_FFTW=0; NVAR=30, mpiifx/OpenMP,
same root Makefile and preserved old object cache/binaries. Build evidence:
`build-final.log`, native smoke results `stellar-tests.log` in that directory.
The cleanup of temporary debug output affects only the smoke executable.

Focused checks passed: native v1 integration and locator; actual BPASS nonzero
photons, time splitting, young-age hold, age/Z range rejection; analytical
group partition with no missing-tail leakage; independent numerical age
integral vs native at relative tolerance 3e-15; nine invalid inputs rejected
(including missing independent binding, production promotion, wrong population,
NaN fraction, invalid SHA and a v1 IMF mismatch). The first exact-node hard-X
comparison exposed ~1.45e-9 relative error in an extremely weak group under the
optimized weighted interpolation. Exact metallicity knots now take their table
values directly, preventing cancellation/contamination; all groups pass the
original tight tolerance. No tolerance relaxation or source-rate adjustment.

Live run `.bpass-native.v0ZwR6/live/run.nml`: same KL16/LC18/SNIa/AGN/DL01
profile with SNRT_STELLAR_SED explicitly set to the new v2 table. Actual STAR
initial mass, age and metallicity feed native interval integration and existing
cell deposition. Four steps complete, exit 0, elapsed 11.01 s. Active sources
per step: 0, 513, 512, 512; the extra source in step 2 is the active AGN.
No chemistry failure, MG nonconvergence or ERROR. All hydro and RT datasets
finite; all NaN_CHK state counts zero. Existing non-cosmological SFRD diagnostic
still prints NaN (also present before this connection), not a NaN evolved state;
this logging issue is not claimed fixed. The numerical high-mass BH seed and
long first two timesteps remain stress-control choices, not a galaxy model.

Final leaf xHII range 0.9426661027300548--0.999999957882685; dust energy positive,
5.842855172006147e-16--2.443067655624263e-14 code energy per volume. Four RT/IR
commits; maximum dust/IR energy balance 6.1550e-10 (native tolerance 1e-9).
Maximum primary dust-ledger relative error 7.1158e-7. HDF5 stores
stellar_sed_enabled=1 and the complete 6449-value table/population identity.
AGB source is loaded; this short run is NOT evidence of terminal AGB events.

Restart `.bpass-native.v0ZwR6/restart/run.nml` from a preserved COPY of step-2
checkpoint completes to step 4, exit 0. Of 94 hydro/SNRT datasets, 91 are exactly
equal, including all SNRT, dust mass and energy; only the three leaf momentum
arrays differ. Maximum absolute difference 1.6543612251060553e-24; maximum
array-normalized difference 1.1254952277483689e-18. An initial bitwise-all-fields
assertion consequently failed; this is recorded as roundoff-level restart
reproduction, not bitwise equality of the entire hydro state. Source-only SHA
change with identical rates in `restart-source-mismatch` is rejected at HDF5
restore, exit 10, before evolution or new output. Logs/inputs are preserved.

Launch policy (effective namelists audited before execution): noutput=1,
aout=2, tout=1e30 (unreached), foutput=2, fbackup=1000000, nstepmax=4. Live
two dumps total 110145532 bytes. Normal restart adds one dump; negative restart
adds none. Including two checkpoint copies: 275363830 bytes across five dump
directories. Reported capacity estimates: 56 MB/dump, 112 MB live plus 168 MB
restart/copies; free space 168--169 TB. `run-reference.sh` in the scratch build
records the exact reference environment for reproduction. No production run,
new gate ladder, external audit, shared-context update, commit or push was done.

### Operator-directed closeout of the selected comparison implementation

After source commit `073daf6` was pushed, the operator objected to repeatedly
extending the final phase. The agreed closeout was then authorized with
"진행해": fix the selected runnable comparison, retain its limitations and
separate further physical improvements/large-run verification from completion.

The [closeout handover](rt_feedback_dust_comparison_closeout_2026-09-07.md)
now fixes the exact local executable/input hashes, effective namelist and
environment, supported single-rank/OpenMP settings, output budget and restart
procedure. Existing BPASS live/restart/rejection logs and native smoke results
were inspected, all listed executable/input SHA256 values were checked, and
the linked libraries resolve in the current workspace. No physics/settings
changed, so the already passed simulation/tests were not repeated. No new
audit, gate framework, binary, output, source-data move or deletion was added.
The unrelated generator edit remains untouched.

**Selected comparison execution implementation: closed.** This does not certify
long-running/distributed production or publication physics. The previously
suggested spectrum-specific transport/mixture implementation is now optional
future physics work, not an open item blocking this closeout. The handover
proposes a separate combined verification task for later authorization; no
follow-on calculation is running or automatically scheduled by this closeout.

### Operator-authorized MPI/GPU/OpenMP extension

New request: "MPI/GPU/OpenMP 확대 ... GPU/OpenMP의 auto dispatch를 필요한 곳에 구현".
This is a separately authorized execution extension, not reopening the closed
physical-input comparison or adding physical models/gate ladders. Work remains
in `/gpfs/kjhan/LRD_JWST`, remote `kjhan0606/LagRamses`, source HEAD
`677a271d8326fe06a374cac24f0800283e334abd` plus the uncommitted changes below.
The pre-existing unrelated generator deletion is preserved and not included.

Implementation:

- Primary SNRT retains its existing GPU/OpenMP auto dispatch and node-local
  MPI rank/device UUID mapping. Dust now shares that mapping for the implicit
  material energy/temperature solve and group emission. CUDA and OpenMP use
  a shared FP64 scalar kernel; native Fortran remains the independent reference.
- `SNRT_DUST_BACKEND=auto|openmp|cuda` defaults to inheriting `SNRT_BACKEND`.
  `SNRT_DUST_GPU_MIN_CELLS=256` is an initial heuristic, not a benchmark-derived
  optimum. Auto considers device availability, owned-cell count and free memory
  divided among GPU sharers; forced CUDA rejects absent/insufficient resources.
  CUDA failures do not cause a CPU replay after device work has begun.
- Material trials publish only on success, preserving the native enclosing
  transaction and MPI failure reduction. IR advection/outer re-emission
  iteration and mechanical feedback remain host-side. No fake CUDA path was
  attached to sparse shared-cell mechanical deposition. No new main namelist
  field or generic parameter-generator modification was needed.
- `mkrun.py` CLI/GUI adds a separate parallel-comparison Run mode with ranks,
  threads and two backend selectors, setup-only/manual `mpiexec`. It uses the
  new executable, keeps the old fixed profile/binary separate, and explicitly
  limits its `I_MPI_FABRICS=shm` launch recipe to one node.

Two directly related native defects were reproduced and corrected:

1. `backup_sink_hdf5` wrote only rank 0's local cloud drift statistics, although
   BH positions use their global sum. Restore replicated that incomplete value.
   Initial MPI2 restart moved BH x from ~0.4988 to ~0.3296 and changed the AGN
   deposition cell (peak gas-energy relative difference ~1). Grid coordinates
   and order were identical, ruling out output ordering. New output reduces
   statistics globally and marks `sink_stat_global_schema=1`; restore contributes
   the sum once, on rank 0, until per-level drift replaces local statistics.
   Legacy multi-rank sink dumps lacking complete statistics now fail explicitly;
   legacy single-rank dumps remain admissible. Old files are not rewritten.
2. Visible GPUs sent `gpu_hydro=.false.` through `godunov_fine_hybrid`; a GPU
   restart crashed there with SIGSEGV. Initialization, hybrid entry and local
   CUDA flux dispatch now all require `gpu_hydro`. This respects the existing
   option and prevents SNRT/dust device availability changing the hydro path.
   The GPU-enabled hydro hybrid's internal crash is not independently repaired
   or qualified by this CPU-hydro profile. Makefile VPATH order is unchanged.

Build and evidence under `.parallel-runtime.luzQV6`:

```sh
make -C .parallel-runtime.luzQV6 -f ../bin/Makefile -j4 \
  HDF5=1 USE_CUDA=1 SNRT=1 DUST_LIVE=1 USE_FFTW=0 EXEC=ramses_parallel_dispatch
```

Final executable `ramses_parallel_dispatch3d`, SHA256
`11755750b5a4e0a0b71f69946fcf2931ef11cf852c20dafabd6ae41724854829`;
`build-dispatch.log` records the build. The earlier binaries and failed runs
are retained. `run-env.sh` records the common physical environment without
launching a job: unchanged KL16/LC18/effective-SSP SNIa, independent BPASS v2,
reference AGN grey transport and DL01 30/70 material. Physical limitations
from the closeout still apply. CUDA visibility was physical GPUs 1,2, with two
OpenMP threads per rank; other GPU processes were not disturbed.

Bounded checks/results:

- `snrt_dust_backend_smoke`: original Fortran vs OpenMP/CUDA v3 and v4,
  including zero-density cells and material-error rollback. Worst relative
  difference 6.186615e-16 (v3 2.958816e-16). Auto with no visible GPU selects
  OpenMP; below-threshold selects OpenMP, threshold zero selects CUDA; forced
  unavailable dust CUDA rejects startup even when primary backend is OpenMP.
  Two ranks sharing one visible GPU both identify `device_sharers=2` and pass.
- Initial four-step `openmp`, `cuda`, `mpi2-auto` runs completed. Compared with
  the old Fortran-material serial baseline, new OpenMP worst array-relative
  difference was 2.12e-14; OpenMP/CUDA primary RT difference 1.61e-7 and dust
  energy 5.52e-8. Serial tests ran concurrently: their ~21 s vs ~11.6 s MPI2
  times are NOT a controlled speedup measurement.
- `mpi2-restart-openmp` preserves the original incorrect-position failure.
  `mpi2-fixed-auto`/`mpi2-fixed-restart-openmp` reproduce its repair before the
  hydro dispatch fix. `mpi2-fixed-restart-cuda` preserves the hybrid SIGSEGV.
  `mpi2-legacy-reject` rejects the incomplete old MPI sink checkpoint, exit 1,
  before evolution/new output. Intermediate binary SHA:
  `18a66611f67b3b765423629e05301c05a011d191f892d9e11d9fea635c3aa78c`.
- Final `mpi2-final-auto`: 4 steps, exit 0; each rank owns 256 cells, selects
  CUDA for primary/material, and has a distinct GPU (`device_sharers=1`).
  CPU hydro stays outside `HYBRID_ENTRY` despite visible GPUs.
- Final `mpi2-final-restart-cuda`: copied step-2 checkpoint to step 4, exit 0;
  91/94 hydro/SNRT datasets exactly equal to uninterrupted execution, including
  all RT and dust arrays. Only momentum differs: max absolute 4.135903e-25,
  array-relative 2.816452e-19. BH position max absolute difference 2.78e-16.
- Final `mpi2-final-restart-openmp`: hides all GPUs, same checkpoint, exit 0;
  max array-relative difference 1.609612e-7 (RT), 9.618529e-8 (dust energy).
  BH position agrees within 2.22e-16. Arrays finite, dust mass/energy positive;
  no MG nonconvergence/ERROR in the three final runs. These are measured
  backend differences, not bitwise cross-backend equivalence.
- mkrun tests: 27 run, 26 pass, one real-display test skipped on the headless
  host; logs in `mkrun-final-tests.log`. No GUI display validation is claimed.

Each launch used a newly created directory and its actual `run.nml` was
reported before execution: nstepmax=4, noutput=1, aout=2/tout=1e30 unreached,
foutput=2, fbackup=1000000. Final three runs have 110147610 bytes each across
two dump directories (including the copied checkpoint in each restart).
Estimated final storage was 336 MB at 56 MB/dump; free GPFS space was 168 TB.
Failed/intermediate outputs remain intact. No production calculation, source
data download, external audit, shared-context write, commit or push was made.

Scope of completion: native dust-material auto dispatch, single-node MPI2
coupled execution and backend-switch restart are exercised. Multi-node scaling,
calibrated GPU speedup thresholds, GPU hydro, CUDA mechanical feedback and
CUDA IR transport are not claimed completed or automatically scheduled.
