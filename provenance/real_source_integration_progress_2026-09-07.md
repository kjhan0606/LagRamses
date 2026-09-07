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
