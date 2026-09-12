# Actual-IC coupled startup and cost pilot

Operator requested the next task after CMB/material integration passed.
This is the first bounded coupled evolution of the existing 128^3 IC, not
a new physics model, calibrated galaxy sample or source-active benchmark.
No new audit is required by the approved audit cadence.

## Effective inputs and physical scope

- Run directory: `/gpfs/kjhan/LRD_JWST/.galaxy-coupled-pilot.ZGTFWW`.
- Effective NML: `run.nml`, SHA256
  `08ff16e8893ff28d7140885062bae3e9ba7170421944a0bb0f8ced60e4d2283a`.
- Pinned executable: `.cosmo-coupling.g8Heia/ramses_cosmo_cmb_v3_3d`, SHA256
  `cfb8b5148e75dd48e4dbd17c2d20dcf47a7d38d1a2f23c3da28a97add64bfcf5`.
- IC: `.galaxy-preenriched.IJ8aeX/ics/level_007`; original MUSIC density,
  velocities and displacements, explicit twelve-element passive files,
  approved Z=1e-10, initially no grains. No new IC generation or resampling.
- 12.5 cMpc, initially z=99; fixed level7, gas+DM, gravity, hydro, kind7
  CHIMES, four C/silicate dust masses and 136-frequency/80-direction live IR.
- Physical CMB `2.727/a`, explicit external material receipt, fixed-group
  a^-3 IR dilution with no spectral redshift, as previously approved.
- No SF, sink/AGN, CR/SGS, Fe/PAH or grain drift. Stellar source and AGN
  environment selectors explicitly unset; yield/population binding remains
  the registered Kroupa effective-SSP reference, not a changed global default.
- Grain growth/sputtering/coagulation/shattering are enabled but an initially
  grain-free, source-free run should not be called a grain-production test.
  The separate dusty manufactured test already exercised nonzero CMB receipt.

## Resource and retention bounds

Submitted Slurm job540015, normal partition: four nodes, eight MPI ranks per
node, four OpenMP threads per rank (128 allocated CPUs), 400 GiB requested
per node. Inter-node MPI is `shm:ofi`, replacing the source environment's
single-node `shm`. Pinning is core/close, OMP dynamic false, 512 MiB thread
stack. This is a bounded 20-minute request, at most42.667 allocated core-hours,
not permission for a parameter ensemble or a long cosmological production run.

The minimum global memory for four dense IR fields is680 GiB. The1,600 GiB
request includes additional room for persistent slot capacity, halos, primary
RT, gas/chemistry arrays and runtime temporaries; peak RSS must be measured,
not assumed from this estimate. grammar-debug has only252 GiB and is therefore
not used for this dense128^3 path. No angular/spectral/IC resolution is reduced.

Stop at nstepmax2, expansion limiter1e-6, no refinement/remapping. Schedule:
noutput1, aout1.1, tout1e100, foutput=fbackup1000000. Expected full dumps: zero.
An uncompressed full dump would be roughly200--300 GiB including persistent
IR; it is not scheduled. Expected new storage below100 MiB (logs/manifests),
shared GPFS free space approximately94 TiB at submission, not a quota.

The wrapper requires native completion, at least64 rank-local CMB commit
records and two IR commits, and rejects known error/nonconvergence messages.
It also requires no output_* directory. Evaluate actual density/NaN,
conservation, thermal/chemistry and phase timing records before interpreting
the result. Separate measured Slurm compute TotalCPU from job allocated CPU;
do not double-count job and step records. No output cleanup is needed unless
raw files were actually generated; inputs and failed-run logs remain preserved.

## Status

Job submitted; completion, physical acceptance, peak memory and timing are
not established by submission. A source-active late/refined interval is still
needed before generalizing this startup cost to galaxy parameter calibration.

First scheduler check: PENDING(Resources), estimated start
2026-09-12T10:15:13 (scheduler-local time, tentative and subject to change).
No physics failure is implied by queueing; no second run or downgraded
resolution was submitted to bypass this resource wait. The submitted batch
will execute both steps and its completion checks without another approval.
