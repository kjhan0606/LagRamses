# CHIMES long-interval numerical repair

User request: resolve the previously failed chemistry interval; do not count
the shorter successful macrostep as repairing it. Work and executions are in
`/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.

## Reproduction and bounded correction

The unchanged input `.effective-coupled.Q9Kw29/retry/run.nml` requests
`courant_factor=.004`; its failing photo interval is exactly
`2166283320464609.2 s` (68.6468 Myr). Diagnostic reproduction is retained at
`.chem-interval.4gVlrR/diagnostic/run.log`. The 157-species BDF/SPGMR solve
repeatedly encounters cancellation-sized negative trace species, including
H2+, and exhausts the existing 64-recovery budget. This is a numerical
failure, not a reason to remove molecules or defer microscopic physics.

A captured actual cell has nH=`7.297525549483983e-4 cm^-3`, 80 directions,
the original molecular shielding/pumping, atomic/molecular spectral banks,
and grain opacity. Its 22,048-byte native input is retained privately as
`.chem-interval.4gVlrR/cell.bin`. A portable run-length IEEE754 hexadecimal
fixture retains its double-precision bits in
`simulation/snrt/tests/fixtures/phase0/chimes_photo_long_interval.hex`.
These are **test inputs**, not retained raw simulation dumps.
Fixture SHA256:
`5710bd69dceb734ca215e60f2d77e150048860b4994a27fa3927636a6b4061af`.
Header is direction count and value count. Each subsequent row is repeat
count plus a 64-bit hex word. Payload order: nH,dt,c_hat,pumping;157 species;
N[9][80];E[9][80];alpha[9][128];shield[2]. There is no phase extension in
this actual failing caller; existing smoke tests cover that API separately.

On this identical input:

| Solve | Result | Measured CPU seconds |
| --- | --- | ---: |
| Original unpreconditioned SPGMR | status51 at t/dt=.3143702324, recovery65, H2+=-4.8148e-35 | 29.1201 |
| Destruction-diagonal left preconditioner | reaches t/dt=1, status0 | 15.6719 |
| Trial restricted to BDF order1 | work limit10000 at t/dt=9.0234e-8; rejected, not adopted | 50.0707 |

Timing is native CPU time on this cell, not a claimed production speedup.
Logs are `replay-unpreconditioned.log`, `replay-preconditioned.log`, and
`replay-order1.log` in the private directory. The replay changed only the
linear preconditioner for the first two rows, not its input or ODE.

`snrt_chimes_photo.cpp` now supplies
`P_ii = 1 + gamma * k_destroy,i` to CVODE's left-preconditioned SPGMR.
The destruction rate sums the **same** atomic shell and molecular
photodissociation coefficients and surviving photons used by the RHS.
Counters and optical depths have identity preconditioning. Off-diagonal
production, secondary-ionization and opacity derivatives are still present
in CVODE's matrix-free Jacobian: they are omitted only from the approximate
preconditioner, not from the equations. The inverse is checked for valid
finite positive denominators and lives in the cell-local solver context.

No species clipping, invented trace floor, species deletion, relaxed
conservation threshold, higher step/retry limit, macrostep reduction or
different radiation source is used. Existing accepted-step positivity and
private-output transaction checks remain. Failure-only messages distinguish
CVODE flags, a recovery limit, and the total work limit; normal recoveries
do not produce production logging. Temporary capture code is **not** in the
production source. No namelist or generator option and no VPATH change.

## Native verification

The portable fixture test calls the actual production photo API, tests
nonnegative finite species/N/E and the existing nucleus/charge/photon/energy
budgets, plus zero-dt identity and invalid-dt rollback. It produces no raw
simulation output.

- Element relative error: `1.738747826517993e-10` (limit1e-8).
- Charge residual per H: `6.882751276668416e-14` (limit1e-8).
- Photon and energy relative residuals: `5.555988861063935e-12` and
  `1.328945225421970e-11` (limits1e-7).
- Native thermal-chemistry smoke: **264 PASS**, no FAIL,
  `SNRT_NATIVE_THERMOCHEMISTRY_OK` (`smoke.log`).
- The checked-in fixture and actual Makefile target also pass:
  `CHIMES_PHOTO_LONG_INTERVAL_OK` (`regression.log`), not just the private
  capture replay. No changes to compiler precision flags or CVODE order.

From the existing isolated `.kind7-dust.PakmEZ` build directory:

```sh
make -f ../bin/Makefile -j1 SNRT=1 HDF5=1 DUST_LIVE=1 DUST_DYNAMICS=1 \
  CHIMES=1 NENER=0 NVAR=199 NVECTOR=32 USE_CUDA=0 USE_FFTW=0 \
  CHIMES_DIR=/gpfs/kjhan/LRD_JWST/.chimes-transition.YdVRvD/chimes \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials-install \
  snrt_chimes_long_interval_test
source ../.effective-coupled.Q9Kw29/environment.sh
./snrt_chimes_long_interval_test \
  ../simulation/snrt/tests/fixtures/phase0/chimes_photo_long_interval.hex
```

## Coupled execution: completed and evaluated

**The reported interval failure is resolved.** Unchanged long input at
`.chem-interval.4gVlrR/long/run.nml`, MPI2/OMP2, four coarse steps, original
star formation, empirical SNIa, Bondi/AGN, RT, CHIMES and C/silicate dust.
All four steps complete (exit0, `Run completed`), including three active
long photo intervals over64 cells each. Four RT transactions and four dust
IR commits pass; no cold-cell rejection, work limit or CVODE failure.
Elapsed run timer920.690s; cooling918.697s. This is a bounded stiff
integration regression, not a recommendation for a69Myr galaxy timestep.

`evaluation-long.json` and empty `evaluation-long.stderr` certify both
step2 and step4 outputs with the **unchanged** evaluator:

| Quantity | Long coupled result |
| --- | --- |
| Maximum element-relative error | 6.351e-16 |
| Maximum normalized charge error | 3.998e-16 |
| Global gas+star+BH mass error, radiated BH mass included | <=2.169e-16 |
| Maximum IR balance residual | 9.2929e-10 |
| Gas temperature at output2 | 50.2536--50.2971K |
| Gas temperature at output4 | 44.7164--44.7544K |
| Actual stars | 64 |
| Cumulative returned stellar mass at output4 | 1.1595904261300422e-7 code |
| Dust mass at output4 | 1.9994373520466815e-6 code |

Gas/grain photon absorption is positive, radiation moments and all157
chemical abundances are finite/nonnegative. Actual accretion and AGN JET
activity occur at steps2--4. No seeded feedback event replaces them.

Final-source short/restart regressions in `short` and `short-restart` both
completed and were evaluated (`evaluation-short.json`, empty stderr):
473 datasets compared,464 bitwise;9 tiny floating differences within the
existing rtol1e-12/atol1e-25 comparison. Hydro/chemistry/dust/RT and all clocks
are bitwise. Element error<=8.48e-16, charge<=5.26e-16; global gas+star+BH
mass including radiated BH mass closes at printed zero residual. Temperature
1604.45--1607.87K. Raw short/restart outputs were removed after evaluation:
414,035,552bytes,8 files; hashes and exact paths in `cleanup-short.json`,
text output metadata retained. Recovery requires rerunning those raw outputs.
The final-source binary is `.kind7-dust.PakmEZ/ramses_chem_fixed3d`, SHA256
`cb37c86878898cfb013b8273e088fbf3a2efd4b41b2c7b47fee403cb01fb9804`.
Later source edits only format/comment the preconditioner and add the test
target, without changing executed arithmetic.
Both live binaries, the native regression executable and final photo source
also have preserved copies in `.chem-interval.4gVlrR`.

The long diagnostic binary is `ramses_chem_preconditioned_final3d`, SHA256
`6b3265162a96195c2e07fd30bf23d985b8e0ce315a52c14c63ce160d106e1af1`.
It includes normal-recovery logging absent in the final binary; identical
preconditioning arithmetic. It passes the original failing interval and
both subsequent intervals over all64 cells and commits each coupled state.
Normal-recovery messages in this diagnostic binary are recovered internal
trials, not published negative abundances or final rejection events.

Both live profiles retain noutput1/aout2/tout1e30/foutput2/fbackup1e6;
~99MiB per dump, two dumps per continuous run. Free space at launch~95TiB.
After evaluation the long raw output was removed:4 files,206,989,104bytes,
exact paths and hashes in `cleanup-long.json`, text metadata retained.
Together with short/restart,621,024,656bytes (~592MiB) were removed. These
raw files require rerunning to recover. Earlier stopped diagnostic/capture
runs created no raw outputs. No active jobs remain. No commit/push performed.

## Scope of the closure

This repairs the demonstrated stiff photo-integration failure without
changing the physical model or admitting previously forbidden model
combinations. Existing photo-then-dark splitting, empirical spectral
closures and separate galaxy calibration/resolution plan remain explicit.
It is not a proof of convergence for every arbitrarily large macrostep.
