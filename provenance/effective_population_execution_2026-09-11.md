# Effective population implementation and coupled closeout

## Step1: mixed v5 + empirical DTD/N100 — PASS

New explicit `population_model='effective_ssp'` (ID2) uses the same v5
single-evolution ordinary source/SED histories, binary_fraction0 as an absent
resolved-binary population, and a baked-in full-initial-SSP empirical DTD.
Only effective mass accounting and empirical event mode are admitted for ID2.
Single-star ID0 still forbids Ia; binary ID1 and strict-WD behavior unchanged.
ID2 is restricted to the user-selected v5 path. No manufactured WD reservoir:
existing effective-SSP checks debit actual remaining stellar mass after all
ordinary channels and verify cumulative prior Ia return. No density added
without particle debit. Restart identity already binds population ID/fraction,
IMF and accounting choices; no new carrier or identity framework required.

Both generators expose ID2 through the shared full editor. The new sidecar
`simulation/snrt/config/fp2_snia_effective_population_runtime_v1.nml` retains
the existing Kroupa Maoz/N100 reference normalization, not an invented Chabrier
conversion. Default IMF remains Chabrier. The live package reweights the same
per-star histories with Kroupa via an explicit private IMF1 history and SED
population2 wrapper; original source data and original wrappers unchanged.

Evidence: `.effective-v5-native-20260911/`, actual-source78648 assertions;
Intel population-contract tests; frontend58 tests, one display-dependent skip.
Live `.effective-population.gQodLz/retry` and `restart-step2`: MPI2/OMP2,
64cells,4steps,86 datasets bitwise identical; positive density/internal energy,
gas+star mass.001code.64 stars at Z=.01, committed age.09484091562Gyr.
Native ordinary return fraction.2323473251193751; empirical Ia events per
initialMsun.0002273111836658904; Ia mass fraction.0003183409915547047.
Measured stellar mass matches both together to1.3605e-16 of initial mass.
Omitting Ia leaves3.1834e-4, so this is a discriminating Ia connection test.
Binary SHA25630927c104aa59631714a34062d3e42bdd5e346432d79bd389c0b11950146ca55.

Preparation failures retained: stale single-star/Chabrier SED wrapper rejected
before evolution; restart initially requested output2 while copying output1.
Corrected inputs use matching IMF/population wrappers and nrestart1 for the
step2 checkpoint. No physical checks or snapshots were changed to pass.
`evaluate.py` compares live masses with the native cumulative source + DTD,
not a separately invented stellar-yield interpolation. Raw cleanup completed:
five HDF files472417600bytes removed, metadata renamed and inputs/binary/logs/
evaluation/cleanup manifest retained. Raw recovery requires a rerun.

## Step2: selected integrated effective profile — configured

C/silicate two-size coadvected grains, CHIMES157 and existing kind7 spectral
coupling; actual stellar sources/v5 + empirical Ia; existing-sink Bondi/AGN
partition reference. NENER0 CPU/OpenMP hydro, HDF5; no CR/MHD/Fe/PAH/drift,
shock-fresh carriers or sublimation in this AGN material layout. These remain
documented other-model limitations, not silently activated combinations.
Growth, sputtering, coagulation/shattering and bounded stellar condensation
use existing operators. No new microscopic closure is claimed.

This is a declared noncosmological effective-model execution domain, not a
calibrated galaxy or permission to extrapolate beyond the source mass/Z/age
and spectrum tables.
The separate science plan is `galaxy_calibration_science_plan.md`.

## Step3: short coupled execution/restart — bounded PASS

Evidence `.effective-coupled.Q9Kw29/short` and `short-restart`.
MPI2/OMP2,64 owned leaf cells, level bounds2--4 with no refined leaves,
NVAR199/NENER0. Final time.00020594code (~.206Myr), maximum stellar committed
age~.154Myr. Ia is admitted but no event is due: its nonzero event/debit
evidence is Step1, NOT this short coupled run. No claim that every source
channel fires simultaneously in this fixture.

Actual stars grow from128 at step2 to256 at step4; stellar returned mass
2.1016243221e-12code. Positive AGN radiation (~2.57e34erg/s) and nonzero jet
accretion occur from the weak existing BH (initial1e-12code). C/silicate
condensation fractions.1,.1,.1; growth/sputtering/coagulation/shattering on.
Grain mass changes from its initial2e-6code; final1.9999976072e-6code. That
combined change is not an isolated calibration of each grain reaction.

Maximum element/charge residuals7.43e-16/4.31e-16; local IR balance7.64e-10.
Gas+star+BH mass with the selected radiative-accretion loss closes to printed
precision. Gas temperature~1604--1608K; density, thermal energy, chemistry and
photon states finite and nonnegative. Local IR balance is not global energy
constancy in a radiating, cooling, accreting system.

Restart:473 hydro/RT/gravity/particle/sink datasets compared;461 bitwise equal.
The12 differences are particle velocities and sink velocity/angular-momentum/
statistics, maximum absolute2.525e-29code, maximum relative-to-peak4.04e-14.
All hydro/chemistry/dust/RT fields, source identities and clocks match exactly.
The evaluator required rtol1e-12/atol1e-25 for nonidentical floating datasets;
it does not report whole-file bitwise equivalence.

Accepted binary SHA256:
`5b445a1db6c871f57af415e553a62b3f15671b49ac9ccd3cba709d022e8f402c`.
Wall times39.929s continuous /29.995s restart. Eight evaluated raw HDF/cooling
files414035552bytes removed; inputs/binary/AGN logs/metrics/text metadata and
cleanup hashes retained. No active jobs remain. No commit/push performed.

### Historical failed long-step attempt (subsequently repaired)

Follow-up user-requested [numerical repair and full long-run evaluation](chimes_long_interval_repair_2026-09-11.md)
now resolves this specific failure. The original interval and four-step
MPI2/OMP2 input complete with destruction-diagonal SPGMR preconditioning,
unchanged rates/tolerances/recovery limits and no abundance clipping.
The failure history below is retained; it is not the current unresolved status.

The initially retained environment restricted RT to level2; all-level IR
correctly rejected it before evolution. `environment.sh` now unsets that
restriction. The first long-step input (courant.004, initial dt~68.65Myr)
then rejected CHIMES photo states with status51 before its first dump.
Counting only consecutive recovery attempts instead of cumulative retries
did not solve it: the trial reached the total integration-work limit
(status4). Both failed jobs were stopped, no raw outputs existed. That
speculative C++ modification was reverted; its diagnostic binary/log remain
as failure evidence. Previous workers' C++ changes were preserved.

Accepted short inputs use the established courant3e-6 and native star-formation
settings. This is not proof of chemical convergence for arbitrary macrosteps,
nor a universal timestep prescription. At that original handoff the~69Myr
photo integration was unresolved and could not be hidden among deferred
microscopic physics. The linked follow-up supplies the missing same-input
repair evidence; the original short run alone did not do so. Neither result
is galaxy calibration or unrestricted macrostep accuracy qualification.

### Reproduction (new run directory required)

The exact evaluated inputs are `short/run.nml`, `short/ic_sink`,
`short-restart/run.nml` and `environment.sh` beneath the evidence root.
Source the environment, run the accepted private binary with MPI2/OMP2 in
a fresh directory containing the copied inputs. For restart, copy the new
step2 `output_00001` to another fresh run directory and use nrestart1.
Audit its output schedule/storage before launching, as required by AGENTS.
Do not attempt to restart from retained `metadata_*`: raw files were removed.

Build recipe (run from an isolated one-level-deep repository build directory,
serial module build, unchanged Makefile VPATH):
```
make -f ../bin/Makefile -j1 SNRT=1 HDF5=1 DUST_LIVE=1 DUST_DYNAMICS=1 \
  CHIMES=1 NENER=0 NVAR=199 NVECTOR=32 USE_CUDA=0 USE_FFTW=0 \
  CHIMES_DIR=/gpfs/kjhan/LRD_JWST/.chimes-transition.YdVRvD/chimes \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials-install \
  EXEC=ramses_effective_dust ramses
```
Local data/libraries and the accepted mixed source package remain required;
this is not a self-contained public data distribution or galaxy calibration.
