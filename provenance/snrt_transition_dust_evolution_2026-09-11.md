# CHIMES transition plus existing grain mass evolution

Status: completed implementation and bounded driver end evaluation.
The initial request below is historical; Fable conditions and results follow.

Operator: implement the next bundle. Project / goal:
`/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`, simulation-ready native
RT/feedback/dust without multiplying instruments or gates.

## Focused plan review request

FIRST answer **Q-GOAL**: does this complete a necessary existing connection
toward the final goal? SECOND **Q-LEAN**: is this minimal, without excess
instrumentation or gates? Then APPROVE / CONDITIONAL / REJECT with concrete
physics/conservation corrections. Read-only; no files changed or jobs run.

The operator approved energy-aware rapid dissociation to atomic CHIMES
above the molecular data ceiling. Kind7 is implemented and native/live/
restart-tested on 10--1e9 K, fixed D03 C/silicate grains. See
`provenance/snrt_hot_transition_plan_2026-09-11.md` final section.

Next connection: enable the existing co-advected two-size C/silicate grain
growth, Tsai--Mathews comparison sputtering, coagulation and shattering in
kind7. Do NOT implement a new dust law or high-T molecular rate extrapolation.
Keep kind5/6 restrictions and all defaults. Kind7 with all rates off must
retain previous behavior. No Fe/PAH, drift, sublimation, condensation or SN
shock injection admission in this bundle. Actual source injection is separate.

Existing `amr_step.jaehyun.f90` calls `dust_mass_advance_level` before the
RT receiver. Existing mass routine collectively stages grain bins, material
energy exchange and CHIMES gas-carrier reconciliation; total elements remain
the conserved gas+solid sum. Existing runtime reads current bins into D03
opacity and uses gas carriers for heat capacity; no frozen opacity cache.
Inspect `dust_mass_runtime.f90`, `snrt_chimes_runtime.f90`,
`snrt_ramses_driver.f90`, `dust_mass_physics.f90`. Important: grain accretion
must reserve molecular nuclei; erosion returns neutral atoms under existing
CHIMES convention without changing surviving ions/molecules or creating E.
No claim to include solid binding energies absent from the existing model.

Implementation: restricted admission change; preserve/repair explicit donor,
energy and opacity wiring as needed; frontend opt-in that leaves old fixed
profile unchanged; existing restart dust parameter identity must distinguish
rate settings. No new per-cell carriers or generic framework. Existing dust
native regression plus one actual mixed cold/hot gas evolution with nonzero
photons and intermediate restart; check elemental, gas+grain material E
exchange locally, updated opacity, no negative states and exact restart.
Keep results/input/log/hash metadata; remove evaluated raw snapshots.

Judge consistency and indispensable issues, not an exhaustive new roadmap.

## Fable disposition and implemented scope

[Focused review](snrt_transition_dust_fable_2026-09-11.txt): CONDITIONAL;
Q-GOAL yes, Q-LEAN yes. Reviewer writes "another session" because it observed
the live working tree: these are this driver's changes, not an ownership
conflict. No reviewer output file write was allowed; its text verdict was
returned successfully and retained here. No follow-up per-helper review.

Accepted the total-H collision density correction, genuinely active dense
cold/hot fixture, explicit existing depletion convention, and lean evidence
retention. The review's phrase "kind5 in H2-bearing cells" is incorrect:
kind5 is hot/atomic and excludes such cells. The correction applies to the
generic grey CHIMES mass path and admitted molecular paths, not admission of
H2 into kind5. It changes neither C/olivine molecular reservations nor PAH
hydrogen chemistry: H is only used as the collision-density parameter here.

The initially proposed neutral-only growth helper was removed: it would
unnecessarily introduce a different donor model. Retain proportional atomic
ion depletion, charge-reconciled electrons, and neutral erosion products.
No new solid binding/recombination energy model is inferred. The preliminary
`mixed` run used that unselected helper; it is NOT final physics evidence.

Kind7 now admits existing four mass/size switches. Condensation, SN shocks,
Fe/PAH, sublimation and drift remain rejected; kind5/kind6 retain fixed
masses. `mkrun.py` exposes a default-off explicit choice and generates all
four switches, with condensation/shocks off. Shared generator/GUI explains
the same restrictions. Existing dust restart identities bind these choices.
No new namelist field, restart version or carrier tensor. Existing RT reads
current bins; no opacity-cache invalidation mechanism was needed.

### Numerical issue exposed by the combined run

The initial dense irradiated run correctly rejected. Independent diagnostic
showed photon survival entries in CVODE's internal state as well as its
returned state were negative (`-1.0107936529880487e-174`, also
`-2.7727785613283782e-272` at normalized t=1). It was NOT merely dense-output
interpolation and was not a negative gas abundance. The prior stop-time fix
alone is insufficient for an exhausted survival ODE.

The mathematically equivalent optical-depth variable now evolves
`d tau/ds = c dt [sum_i(n_i sigma_i) + alpha_dust]`, with survival
`exp(-tau)` used in every atomic/molecular/grain capture rate. Species and
all capture/energy counters still evolve together. Trial tau is constrained;
accepted negative states still reject. No published negative clipping.
Native regression showed 1e-9 local tolerance was insufficient for the
nonlinear survival versus integrated ledger: tightened to 1e-11, with
correspondingly tighter absolute tolerances. Acceptance budgets unchanged.

That exposed a separate storage boundary: a positive FP64 photon survivor
can round to FP32 N=0 while leaving FP64 energy. Existing collective commit
correctly rejects that malformed pair. The live adapter now rounds only
subnormal FP32 survivor packets together at fixed E/N, requiring the sum
of absolute energy changes <=64 FP64 epsilons of incoming cell radiation.
Significant-tail inputs still reject unchanged. This is quantified storage
roundoff, NOT additional absorption or gas heat. No general photon floor,
untracked physical dissociation cost or new per-step grain diagnostics.

## Evaluation inputs

Root: `/gpfs/kjhan/LRD_JWST/.chimes-dynamic.yTEzYs`. All native source
changes use the existing CPU/OpenMP build and preserve Makefile VPATH.
CHIMES ABI6/data remain from the preceding transition bundle. Existing
native mass smoke (growth, sputtering, size exchange, material energy,
element reservation) passes; frontend 49 tests, 48 pass plus display skip.
Native thermochemistry now includes dense photon exhaustion and paired
storage-tail acceptance/rejection cases; final result: 254 PASS.

The dense fixture uses 4^3 hydro, NVAR187/NENER1, MPI2/OMP2, advective CR,
no self-gravity/cosmology, stars/AGN/sinks inactive. Carbon and olivine both
have two nonzero size bins. Cold dense half has n_H~2170 cm^-3; hot diffuse
half starts above the molecular ceiling. Four steps with mass evolution,
then nonuniform LW/HI/hard photons inserted into a private midpoint copy
for four resumed steps and an intermediate restart. This is not an actual
stellar/AGN source, galaxy calibration or hydrodynamic convergence study.

Schedules audited before launch: noutput1, aout2/tout1e30 outside the test,
foutput2, fbackup1e6, measured dump ~9.8 MiB, available storage ~100 TiB.
Diagnostic one-step run stopped at step3 and produced no additional dump.
Failed trials and their logs are retained as diagnostics, not successful
evidence; raw output is removed only after the final evaluation.

## Completed end evaluation

- Final native `complete-native.log`: **254 PASS**, including all prior
  atomic/molecular/transition tests and three new dense-exhaustion/storage
  cases. `mass-final.log` passes the existing mass/material/size tests.
  `frontend-complete.log`:49 tests,48 pass,1 optional display skip.
- `dense` four-step integration and `dark-final-restart` two-step continuation:
  **394 physics datasets bitwise identical**. The count is eight fewer than
  the prior 402 because Poisson/gravity is intentionally inactive. This also
  checks that the photo solver replacement preserves the no-photon result.
- `radiation` four resumed steps and `radiation-restart` two-step continuation:
  **394 datasets bitwise identical**. Both complete successfully; all states
  finite and nonnegative. Final maximum element residual2.96e-16, charge
  residual2.16e-16. Maximum IR commit balance1.2372e-10.
- Starting dust/rho=.002 in both regions. At dense step4, cold dust/rho is
  .00229604852 and hot .00199516322. At irradiated step6 these become
  .00247638083 and .00199500129. Both C and silicate bins evolve; cold
  grains grow and hot grains erode. Final temperatures51.168--552289.157 K.
- Cold step4 total n_H=2173.474 cm^-3 versus molecular-depleted970.731:
  the previous erroneous density would select shattering instead of
  coagulation here. Step6 depleted value694.581, still on the wrong side
  of the threshold without the correction. Native tests independently
  cover both size-exchange directions and the material exchange budget.
- Primary photon energy decrease vs gas plus grain capture: relative
  error**1.68111e-8**. Grain capture is nonzero every step and changes from
  2.62396e68 to2.70352e65,2.72570e65,2.74344e65 eV. This variation reflects
  both evolving radiation and opacity, not a controlled opacity-only
  measurement. Driver inspection confirms opacity is rebuilt from new bins.
  These are whole artificial-box totals, not luminosities per cell.
- No global gas+grain+radiation closure claim: optically thin CHIMES cooling
  leaves the modeled system without an explicit transported loss reservoir.
  Material exchange conserves represented energy natively; no new global
  cooling ledger was added just to manufacture such a claim.

Final executable `.chimes-band-live.PmDxvQ/ramses_dynamic_complete3d` SHA256
`0329c433a9c616a9bb7a475966d22f974fd37bf9b2dab59fc48010c5f64ec58e`.
Dark original uses the same mass code before the photo-only numerical repair;
final no-photon replay is exact. ABI6 library hash remains
`c332140e4447be0916f4602172df5f6c7b133e3a5178dfeba5586fb181fa8a78`.
The final four-step radiation run took173.4 s; this is a stiff, extreme
density-contrast correctness fixture, not a production throughput benchmark.

Source diff whitespace check passes (vendor patch context excluded; upstream
reverse-apply check is retained from the preceding bundle). Effective nmls,
environment, build logs/hashes, seed.py, evaluate.py and results.json remain
in the private evidence root. `cleanup_manifest.json` enumerates exact raw
targets and HDF5 hashes, including failed/unused diagnostic copies. The raw
outputs have no retained backup and require rerunning these fixtures to recover.
Cleanup complete:18 exact directories,183,157,514 apparent bytes (~174.7 MiB).
No previous bundle, shared physical bank or production output was removed.

No new default, physical model, mandatory audit or implementation queue was
created. Rate-off mass settings retain their previous meaning; the equivalent
photo ODE repair can change floating-point values for nonzero radiation, so
bitwise agreement with older illuminated binaries is not claimed. No commit
or push was requested/performed in this continuation.
