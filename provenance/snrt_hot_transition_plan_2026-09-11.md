# High-temperature continuation: focused physical design review

**Current status:** operator approved recommendation B ("권고안대로 진행해").
Implementation and bounded live evaluation are complete; see the final
section. The original review request and preapproval checkpoint below are
retained as history, not a pending approval or current implementation state.

Project / goal: /gpfs/kjhan/LRD_JWST, kjhan0606/LagRamses. Simulation-ready
RT/feedback/dust, not proliferating gates or Python instruments. Operator
requested commit/push (done:9fcab72) and completion beyond the cold-only
10--10^4.98 K admission. Interpret high-temperature coverage within actual
nonrelativistic data/EOS validity (existing atomic network through1e9 K),
not literally infinite temperature. Existing fixed-grain restrictions stay.

Read-only focused scientific design review. FIRST answer Q-GOAL (does the
design deliver the requested goal honestly?), SECOND Q-LEAN (minimum
necessary implementation, no extra framework). Then give a definite
recommendation on the model choice below and essential corrections.
Do not edit files, send messages or launch jobs. Do not spawn subagents.

Already completed: cold157-species molecular + competing grain spectral
photo/dark live receiver, MPI2/OMP2 restart402 datasets exact,191 native
tests. Hot atomic (dust-free) separate mode covers T>1e5 through1e9 K.
The cold mode deliberately rejects thermal trials beyond the actual table
max. Main molecular-cooling logT grid1..4.98. H2-He k0 is unphysical at
logT9 (logk0=7.18), critical-density log=-62468.96. Upstream
zero_molecular_abundances() removes molecular carriers at entry> Tmol;
it does not account for dissociation binding energy. Raising Tmol is NOT
an acceptable fix. A returned successful binary is not scientific validation.

Need determine the leanest defensible high-T extension:

A. Retain a full time-dependent molecular network at high T. This requires
source-supported replacement high-T molecular collision/line cooling rates,
not arbitrary clamping. Is there a realistically complete data route within
the current species set? Identify exact minimal inputs; don't invent fits.

B. Explicit approximate atomic high-T transition, retaining the cold network
below its domain edge and using CHIMES atomic NEQ above. Before any topology
change, conserve each nucleus, preserve net charge, account for molecular
dissociation/ionization binding cost and adjust thermal capacity/temperature.
Use actual thermochemical data, not a generic energy-per-molecule. Follow
thermal crossing as an event/substep, rather than retaining molecular fits
outside their domain in internal CVODE trials. No newly imposed equilibrium
metal ionization. Distinguish rapid-dissociation approximation from a true
kinetic shock model, particularly molecular-rich newly heated low-density gas.
Does this meet the request as an explicit model, or is it a material physical
choice that must be offered to the operator before enabling? Can energy-aware
partial dissociation at the threshold be made consistent without an ad hoc
temperature clamp or a new chemical-energy field throughout the old network?

C. If neither can honestly be completed with existing data, identify the
specific physical choice/data blocking admission, rather than proposing
many new gates or declaring the cold mode to be the requested completion.

Inspect patch/lagRamses/snrt_chimes{.f90,_bridge.c,_photo.cpp,_runtime.f90},
snrt_ramses_driver.f90, and upstream local
.medium-pah-charge.wBiAfM/chimes/src/{chimes.c,chimes_cooling.c,update_rates.c}.
Existing review/decision: provenance/snrt_chimes_molecular_coupling_plan_2026-09-11.md
and ..._implementation_2026-09-11.md. Give file/line evidence for substantive
claims. Prioritize a workable physical decision over implementation cosmetics.

Literature context found by driver: Richings et al2014 original network
https://arxiv.org/abs/1401.4719; recent hybrid-CHIMES
https://academic.oup.com/mnras/article/543/2/891/8241358 explicitly has finite
tables and omits electron-electron bremsstrahlung relevant around>1e9 K;
its high-T equilibrium tables are not full time-dependent molecular rates.

## Driver's independent findings after the request

- Upstream HEAD checked against origin is still
  `a58e5c0311993b51abc63d84fff2958a0104f6d6`; this is not just an obsolete
  local checkout lacking a newer high-T molecular implementation.
- [Official CHIMES examples](https://richings.bitbucket.io/chimes/user_guide/ChimesDriver/Examples.html)
  explicitly retain Tmol_K=1e5 and remove molecules above it, even for
  their temperature grid extending to1e9 K. That grid is not evidence of
  a complete high-temperature molecular kinetic network.
- [ATcT thermochemistry](https://atct.anl.gov/Thermochemical%20Data/version%201.130/)
  offers source-backed0 K energies, including molecular ions/isomers;
  it is a route to binding-energy accounting, NOT collisional rate data.
  [NIST ion thermochemistry](https://webbook.nist.gov/chemistry/ion/)
  explains why mixing0 K ionization energies and298 K formation enthalpies
  needs care. Do not use a generic4.48 eV for all twenty molecules.
- Native diagnostic `.chimes-hot-transition.ddwgtK/topology_probe.f90`
  linked to the current actual bridge/library reproduces the legacy
  topology energy defect: pureH with n(H2)/nH=.1, n(HI)/nH=.8,
  T=101000 K, nH=1 cm^-3, no photons, dust or CR, dt=1e-12 s returns
  status0, H2=0 and unchangedT. Thermal energy increases by
  **0.11111111111111110** solely because particle number increases.
  See `topology.log` and `build.log` in that directory. This is not a
  RAMSES calculation and produces no raw simulation snapshots.
  Existing cold/hot spectral admission avoids this molecular input case;
  blindly relaxing those protections would expose it to the new live path.
- Upstream chooses mol_flag once at entry (`chimes.c:227`), while
  `rate_equations.c:549,584` updates trialT and evaluates rates at each RHS.
  A valid transition must therefore handle internal thermal crossings,
  not just inspect the temperature returned after the whole hydro step.

Operator clarification was requested: an energy-aware rapid-dissociation
approximation changes shock kinetics relative to a fully time-dependent
high-T molecular network. Neither the approximation nor a literal
unbounded-temperature completion is silently activated by this record.

## Review received and driver disposition

[Fable review](snrt_hot_transition_fable_2026-09-11.txt) completed successfully
with no CLI error. It recommends B as an explicit opt-in physical model,
not a silent relaxation of kind6. The driver agrees with the need for a
source-backed reaction-energy debit, charge-preserving atomic products,
temperature recomputation, explicit model/restart identity and handling of
thermal crossings. A successful cold-only run is not feedback-heated-gas
qualification; the prior completion was explicitly bounded.

Do NOT adopt the review's stronger unsupported generalizations:

- Its dissociation/cooling ratio is not established for arbitrary density,
  metallicity, charge fractions, grain inventory or radiation fields. A
  finite-step low-density molecular shock remains a distinct physical case.
- Its statement that the He critical density already underflows at logT6
  is incorrect in FP64: log10(ncrit)=-57.27 is representable. The retained
  table becomes genuinely unrepresentable further up; rate validity is a
  separate issue from floating-point range.
- A status50 caused by an internal CVODE trial does not by itself locate a
  physical accepted-solution crossing. Bisection solely on that status may
  mistake trial overshoot for an event. Any eventual transition implementation
  must establish the actual crossing state and must not atomize earlier at
  a lower-temperature entry merely to make the solve pass.
- Identical endpoints for differently timed entry/photo/dark events are
  not generally a valid test: their remaining chemical evolution differs.
  Check conservation at each event and timestep convergence instead.
- The review did not inspect the HDF5 values directly. Its assertion that
  no high-T database exists anywhere is broader than this investigation:
  no complete compatible source was established here, not a universal proof
  that no relevant data exist.

The operator's model-choice question is still unanswered at this checkpoint.
No unapproved instantaneous conversion, table clamping, Tmol increase,
new runtime selector or broadened production claim was implemented. The
completed cold-coupling commit9fcab72 is pushed and verified on origin/main;
these new design/evidence records are local. No RAMSES job was launched,
so no additional simulation raw-output cleanup is due.

## Approved implementation and live evaluation

Operator approval selects B, not high-T molecular kinetic extrapolation.
Implemented opt-in `chimes_transition_d03_maxent128_fs2010_v1` (kind7):

- All twenty molecular species use source-backed atomization energies
  from [ATcT v1.130](https://atct.anl.gov/Thermochemical%20Data/version%201.130/),
  **0 K**, not mixed 0/298 K enthalpies. Original IDs and transcription are
  in `simulation/snrt/data/chimes_atomization_atct1130.json`; compiled values
  in `patch/lagRamses/snrt_chimes_atomization.h`. H3+ uses para ground state,
  CH2 triplet and C2 singlet. Fragment charge is assigned to the lowest
  ionization-cost constituent atom; this is an explicit channel assumption.
- Preserve each nucleus and charge/electrons, subtract binding energy from
  thermal E, then recompute T with the existing translational heat capacity.
  A failure leaves the caller's transaction unpublished. No arbitrary
  constant per molecule and no thermal energy creation on topology change.
- Entry, photo heating and accepted CVODE dark upward roots all connect to
  the helper. ABI6 adds optional hooks without changing upstream structures.
  `CVodeSetStopTime` prevents roots beyond the requested interval. Boundary
  rate continuation is only for off-domain solver trials, not accepted
  high-T molecular evolution or a physical table extrapolation.
- Complete the remainder with atomic NEQ; cold molecular formation may
  resume next call. Finite 10--1e9 K domain, first-order splitting, instantaneous
  rather than kinetic dissociation. Existing kind5/kind6 remain unchanged.
- Actual RAMSES uold chemistry/thermal state, scattering-only transport,
  competing gas/grain primary absorption and material/IR commit are wired.
  Fixed carbon/silicate two-size grain comparison only: no Fe/PAH, mass
  processing, relative motion or sublimation; not high-T grain survival.
- Native restart14 / chemical identity7 / tested HDF5 format64 bind the new
  selector and ATcT SHA256. Carrier width unchanged. `mkrun.py` and generic
  generator/GUI document and export the opt-in. No new namelist variable.

### Evidence retained

Private root `/gpfs/kjhan/LRD_JWST/.chimes-transition.YdVRvD` retains
inputs, ABI6 build, original ATcT HTML, native/frontend/metadata logs,
compact `results.json`, evaluation and fixture scripts, and cleanup manifest.
These scripts prepare/measure the test; runtime physics is native C/Fortran.

- Native thermochemistry: **251 PASS**, including twenty molecular
  energy/nucleus/charge closures, unaffordable conversion rollback, finite
  domain, entry/photo/dark crossings, and old-mode regressions. Dark root
  located at t=106164.497829073 s, T=95499.2627963112 K. The earlier photo
  fixture did not supply enough photons to cross; its corrected exposure
  crosses to 433756.6483 K with one event and explicit energy closure.
  `native3.log` is the final result; earlier logs are retained diagnostics.
- Restart/identity metadata smoke passes, including altered identity refusal.
- Frontends: 49 tests, 48 passed and one display-dependent skip.
- Actual hydro/CHIMES/D03/advective-CR run: 4^3 cells, NVAR187, NENER1,
  MPI2/OMP2, no CUDA. Cold/hot molecular initial regions, four steps,
  32 atomization events. Stars, sinks and AGN intentionally inactive.
- Controlled copied-checkpoint heat pulse in two cold cells plus nonuniform
  LW/HI/hard photons, four resumed steps: two additional atomizations,
  nonzero gas and grain absorption, finite positive chemistry, and unchanged
  integrated grain masses. This validates the thermal receiver, NOT an
  actual SN/AGN source calculation or galaxy calibration.
- Dark and irradiated intermediate restarts each reproduce **402 datasets
  bitwise** (hydro, SNRT, gravity, particles). The irradiated final gas
  spans 1622.085--556153.655 K. Primary photon energy budget relative error
  **1.25423e-8**, maximum IR commit imbalance **1.0437e-13**.
- Effective inputs are `mixed/physical.nml`, `mixed-restart/physical.nml`,
  `irradiated/physical.nml`, `irradiated-restart/physical.nml`. All schedules
  audited before launch: noutput1, aout2/tout1e30 outside interval, foutput2,
  fbackup1000000. Each measured dump about9.8 MiB, free space101 TiB.
  Raw snapshots are removed after completed evaluation per operator policy;
  their hashes/byte counts and exact targets are retained in the manifest.

Build identity: binary `.chimes-band-live.PmDxvQ/ramses_transition3d`, SHA256
`f221ecb9fd1a1a1cd7b85331fc219f6d5cc20b3592f105c092b0acd0b33ce179`;
ABI6 `chimes/libchimes.so`, SHA256
`c332140e4447be0916f4602172df5f6c7b133e3a5178dfeba5586fb181fa8a78`.
ATcT JSON SHA256
`0bfa808eda56b7ea41b1a9e083e5379106caa2960f4b8a242f6e2e902606a045`.
Atomic/molecular/main data identities remain those of the approved cold
comparison. Makefile VPATH unchanged; the new header is an explicit build
dependency. Durable upstream receiver patch passes reverse-apply check
against the tested ABI6 tree (upstream context whitespace preserved).

No additional model, automatic audit, production run or new completion gate
is created. This completes the approved B connection and bounded integration;
it does not qualify full high-T molecular kinetics, relativistic EOS, dust
destruction/survival, cosmological calibration or every feedback source.
