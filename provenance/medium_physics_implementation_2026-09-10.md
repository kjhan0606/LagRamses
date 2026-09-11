# Medium-term physics implementation

Operator: implement the eight medium-term physics groups listed after the
bounded gas-MHD completion. This is not authorization to invent missing
stellar data or to claim universal production qualification.

## Scope and ordering

Continuation2026-09-11: the preapproved remaining bundles1--5 are tracked in
[the execution record](remaining_physics_bundle_execution_2026-09-11.md).
Bundle1's kind7 C/silicate condensation, shocks, sublimation and relative
dynamics connection has [completed native/live/restart evaluation](kind7_dust_process_implementation_2026-09-11.md).
This supersedes older limitations for that bounded connection below, not
the outstanding PAH/Fe or stellar source-data limitations.

Operator continuation (2026-09-10): complete group 7 and then all remaining
numbered medium-term groups; implementation is preapproved. Continue the
existing work streams below, without per-helper approval/audit gates. First
repair the exposed cold absolute-emission material bracket using the actual
DL01/Fe material model and band Planck law, not a temperature floor or a
shortened timestep. Remaining source-dependent physics must use matched
physical inputs; approval does not supply absent yields or binary histories.

1. Native dust/radiation improvements: metallic-Fe mass evolution, followed by
   the data-supported PAH/spectral/angular extensions and a collision-aware
   two-size/multibin decision. Preserve currently approved defaults.
2. Stellar population/source improvements: common-population SED and returns,
   microscopic Ia histories, missing mass/Z intervals, time-dependent wind
   composition/velocity and decay, and alternative massive-star/PPISN/PISN
   sources. Use matched, traceable input data; do not relabel independent
   BPASS and LC18 populations as one population.

These are work streams, not additional approval gates. Dust--magnetic-field
coupling, anisotropic CR transport, and nonideal MHD remain out of scope.

### Original numbered scope: current disposition

| Group | Current disposition (not completion claims) |
| --- | --- |
| 1: same-population SED/feedback | Bounded [PARSEC high-mass common-population connection](parsec_common_radiation_implementation_2026-09-10.md) passes native Intel/GNU and MPI2 exact restart: same14--600Msun nodes, IMF cells, Z and terminal ages for Q/E and eleven-element feedback. [Five-Z extension](parsec_metallicity_extension_implementation_2026-09-10.md) now adds .017/.02/.03 to the same matched source and passes native/MPI2 exact restart. [Eleven-Z printed-precision/RATE extension](parsec_low_z_precision_implementation_2026-09-10.md) now reaches1e-11--.03 on495 matched nodes; GNU native and MPI2 exact restart PASS after repairing Intel FTZ/DAZ; evaluated raw outputs deleted. Explicit Q5/Planck atmospheric closure and full-track Planck missing tails; NOT a full SSP. Below14Msun radiation/returns remain unmatched. BPASS and legacy LC18/Monash are still independent populations. |
| 2: microscopic Ia population | Still incomplete. Keep effective DTD/N100 available. The targeted COSMIC grid has illustrative primary-IMF/system weights, but is not a matched individual-star population or a near-Chandrasekhar/N100 event sample; revised retention requires coupled binary re-evolution. |
| 3: missing stellar mass/Z intervals | Still incomplete. Existing low-Z/AGB7 and solar 9--13 extensions remain. [PARSEC v4](parsec_pair_feedback_implementation_2026-09-10.md) supplies an opt-in nonrotating14--600Msun five-fate alternative, now extended to [Z=.008/.014/.017/.02/.03](parsec_metallicity_extension_implementation_2026-09-10.md) with matched photons, unchanged original-node payloads and exact native/live source binding. [Eleven-Z comparison](parsec_low_z_precision_implementation_2026-09-10.md) adds selected low-Z branches with gross yields unchanged and printed-bound baryonic reconstruction; excluded author1e-6/.001 are explicitly neighboring-Z mixtures, not their original fates. Z=0, unmatched lower masses and other populations remain incomplete. |
| 4: phase winds/net/isotope decay | [PARSEC phase-wind energy](parsec_phase_wind_implementation_2026-09-10.md) connects actual M/L/Teff/XH to a named cool/hot/H-poor escape-speed comparison, with positive eleven-element histories and exact endpoints; Intel/GNU and MPI2 exact restart pass. [LC18 prompt decay](lc18_prompt_decay_implementation_2026-09-10.md) now corrects fast Ni/Co-to-Fe at source release, with unchanged native mass/energy/AGB/Ia budgets and default. This is a half-life<=100yr coarse-grained model, not full/live radioactive evolution. Phase classification/LBV/bi-stability and long-lived isotope/dust-lattice evolution remain limited/incomplete. |
| 5: alternate massive-star/P(P)ISN | Bounded [PARSEC five-fate native implementation](parsec_pair_feedback_implementation_2026-09-10.md) passes physical-source tests and MPI2/OMP2 thermal+CR evolution/exact restart. Pair ejecta/energy have their own channel, not a second remnant. Other Z/rotation populations and resolved pulse chronology are not supplied; this is not a general population closure. |
| 6: collision-aware size-distribution choice | Bounded native comparison/decision complete: [collision evidence](dust_collision_comparison_2026-09-10.md). Coagulation/fragmentation and 16/32/64-bin comparisons implemented. Live two-size retained; no equivalence or live multibin claim. Diffuse boundary fate remains a limitation. |
| 7: PAH/Fe | Fe kinetics and cold radiation are implemented/tested. Fixed-H PAH charge/IR/CHIMES coupling and restart are tested. Optional H0--13 comparison connects H loss/attachment, state-dependent masses and binding energy through hydro. [H2 vacancy-refilling capture](pah_h2_implementation_2026-09-10.md) now adds a cation-only, M13-bound-rate subchannel with finite CHIMES donor and native/live conservation/exact restart. Carbon-skeleton destruction, H2 formation/single-vacancy/superhydrogenation pathways, higher charge states and general hard-photon survival remain incomplete. |
| 8: spectral/angular RT | Angular refinement implemented: 80/320/720 choices and [angular evidence](snrt_angular_implementation_2026-09-10.md). Native opt-in H/He intragroup N/E reconstruction, evolving Verner absorption, actual energy into chemistry and MPI2 exact restart are [implemented/tested](snrt_band_implementation_2026-09-10.md). BPASS Q/E reaches native age/Z integration and injection: [stellar energy evidence](snrt_stellar_energy_implementation_2026-09-10.md). Optional [node-resolved FS2010 secondary electrons](snrt_secondary_spectrum_implementation_2026-09-10.md) now pass native physical tests and MPI2 exact restart. Fixed/v2 defaults retained. This does NOT resolve group1's population mismatch; [D03 grain128-node absorption/scattering/heat](snrt_d03_spectrum_implementation_2026-09-10.md) now passes native and MPI2 exact restart for static four-bin grains. [Static Fe six-bin spectral connection](snrt_fe_spectrum_implementation_2026-09-10.md) now passes native and MPI2 exact restart with a synthetic sub-eV source; CHIMES-free static layout and missing Fe aggregate face flux are fixed. Retained <=4eV/300K limits exclude general stellar/AGN Fe heating. CHIMES reaction-table spectral coupling, general Fe photoelectron/charge closure, PAH and relative-motion spectral extensions remain incomplete. |

Group 8 continuation: the [CHIMES atomic spectral building block](snrt_chimes_spectrum_implementation_2026-09-10.md)
now supplies native directional N/E moments for 311 reactions / 682 partial
shells, with Fortran bindings and passing existing chemistry smoke tests.
The [conservative photo operator and hot-cell CHIMES connection](snrt_chimes_photo_implementation_2026-09-11.md)
now evolve finite photons, atomic species and shell-wise secondary energy
with CVODE and connect to actual nonradiative chemistry/cooling. Native142
checks and first-order dt convergence pass; cold molecular, grain competition
and live driver/restart wiring remain unfinished. No runtime guard was removed.

**2026-09-11 update (supersedes the preceding native-only checkpoint):**
the [cold molecular connection](snrt_chimes_molecular_coupling_implementation_2026-09-11.md),
[energy-aware hot continuation](snrt_hot_transition_plan_2026-09-11.md), and
[existing grain-mass evolution connection](snrt_transition_dust_evolution_2026-09-11.md)
now include native/live MPI2/OMP2 and restart evidence. CHIMES atomic/molecular
reactions and competing C/silicate grain spectra are no longer an unconnected
building block. Kind7 can use existing growth/sputtering/size exchange, with
explicit rapid-dissociation and depletion conventions. Fe/PAH/drift spectra,
condensation/SN shocks and sublimation in this particular path remain excluded;
other numbered population/data constraints are not silently marked complete.

No remaining group has been silently reclassified as complete, rejected, or
awaiting another operator approval. The full requested implementation is
still unfinished. "Data-dependent" is not a claim that all public sources
are unavailable: ingestion, model selection and native coupling still have
to be done. Sources checked in this continuation include
[BPASS outputs](https://bpass.auckland.ac.nz/9.html),
[super-AGB sources](https://www.astro.ulb.ac.be/~siess/pmwiki/pmwiki.php/StellarModels/SAGB),
[HW02 helium-core yields](https://2sn.org/DATA/HW01/), and the PAH sources below.

## PAH carbon-fragmentation disposition (2026-09-10)

[Focused Fable review](pah_fragmentation_fable_2026-09-10.txt) does not
approve immediate live carbon destruction. Driver agrees with the concrete
gaps: fragment chemistry (C2H2/C2H absent from CHIMES157), daughter optical
coverage below C22, and reaction enthalpies distinct from activation barriers.
[Acquired sources](pah_fragmentation_sources_2026-09-10.md) do not yet close
these gaps. No zero rates, inert C22 endpoint, or invented atomization was
added. Source selection remains our work, not a new operator approval gate.
Group 7 stays unfinished. Proceed with the already-approved group-6 native
collision-aware two-size/multibin comparison; no new live carrier tensor.

## First implementation: metallic-Fe mass evolution

Deliver a selectable native neutral-grain comparison using the existing two
Fe mass carriers, gas-phase Fe depletion, CHIMES, and dust enthalpy. Compare
with the currently fixed-Fe baseline. Growth must consume gas Fe only after
olivine's reservation; erosion must return the same element to gas. No
renaming C/silicate physics as an independently validated Fe sputtering law.
Select the rate law only after checking a primary source. Unknown Fe rates
or new required material inputs are a scientific dependency, not zero rates.

Include the namelist generators, restart parameter identity, conservation,
zero-rate identity, and a small native integration in this implementation.
If relative dust motion is supported, transfer donor momentum and account
for kinetic energy; otherwise reject that combination explicitly, never
silently evolve only its masses. Do not relax the existing Fe optics bounds.

### Plan-review disposition

[Fable's single plan review](medium_physics_fable_plan_2026-09-10.txt) found
Q-GOAL conditionally suitable and Q-LEAN acceptable. Adopted its substantive
condition: ship an Fe-specific erosion law with growth, not permanent
growth-only Fe. Acquired the Choban 2026 Table 3 Fe fit to Nozawa 2006;
see [pinned sources and assumptions](medium_fe_sources_2026-09-10.md).
The existing mass smoke was extended; no new test framework was added.
Both absolute donor momenta (Fe phases 5/6) and collective staged commit are
wired. Unresolved SN shocks reject with kinetics enabled. Existing global
growth/sputtering switches remain available for controlled process isolation.

Implementation adds `dust_fe_kinetics` (default false) and
`dust_fe_sticking` (default zero; explicit comparison probability when on).
The runtime uses total gas hydrogen nuclei, including molecular hydrogen,
for sputtering, while keeping molecular reservations for C/olivine growth.
Fe has no molecular carrier in the selected CHIMES network. Restart keeps
the old Fe attribute byte-for-byte when off and appends the exact kinetic
coefficients/domain/conventions when on. Both CLI and shared GUI are updated.

### Scope still remaining

This Fe feature is **one part of group 7**, not completion of all eight
approved medium-term groups. Same-population SED, microscopic Ia, stellar
mass/Z gaps, phase-resolved yields/winds/decay, alternative high-mass/PISN
inputs, collision-aware multibin comparison, remaining PAH/Fe high-energy
physics, and spectral/angular/cross-group RT are not marked complete here.
Do not turn the source constraints below into invented zero yields or new
mandatory infrastructure gates.

### Driver end evaluation: bounded Fe feature complete

Native implementation, generator wiring, physical-source binding and short
integration/restart are complete. This is a **bounded comparison PASS**, not
publication qualification or completion of group 7/all eight groups.

- Build: `.medium-fe-build.kl60Z7/ramses_fe_kinetics3d`, Intel MPI/ifx,
  NVAR338/NENER1, SNRT/DUST_LIVE/DUST_DYNAMICS/CHIMES/Fe/PAH/HDF5,
  CPU/OpenMP hydro. SHA256
  `e432875a5a89ebd376fdba3a5f498118363240fe5ef925abd90a9bfd91b44514`.
  A subsequent source edit clarifies a comment only; executable behavior is
  unchanged. VPATH order was not modified. The clean parallel attempt hit
  existing incomplete module-order dependencies; serial builds were used.
- Existing `dust_mass_smoke` passes Intel and GNU. GNU used bounds checks
  and invalid/zero/overflow traps. Fe cases cover polynomial units/domain,
  seed-only analytic growth, hot erosion, shared reservoir, timestep split,
  scalar/momentum parity, kinetic-energy heat, zero-rate identity, atomic
  failure, and incompatible unresolved SN shocks. No new smoke framework.
- Existing GUI/wizard suite: 47 tests run, 46 pass, one display test skipped.
  Default/off, enabled options, nonfinite/invalid probabilities and
  unsupported combinations are covered.
- Final effective input: `.medium-fe-live.S0OD4h/final/physical.nml`;
  MPI2/OMP2, noncosmological periodic 8^3 hydro, CHIMES, live IR,
  C/silicate/Fe/neutral PAH and relative dust motion, no star formation or
  AGN. `t_star=eps_star=0`; density=10 code units, initial gas pressure=1e-5,
  Fe sticking=.3, C/silicate sticking=0 to isolate Fe growth, thermal
  sputtering enabled, size collisions disabled. NML/environment retained.
- Two steps, final code time `3.43234212e-11`. The final Fe mass grows by
  `4.116040841495305e-12` relative to its initial value; this tiny signal
  verifies wiring, **not** a macroscopic growth calibration. Mean total Fe
  `0.03255813953488373` and total metal `0.10099999999999999` stay fixed;
  gas Fe is `0.0062790697674007375`. Dust aggregate/species closure residual
  is at most `3.469446951953614e-18` (code mass-density units).
- Final restart input: `.medium-fe-live.S0OD4h/final/restart/physical.nml`.
  From output1, the same final binary reaches output2. All **1063 datasets**
  are bytewise-equal numerically to the uninterrupted output; no floating
  dataset is nonfinite. The full Fe physics identity attribute also matches.
  Whole-file hashes differ because restart/header provenance differs.
- Final outputs were 60,045,368 bytes per HDF5 dump. Evaluation is complete;
  [cleanup](medium_fe_raw_cleanup_2026-09-10.md) removed 18 raw files
  (560,331,792 bytes), retaining inputs/logs/build identity and these results.

### Exposed existing limitation (not silently waived)

**Superseded below by the cold-continuation implementation and integration.**
The following paragraph records the original reproduction, not the current
status of the 5 K bracket.

Longer cold, dark PAH/Fe integration reaches below the supplied IR/material
temperature range (5 K minimum in this contract; the separate 10 K background
parameter is not the PAH absolute-emission lower bound). The analytic Fe material solve rejects its lower energy
bracket. It occurred with Fe kinetics both on and off, and with/without
relative motion. The initial Fe mixture was independently checked to be
20 K, not an invalid initial temperature. Shortening the fixture timestep
keeps two steps in the supported range; it does not repair long-time cold
evolution. The earlier runs and rejection logs remain in
`.medium-fe-run.wLeMv4`, `.medium-fe-off.4hsBmr`, and
`.medium-fe-live.S0OD4h/{coadvected,short}`. An exploratory interpolation
change did not fix the rejection and was reverted: Fe's actual IR callback
uses analytic mixture enthalpy. No radiation tolerance/domain was relaxed.

Current limitations remain neutral/unfocused grain rates, a fixed low-Z
projectile mixture, representative (not evolving) radii, no adsorption
latent heat/nonthermal erosion, and the original cold/electric-only optical
scope. Extending PAH/Fe low-temperature material coverage belongs with the
remaining group-7 physical work, not a newly invented generic gate.

### Group 7 continuation: absolute cold radiation (2026-09-10)

Implemented native bulk C/olivine/Fe enthalpy below 5 K, with zero-point
energy excluded. The DL01 low-temperature Debye integrals give a T^3 carbon
term and T^3+T^4 silicate terms; each material is joined to its existing 5 K
knot. Fe retains the existing HD17 electronic+Debye law and JANAF anchor.
The existing 5--3000 K C/olivine interpolation and non-PAH net-bath receiver
are unchanged. Sources: [DL01](https://arxiv.org/html/astro-ph/0011318),
[HD17](https://arxiv.org/abs/1611.08607).

The PAH shared absolute-field callback now supplies actual IR quadrature
photon energies to the native Fe/bulk receiver. Below its first knot it
solves U(T)+dt*P(T)-Q(T)=Uold+dt*H on [0,5 K], with each band's exact Planck
ratio, rather than flooring T or extrapolating a grey T^4 luminosity.
The same analytic material inverse is wired through mass exchange, phase
enthalpy advection and the next RT call, including the Fe-free PAH host.
No new namelist parameter: existing PAH selection enables this correction;
ordinary non-PAH paths retain their bath semantics. Fe material and PAH
restart identities are version 2; incompatible old identities reject.

Evidence (reuse of existing native smoke programs; no new framework):

- Intel native mass and radiation/backend suites pass. GNU mass suite
  passes with bounds checks and invalid/zero/overflow traps. New cases cover
  bulk U(0)=0, cold inverse, Fe and Fe-free absolute cooling, gas exchange,
  exact per-band Planck comparison and zero-energy equilibrium.
- Physical-opacity dark-cell cases cool from 4 K to
  0.10290361146, 0.10353417472, 0.10189679145 K (Fe, gas-coupled Fe,
  Fe-free host). This tests the *existing finite spectral quadrature*, not
  unrestricted centimetre-opacity accuracy or the presence of a CMB.
- Clean forced serial rebuild: `.medium-cold.LVkaiC/build-native.log`;
  `ramses_medium_cold3d`, SHA256
  `6670007bb799782bd36679d61336f10b38472886f5de35d9dc0f27a3cec0b430`.
  Intel MPI/ifx, NVAR338/NENER1, CPU/OpenMP, SNRT, CHIMES, live PAH/Fe and
  dust relative dynamics. VPATH order was not changed.
- Effective input `.medium-cold.LVkaiC/live.S2wuGA/physical.nml`, MPI2/OMP2,
  periodic noncosmological 8^3, no SF/AGN, same Fe growth isolation as above.
  Courant factor 5e-5 is the previously failing value, NOT the 1e-10
  workaround. Four steps finish in 74.67 s, time=3.43234212323913e-5 code,
  one million times the previous two-step cumulative time.
- Final independent low-T material inversion gives T=2.7518163--2.7525997 K
  across all 512 cells. Total Fe mean=0.03255813953488373 and total metal
  mean=0.10099999999999999 are unchanged. Fe grain growth=4.093978745e-6
  relative; dust aggregate closure <=6.94e-18 code density.
- `.medium-cold.LVkaiC/restart.QmAkTN/physical.nml` restarts output1 (step2)
  to step4 in 37.16 s. All 1063 datasets exactly match the uninterrupted
  output2; all floating datasets are finite. Fe, PAH and phase-dynamics
  identity attributes match. Whole-file hashes differ due to headers.
- Evaluated raw output cleanup: [exact files and hashes](medium_cold_raw_cleanup_2026-09-10.md),
  249,036,352 bytes removed; inputs, logs and binaries retained.
- Later source changes only replace literal 5 K with the identical named
  first knot, clarify an identity comment, and reject an ambiguous *unused*
  alternate bath argument in the extended absolute API. The final native
  smoke run is retained separately; these do not change this live call's
  inputs or numerical branch.
- Final-source executable is `.medium-cold.LVkaiC/ramses_medium_cold_final3d`;
  successful build retained in `build-final-native.log`. The full evolution
  and restart measurements above refer specifically to the hashed earlier
  executable, not an unperformed second evolution with this final build.

Driver evaluation: the reproduced 5 K failure is resolved in the bounded
live comparison, without timestep shrinkage or tolerance relaxation. This
is NOT completion of group 7 or of all eight groups. PAH charge-dependent
optics, H/C chemistry, destruction/photoelectron ledgers and hard-photon
coupling remain to implement; the existing neutral/fixed-H guard stays.

Source follow-up for that implementation: Montillaud et al. 2013 provides
charge/H-loss/addition prescriptions but explicitly does not model carbon
skeleton loss. Its approximations for missing molecular rates must remain
named approximations, not measured values. Malloci's database supplies
normal-coronene neutral/cation cross sections and vibrational modes;
those cannot simply be relabelled for every hydrogenation state. See
[reaction model](https://arxiv.org/html/1301.6507v1),
[coronene UV data](https://astrochemistry.oa-cagliari.inaf.it/database/coronene/coronene_sigma.html),
[coronene vibrational data](https://astrochemistry.oa-cagliari.inaf.it/database/coronene/coronene_vib.html).
All remaining numbered groups retain the operator's preapproval. No generic
infrastructure gates or user approval waits were introduced by this repair.

### Group 7 continuation: native PAH charge and gas-electron receiver

At the component milestone, implemented in the existing Fortran PAH module and C/Fortran CHIMES adapter;
no Python runtime, new test framework, extra approval gate or namelist option.
This is **local receiver implementation, not live charged-PAH completion**.
The live `pah_neutral_absolute_v1` model and its 4 eV guard remain unchanged.

- Neutral/cation event populations preserve molecule count and pre-existing
  vibrational energy during photoionization. Actual absorbed photon counts
  split into nonionizing vibrational captures, stored ionization energy and
  gas photoelectron heat. Exhausting the frozen neutral inventory rejects
  the entire trial, rather than clipping the photon budget.
- The C24H12 comparison uses Montillaud et al. (2013), eq. 1/table 4:
  IP=7.02 eV and alpha=1e-5 sqrt(300/T) cm3/s. Recombination is a model rate,
  not a measured coronene coefficient. The implemented admission is photons
  <=13.6 eV and gas T=10--10000 K; this is not universal survival validity.
- Electron/cation recombination counts use the analytic bimolecular solution,
  including electron depletion and the equal-abundance limit. Captures add
  IP plus kT to the neutral vibrational distribution. The kT electron energy
  is an explicit Maxwellian/sigma-v proportional to E^-1/2 assumption,
  not the bulk mean 3kT/2 or a measured photoelectron spectrum. Number and
  energy are preserved by interpolation. Above-grid events reject except
  roundoff-scale tails, whose molecules/electrons remain unreacted.
- `pah_charged_absorbed_step` composes photoionization, charge-resolved
  vibrational absorption/IR emission and recombination into one atomic local
  transaction. First-order splitting and frozen gas temperature are explicit;
  the enclosing gas solve must account for both heat and electron count.
  Inputs are already-debited, charge-partitioned captures. Stationary grains
  only; moving-grain work/Doppler physics is not silently inferred.
- `pah_charge_optics` reads the original neutral **or ionized** Draine table
  and rejects the wrong charge header. These are LD01 phenomenological PAH
  optical data, not molecule-specific Malloci spectra. Both charge states use
  the existing DL01 mode prescription: a named approximation, not a full
  reproduction of Montillaud's molecular model. H-state evolution, anions,
  dications and carbon-skeleton destruction are not provided by this receiver.

Sources: [Montillaud et al.](https://arxiv.org/html/1301.6507v1),
[Draine original optical tables](https://www.astro.princeton.edu/~draine/dust/dust.diel.html).
Original ionized data retained under `.medium-pah-charge.wBiAfM/`:
`PAHion_30.gz` SHA256
`0a7b7880ad5916daec1a76d14310933d09b46b79b30d84624f85c9322587b836`;
decompressed `PAHion_30` SHA256
`f0ef0128504d8d3d9e57d4b4443792bc601fcbcab5cf8bf9542cd522a9063d05`.
Direct source: `https://www.astro.princeton.edu/~draine/dust/diel/PAHion_30.gz`.

CHIMES previously reconstructs electrons from **gas-only** ionic charge.
The new signed-solid-charge entry points reconcile gas+solid neutrality and
retain it through the actual chemistry solve. The reproducible external
patch is now receiver ABI5 with one optional constraint callback; structures
and the 157-species network layout are unchanged. Its cell-specific context
is immutable during each chemistry substep and private to each OpenMP cell.
It is not itself a charging-rate RHS. Old ABI4 libraries retain the original
zero-charge path and explicitly reject nonzero solid charge (status 4).
Optional callback discovery uses dynamic lookup; the initial weak-data-symbol
approach failed an actual ABI4 load test and was replaced, not waived.

Built a separate library from pinned CHIMES
`a58e5c0311993b51abc63d84fff2958a0104f6d6` in
`.medium-pah-charge.wBiAfM/chimes`; original ABI4 library was not overwritten.
ABI5 library SHA256:
`aa9427c4f558f2a9179deff03dca88d89243ff94a8e639bb17a8ee19d32972c4`.
The source patch reverse-apply check matches that source tree.

Evidence, reusing existing native smoke programs:

- Intel mass suite and GNU mass suite (bounds checks and invalid/zero/overflow
  traps) pass, including original neutral+ion optical data. Coupled local
  relative energy residual is 3.74e-14 (Intel), 3.70e-14 (GNU). In the
  physical-optics fixture, gas heat=1.83721581393161e-23 erg/cm3,
  ionization-energy change=4.32682252431721e-23 erg/cm3, and IR emission over
  the 1 s step=1.93912161104615e-22 erg/cm3. These are one-cell checks,
  not an astrophysical calibration or a spatial charged-PAH evolution run.
- Native CHIMES suite passes with ABI5 and the preserved ABI4 library.
  Signed solid-charge states preserve total charge and nuclei; four different
  charge values reproduce serial solutions exactly under OpenMP. Zero charge
  exactly follows the original entry point. Real PAH-produced photoelectrons
  and their heat are accepted by the actual CHIMES thermal/chemical solver.
- Early and late failures preserve outputs; high photon energy, missing
  neutral inventory and resolved excitation-grid overflow reject.
- Logs: `.medium-pah-charge.wBiAfM/{build-mass.log,mass.log,`
  `build-chemistry.log,chemistry-abi5.log,chemistry-abi4.log,`
  `gnu/build-ordered.log,gnu/run.log}`. An initial GNU command's module-order
  failure was corrected; the ordered build and run are the successful ones.
- Full native hydro/CPU/OpenMP link also passes: NVAR338/NENER1,
  SNRT/DUST_LIVE/DUST_DYNAMICS/Fe/PAH/CHIMES/HDF5, existing VPATH unchanged.
  `.medium-pah-charge.wBiAfM/ramses_pah_charge_components3d`, SHA256
  `4239e2d35935851de3294106b0f40111816a540c376a82cfe1d6d20a8f289f9e`;
  log `build-native.log`. This binary was built, **not evolved** in a new
  RAMSES run. No raw simulation dumps were generated in this continuation.
  The original-table local test uses `SNRT_PAH_ION_TABLE` pointing to the
  retained original ionized file; this is not a new live runtime option.

Driver assessment: bounded native local receiver and gas charge interface
PASS. Full charged-PAH spatial evolution is not yet qualified or enabled.

Remaining integration within the already approved group 7: persistent
charge-resolved carriers and their H/C/energy accounting, injection/advection/
restart identity, charge-dependent shared opacity and actual primary capture
partition, coupled gas temperature/electron update, and consistent exclusion
of explicit PAHs from any duplicate empirical grain-recombination term.
H-loss/addition, destruction and broader charge/hard-photon physics remain
required before relaxing live survival limits. No production charge switch
has been exposed prematurely; mkrun/generator namelists therefore did not
change in this continuation. Do not describe this as completing group 7 or
any of the other seven medium-term groups.

Group-8 wiring inspection: the live primary moving-scatter call preserves
per-group N while carrying an independent energy moment, and the IR entry
is energy-only (its auxiliary scalar is not a physical photon count).
Reassigning an entire group by its mean E/N would hide photons on either
side of a spectral boundary and is not a valid completion of the planned
mixed-spectrum redistribution. No such shortcut or unused remapping helper
was added. A physical within-group spectral reconstruction and absorption
coefficient update remain actual implementation work, not a new audit gate.

### Group 7 continuation: fixed-H charged PAH live wiring

The former local-only receiver is now connected to hydro, primary RT, shared
IR, gas heat/electrons, CHIMES and HDF5 restart as the separately selectable
`pah_charge_fixed_h_v1` comparison. The unchanged neutral option still has
its 4 eV bound. This is **bounded runtime completion, not completion of group
7 or all eight groups**. H loss/addition, destruction, anions/dications and
general PDR/hard-source survival remain actual implementation work.

- `DUST_PAH_CHARGE=1` requires `DUST_PAH=1` and adds 128 cation carriers;
  hydro NVAR443/NENER1 with CHIMES and no Fe/relative motion. The 256 states
  participate in carrier normalization, H/C reservation, stellar removal,
  neutral source injection, coarse/fine/MPI transport and checkpoint payload.
  No VPATH reorder. Both namelist generators and GUI expose the same option.
- Cation mass carries excitation plus IP, avoiding duplicate stored energy.
  Primary/IR opacities use neutral and ionized original tables. Primary
  capture partitions retain the full-step initial charge mixture across IR
  substeps. Gas photoelectron heat, recombination energy and electron number
  are staged with the same atomic commit; electron changes update gas heat
  capacity. Gas chemistry and PAH charging remain explicitly first-order
  split. This is not a converged implicit unified chemical network.
- CHIMES requires ABI5 for the new option and includes solid charge in
  initialization, transport reconciliation, evolution and restart checks.
  Explicit PAHs are excluded from empirical CHIMES grain charging terms.
  Live closure revealed an m_H convention issue before final acceptance:
  dividing PAH mass by 300 is not the molecular count in CHIMES carrier
  units (1.66e-24 g versus the PAH proton-mass convention). The implementation
  now converts by `atomic_mh/dust_pah_molecule_g`, and the native smoke checks
  this conversion. The superseded `evolution/` run is not acceptance evidence.
- The model fixes C24H12, photons <=13.6 eV and gas 10--10000 K where PAHs
  exist; no cosmology, Fe or relative dust motion. Both charge states share
  the declared DL01 vibrational-mode approximation. Newly generated README
  and UI state the restrictions, rather than promising generic survival.

Driver end evaluation: **PASS for this bounded live comparison**.
Evidence directory `.pah-charge-live.9ReUG2/`:

- Final executable `ramses_charged_pah_fixed3d`, Intel MPI/ifx, CPU/OpenMP,
  SHA256 `567a8f88cd14f4a015e4ef04431fea3fd57f29c5dff151990ad914fa1a8cb26a`.
  Native build logs retained; ABI5 dependency is the pinned library above.
- Existing physical-optics backend tests: neutral relative energy residual
  5.018e-13; charged mixed/pure-PAH capture, gas heat, electron count, IR and
  hard-photon rollback pass with relative energy residual 6.312e-14. The
  unit-convention regression passes in `charged-units.log`. Existing mass
  and native thermochemistry suites also pass; no new test framework.
- GUI/wizard: 47 tests run, 46 pass, one display test skipped (`gui-final.log`).
- `fixed.9jxVyS/physical.nml`: periodic noncosmo 8^3, MPI2/OMP2, four steps,
  no stars/AGN, no incident radiation, active initially nonuniform PAH
  neutral/cation populations and gas velocity. It tests advection,
  recombination, gas electrons, IR and next-step CHIMES reconciliation.
  This run does **not** test a live FUV source; physical FUV capture partition
  is covered by the native mixed receiver test above.
- At steps 2/4, maximum gas+solid charge residual/electron count is
  7.213e-16 / 1.447e-15; gas+solid elemental residual is
  2.168e-16 / 4.337e-16. PAH mean mass remains 2.25e-4 code density;
  relative mass residual <=6.67e-16. Mean cation mass decreases from initial
  7.5e-5 to 2.472094492364652e-5 and 1.493921270908681e-5. All populations
  remain finite/nonnegative. Gas T is about 139--141 K in the saved states.
  Shared IR commit relative balance <=2.544e-16. This is not a claim that
  hydrodynamic E alone is conserved in a cooling run.
- `restart.TALRDn/physical.nml`: restart from step2 to step4; all 1329 hydro,
  3 RT, 12 gravity and 18 AMR datasets match continuous execution exactly,
  as do the PAH identity and step counter. Retained compact results:
  `physical-check.txt`, `restart-check.txt`. Wall times were 72.694 s
  continuous and 39.288 s restart. Launch output policies were audited
  before each run; raw dump cleanup is recorded separately.

This evaluation does not establish charged-PAH CUDA kernels, a mixed-refinement
charged run, stellar-source evolution with this option, timestep convergence
or general astrophysical survival. No new approval wait or per-helper audit
was introduced. Remaining original groups retain their table dispositions.

## H-state PAH continuation

The new optional `pah_hydrogen_m13_dl01_v1` comparison is implemented in
native Fortran and connected to the live hydro driver, not only an analysis
script. It carries 14 H states, two charges and 128 excitation states
(3584 masses; hydro NVAR3771). Existing neutral/fixed-H choices are retained.
The generator and GUI expose the new choice and its restrictive domain.

The H-loss/attachment prescription comes from
[Montillaud et al. 2013](https://arxiv.org/html/1301.6507v1), with generic
DL01 harmonic modes instead of molecule-specific spectra. H loss competes
with photon excitation and IR cooling in the same backward-Euler system;
attachment consumes finite gas HI, and photoelectron/recombination coupling
shares the gas electron ledger. Molecular binding energies and per-H mass
are included in the same transaction. Fully dehydrogenated C24 remains a
solid skeleton. Normal-H optics/cooling and constant IP across H states are
declared approximations. No carbon fragmentation, H2 chemistry, negative or
higher positive charges, or general hard-source survival is implied.

Evidence: `.pah-hydrogen.6fDwRw/` (Intel MPI/ifx and independent GNU bounds,
invalid/zero/overflow-trap component tests):

- `mass-final.log`, `gnu/run-final.log`: analytic competing loss/cooling,
  finite-H donor, binding/excitation/IR conservation, charge, overflow rollback
  and H0 skeleton retention pass. DOS grid refinement 1 -> 0.5 meV changes
  sampled 10--64 eV H12 rates by at most 0.231%; this is not full physical
  excitation-grid convergence. GNU bounds checking exposed a bond-array
  slice mismatch; corrected, then rerun successfully.
- Original physical neutral/ion optics tests: 12 eV absorption releases
  1.99852905352283e-11 H atoms/cm3 in the chosen native case, with energy
  residual at roundoff. Coupled H/charge attachment and recombination energy
  residual is 4.77e-17 in the final Intel test. Existing neutral/fixed-H
  backend tests pass in `neutral-final.log`/`charge-final.log`.
- `gui.log`: 47 tests, 46 pass and one display test skipped.
- Initial `live/` and equivalent `fast/` runs use noncosmo 8^3 periodic gas,
  MPI2/OMP2, four steps, nonuniform H12-neutral/H11-cation populations,
  advection, CHIMES and IR. No stars/AGN or incident FUV. H and charge evolve;
  local physical-optics tests, not this live fixture, cover FUV absorption.
- The general solver's steps 2/4 have maximum gas+solid nuclear residual
  4.337e-16 / 3.601e-16, charge/electron residual 4.171e-15 / 8.654e-16,
  and mean skeleton-count relative residual <=5.56e-16. Mean PAH mass changes
  from 2.25e-4 to 2.253693286917134e-4 / 2.253100766788167e-4 code density,
  as it must when H exchanges with gas. All populations are finite and
  nonnegative. The optimized zero-absorption triangular branch matches the
  general solver's step2 hydro state exactly.

The initial restart **failed**: the enlarged physical identity exceeded the
legacy HDF5 attribute-size limit. The prior writer printed HDF5 diagnostics
but allowed the run to finish without that attribute. Those initial dumps
are therefore NOT accepted restart evidence. The H-state identity now uses
a distributed dataset with immediate checked readback; old neutral/fixed-H
attribute formats stay unchanged. The reader requires the appropriate
representation, exact length and coefficients, with collective error handling.
Fresh integration/restart evaluation of this correction: **PASS for the
bounded H-state comparison**, not a complete PAH survival model.

- Final binary `ramses_pah_h_checked3d` SHA256
  `ab678292556155052c7e889d8a08deacedb63c7d4f963c508748a9dfc4e6d1e6`;
  `build-checked-final.log`. CPU/OpenMP, Intel MPI/ifx, pinned CHIMES ABI5.
- `checked.2XffQ5/physical.nml`, MPI2/OMP2, noncosmo 4^3, four steps;
  `replay.ITn60C/physical.nml`, steps2 -> 4. This smaller fixture isolates
  the checkpoint correction; it does not replace the preceding 8^3 physics
  comparison or establish spatial convergence. Effective output policy:
  noutput1, aout2, tout1e30, foutput2, fbackup1e6, two ~15 MiB dumps.
- All 7542 hydro, 2 RT, 8 gravity, 12 AMR datasets and the 11532-coefficient
  PAH identity dataset match continuous/restart execution exactly. No HDF5
  diagnostics, MG failures or nonzero NaN counters. IR relative balance
  <=1.223e-16. Nuclear residual <=4.802e-16, charge/electron residual
  <=2.589e-15 and mean skeleton-count residual <=3.34e-16. Mean H per
  skeleton changes from initial 11.665924 to 12.157814 / 12.078899; the
  initially H11 cations attach gas H, and excited grains can lose H again.
- `evaluation.txt` and the one-off `evaluate.py` retain the compact binary
  comparison. Wall times: 56.362 s continuous, 35.335 s restart.
- The optimized/general solvers' original 8^3 step2 and step4 states match
  exactly in all 11313 hydro, 3 RT, 12 gravity and 18 AMR datasets. Wall
  times 395.585/404.375 s are not a controlled throughput benchmark and
  do not imply a substantial runtime speedup.

No new default, tolerance relaxation, per-helper audit or approval wait.
This does not establish PAH-H CUDA kernels, stellar-source evolution,
mixed-refinement evolution or physical time/spatial/excitation convergence.
Raw cleanup is recorded in `pah_hydrogen_raw_cleanup_2026-09-10.md`.

### H2 vacancy-refilling follow-up (implemented)

`pah_h2_rehydrogenation_v1` extends the existing H-state comparison without
adding carriers, changing CHIMES ABI, or changing defaults. Cation H0--10
captures molecular H2 to refill two vacant sites at the M13 bound rate. It
is not a bound on total H2 effects, and does not include H2 formation,
single-vacancy abstraction or molecular superhydrogenation. Both generators
are updated. Native conservation/parity/refinement and one MPI2 live/restart
comparison pass; [scope and evidence](pah_h2_implementation_2026-09-10.md).
This does not resolve or replace the following carbon-fragmentation work.

### Active carbon-fragmentation follow-up (not implemented)

Source investigation continues; this is not a new approval gate or a claim
that public rates do not exist. [Lange, Dominik & Tielens 2025](https://doi.org/10.1051/0004-6361/202347722)
provides H/H2/C2H2 competing channels and an X-ray energy-partition model.
Its pyrene-based fragmentation extrapolation is explicitly uncertain and
can overestimate carbon-backbone breakdown. Its cascade stops at backbone
damage; it does not supply the subsequent daughter population chemistry
needed for a conservative live galaxy calculation. Do not turn its damage
probability into immediate atomization of the whole C24 molecule, and do
not deposit the entire X-ray energy as PAH vibration.

[Murga et al. 2020](https://arxiv.org/html/2007.06568v1) tracks changing H/C
and distinguishes C2H2 release (enough H present) from C2 release after
dehydrogenation. [Murga et al. 2022](https://doi.org/10.1093/mnras/stab3061)
extends size/H/charge evolution and compares carbon-loss assumptions. These
are useful physical prescriptions, but copying rates alone does not connect
their fragments to this code's gas chemistry, radiation and energy budget.
The current CHIMES157 receiver includes C2 but not C2H2. Required next
implementation is a consistent daughter/fragment inventory and its energy,
charge and gas-chemistry coupling, alongside the selected fragmentation law.
These tasks remain in approved group7; they are not silently deferred to
long-term work or marked finished by the H-state tests.

## Completion discipline

Use the operator-approved selective Fable policy in AGENTS.md: no automatic
per-bundle review; call for substantive new physics/coupling or an unresolved
scientific design question, asking Q-GOAL then Q-LEAN. Focused existing native
tests and driver end evaluation remain. No new test framework, dense
synthetic atlas, or per-helper review. A bounded comparison is not a new
default or a complete physical model. Report completed runtime capability
separately from data-blocked and not-yet-implemented items. Remove evaluated
raw simulation dumps per AGENTS.md, retaining compact evidence and inputs.

## Initial source constraints

- LC18 phase rows do not supply all eleven surface elements at every phase.
- Huscher AGB tables omit S/Ca/Fe and full stellar lifetimes; they cannot
  silently replace the complete eleven-element stellar source.
- Existing COSMIC histories are unweighted individual binaries and do not
  establish an SSP Ia DTD or the instantaneous explosion ejecta.
- BPASS's staged IMF/mass bounds and binary histories do not match the
  currently selected feedback population.
- Current Fe material/optics treatment is a cold, electric-only comparison;
  charged grains, magnetic absorption and unrestricted hard photons are not
  implied by adding a mass exchange operator.
