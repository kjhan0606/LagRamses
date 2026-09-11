# Single-photon complete-atomization yield comparison: proposed plan

Status: native helper and mixed receiver implemented and tested; NOT a
declaration that the production channel/integrated run is complete.
User approval: "근사모형, 명시적
비교로 진행". The driver's design choice excludes ensemble pooling.
Checkout verified `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`;
shared context and AGENTS read. The subsequent user directive authorizes
`dust_pah_atomization.f90`, the atomization smoke section, and the actual
`dust_pah_mixed.f90` kernel integration. Main owns the outer receiver, source,
selector, live-model identity and IR population-loss admission. No RAMSES/MPI
jobs or commits. The one selective Fable review is recorded in
`pah_atomization_fable_2026-09-11.md`; it was not repeated.

## Model and scope

Keep the existing C24 population `(128,0:13,0:1)`; remove a molecule only
when ONE absorbed photon plus THAT molecule's own excitation can pay for
complete atomization. Yield is one above the energy threshold, zero below.
No energy sharing between molecules, budget allocator, new carbon-size
axis, C22 optical surrogate, acetylene/vinylidene substitution, or thermal
atomization inferred from a cell-averaged energy reservoir.

This is an energy-allowed unit-yield comparison, not measured kinetics or
a prediction of real X-ray branching. It assumes complete retention of the
absorbed photon's available energy, no escaping photoelectron or fluorescence,
and the charge-preserving atomic endpoint below. It bounds the yield of
this specified single-photon endpoint, not all cheaper carbon-fragment loss.

## Reference energy: sourced anchor, explicit thermal-model conversion

Use neutral ground-state H and C atoms as zero chemical energy. NIST gives
gas-phase C24H12 formation enthalpy at 298 K `295 +/- 11 kJ/mol`:
[NIST coronene / Roux et al. 2008](https://webbook.nist.gov/cgi/cbook.cgi?ID=C191071&Mask=2F&Units=SI).
Do not relabel that value as a measured 0 K atomization energy.

With the existing 102 DL01 modes in `dust_stochastic.f90:dust_pah_modes`,
the native harmonic formula gives `Uvib(298.15 K)=0.3067980323 eV`.
The nonlinear ideal-gas thermal enthalpy increment is `Uvib+4 kT`, or
`39.51733816 kJ/mol`. The elemental increments are graphite `1.050` and
H2 `8.468 kJ/mol` from [NIST/CODATA reference states](https://cccbdb.nist.gov/refstate.asp).
Thus, using kJ/mol throughout the formation calculation:

```text
Hf0(PAH) = 295 - 39.51733816 + 24*1.050 + 6*8.468
         = 331.49066184
A12 = [24*711.396 + 12*216.034 - Hf0(PAH)] / 96.48533212
    = 200.38715640 eV per C24H12 molecule
```

The atomic values are already pinned in `snrt_chimes_atomization.h` and
`simulation/snrt/data/chimes_atomization_atct1130.json`. The input enthalpy
uncertainty contributes about 0.114 eV; it does not include the generic-mode
thermal correction or the existing H-state approximation. Keep these model
choices in the restart identity. Keep the existing PAH IP, 7.02 eV, consistent
throughout; this plan does not substitute a different experimental IP.

## Individual event and exact material ledger

For state s=(j,h,q), define (all quantities here in eV):

```text
E_PAH(s) = u(j) + B(h) + q*I_PAH - A12
I_C      = (1797.849 - 711.396)/96.48533212 = 11.26029186
D(h,q)   = A12 - B(h) + q*(I_C-I_PAH)
Y(s,eps) = 1 if eps + u(j) >= D(h,q), otherwise 0
```

`B(12)=0`, `B(0)=48`, `B(13)=-3.2` are the existing H model's offsets.
`eps` is the individual absorbed photon energy available in the material
frame, not gas heat, a cell energy sum, or a nominal reference energy that
ignores transport energy correction/mechanical work.

For every destructive event:

```text
C24Hh  + photon -> 24 C      + h H
C24Hh+ + photon -> 23 C + C+ + h H
Qexcess = eps + u(j) - D(h,q) >= 0
```

Destroy one parent, consume one photon, add `(24-q) C`, `q C+`, `h H`,
and change no electron or H2 count. Transfer Qexcess to gas thermal energy
with an explicit immediate-local-accommodation assumption for product kinetic
energy; do not both retain fragment KE and add the same energy as heat.
For a cold neutral H12 parent and a 250 eV photon, Qexcess=49.61284360 eV;
for its cation, Qexcess=45.37255174 eV.

This obeys `eps = delta(E_PAH + E_gas,chemical) + Qexcess` event by event.
The 64 eV vibration grid need not be extended: a destructive event exits to
gas directly, never through a fictitious >64 eV PAH bin. Subthreshold photons
do not destroy a grain via this channel, however large their aggregate energy.
Normal earlier excitation of the SAME surviving grain remains represented by u.

## Native event integration contract for main

The charged `dust_pah_mixed.f90:material` currently derives capture counts
from energy divided by a fixed group energy; `pah_coronene_photoionize`
rejects occupied photon energies above 13.6 eV. `pah_live_prepare` also masks
hard primary opacity. Therefore this cannot be activated by changing only
the 64 eV grid or removing the 13.6 eV check.

Main's absorption interface must provide accepted PAH events and their
actual energies, preserving charge and the chosen spectral resolution.
An energy/count ratio is an event energy only under an explicitly declared
monochromatic group/ray closure. Two moments alone cannot determine a hard
tail across a nonlinear threshold; never merge soft/hard streams first.
Hard-band absorption coefficients must be supported by the selected optical
source/domain; do not invent a continuation or confuse scattering with absorption.

Route destructive captures BEFORE the old photoionization/vibrational call.
Remove both their number and full energy from that call; its other captures
continue normally. Derive rates from per-grain absorption, not from an
independent guessed Arrhenius rate. In a frozen-state substep with the current
H-independent cross sections, distribute charge-resolved captures by
`a(g,s)=captures(g,q)*p(s)/sum_q(p)`, then apply Y separately to each state.
This is event weighting, not energy pooling. Accumulate all groups before
checking `sum_g(a*Y) <= p(s)`; never sort grains by cheapest energy cost.
Recompute depletion/opacity with the existing trial/substep machinery. A
frozen-opacity trial that spends photons on exhausted parents must be retried
with smaller transport/material steps, not silently clipped after absorption.
The split requires timestep convergence, not a claim of exact finite-step
competing-event statistics. Each accepted event is debited exactly once.

Important bounded admission: the old receiver cannot process a non-destructive
hard event merely because this new channel returned Y=0. Initially admit only
supported soft events plus hard events guaranteed to atomize their addressed
states. A cold-state sufficient hard threshold for ALL current H/charge states
is about 207.82744826 eV. Intermediate occupied hard bands remain unsupported
until main supplies their conservative non-destructive processing. Do not
silently make those bands transparent, clamp excitation, or thermalize them
through an invented fallback. This is the specific integration decision still
needed for a broad hard spectrum, not another request for a source-gap review.

Native pure helper API (module `dust_pah_atomization_physics`, source file
`dust_pah_atomization.f90`; implemented, no persistent event tensor):

```fortran
pah_atomization_bond(hydrogen) ! erg, existing H-state offset
pah_atomization_threshold(excitation, bond, charge) ! erg, minimum photon
call pah_atomization_batch(excitation, bond, hydrogen, charge, &
    photon_energy, captured, old, next, destroyed, receipt, ierr)
```

State vectors have shape `(ns)`; photon energy, captures and destroyed counts
have shape `(ng,ns)`. Hydrogen and charge are integers. Energies use erg/event;
counts use cm^-3. `next`, `destroyed`, and `receipt` are unchanged on failure,
overwritten on success. Status 0=success, 1=invalid/nonfinite, 2=parent
exhaustion; no clipping or rescaling occurs. The helper has no cell energy
budget, rate fit or optical-coverage override.

`type(pah_atomization_receipt)` exposes `hydrogen`, `carbon` (neutral),
`carbon_ion`, `photon_number`, `photon_energy`, `gas_heat`,
`pah_level_removed`, `binding_increase`, and `gas_ionization_increase`.
Its energy identity is `photon_energy = gas_heat - pah_level_removed +
binding_increase + gas_ionization_increase`. Public erg constants are
`pah_atomization_a12`, `pah_atomization_ip`, `pah_atomization_carbon_ip`.

The retained neutral `.dust-pah-native.LWUJVc/PAHneu_30.dat` and ionized
`.medium-pah-charge.wBiAfM/PAHion_30` tables both have 30x1201 data rows,
lambda=0.001--1000 micron, or 0.00123984--1239.841984 eV. In the actual
nine-group monochromatic closure, ONLY group 8 at 869.634149 eV is a supported
hard event. Its full 500--2000 eV parent interval is not entirely covered:
this comparison is monochromatic, not an asserted broadband hard tail.
250 eV examples above are helper algebra tests, not the live group selection.
All other occupied unsupported groups remain rejected.

## Derived binding and source convention: no fake heat

Derive `Ebind_PAH=-A12*N_PAH` from populations, with no new advected scalar.
Leave the stochastic levels nonnegative relative to their existing references.
In `dust_pah_mixed`'s initial-reference material residual add
`A12*(Nold-Nnext) + I_C*delta_Cplus`; vibration, H bond and PAH IP changes
are already in the population energy. H2's existing chemical correction remains.
No second debit of the full D is allowed after this residual accounts for it.

`snrt_chimes_atomize` supplies the atomic/ionic energy convention, not the
reaction operator for PAHs: do not invoke it on unrelated gas molecules.
Insert the explicit C/C+/H products in CHIMES, preserve charge, and recompute
particle heat capacity. Prevent subsequent dust-depletion reconciliation from
redistributing these products over arbitrary carbon ions/molecules. Transfer
the removed `(288+h)*mH` phase mass and donor momentum; mechanical mixing heat
is separate. Keep graphite latent energy and PAH binding in one reference,
without charging graphite latent heat again for PAH-to-gas carbon.

Selected by main: atomic gross-ejecta donors with locally retained formation
energy. The source gas correction is `+A12*Ninject-Uvib_inject`, only for
actual opt-in source condensation. The increase in gas heat is balanced by
the newly injected negative PAH binding energy. This supersedes the original
precondensed-source recommendation and the corresponding Fable receipt
condition: NO formation-escape ledger is also applied. The precondensed
convention would be a different source model, not an additional correction.
Use the same reference/receipt in both comparison branches and on initialization;
introducing the reference itself must not heat an existing cell.

## Native verification and pending integration

Private build: `.pah-atomization-native.MFwJp8` (not any shared PAH build
directory). GNU 13.2, `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`.
Independent optimized private build: `.pah-atomization-opt.oZILRY`, `-O2`
with the same bounds/FPE checks. Both full smoke executions exited 0 after
main's module-only rename to `dust_pah_atomization_physics`; no namespace
collision remains in these private builds. Logs are `mixed_checked.log`
and `mixed_optimized.log` in the respective directories. The helper and
five mixed-case output values agree at printed precision across builds.
The complete existing `dust_mass_smoke` passed with the new helper tests:
native-mode reconstruction of A12, both charges/H0/H12/H13 at adjacent-FP64
threshold energies, individual hot-grain assistance, no soft-event pooling,
one hard event/one parent, nuclei/charge/material receipts, atomic-source
formation/destruction closure, nonfinite and overflow rejection, exact
exhaustion and bitwise failure rollback.

The actual mixed test is opt-in with `SNRT_PAH_ATOMIZATION_SMOKE=1` plus
the two retained optical-table environment variables. It exercises group 8,
partial and complete loss, both charge endpoints, shared native IR balance,
unsupported-group rejection and a late-cell failure after local receipts
were staged. The full checked smoke run, including existing original-optics
H/H2/catalytic tests, passes. `mixed_checked.log` records energy residuals
relative to injected primary energy: partial mixed-charge loss 2.97e-16;
complete neutral and cation loss 0 at printed precision; soft-only 1.49e-13;
soft+hard 1.48e-16. The soft-only/soft+hard cases check that residual captures
still reach the legacy excitation receiver without a duplicate hard debit.

The test uses a finite 1 s step. A preliminary 1e-12 s fixture formed a
minute H13 cation population via reattachment and was rejected inside the
existing `pah_hydrogen_absorbed_step`; its checks were NOT relaxed or edited.
Total H (including bound H and H2), rather than frozen gas H, is checked
because surviving PAHs can reattach the atomization products.

The fixture sets in-memory contract axes from the compiled D03 data and
loads the retained original PAH tables. This qualifies the actual native
mixed/IR kernel; it is not an outer CHIMES/source/RAMSES integration run.
Neither MPI nor RAMSES was launched. Existing dirty changes outside the
owned atomization sections were preserved. No Makefile changes were made
by the helper/mixed owner; main owns the expanded smoke link dependencies.

The mixed public tail is now `...,gas_atomic_h,gas_molecular_h2,gas_atomic_c,
gas_carbon_ion`, all gas vectors `(ncell)` in cm^-3. The C/C+ arrays are
required only by the new mode. Captures are state/charge partitioned before
the hard events are removed; soft captures continue through the existing
H2 catalytic step. C/H particle capacity, local excess heat and the A12/IC
material residual use the same scratch transaction. The IR call permits
population loss only for this named, stationary, no-Fe endpoint.

Driver implementation choice (still pending integration): use the alternative
**atomic donor with locally retained formation energy** above, explicitly
for this opt-in model. The source already provides atomic gross elemental
ejecta and applies a phenomenological PAH condensation fraction. The signed
correction `+A12*Ninject-Uvib_inject`, together with the derived negative
PAH binding energy, closes that atomic-source energy budget without a new
escaped-energy carrier or an unconsumed diagnostic receipt. It is not
measured condensation radiative transfer. Existing cells receive no heat
merely from selecting the reference; only actual new condensation releases
binding energy. Default source handling is unchanged. This supersedes the
precondensed-ejecta recommendation for implementation, not its validity as
an alternative convention. Do not apply both conventions.

Native tests: threshold below/equal/above; neutral/cation and H0/H12/H13;
hot-state individual threshold; many subthreshold photons give zero direct
atomization; one hard photon cannot destroy two grains; gas H/C/charge and
full energy closure; no duplicate heating/photoionization; exhaustion retry
and failure rollback; source-reference closure; unsupported-band rejection;
old-option numerical regression. Main owns integration and the eventual
integrated dust run. Existing MHD work is outside this comparison's scope.

## Integrated Intel build and native regression

Full CPU/HDF5/CHIMES build with NENER1/NVAR3773 succeeded in
`.remaining-dust.d6oJ1d`; frozen executable SHA256
`b3c58d44bbba17cd84ecb41ad0390854ac255f8073c17eb7588a7a0c59591fc9`.
Both existing mass/backend smokes and the atomization opt-in passed; logs
`build-smoke-rounding-final.log`, `dust-mass-rounding-final-smoke.log`,
`dust-backend-rounding-final-smoke.log` retain the results.

The initial mass-smoke link lacked its newly consumed optical/mixed objects;
only target dependencies/link inputs were repaired, without changing VPATH.
Intel revealed a1-ULP independently evaluated H4 bond difference: the test
now measures ULPs (maximum1, allowed4), not exact-bit algebra equivalence.
An exact-zero invalid-ground fixture now uses exported A12, not an
independently rounded reconstruction. Photon threshold, exhaustion and
overflow acceptance were NOT relaxed.

Actual source injection initially divided A12 by the legacy fixed-H molecule
mass, under-crediting formation heat by~0.755%. It now uses the exact same
H12 state mass as the injected excitation and population. Native actual
injection-mass/binding cycle closes with zero measured energy residual.
Five mixed atomization cases have maximum relative energy residual2.182e-13;
source-default parity and late-cell rollback remain tested. These are native
implementation results; actual MPI/CHIMES evolution is evaluated separately.

Live integration corrections (preserved initial/retry01 failures in
`.pah-atomization-live.Bs73kx`): the material interface initially compared
the supplied gas heat-capacity coefficient against exact SI kB, rejecting
the pinned CHIMES rounded BOLTZMANNCGS. The receiver now consistently uses
the caller's already positive/finite EOS coefficient. The native actual
CHIMES-coefficient case closes energy to2.97e-16 relative. No thermal
acceptance tolerance was increased.

Next, a first-step failure was independently reproduced in both catalytic
and atomization selectors by `.pah-live-native-diagnose.8b02hQ`, at about
153K with tiny IR captures. H13 abstraction's subtractive survivor update
produced a negative subnormal population (-8.9199e-319) under the optimized
no-FTZ build. The equivalent update `(1-f)*old_H13`, with explicitly finite
0<=f<=1, preserves nonnegativity without clipping. Both native reproductions
pass after this change. Reaction rates, donor/event counts, chemical-energy
receipts and existing conservation tolerances are unchanged. The corrected
actual live/restart execution remains a separate required comparison.

## Final bounded MPI live/restart result (2026-09-11)

That comparison is now complete in
`/gpfs/kjhan/LRD_JWST/.pah-atomization-live.Bs73kx`.
Final binary SHA256:
`c923fd98f3c9c4adccdbbfbc6d8d15fcb8c7b48837f982535509f367d908e305`;
released from `.pah-h13-survival-binary.3ZICPB/ramses_dust_comparison3d`.
Both controls use MPI2/OMP2,64 cells,4 steps and a successful step2 checkpoint
restart to step4; NVAR3773/CHIMES1/DUST_IRON1/DUST_PAHH1, runtime Fe/drift off.
Initial neutral/cation H13 and H2 are nonzero. All grain condensation is zero.
Native SF forms256 stars with initial mass1.995641472410421e-6 code. Ordinary
material feedback is the declared synthetic fixture, not a physical stellar
population claim. Radiation is only fixed group8 at869.6341490248457eV.

The `retry02` rate1e36 photons/s/initial Msun is a weak, non-discriminating
activity control, NOT a radiation/atomization PASS. Its execution and
conservation/restart checks pass. The otherwise identical `flux42` rate1e42
control has measurable excess PAH destruction. Both have zero FINAL stored
primary photons: this is not evidence of absent source or underflow. Actual
native source activation, group8-only SED binding, and the differential
destruction/carbon receipts establish positive radiation-dependent activity.

At step2, strong-minus-weak destruction is2.861607536273037e-13 code PAH-number
density,3.828194970747485e-7 of the initial population. Neutral C and C+ gains
are6.820325895360854e-12 and4.753219006780109e-14; their sum equals24 times
destroyed PAH number to2.3685e-10 relative. At step4, differential destruction
is1.7181553638102171e-12,2.2985100644750014e-6 of initial PAH number. ALL gas
carbon nuclei gain4.1235728734110096e-11, matching24 times destruction to
6.4626e-11 relative. Subsequent CHIMES chemistry moves some returned C/C+
into CO/CH species; a C/C+-only cumulative closure would be incorrect.
The initial analysis attempt making that assertion is retained separately
as `flux42/evaluation_endpoint_attempt.txt`; final evaluation includes every
gas carbon-bearing species and retains the direct endpoint check at step2.
Stellar birth/current mass, birth epoch, type and identity match the weak
control exactly, as do gas density, total elemental carbon and inert silicate.
This separates the positive response from SF/transport; no skeleton-constant
assertion or fabricated thermal reservoir is used.

Maximum strong-run relative residuals: nuclei8.674e-16, charge4.248e-15,
accepted live IR energy6.2422e-16 (restart4.0572e-16). Gas plus stellar mass
is10 code exactly. IR energy is the accepted local coupled balance, not
global energy constancy in a cooling/feedback run. Material identity has
11587 entries. All7591 checked restart datasets match exactly: hydro7546,
snrt3,gravity8,amr12,header1,particles17,coarse3,domain1; step/time/NVAR
attributes also match. Strong live/restart wall times106.685/77.006s;
weak live/restart118.575/86.465s. No further experiments or source changes
were needed after this bounded closure evaluation.

Compact evidence: package `EVIDENCE.md`, each control's `evaluate.py` and
`evaluation.txt`, effective NML/env/SED, binaries, logs and output text
metadata. Authorized cleanup completed for BOTH evaluated controls:
16 HDF5/cooling files,116816128bytes (~111.4MiB), removed after paired
evaluation and restart. Eight `metadata_00001/metadata_00002` directories
retain the output text records. `cleanup_manifest.json` records exact raw
paths/sizes/SHA256; `cleanup.log` confirms completion. Raw files are not
retained; reproducible inputs/binaries remain. Original failed attempts are
untouched. This is a named mono-group comparison, not broadband PAH survival,
an actual AGN/physical stellar SED, or a production-science qualification.
