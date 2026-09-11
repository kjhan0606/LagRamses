# Approved remaining physics bundles: execution record

User preapproves bundles1--5 without routine intermediate approval pauses.
Workspace `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
This is the existing scope, not five additional gates. Driver evaluates
completed work; one substantive design review may cover a physical extension.

## Current result: approved named comparisons implemented and evaluated

The explicit-approximation scope is complete: PAH catalytic H2 and
single-photon atomization, Fe UV/thermal limits, mixed/truncated low-mass
sources, frozen-contact effective SNIa, wind bi-stability velocity and
transparent Al26/Fe60 material. Defaults are unchanged. Live/restart evidence
is recorded below where applicable; wind has actual-grid/native source/SED
evidence and Fe has bounded live (not restart) evidence.

This does NOT complete a daughter-PAH/higher-charge network, keV metallic-Fe
electron cascades, one matched full-mass stellar population, microscopic
binary ownership, LBV/full rotation/pulse chronology or radioactive heating.
Those are documented physical exclusions, not additional gates added to this
approved comparison implementation. Historical source assessments below
explain why the explicitly selected approximations were used.

## Bundle1: kind7 dust processes, completed bounded connection

Final status: [driver evaluation PASS](kind7_dust_process_implementation_2026-09-11.md).
The progression below retains the encountered failures, not current blockers.

See [approved design and Fable disposition](kind7_dust_process_bundle_plan_2026-09-11.md).
Native build and264 chemistry assertions pass. Existing native dust mass,
sublimation, shock, phase transport/drag and rollback tests pass. Frontends:
49 passes, one display-dependent skip. These native tests alone do not
constitute a live pass; the final live/restart evaluation is linked above.
Private evidence `.kind7-dust.PakmEZ/`; final raw cleanup is recorded below.

Initial fixture startup inherited an incompatible AGN reference selector;
corrected to `legacy` with AGN inactive. Its warm2500K grain IC then exceeded
the inherited300K material grid, correctly rejecting in the coadvected mass
step and relative enthalpy face pack. Preserve both failed logs. The current
`coadvected-hot`/`relative-hot` inputs use the existing3000K material contract
from the completed sublimation work, SHA256
`8bc8c51d07118b2cd63121e4e379f8c73ba01c836f2cefe02a1d975e6fdf706c`.

The coadvected hot profile completed two steps, with actual stellar radiation
in step2. Global gas+stellar mass residual is zero at printed FP64 precision;
maximum element/charge relative residuals8.471e-16/4.270e-16; maximum IR
balance1.662e-12.64 stars formed in step1 and returned3.6367041692e-9 code mass
by step2. Compact evaluation is `coadvected-evaluation.json` in that directory.
Relative hot integration exposed a separate cold-domain connection: its
absolute C/silicate IR callback omitted the already implemented Planck
low-temperature branch and still inherited net-bath admission. This was
repaired consistently in material, phase enthalpy and split-sublimation paths,
without adding a bath heat source or relaxing energy acceptance.

The exact native source re-evaluation (`source_check.f90`/`.log`) uses the
actual Chabrier IMF (`imf_id=2`), not the initially assumed Kroupa selection.
At Z=.02 and age0--5.148513184858701Myr, wind returns0.1212097691Msun and
5.851039824e-4Msun of condensed grains per initialMsun. Ordinary SNII is
still zero; the pair channel (PPISN at this node) supplies1.136851765e44erg
per initialMsun. That channel enters the existing SNII+pair shock-energy
sum, so no extra SNII-only live gate is required. This native calculation
verifies actual source input, not independently measured live destruction.

Relative cold continuation reached step1, then exposed the strict one-way
sublimation momentum check at an unchanged mass. Merely expressing cgs/code
conversion as a survival ratio is not enough under the production optimizer:
an isolated ifx-O3 test of40000 no-erosion conversions gives1250 apparent mass
increases, versus zero identity failures with an explicit equality branch.
The current code preserves exact unchanged mass before rescaling nonzero
erosion and retains strict no-growth phase transfer. Diagnostic messages now
identify phase, atomic-budget, chemical-reconciliation and energy failures.
The latest live pair is `relative-identity`/`relative-identity-resume`, binary
SHA256 `b3395e1d5585f2bf1de887902d30e2c72295144b874104754ee0d695ec4398c3`.
Both completed and pass conservation;442 datasets and physics/clock
attributes match bitwise at both steps. All9 evaluated raw directories
(about566MiB) were removed; `cleanup.json` and text metadata retain the record.

## Bundles2--4: pre-approval source assessment (historical)

2. PAH carbon loss lacks a complete daughter thermochemistry/optics network:
   C2H2 is absent from CHIMES157 and C22->C20 leaves the current optical table.
   Primary-paper product energies are being checked, rather than treating
   activation barriers as enthalpies. Fe hard-photon absorption still lacks
   a selected metallic-Fe electron/charge/energy partition. Existing FUV PAH
   charge/H and cold Fe mechanics must not be advertised as hard-source
   survival. [Prior source assessment](pah_fragmentation_sources_2026-09-10.md).
   Source check correction: [Hensley & Draine2017](https://arxiv.org/abs/1611.08607)
   does provide metallic-Fe UV yields, size/charge-dependent thresholds,
   photodetachment and collisional charging with a demonstrated13.6eV cutoff.
   UV charge/heat wiring is therefore implementable, not an absence of all
   metallic-Fe data. Its discussion does not provide a complete keV Auger/
   secondary/fluorescence cascade; extending that UV model to all nine photon
   groups would still invent the missing hard-photon partition.
   Native integration review: a finite-time nearest-neighbour charge solver
   can use event counts to update gas electrons, grain charge energy and
   sensible heat consistently. An equilibrium charge reset cannot silently
   supply/remove electrons or energy. For the existing10/100nm grains, the
   HD2017 electron-autoionization and13.6eV photon bounds imply174/9237
   charge states, about75kB per cell if every population is a passive FP64
   carrier. This is not a production layout decision. Ion collisions can
   cross the photon-only upper bound, so it is not a reflecting ion boundary.
   A compact charge representation and neutralized-ion desorption/energy
   accommodation need a declared physical closure. Negative grains may
   emit ionizing electrons even under a UV cutoff; all-heat deposition is
   not universally valid. Existing Fe Qabs alone does not determine the
   attenuation length needed by the yield enhancement: use the same
   source dielectric model, not an inferred unique refractive index.
   No disconnected charge helper or additional hydro state has been added.
3. Existing PARSEC archives include2--12Msun tracks and Q data, but published
   integrated ejecta begin at14Msun. The lower tracks end before complete
   TP-AGB; below2Msun and exactZ=0 are absent. Pre-terminal coverage is possible
   only with an explicit track-end validity condition, not a fabricated death
   or remnant. Existing KL16/Fishlock/Sukhbold alternatives are different
   populations, not a same-evolution completion.
4. Existing COSMIC donor/WD genealogy is useful but its0.15Msun He-triggered
   CO-WD destruction is not a justified N100 event. Saved histories lack
   incident transfer, coupled retention/outflow/ignition and channel-matched
   ejecta. Enlarging that grid or normalizing its events cannot repair this.
   Preserve the already approved empirical DTD/N100 path unchanged.

Operator answered explicitly: proceed with approximation models as named
comparisons (2026-09-11). This resolves the model-choice pause; bundles2--5
remain preapproved without routine intermediate approval. Preserve defaults,
select the mixed/approximate models explicitly, bind their identities and
document domains and omitted physics. Approximation approval does not turn
an activation barrier into a product enthalpy, or a mixed stellar grid into
one self-consistent evolutionary population. Source-supported work and
declared closures now proceed together.

## Bundle5: source-supported implementation rationale (historical)

LC18/NUBASE long-lived inventories are present: this is missing native
source/transport/decay wiring, not missing nuclear data. Revisit the
[existing radioactive material design](stellar_radioactive_inventory_plan_2026-09-10.md)
for the now explicitly requested long-lived evolution, rather than repeating
the already completed prompt Ni/Co projection. Source-time convolution and
gas elemental changes are required; another diagnostic isotope helper is not
completion. Dust-lattice transmutation and radioactive energy deposition
remain distinct physical couplings, not implicit behavior of an abundance
option.

[Selective Fable review](remaining_radioactive_fable_2026-09-11.txt) finds
the abundance-only model scientifically admissible but questions its value
with isotope-sensitive dust/MeV consumers excluded. It recommends reducing
any live continuation to Al26/Fe60 on the existing prompt-projected LC18
source, not undoing the prompt Ni/Co correction. Its few-node abundance
ratios are review estimates, not a demonstrated universal sensitivity bound.
The driver does not equate small abundance changes with mathematically zero
effect, or treat this recommendation as cancellation of the user's bundle5.
The reduced two-parent live isotope implementation is now in progress;
completion requires source convolution, gas transport/decay and integration
evidence, not merely the presence of its modules.

## Approved comparison execution

- [PAH catalytic H2](pah_catalytic_comparison_implementation_2026-09-11.md):
  native physics and MPI2/OMP2 restart evaluation complete; four evaluated
  raw output directories removed, compact evidence retained.
- [Four-plan selective review](remaining_comparison_fable_2026-09-11.md):
  SNIa and wind PASS, Fe and mixed low-mass CONDITIONAL. Driver disposition
  records the necessary accounting clarifications and rejects a universal
  sub-9-eV electron bound inferred by the reviewer.
- [Mixed low-mass connection](parsec_mixed_lowmass_evidence_2026-09-11.md):
  actual979-node source evaluation plus MPI2/OMP2 live/restart PASS. Raw
  cleanup complete; mixed/truncated source limitations remain explicit.
- Frozen-contact effective SNIa: MPI2 live event and restart PASS in
  `.snia-live.hUBFAz`,83 final datasets exact. Actual Ia-only interval mass,
  Fe and energy receipts agree with the selected N100 event source; the
  effective-SSP approximation and exactZ=.01 domain are unchanged.
- Al26/Fe60: MPI2/OMP2 live source/transport/decay and restart PASS in
  `.radioactive-live.WwImg6`,83 final datasets exact. Actual age-convolved
  source signal agrees within measured subsequent stellar uptake and
  distinguishes the incorrect fresh-at-step-end treatment. Three evaluated
  raw directories and the restart symlink removed; compact evidence retained.
- Fe UV/thermal limits: both actual MPI2/OMP2 three-step comparisons PASS
  in `.fe-photon-live.z1eGRu`. Positive source/absorption, element/charge
  closure and admitted material temperatures are verified. Thermal Fe
  attribution is inferred from unchanged opacity weights, not a saved
  Fe-only final-energy field; UV collision heat is not photon heating.
- PAH single-photon atomization: actual MPI2/OMP2 weak/strong comparison
  and restart PASS in `.pah-atomization-live.Bs73kx/flux42`. Relative to the
  weak source, final PAH loss is2.29851e-6 of initial number; all gas carbon
  nuclei recover24 per destroyed grain to6.46e-11 relative. Identical star
  masses/birth epochs separate the effect from stellar uptake.7591 restart
  datasets match. Final zero primary fields alone do not imply no source:
  the paired comparison demonstrates measurable material deposition.
  Both evaluated PAH controls'16 raw HDF5/cooling files (116816128bytes)
  were removed after comparison; eight metadata directories plus inputs,
  binaries, logs, metrics and hashes remain. Fe's two evaluated raw outputs
  (~30.1MB) were likewise removed. Recovery requires rerunning retained
  inputs; no raw backups are implied by the compact evidence.
  None is silently enabled by approximation approval. Final frontend
  suite:57 passes/one display-dependent skip;
  native configuration:48 assertions pass (`.frontend-final.NGzCu0`).

Wind tracks contain continuous M/L/Teff/abundances, not LBV eruption histories.
Woosley2017 Table6 does contain actual two-pulse times/masses/energies for
T125A/B; complete matched eleven-element pulse yields are not supplied there.
Do not state that all pulse chronology is absent, or partition PARSEC totals
arbitrarily among these different stellar models.
