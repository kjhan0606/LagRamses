# Approved Fe photon comparisons

These are explicitly approximate comparisons, not full metallic-Fe
photoelectric, magnetic, or keV/Auger physics. No stellar-source changes.

* `fe_electric_compare_v1` remains unchanged: primary absorption <=4 eV.
* `fe_thermal_limit_v1`: fixed-group primary absorption through 10000 eV,
  all accepted energy retained by the existing solid sensible reservoir.
  This is an upper-retention limit, not a measured photoelectric yield.
* `fe_uv_cycle_v1`: fixed-group <=13.6 eV, HD2017 photoemission and
  photodetachment with the existing composite Fe dielectric, classical OML
  collisions and their matching kinetic-energy moments. Equilibrium charge
  populations are worker-local scratch, not persistent hydro fields.

Both new comparisons retain the existing <=300 K Fe material restriction,
static/coadvected grains, and exclude PAH/relative motion and kind7. UV
requires CHIMES157. Setup/HDF identity plumbing is owned by the main agent.

UV uses stationary adjacent-edge catalytic cycles, so no untracked net
electron or grain-charge creation occurs. H+/C+ neutralization returns a
ground-state atom with full accommodation, outgoing mean kinetic energy
2 k Tdust, no adsorption residence and no extra recombination photon. This
is a declared surface-outcome approximation. Gas ionization energy is
debited by the actual H+/C+ -> H/C state change, and the released
13.598440002498787/11.260291860855007 eV per ion is credited to Fe sensible
heat. These are the same ATcT v1.130 H/C formation-enthalpy differences in
`snrt_chimes_atomization.h`, not aliases of the rounded FS thresholds.
It is NOT also added to gas
thermal heat. For an ion/electron capture cycle, the explicit balance is
ion KE + electron KE + ion IP = returned atom 2 k Tdust + Fe heat.
For a photoemission/electron-capture cycle, photon energy + captured
electron KE = escaping electron KE + Fe heat; grain F cancels between the
two opposite edges. The escaping energy is then split into gas heat,
actual secondary ionization binding changes and an excitation ledger.
The FS callback supplies unchanged target-limited ENERGY fractions. Fe
secondary event counts divide by the same H binding cost above; no offset
heat is inserted to conceal a binding mismatch. The existing CHIMES helium
secondary costs, 24.59/54.42 eV, are retained; there are no Fe/He ion-impact
channels in this comparison. The independent CHIMES gas solver is unchanged.
These are two complete catalytic cycles, not one photoemission plus ion
capture (both of those alone increase Z).

Escaping-electron distributions,
not merely their mean, feed the existing FS2010 target-limited partition.

The actual split is transport/C-sil absorption -> CHIMES -> Fe UV photon
and catalytic step -> existing mixture IR/exchange -> collective commit.
Only UV removes Fe absorption from the earlier transport opacity. Fe bulk
absorption and photodetachment then compete in one photon debit. Existing
Fe scattering/IR/material remain. Neutral-H-only Fe collision area in IR
avoids double counting the explicitly processed H+ impacts.

Gas/photon reservoirs evolve for finite time, with rates recomputed as
they change; equilibrium is used only to eliminate the trace grain charge
distribution. The comparison requires omitted absolute charge/electron
inventory <=1e-3, omitted absolute charge energy/material-plus-gas energy
<=1e-3, and charge relaxation time <=1e-2 of the integration/forcing time.
No grain-charge projection or arbitrary thermal correction is permitted.
Finite/positive state, elemental, charge, photon and energy checks precede
any publication. Signed collision receipts are not photon absorption.
The relaxation test uses a charge-diffusion estimate, 2 max(Var Z,1) /
sum f_Z (Jup+Jdown), not a certified slowest eigenmode. Gas T is restricted
to (0,10000] K. Explicit substeps limit each source forcing fraction to
0.05, cap at 10000, and freeze the pre-source dust T for this first-order
split; final material temperature is checked by the existing implicit IR.
The source receipt requires the version-4 exchanging IR contract. Fixed
group photon storage retains the pre-existing FP32 transport precision.

Aggregate Fable review correction: the blanket '<9 eV photoelectrons'
statement is false for negative grains. At a=10 nm, Z=-10, the fourth
reference photon (12.2954112567 eV) has IP=3.1089942468 eV and detaches an
electron of 9.1864170099 eV. At Z=-50 it is 14.9462752 eV. Population
weights can be tiny but are not assumed zero; the existing FS callback
is evaluated at the actual emitted energies/distribution quadrature.
The stationary high-Z outward current must be <=1e-12 of total charging
activity; a non-negligible tail rejects instead of reflecting it silently.

Native checks reuse the project's standalone smoke style: thermal-limit
identity/retention, threshold/negative-charge photo kernels, OML moments,
UV actual photon/species/heat exchange, trace/relaxation rejection and
atomic rollback. Native/private builds only; main owns live MPI evidence.

Main's aggregate review: Q-GOAL first (bounded approximate live coupling
without full-physics claims), Q-LEAN second (one native kernel plus existing
driver/IR transaction, no charge carriers or new audit framework).

## Implemented interface and evidence

Existing namelist field `dust_iron_model` selects the two new models; no
new namelist parameters/passives. CPU OpenMP material backend, fixed groups,
existing 10/100 nm six-material optics. UV additionally requires
`dust_cooling='chimes_neq_v1'` and the version-4 exchanging IR contract.
UV live speed is cm/s (`snrt_c_cgs*reduced_c`), unlike the dimensionless
speed fraction in external CHIMES controls. All persistent cell/RT/IR
publication remains inside the existing collective transaction.

`dust_mass_physics`: `dust_fe_photon_identity_n=12`, function
`dust_fe_photon_identity()` returns schema, model ID (0/1/2/3), primary
ceiling, Tdust ceiling, trace-charge bound, trace-energy bound, relaxation
ratio, then the five OML/accommodation/ground-return/cycle/fixed-energy
convention identifiers. Main appends this only for new selectors.

`dust_iron_photons`: `fe_photon_data_identity_n=95`, function
`fe_photon_data_identity()` has schema 2. Fields 1--35 bind Coulomb/work/
escape-length constants, nine reference energies, eighteen attenuation
lengths, source fraction, tail bound, gas-T ceiling and step cap. Fields
36--39 are zlo/zhi; 40--53 HD coefficients; 54--57 OML coefficients;
58--61 physical constants; 62--67 numerical choices; 68--83 quadrature;
84--88 ATcT inputs/conversion; 89--92 gas binding costs; 93--95 UV ceiling,
light speed and gas heat-capacity factor. Main appends it for UV only.
All entries are the parameters used by the kernel, not separately rounded
metadata. `fe_gas_ip_ev(4)` publicly exposes the gas reference costs.

Fe-only Makefile additions: pure `dust_iron_photons.o` in the base module
list (no CHIMES link dependency), dependencies on dust mass/Fe optics/data,
CHIMES live adapter and backup-HDF ordering, and standalone
`dust_iron_photons_smoke`. Its CHIMES-enabled form uses the actual pinned
FS callback; disabled form uses an explicitly all-heat algebra fixture.
RADIOACTIVE and VPATH were preserved.

Private evidence directory: `/tmp/fe-photon-native.x32KBD`.
Exact reproducer: `bash /tmp/fe-photon-native.x32KBD/check_fe.sh`.
Logs: `check-fe.log`, `native-smoke.log`, `native-fs-smoke.log`, and the
per-module `*-compile.log`. Tests include hard full-retention receiver,
H/C chemical energy closure, negative-charge energetic electrons, exact
serial/two-worker agreement, rejection with unchanged outputs, and actual
nonzero FS excitation. Live adapter/driver compilation uses existing
read-only dependency modules; the no-dust NVAR21/NENER1/RADIOACTIVE compile
checks the disabled source path, NOT a fresh NVAR21 ABI build. Main owns
full-build/MPI/restart validation; none was launched here.

## Final bounded live evidence (2026-09-11)

The subsequent main-authorized live tests used one PRIVATE copy of the
combined CPU/CHIMES1/HDF5/NENER1/NVAR3773 executable, SHA256
`b3c58d44bbba17cd84ecb41ad0390854ac255f8073c17eb7588a7a0c59591fc9`,
at `/gpfs/kjhan/LRD_JWST/.fe-photon-live.z1eGRu/ramses_fe_frozen3d`.
No shared-source changes or extra simulations were needed for this evaluation.

Exact launches, from `/gpfs/kjhan/LRD_JWST/.fe-photon-live.z1eGRu`:

```sh
bash ./run-case.sh thermal-retry
bash ./run-case.sh uv-retry
python3 evaluate.py > evaluation.json
```

The runner verifies the frozen hash, sources each private `environment.sh`,
then invokes `mpirun -np 2 ../ramses_fe_frozen3d physical.nml` with OMP2.
Both runs completed three steps, ~67.523 yr, 64 uniform level-2 leaves
(four grids/rank) and 192 actual stellar particles. Wall times including
MPI startup were 10.093 s thermal and 13.490 s UV. Rank-1 active emitting
sources were 0,32,64 over the three intervals. Spectra are deliberately
synthetic constant per-initial-Msun reference-group controls: group 9 at
4023.594574 eV for thermal, group 4 at 12.295411 eV for UV. These are not
PARSEC/matched physical SSP spectra, AGN, kind7, restart or AMR evidence.

The initial `thermal/run.log` records a FAILED namelist admission, despite
exit status zero; it produced no dump and is retained unchanged. Its
all-feedback-off request contradicted `user_selected_model_v1`. The
authorized retry restores wind/AGB/SNII/SNIa and the pre-existing matched
binary effective-SSP SNIa contract; it does not retag a source or change
Fe physics. Condensation remains zero. Native mass/energy feedback is
enabled, not assumed absent: final `mp0-mass` is exactly zero for every
star at this precision, so no measurable mass return occurred in these
short intervals. SNIa's 40-Myr minimum delay is well beyond the test age.
The earlier claim that all channels were disabled applies only to the
rejected inputs. No claim of an exercised wind/SN feedback receipt is made.

| Measured final-state/transaction check | Thermal | UV |
| --- | ---: | ---: |
| Maximum element relative discrepancy, gas species + C/sil/Fe solids | 6.425e-16 | 1.028e-15 |
| Maximum charge relative discrepancy | 1.378e-16 | 2.786e-16 |
| Gas + stars mass relative error against initial box mass | 1.826e-16 | 1.826e-16 |
| Common dust temperature (K) | 10.00645 | 10.09350 |
| Gas temperature (K) | 105.71627 | 105.71624 |
| Largest logged IR balance residual | 3.255e-13 | 9.8312e-10 |

All 3773 stored hydro fields and the primary/IR payload are finite; species,
photons and IR energies are nonnegative, thermal energy positive, all three
RT/IR transactions committed, and chemistry failure counts are zero. The
IR residual is the existing solver's local transaction check, NOT an
independent whole-history gas/chemical/radiation energy closure measurement.
Final solid/IR energies and exact nucleus/charge bookkeeping are retained
in `evaluation.json`; no additional energy tolerance was changed.

Fe attribution is deliberately distinguished between the cases:

* Thermal: 3.400329e44 photons emitted (birth-epoch/initial-mass integration),
  3.397594e44 remaining, and 2.197094e44 eV logged gas absorption. The
  remaining photon-budget debit is approximately 8.806644e44 eV. The
  unchanged six-material fixed-group absorption coefficients put
  17.4338978% into metallic Fe, implying approximately 1.53534e44 eV
  retained by Fe before common-mixture IR/exchange. This is an inferred
  photon-residual attribution subject to existing FP32 transport and log
  precision, NOT a directly saved per-Fe capture counter or final Fe-only
  sensible energy. Both Fe radii have positive group-9 absorption opacity
  (243.1184 and 239.9938 cm2/g). Scattering is not included as absorption.
* UV: the directly logged Fe photon receipts are 0, 4.9462044725e44 and
  1.4246217339e45 eV, totaling 1.91924218115e45 eV. The associated catalytic
  gas-thermal/solid/excitation ledger sums are -2.59171915671e49,
  +1.3678305287e52, and 0 eV. Solid heat includes ion/electron captures and
  released gas binding energy, including in the first dark interval;
  it must NOT be equated to photon energy. The equilibrium-trace and
  finite-time budget checks succeeded, but their maxima are not recorded
  in this live log. This live case does not demonstrate nonzero secondary
  excitation; that evidence remains the existing native FS smoke.

The native evidence has been preserved on GPFS at
`/gpfs/kjhan/LRD_JWST/.fe-photon-native.PFk9fW`, including historical logs.
Its `check_fe.sh` working path was updated to GPFS; the original `/tmp`
tree was not deleted. Earlier temporary-path instructions are historical.

Effective retry NMLs/environments/SEDs, failed and completed logs, frozen
binary, launch identities, snapshot metadata, compact `evaluation.json`,
and cleanup manifest remain in `.fe-photon-live.z1eGRu`. Each completed
run produced one ~15.05-MB raw output directory; after evaluation only
those two raw output directories are removed per project retention policy.
No unrelated outputs or checkpoints are touched. This completes the two
bounded Fe comparison live cases, not full metallic-Fe hard-photon physics.
