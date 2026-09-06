# P4-A: RAMSES hydro snapshot ingestion

## Scope

P4-A makes a RAMSES hydro state usable by the static S_N radiation solver. It does not start the dual-AGN science analysis and it does not depend on a RAMSES `rt` output. Radiation sources are a separate, explicit catalogue.

## Canonical staging format

`snrt_core.snapshot` defines a versioned HDF5 input file. The file has proper-cgs fields:

- `/gas/hydrogen_number_density_cm3` and `/gas/helium_number_density_cm3`: nuclei number density.
- `/gas/temperature_k`: gas temperature.
- `/gas/dust_relative_abundance`: dimensionless multiplier for the supplied dust cross section.
- `/gas` attribute `dust_relative_abundance_origin`: `direct` or
  `metallicity_solar_times_dust_to_metal`.
- `/gas/velocity_cm_s`: optional three-component proper velocity.
- `/gas/metallicity_solar` and `/gas/dust_to_metal`: optional composition fields.
- `/ionization/x_hii`, `/ionization/x_heii`, `/ionization/x_heiii`: initial primordial number fractions.
- `/ionization/x_h2`: optional molecular fraction.
- `/grid/cell_level`: optional AMR source level for each static cell.
- `/sources/cell_index` and `/sources/photon_luminosity_s`: optional source positions and photon-number luminosities per RT group.
- `cell_width_cm` and `grid/left_edge_cm`: proper grid geometry.

Format version 3 is the current writer format. The reader remains compatible
with version 1 and 2 files when their dust state is unambiguous. A legacy file
with non-zero dust but no `dust_relative_abundance_origin` attribute is rejected
because its abundance origin cannot be reconstructed; zero-dust legacy files
retain the `direct` default. New derived staging records the exact product
`metallicity_solar*dust_to_metal` and validates it on read.

The source catalogue is deliberately not inferred from gas cells, particle masses, BH accretion rates, or a RAMSES `rt` output. The science workflow must supply an audited stellar/BH SED-to-group luminosity conversion.

## Current snapshot audit

The metadata at `/gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016/info_00016.txt` loads as a `yt` `RAMSESDataset` with a 1024-cubed root domain and redshift 3.7993. Its field index cannot be built when particle files are present: the particle header gives a two-value `nstar_tot` record where the stock `yt` frontend expects one value. A project-local hydro-only view, containing copied text metadata and links to only `amr`/`hydro` rank files, builds successfully. The available hydro fields are density, pressure, three velocities, and magnetic left/right fields; this output has no cooling fields.

This is a particle/frontend compatibility boundary, not evidence of bad hydro data. `stage_ramses_hydro_only` creates the safe view below a caller-specified project scratch subdirectory and uses the standard AMR/hydro reader. No source snapshot file is changed. An audited particle frontend patch remains a separate future task if particle catalogues must be read through `yt`.

### Registered transitional checkpoint

The stopped `lageunha` checkpoint at
`/gpfs/kjhan/Phase0_validation/10cMpc_h_z0_epsstar002/output_00011` is a
complete native binary package and is now the registered comparison state. Its
native `hydro_*.out*` files are the standard primitive-output route: density,
three velocities, pressure, then density-normalized passive scalars. This is
different from the raw-conservative HDF5 route documented below. The
checkpoint-specific map is
[`config/feedback_transition_phase0_native_field_map.json`](config/feedback_transition_phase0_native_field_map.json).

For this output, the core internal layout is `uold_1` density,
`uold_2..uold_4` momentum, `uold_5` total energy, `uold_6` total metal, and
`uold_7..uold_17` candidate Phase 0 element fields. The native serialized
values are respectively density, velocity, pressure, metallicity fraction,
and scalar fractions. A converter must therefore not apply the HDF5
conservative restriction/division rule to the native payload. The element
names behind `scalar_01..scalar_11` remain an inference until a self-describing
export or exact pre-selector source snapshot is available.

The native preflight ledger is
[`data/feedback_transition_phase0_output_00011_native_audit.json`](data/feedback_transition_phase0_output_00011_native_audit.json).
It confirms all 32 rank files for `amr`, `hydro`, `grav`, and `part`, but no
`rt`/`sink` component or force-source ledger. It is therefore eligible for
hydro/field-map and rasterizer development, not for a production source,
dust, or live-feedback comparison.

### Native yt field-order gate

The native `hydro_file_descriptor.txt` is the old unversioned descriptor
format. yt 4.4.2 does not parse that format and, for `NVAR=17`, falls back to
an MHD field layout. If the adapter accepts that fallback, its `Pressure` field
actually points at a later passive scalar and the derived temperature becomes
nonphysical. The hydro-only adapter therefore accepts an explicit ordered
`hydro_fields_in_file` list; the checkpoint-specific order is
`Density, x-velocity, y-velocity, z-velocity, Pressure, Metallicity,
scalar_01..scalar_11`. The temporary hydro view remains source-preserving.

The reproducible preflight driver is
[`tools/p0_stage_native_hydro_probe.py`](tools/p0_stage_native_hydro_probe.py).
It produced the small canonical interface artifact
[`data/feedback_transition_phase0_output_00011_hydro_probe.h5`](data/feedback_transition_phase0_output_00011_hydro_probe.h5)
and manifest
[`data/feedback_transition_phase0_output_00011_hydro_probe.json`](data/feedback_transition_phase0_output_00011_hydro_probe.json).
For the stopped state, the available atlas starts at `a=0.1639`, later than
the checkpoint `a=0.14854`; therefore this probe uses explicitly recorded
`mu=0.6` pressure conversion. It is a wiring/units preflight only, with
production gate closed for missing solar-normalized metals, dust, H₂, and
source luminosity fields.

### Native particle audit

The native `part_00011.out*` payload is also a compact cuRAMSES particle
binary, not a stock yt particle stream. Its header stores `nstar_tot` as an
8-byte `LONGINT`, while the yt 4.4.2 default handler expects a 4-byte record;
therefore the full RAMSES frontend cannot yet be used to read the particle
catalogue. The source-preserving audit driver is
[`tools/audit_native_particle_output.py`](tools/audit_native_particle_output.py),
and its result is
[`data/feedback_transition_phase0_output_00011_particle_audit.json`](data/feedback_transition_phase0_output_00011_particle_audit.json).

The audit streamed the type record without copying the 15.7 GB particle
payload and verified all 32 rank files: 134,217,728 DM particles, 42,342
stellar particles, and no sinks. Header totals and type counts agree. Stellar
positions, ages, metallicities, and source luminosities have not yet been
decoded into a photon-emitting source catalogue; the SNRT production source
gate therefore remains closed.

The dedicated metadata hand-off is now available from
[`tools/read_native_stellar_catalogue.py`](tools/read_native_stellar_catalogue.py):
[`data/feedback_transition_phase0_output_00011_stellar_catalogue.csv`](data/feedback_transition_phase0_output_00011_stellar_catalogue.csv)
and its manifest
[`data/feedback_transition_phase0_output_00011_stellar_catalogue.json`](data/feedback_transition_phase0_output_00011_stellar_catalogue.json).
It contains all 42,342 stars with normalized positions, current and initial
mass, native birth fields, and a proper-time age range of `0.002--750.1 Myr`.
The reader records `birth_epoch` as conformal time rather than incorrectly
using it as a scale factor; no SED or grouped photon luminosity is attached.

## Adapter contract

`stage_ramses_with_yt` and `stage_ramses_hydro_only` require an explicit `RamsesFieldMap`. The HDF5 route additionally requires the versioned map at `config/p4_hdf5_field_map_pilot.json`: every raw `uold_N` field declares its conservative quantity, unit, and volume or mass-weighted restriction. It never guesses field names, units, or thermal composition. A constant is allowed only with a reason and is marked as a pilot fallback in metadata. The adapter converts density to H/He number densities with the declared hydrogen mass fraction (default 0.76), then writes the canonical HDF5 file only after complete AMR leaf coverage has been assembled.

The HDF5 driver is `tools/p4_stage_hdf5_level15.py`. It scans all populated AMR
levels, deposits only `son_flag = 0` leaves, rejects incomplete coverage, and
records the density mass-balance error. `--preflight-only` validates all mapped
groups and dataset lengths without reading the field payload or writing output.
`--require-production-contract` is a hard gate: velocity, thermal state,
metallicity, dust-to-metal, H/He initial ionization, H₂, and source fields must
be mapped to snapshot datasets rather than constants. The pilot map intentionally
does not pass this gate because the current HDF5 export has no audited dust or
non-equilibrium chemistry fields.

## output_00017 audit

The exact producing build is recorded in the output sidecars as
`cuRAMSES-kjhan`, patch `Horizon5-master-2`, commit
`c59497a54b06eb655416dd04b4e02945fe8421eb`, with `NVAR=11`,
`gamma=1.6666667`, and `nlevelmax_file=15`. Its HDF5 writer stores raw
conservative `uold_N` arrays, not primitive pressure or velocity fields. The
active field map therefore uses:

- `uold_1`: density;
- `uold_2..uold_4`: momentum density, restricted conservatively and divided by
  the deposited density for velocity;
- `uold_5`: total-energy density, converted to thermal pressure with
  `P=(gamma-1)(E-|rho v|^2/(2 rho))`;
- `uold_6`: metal mass density (`imetal=6`), converted with
  `Z/Z_sun=(uold_6/uold_1)/0.02`;
- `uold_7`: delayed-cooling reservoir (`idelay=7`), not used as metallicity;
- `uold_8..uold_10`: three element-like slots from the yield-table path, kept
  out of the RT composition contract until their exact checkpoint mapping is
  independently certified; `uold_11` is zero in the inspected output.

The real 866,729,878,508-byte HDF5 payload remains an external source asset;
it is not copied or hashed as part of this code/data-contract step. The small
`info_00017.txt`, hydro descriptor, and build sidecars are available under the
`/gpfs` project tree. Both thermal atlas files were copied there and verified
against their source SHA256 values. A real-value interface staging run now
exists; a production interpretation and `--require-production-contract` run
remain blocked by the unresolved dust and initial chemistry fields, not by the
AMR/HDF5 layout.

Every future production staging run should record these items alongside the HDF5 file:

- snapshot path and output number
- bounding box, selected AMR level, and proper cell width
- field map and each source-field unit
- hydrogen fraction, dust normalization, and initial-ionization recipe
- source-catalogue construction and SED group definitions

## Local dependencies and check

Install optional staging dependencies only in `snrt/.venv`:

```bash
./.venv/bin/pip install -r requirements-ramses.txt
PYTHONPATH=. .venv/bin/python tests/p4_ingestion.py
PYTHONPATH=. .venv/bin/python tests/p4_hdf5_staging.py
PYTHONPATH=. .venv/bin/python tests/native_stellar_catalogue.py
```

The checks validate the v3 cgs H/He conversion, dust-origin attribute, geometry, source-group
catalogue, HDF5 round-trip, explicit field map, AMR coverage, and density mass
balance. They do not touch the production RAMSES snapshot.
