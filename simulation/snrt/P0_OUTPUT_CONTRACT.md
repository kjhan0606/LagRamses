# P0 Output Contract: lagRamses to S_N Input

Status: Draft v0.3 — output_00017 staging validated, production composition pending
Date: 2026-08-31

The canonical static input writer is now format version 2. The implementation
and synthetic contract test live in
[`snrt_core/snapshot.py`](snrt_core/snapshot.py) and
[`tests/p4_hdf5_staging.py`](tests/p4_hdf5_staging.py). This closes the file
layout and conservative-leaf bookkeeping portion of the contract; it does not
identify missing physical fields in an existing checkpoint.

## 1. Confirmed checkpoint structure

The inspected `cuRAMSES-kjhan` HDF5 checkpoint writer for the
`Horizon5-master-2` build provides the active output contract. It writes raw
conservative variables. The following groups are confirmed.

| Checkpoint group | Confirmed content | S_N use |
| --- | --- | --- |
| `/amr/level_*` | grid centers, refinement flags, CPU map | reconstruct leaf-cell geometry and rasterize a nested Cartesian domain |
| `/hydro/level_*` | one raw conservative dataset per `uold` variable | gas density, momentum, total energy, metallicity/passive scalars after variable-map decoding |
| `/gravity/level_*` | potential and acceleration | optional local-RHD boundary and diagnostic input |
| `/particles` | position, velocity, mass, ID, particle type, birth epoch, and metallicity when enabled | stellar-source catalogue and stellar SED assignment |
| `/sinks` | ID, mass, position, velocity, sink time, BH growth ledgers, stored feedback energy, angular momentum, spin, and efficiency | AGN source catalogue and AGN SED normalization |

The normal output driver also writes separate hydro, particle, sink, and RT
checkpoint files. The HDF5 chemistry/RT field content has not yet been
established as part of this contract and must not be assumed by the converter.

## 2. Converter input schema

The converter will produce a self-contained HDF5 input file for one static
S_N domain. Every field is in CGS units and every cell represents a physical
volume, not a RAMSES grid index.

### Required cell fields

| S_N field | Source | Status |
| --- | --- | --- |
| `rho` | decoded `uold_1` | confirmed raw availability |
| `velocity[3]` | `uold_2..uold_4` momentum densities divided by deposited density | implemented and checked against the raw-conservative writer |
| `temperature` | `uold_5` total-energy density reduced to pressure, then inverted with the thermal atlas | implemented with HDF5 header `gamma=1.6666667` |
| `metallicity` | `uold_6` metal mass density divided by `uold_1` and `0.02` | implemented; `uold_6` is confirmed by `imetal=6` in the producing build |
| `x_HI`, `x_HII`, `x_HeI`, `x_HeII`, `x_HeIII`, `x_H2` | RT/chemistry checkpoint | required but unconfirmed in HDF5 |
| `dust_to_metal` | explicit snapshot field, or named subgrid prescription | required; pilot map uses a recorded non-production placeholder |
| `cell_mask` and `level` | AMR grid center and refinement flags | leaf coverage is checked; `cell_level` is written in v2 |

The HDF5 adapter takes the explicit field map in
[`config/p4_hdf5_field_map_pilot.json`](config/p4_hdf5_field_map_pilot.json).
Density, momentum density, total energy density, and metal density are
volume-conservative. Primitive fields derived from those conserved quantities
are calculated only after restriction. Thus the current output is a validated
coeval interface snapshot, not a physical dust/chemistry production snapshot.

### Required source fields

| S_N field | Source | Status |
| --- | --- | --- |
| stellar position, mass, age, metallicity | `/particles` | confirmed, subject to particle-type decoding |
| stellar luminosity per group | external stellar SED table | converter responsibility |
| AGN position, velocity, BH mass, spin | `/sinks` | confirmed |
| AGN luminosity per group | accretion history plus AGN SED table | requires luminosity convention below |
| AGN emission time and source interval | output time plus BH accumulator interval | required metadata |

## 3. AGN luminosity convention

The sink checkpoint records `dMsmbh`, `dMBH_coarse`, `dMEd_coarse`, and
`eps_sink`. These quantities are sufficient to identify the intended
accretion state, but they do not by themselves define a safe luminosity unless
the accumulation interval and radiative-efficiency mass convention are
recorded with the snapshot.

The converter must therefore receive an explicit `mdot_BH` source field or a
sidecar source ledger containing:

1. the elapsed physical interval for each accumulator;
2. whether the recorded increment is inflowing gas mass or retained BH mass;
3. the radiative efficiency used for that interval; and
4. the intrinsic AGN SED family and normalization convention.

Until this ledger exists, S_N test runs may use controlled analytic AGN
luminosities but must not infer a scientific AGN luminosity from `dMsmbh`
alone.

## 4. Rasterization rules

1. Select a physical cube centered on the primary BH or on the dual-AGN
   barycenter.
2. Retain only AMR leaf cells and deposit mapped quantities onto the static
   Cartesian mesh using the map's declared averaging rule.
3. Preserve density, momentum, total energy, and mapped metal mass during
   restriction; derive velocity, thermal pressure, and `Z/Z_sun` afterwards.
4. Refined AMR cells overwrite their parent contribution; no parent and child
   values may be counted twice.
5. Record the input snapshot ID, expansion factor, physical time, converter
   version, and every source-table hash in the S_N input metadata.

## 5. Remaining production blockers

1. Certify the remaining `uold_8..uold_11` element/passive-scalar mapping for
   the exact checkpoint and determine whether any RT chemistry fractions are
   stored; add an explicit exporter if not.
2. Establish the active dust model and whether dust is a stored field or a
   metallicity-derived subgrid quantity.
3. Define the sink accumulator interval and the retained-mass versus inflow-
   mass convention used to construct AGN luminosity.
4. Identify the star-particle type and age convention used by the active
   output.

No lagRamses source modification is authorized by this document. The current
SNRT converter remains limited to the explicit pilot map until these items are
resolved and the production-contract gate passes.
