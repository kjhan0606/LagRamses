# P0 Output Contract: lagRamses to S_N Input

Status: Draft v0.1  
Date: 2026-08-28

## 1. Confirmed checkpoint structure

The inspected `patch/lagRamses` HDF5 checkpoint writer provides a practical
starting point for the S_N converter. The following datasets are confirmed.

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
| `velocity[3]` | decoded momentum conservatives divided by density | confirmed raw availability |
| `temperature` | total energy, velocity, EOS, and cooling convention | requires exact `uold` variable map |
| `metallicity` | passive scalar or dedicated hydro variable | requires exact `uold` variable map |
| `x_HI`, `x_HII`, `x_HeI`, `x_HeII`, `x_HeIII`, `x_H2` | RT/chemistry checkpoint | required but unconfirmed in HDF5 |
| `dust_to_metal` | explicit snapshot field, or named subgrid prescription | required; no silent default |
| `cell_mask` and `level` | AMR grid center and refinement flags | confirmed raw availability |

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
2. Retain only AMR leaf cells and deposit their volume-weighted conservative
   quantities onto the static nested S_N blocks.
3. Preserve mass, each passive scalar mass, and total gas energy during every
   restriction operation.
4. Refined AMR cells overwrite their parent contribution; no parent and child
   values may be counted twice.
5. Record the input snapshot ID, expansion factor, physical time, converter
   version, and every source-table hash in the S_N input metadata.

## 5. P0 blockers to close before implementation

1. Decode the exact hydro `uold` variable ordering for the active HR5 build.
2. Establish whether HDF5 outputs all RT photon groups and non-equilibrium
   chemistry fractions; add an explicit exporter if not.
3. Establish the active dust model and whether dust is a stored field or a
   metallicity-derived subgrid quantity.
4. Define the sink accumulator interval and the retained-mass versus inflow-
   mass convention used to construct AGN luminosity.
5. Identify the star-particle type and age convention used by the active
   output.

No lagRamses source modification is authorized by this document. The next P0
inspection resolves these five items and produces a minimal, versioned input
schema for the converter.
