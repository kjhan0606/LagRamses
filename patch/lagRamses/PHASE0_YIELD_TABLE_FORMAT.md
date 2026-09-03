# Phase 0 yield-table format

## Purpose

This is the canonical interchange format for the Fortran stellar-enrichment
modules. It is an ASCII format for provenance and development, and can be
converted into a read-only Fortran module for production runs.

## One data row

Each non-comment row contains the following whitespace-separated fields:

```text
channel initial_mass birth_metallicity age_yr returned_mass remnant_mass energy
momentum_x momentum_y momentum_z ejecta_H ejecta_He ejecta_C ejecta_N ejecta_O
ejecta_Ne ejecta_Mg ejecta_Si ejecta_S ejecta_Ca ejecta_Fe net_H net_He net_C
net_N net_O net_Ne net_Mg net_Si net_S net_Ca net_Fe
```

Comment lines begin with `#`. The units are:

```text
initial_mass, returned_mass, remnant_mass, ejecta_* : Msun per initial star
birth_metallicity                                      : mass fraction
age_yr                                                 : yr on disk; converted once to age_gyr in memory
energy                                                 : erg per initial star
momentum_*                                             : g cm/s per initial star
net_*                                                  : Msun per initial star
```

The element order is fixed as:

```text
H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe
```

Channel identifiers are:

```text
1 = massive-star wind
2 = AGB wind
3 = SNII
4 = SNIa
5 = PISN
```

## Quantity definitions

`ejecta_*` is the actual material returned to the gas and must be
non-negative. `net_*` is the newly produced or consumed material relative to
the initial composition and may be negative. Only `ejecta_*` is a gas-mass
source term.

The table values are cumulative through the stated stellar age. A source
interval is represented explicitly by `previous_age_gyr` and
`current_age_gyr`, and is obtained by subtracting the value at the previous
age from the value at the current age. The conversion from on-disk `age_yr`
to in-memory `age_gyr` occurs exactly once in the table loader (and the
embedded generator emits the same in-memory unit). A table that contains
only a final lifetime-integrated yield is not a time-dependent table.

## Source-cell mass assignment

`stellar_yield_table_t%mass_assignment_mode` is explicit runtime metadata.
The legacy/default `linear` mode is suitable only for quantities whose source
model authorizes continuous mass interpolation. The
`piecewise_constant_source_node_mass_cell` mode assigns an interior query to
the left source node using the half-open convention `[m_i, m_{i+1})` and uses
the exact node at the final edge. This prevents a terminal fate, remnant, or
event outcome from being fabricated by interpolating across source nodes.
The mode does not authorize interpolation in metallicity, rotation, engine
branch, or age; an approved fate resolver must supply those axes as exact
source-node coordinates and carry the source provenance key.

## Legacy-table conversion rules

The legacy AGB and massive-wind headers contain useful elemental net yields,
wind masses, and metallicity grids. They do not by themselves provide a
resolved release history. They may be converted to this format only after an
explicit release-time model or a source table with an age axis is supplied.
The converter must never silently assume that the full AGB yield is released
at birth or uniformly in time.

The legacy SN table is also a reference implementation, not the final
production calibration. Its mass and metallicity interpolation must be
reviewed together with the selected yield literature before being used for a
publication run.

## Embedded production mode

The generator
`generate_stellar_yield_module.py` accepts this canonical ASCII format and
produces a Fortran module containing read-only arrays plus a loader for
`stellar_yield_table_t`. The generated module records the input path and
SHA256 digest. The external reader remains available for development and
table comparison.
