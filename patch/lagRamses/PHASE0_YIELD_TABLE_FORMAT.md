# Phase 0 yield-table format

## Purpose

This is the canonical interchange format for the Fortran stellar-enrichment
modules. It is an ASCII format for provenance and development, and can be
converted into a read-only Fortran module for production runs.

## One data row

Each non-comment row contains the following whitespace-separated fields:

```text
channel initial_mass birth_metallicity age returned_mass remnant_mass energy
momentum_x momentum_y momentum_z ejecta_H ejecta_He ejecta_C ejecta_N ejecta_O
ejecta_Ne ejecta_Mg ejecta_Si ejecta_S ejecta_Ca ejecta_Fe net_H net_He net_C
net_N net_O net_Ne net_Mg net_Si net_S net_Ca net_Fe
```

Comment lines begin with `#`. The units are:

```text
initial_mass, returned_mass, remnant_mass, ejecta_* : Msun per initial star
birth_metallicity                                      : mass fraction
age                                                    : yr
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

The table values are cumulative through the stated stellar age. A timestep
source is obtained by subtracting the value at the previous age from the
value at the current age. A table that contains only a final lifetime-
integrated yield is not a time-dependent table.

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

