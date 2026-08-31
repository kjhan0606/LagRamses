# Offline thermal atlas for runtime interpolation

## Division of work

The Grackle generator is an offline preprocessing program. It runs once per scale factor and produces immutable subtables. The simulation/runtime reads one HDF5 thermal atlas and interpolates it in this order:

1. bracket the current simulation scale factor `a`
2. trilinearly interpolate each bracketing subtable in `log n_H`, `log Z/Zsun`, and `log T`
3. linearly interpolate the two results in `a`

The atlas fields are `mu(a,n_H,Z,T)` and signed `net_rate(a,n_H,Z,T)`. `mu` enters the pressure-to-temperature inversion. `net_rate` is a background UVB/metal thermal term. It is not a replacement for local S_N photo-heating.

## Table schedule

The initial scale-factor schedule is [p4_thermal_atlas_scale_factors.txt](config/p4_thermal_atlas_scale_factors.txt). It has tighter coverage around `z>=14`, the scientific target, and includes the current `z=3.799` P4 snapshot only as an ingestion validation point. Every actual hydro output used by the production simulation should be inserted into this schedule before the offline atlas is frozen.

## Required future extension

The current hydro output has no metallicity field. Runtime staging therefore uses `Z=10^-6 Zsun` explicitly, not an inferred metallicity. Production output must write metallicity, and the atlas must receive that cell-wise field. Shielded and unshielded UVB tables should remain distinct atlas families rather than being mixed by an undocumented interpolation parameter.
