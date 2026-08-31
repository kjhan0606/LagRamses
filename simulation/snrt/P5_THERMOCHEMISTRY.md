# P5 static thermochemistry plan

P5 couples the P2 local photo-heating rate to the offline Grackle background
thermal atlas. The table convention is explicit: `net_rate > 0` heats and
`net_rate < 0` cools, both in `erg cm^-3 s^-1`. Local S_N photo-heating is not
contained in the atlas and is added exactly once.

`snrt_core.jax_thermal_atlas` holds the 4-D `(a, n_H, Z, T)` atlas as JAX
arrays and performs edge-clamped multilinear lookup inside XLA. The current
atlas has 1,228,800 entries and occupies about 5 MB as float32, so it is small
enough to replicate on a TPU device.

The pending operator will use a fixed-count subcycle within each transport
step: update photon/ionization state, then solve the backward-Euler energy
equation with the temperature-dependent signed Grackle rate using 24 fixed
bisection iterations. The residual is `u(T_new)-u_old-dt*(photoheat+R(T_new))`.
The table temperature bounds are used only if that residual has no in-table
root. Fixed subcycle and bisection counts preserve static TPU control flow.
Subcycle convergence must still be measured against a larger fixed count before
use.

The current P4 atlas ends at `a=0.20836524`, while output 00017 is at
`a=0.208497764676753`; runtime lookup therefore clamps at the high-a edge.
This is acceptable only for the current short numerical pilot. A production
run must build an exact-time or bracketing Grackle subtable before P5 results
are interpreted.
