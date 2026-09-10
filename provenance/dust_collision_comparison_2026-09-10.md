# Collision-aware grain-size decision: native reference complete

Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
Approved medium-term group 6: compare two-size and collision-aware multibin
physics and decide whether replacing the live representation is justified.
This is **not** a selectable live multibin hydro implementation.

## Implementation and assumptions

`dust_mass_physics.f90` now has an actual native coagulation/shattering
operator, not just the earlier imposed radius shift. The existing smoke
executes it. No Python framework, NML change, new passives, VPATH change,
RAMSES launch or snapshot was needed; current production defaults remain.

Physical model: [Hirashita & Aoyama 2019, sections 2.4/2.5](https://arxiv.org/html/1810.07962v1).
Compact spherical same-material grains, geometric collision cross-section,
their empirical turbulent speeds, eight fixed Gauss-Legendre relative angles.
Average the collision-weighted yields, not yield at average speed. Unordered
pairs include both targets; self pairs have the 1/2 event factor. Coagulation
sticks without an extra threshold. Shattering uses material disruption
strength, a remnant, and an analytic power-law fragment mass spectrum.

Numerics: fixed-pivot number/mass-conserving single daughters; analytic
fragment mass integrals in geometric-pivot cells (fragment number is only a
discretization approximation). Explicit Euler limits net donor depletion to
2%, then 1% in the comparison. Same-bin remnants are included in that net
rate, avoiding fictitious gross-collision stiffness. Positivity is not
enforced by clipping; no mass renormalization. Inputs and product budgets
are checked, updates staged, failures return unchanged input state.

Mass outside the radius grid is a separate inert numerical reservoir. It
does not collide again, become gas, or acquire PAH chemistry. No opacity,
gas heating, kinetic energy or dust-drift coupling is inferred. This is not
an exact reproduction of the paper's discretization or random-angle solver.

## Paired physical results

Initial MRN a^-3.5 on 0.005--0.25 micron, dust/gas=.01, gas mass per H=1.4 mp.
Each material is an isolated comparison with that total dust abundance,
not simultaneous C and silicate each doubling an actual mixture's dust mass.
Grid: 16/32/64 pivots from 0.0003 to 10 micron. Carbon density 2.2, silicate
3.3 g/cm3. Dense: nH=1000/cm3, T=100 K, Mach1, 3 Myr. Diffuse: nH=.3/cm3,
T=8000 K, Mach3, 30 Myr. Small means a<0.03 micron; all mass fractions below
are normalized to **initial** dust mass, not surviving dust mass. Area is
geometric projected area per initial dust mass, not wavelength opacity.

| Material / regime | Two-size small fraction | 64-bin small fraction | Two-size area / 64-bin area | Mass below minimum radius |
| --- | ---: | ---: | ---: | ---: |
| C / dense | .05171 | .12334 | 1.1165 | 0 |
| Silicate / dense | .06998 | .13690 | 1.2146 | 0 |
| C / diffuse | .35103 | .61307 | .16517 | .28792 |
| Silicate / diffuse | .31748 | .54217 | .27477 | .07763 |

The initial true small fraction is .23875 versus the 64-pivot projection
.22698; the existing two-size radii already overrepresent the initial MRN
area by a factor 1.957. Therefore these are **closure AND representation**
differences, not pure collision-rate errors. The initial 64-pivot area is
within .183% of the continuous MRN value. At the endpoint the 32/64-bin
area changes are .53--2.01%; small fractions change by .0030--.0049 absolute.
Halving the depletion fraction changes area/small observables by <=.254%
across the matrix. No monotonic convergence of every moment was assumed.

## Driver end decision

[One Fable plan review](dust_collision_fable_2026-09-10.txt) approved the
comparison, asking for an explicit criterion and a smaller matrix. The
[driver response](dust_collision_plan_2026-09-10.md) specifies 10% area and
.05 absolute small-fraction differences as a bounded adequacy criterion,
not a universal physical acceptance threshold. We retained both cheap
material refinements because disruption strength changes the cascade.
The driver performs this end evaluation under the operator's latest role.

1. **Do not claim two-size and multibin equivalence.** Dense comparisons
   exceed that criterion; the difference is larger than measured time/grid
   changes. Continuous size information matters for surface-dependent work.
2. **Retain the live two-size default.** HA19 is a different approximate
   closure, not an observed truth. In diffuse gas 7.8--28.8% is in a
   noninteracting numerical reservoir; that prevents physical qualification
   of the long-time shattering prediction even though total mass closes.
   Do not silently discard it or reinterpret it as gas conversion.
3. A multibin science option has demonstrated value for size/area-sensitive
   predictions, but this comparison does not authorize calling an unwired
   mass spectrum a complete live RT model. Size-dependent optical/thermal
   coupling and a physical sub-grid-size fate would be needed to replace the
   live representation. This is a recorded limitation, not another invented
   gate or automatic hydro-variable expansion in this bundle.

Group 6's bounded comparison/decision is complete. Groups 1--5, remaining
group-7 PAH carbon chemistry, and group 8 are **not** marked complete.

## Native evidence

Existing `.pah-hydrogen.6fDwRw/dust_mass_smoke`, Intel MPI/ifx, passes full
native suite: `build-collision-verified.log`, `test-collision-verified.log`.
SHA256 `9d21f9e788e7694105db9da83e78b55675b43253bbe3b549082551edf928a3b9`.
No new RAMSES evolution or MPI/CUDA claim is made for this reference.

GNU full suite also passes in `.pah-hydrogen.6fDwRw/gnu/`:
`build-collision-verified.log`, `run-collision-verified.log`, with
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Executable SHA256 `7a62432b6d4e5ab64ab73b5fc8de9cc12cfe55e6578888959594200979634521`.

Mass plus explicit boundary reservoir closes to 4.94e-12 relative or better
across both compilers (Intel maximum 3.37e-12).
Constant-kernel analytic total-number relative errors halve from .000999
to .000500 on timestep refinement. The equal-speed isotropic collision
mean differs from its analytic 4v/3 by .0002534 with eight angle nodes.
Zero time/rate identity, positive states, invalid-kernel rollback,
corrupt product mass-budget rejection, analytic fragment underflow and full
coagulation upper-overflow are covered. A first upper-overflow analytic
test used too coarse an Euler step for its stated tolerance; halving that
test's step resolves it without loosening tolerance (failed logs retained).

Source hashes:

- dust_mass_physics.f90: `cf89c8371fd105c49a418661f3c9bce947b85224b5d11dc90bcf0a6739303009`
- dust_mass_smoke.f90: `db3a0739938ac1554f4cae93e4af8d334b235f34f001a52667e5ed0f43fd036d`

Only compact native logs and binaries were generated. No raw simulation
outputs were created, so no output cleanup was required for this bundle.
