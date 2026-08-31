# Phase 0: Fortran-native stellar enrichment design

## Status

Phase 0 establishes the production interface for stellar mass return and
element-by-element enrichment in lagRamses. It does not yet change the
hydrodynamics, cooling, RT, dust, or star-particle evolution code.

The C libraries under `Eunha.A1/LagEunha/src/AstroModules` are reference
implementations and table sources only. Production runtime code will be
Fortran and will live in `patch/lagRamses/`.

## Runtime feedback selection

The compiled executable retains both feedback implementations.  The runtime
selection is made in `&STELLAR_ENRICHMENT_PARAMS`:

```text
feedback_mode='channel_resolved'  ! default
feedback_mode='legacy'            ! historical reproduction only
```

`channel_resolved` consumes the explicit per-channel mass and energy ledger.
Only the SNII returned mass loads the delayed-cooling reservoir; winds, AGB,
SNIa, and PISN do not.  Energy is read in erg from the normalized yield-table
contract and converted once using the RAMSES energy unit.  This path therefore
does not use the legacy `ESN2/NSN2` normalization.

`legacy` preserves the historical `feedback.kjhan3.f90` path, including its
total-`mloss` delayed-cooling update and legacy SN energy expressions.  It is
provided for controlled reproduction and must not be interpreted as the
corrected default physics.

## Scientific scope

The source model will support four stellar channels:

1. Massive-star winds.
2. AGB winds and planetary-nebula ejecta.
3. Core-collapse supernovae (SNII).
4. Type-Ia supernovae (SNIa).

Pair-instability supernovae (PISN) remain an optional channel and must not be
enabled implicitly. Dust production, dust destruction, and RT coupling are
subsequent phases, but this phase must expose the elemental source terms they
will require.

## Common element convention

The production interface uses eleven conserved physical species:

```text
H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe
```

`Z_TOTAL` is a diagnostic quantity derived from the physical species. It is
not an additional independent conserved field. This avoids double counting
metal mass when both individual elements and total metallicity are present.

The source interface must distinguish:

- `ejected_mass(k)`: actual mass of element `k` returned to the gas; it is the
  quantity deposited into gas cells.
- `net_yield(k)`: newly produced or consumed mass relative to the initial
  stellar composition; it may be negative and is diagnostic only.
- `returned_mass`: sum of the actual ejected elemental mass.
- `remnant_mass`: mass retained by the stellar remnant.

## Table contract

The Fortran table reader must provide a single normalized interface regardless
of the original table format. Every channel must document:

- Initial stellar mass range.
- Birth metallicity range and interpolation coordinates.
- Stellar age or delay-time range.
- Actual ejecta and net-yield definitions.
- Remnant mass prescription.
- Energy and momentum source units, when available.
- Extrapolation and out-of-range behavior.
- Literature source and table version.

AGB tables that only provide lifetime-integrated yields are not sufficient for
an exact time-dependent source. They must be combined with an explicit mass
loss history or a declared release-time prescription. The release-time
approximation must be kept separate from the elemental yield interpolation.

## Runtime source contract

Star particles represent single-age stellar populations. The source engine
must return cumulative quantities at two ages and deposit only their
difference:

```text
source = cumulative(age + dt, birth_Z, initial_mass)
       - cumulative(age,      birth_Z, initial_mass)
```

This rule applies independently to returned mass, every elemental ejecta
component, energy, and momentum. It prevents repeated injection when a star
particle is visited on multiple timesteps.

The production Fortran API should have the following logical layers:

```text
stellar_yield_tables
    load tables once per MPI rank
    interpolate cumulative yields and instantaneous rates

stellar_enrichment
    evaluate channel sources for one star particle and one timestep
    combine winds, AGB, SNII, SNIa, and optional PISN

stellar_deposition
    distribute source terms to neighboring AMR gas cells
    enforce positivity and local/global conservation
```

The exact Fortran procedure names are deferred until the existing RAMSES
particle and feedback interfaces are mapped in Phase 0 implementation work.

## Required star-particle state

The source calculation requires, at minimum:

- Formation time.
- Initial stellar-particle mass.
- Current stellar-particle mass or an explicitly documented fixed-mass
  convention.
- Birth total metallicity.
- Birth elemental abundances if net yields are reconstructed from initial
  composition.
- IMF and population mode, if these are not global parameters.

The mass ledger must satisfy:

```text
initial_mass = living_mass + remnant_mass + returned_mass
```

The implementation must choose one particle-mass convention and apply it
consistently to gravity, feedback, restart files, and diagnostics. A source
term must not be added to gas while silently leaving an incompatible stellar
mass ledger.

## Gas-cell deposition contract

For cell `i` with volume `V_i` and deposition weight `w_i`, the elemental
update is conceptually:

```text
delta_rho_k(i) = w_i * delta_M_ejected(k) / V_i
rho_k(i)       = rho_k(i) + delta_rho_k(i)
```

The same source event must update total gas mass, thermal energy, and
momentum using the existing RAMSES feedback conventions. The deposition layer
must not directly overwrite metallicity; it must update conserved densities
and derive abundances afterward.

AMR and MPI requirements:

- Use the existing particle-to-cell neighbor selection and boundary exchange.
- Accumulate all source channels before applying the cell update.
- Preserve conservation across MPI domains and refinement levels.
- Apply positivity limits without creating or destroying unaccounted mass.
- Make the result deterministic under restart and domain decomposition, within
  the declared floating-point reduction tolerance.

## Element-field decision

Phase 0 will define the field layout before implementation. The preferred
layout is physical elemental mass fractions or conserved elemental densities
for the eleven species above, with total metallicity derived as a diagnostic.
If the first production run cannot afford all fields, the reduced layout must
be explicitly recorded and must retain at least H, He, C, O, Mg, Si, and Fe for
the planned dust and RT work.

## Validation gates for the next phase

The Phase 1 implementation may start only after these contracts are fixed:

1. A zero-age/zero-timestep source returns exactly zero.
2. Time-integrating the source reproduces the cumulative table result.
3. Returned mass plus remnant and living mass closes the stellar mass budget.
4. Elemental ejecta sum to the declared returned mass within tolerance.
5. Net yields are never used as actual gas mass ejecta.
6. A single star particle gives the same integrated result independent of
   timestep subdivision.
7. MPI, AMR-level, and restart tests preserve the global mass and elemental
   budgets.

## Explicitly out of scope for Phase 0

- Dust grain growth, shattering, coagulation, or sputtering.
- RT photon source coupling.
- Recalibration of star formation or feedback efficiencies.
- Production simulation submission.
- Replacing the existing scalar-metallicity path before the new path has a
  conservation test.
