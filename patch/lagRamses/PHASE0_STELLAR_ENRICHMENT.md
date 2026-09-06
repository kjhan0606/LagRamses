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
It also requires an external canonical table named by `PHASE0_YIELD_TABLE`;
an unset, empty, truncated, or nonexistent path is a startup error and never
selects the compiled synthetic fixture. The IMF, population model, and all
five channel mass windows are read from the same namelist, for example:

```text
imf_id=2                         ! Chabrier (default when omitted)
population_model='single_star_ssp'
yield_source_basis='per_star_cumulative'
imf_mass_min_msun=0.08
imf_mass_max_msun=120.0
binary_fraction=0.0
channel_mass_min_msun=0.8,1.0,8.0,3.0,140.0
channel_mass_max_msun=120.0,8.0,40.0,8.0,260.0
```

Except for `imf_id`, these population-contract fields are mandatory for
`channel_resolved`; they must not be inherited from compiled defaults.
An omitted `imf_id` selects Chabrier on every read, never the previous
namelist's choice. Explicit invalid IDs are errors. The current production
integrator accepts only `per_star_cumulative`, because it applies one explicit
IMF convolution over the declared IMF mass support.  Per-event and already
SSP-integrated cumulative/rate bases are recognized but rejected until their
separate evaluators exist, preventing a second IMF convolution.  The effective
values and external table path are written to each RAMSES info file. Queries
outside the table's mass, birth-metallicity, or age domain fail instead of
clamping to an endpoint.

The implemented IMF shapes have a lower supported mass of `0.08 Msun`.
Channel-resolved startup rejects a configured IMF lower bound below that
value; the upper bound and every enabled wind/AGB/SNII channel window must be
explicit and internally nested.

IMF choices in `stellar_enrichment_params` (existing IDs preserved):

| imf_id | Model |
|---|---|
| 2 | Chabrier — default |
| 1 | Kroupa |
| 0 | Salpeter |
| 4 | Miller–Scalo, continuous three-power-law approximation |
| 3 | Historical Pop III — separate, not an ordinary-population default |

The existing Chabrier form is the individual-star lognormal
(`mc=0.079 Msun`, `sigma_log10=0.69`) joined continuously to
`dN/dm ∝ m^-2.3` above 1 Msun. Miller–Scalo uses `dN/dm` slopes
-1.4, -2.5 and -3.3, with breaks at 1 and 10 Msun and continuity coefficients
1, 1 and `10^0.8`. This is the piecewise approximation also documented by
[BASTA](https://basta.readthedocs.io/en/devel/_modules/priors.html#millerscalo1979),
not the original lognormal parametrization. Its end slopes are explicitly
continued over the configured support, including 0.08–0.1 and 100–120 Msun
when using this project's usual 0.08–120 interval. This is a modeling
convention, not additional stellar-yield support.

All choices normalize `integral(m*phi(m) dm)=1` over the configured mass
support; per-star yields receive exactly one IMF convolution. Changing the
IMF does not automatically reweight an already IMF-integrated SED or DTD
table. In particular, the approved N100/Maoz SNIa contract remains on its
explicit Kroupa basis: a different runtime IMF needs a matching approved
conversion, and the existing identity mismatch check must continue to reject
an inconsistent combination.

`binary_ssp`, SNIa, and PISN settings are parsed but are rejected by production
initialization at this gate.  This is a temporary fail-closed admission control,
not a reduced science scope.  Binary-population normalization and fate
ownership, the SNIa DTD/event model, and the PISN/PPISN population/core-mass
decision are mandatory active feedback gates.  A later reviewed PISN gate may
approve an explicitly disabled science configuration; it may not be skipped or
activated through the generic mass interpolator.

Only the SNII returned mass loads the delayed-cooling reservoir; winds, AGB,
SNIa, and PISN do not.  Energy is read in erg from the normalized yield-table
contract and converted once using the RAMSES energy unit.  This path therefore
does not use the legacy `ESN2/NSN2` normalization.  The separate legacy Sedov
call is also disabled in this mode, preventing a second SNII injection; any
mechanical source is taken from the channel-resolved momentum ledger.

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

Pair-instability supernovae (PISN) remain conditional on the selected stellar
population and must not be enabled implicitly.  Its explicit eligibility and
disable/enable decision is nevertheless a required gate. Dust production, dust
destruction, and RT coupling are subsequent phases, but this phase must expose
the elemental source terms they will require.

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
must return cumulative quantities at two explicitly named ages and deposit
only their difference:

```text
source = cumulative(current_age_gyr,  birth_Z, initial_mass)
       - cumulative(previous_age_gyr, birth_Z, initial_mass)
```

This rule applies independently to returned mass, every elemental ejecta
component, energy, and momentum. `indtab` is the persisted code-time marker
for `previous_age`; a pending interval is opened before evaluation and is
committed only after all gas-cell and particle updates succeed. A repeated
or restarted committed age is therefore an exact no-op, while a failed
interval remains retryable. This prevents repeated injection when a star
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

## HDF5 restart state

For `star` or `sink` runs, the `/particles` HDF5 group uses
`stellar_state_schema_version=1` and stores the active-particle datasets
`birth_epoch`, `tpp`, `mp0`, and `indtab` (plus `metallicity` when enabled).
The datasets use the same packed active-particle order as position, mass, and
particle identity. `tpp` is the formation-time marker used by the RAMSES
runtime, `mp0` is the initial stellar-particle mass, and `indtab` is the last
successfully committed cumulative-release cursor.

The HDF5 restore path requires the schema marker and all three release-state
datasets. A missing or unknown version is a hard restart error; it must not
silently retain `init_part` zero values. This protects both continuation
determinism and the no-double-counting invariant. Binary particle backups
already carry the same three fields and remain a separate compatibility path.

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
4. Tracked elemental ejecta do not exceed returned mass; the non-negative
   residual `returned_mass-sum(tracked)` is retained in gas mass and, together
   with tracked C--Fe, is deposited into the generic-metal field.
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
