# Native stellar endpoint transaction

The operator instructed the driver to decide directly. Driver review adopts
a synchronized endpoint transaction for stellar return and birth. The
existing Astra reviewer independently recommended this ordering; that advice
is not a Fable verdict. Q-GOAL: this connects actual stellar mass, energy and
element exchange to evolving native QUMOND gravity. Q-LEAN: reuse the existing
channel-resolved feedback, cooling and star-formation routines, with a narrow
scheduling/roster interface. No new disc survey or force stencil is required.

The full goal still includes the Lambda-CDM/QUMOND comparison, stellar/AGN
feedback, enrichment and cooling/heating without dust or ray tracing, and
subsequent BIMOND/FDM. This stellar implementation does not close sink/AGN
coupling, MPI, regridding, restart, cosmology/ICs or galaxy calibration.

## Physical and numerical ordering

`RAR_QUMOND_STELLAR_SOURCES=1` requires the existing native common-time mode,
PIC and stars, the channel-resolved stellar build and feedback selection,
and `f_w=0`. The old gas/particle roster first completes gas transport and
particle KDK, including the synchronized endpoint gravity half-kicks.

At that endpoint, repartition particles onto their actual levels. Initialize
`unew` from synchronized `uold` on every level before any return. Apply the
existing SSP source once per level using current positions and velocities;
reverse-exchange contributions and commit physical leaves directly. Do not
call `set_uold` here because it would add another gravitational half-source.
Reject deposition into covered cells. Apply original cooling once, then form
stars on physical leaves, retaining the existing fine-to-coarse ordering.
Suppress the original recursive stellar/cooling calls only in this mode.

The selected isotropic source transfers stellar mass with its synchronized
bulk momentum and kinetic energy plus the existing channel energy/element
budget. Surviving stars keep position and velocity. The existing SSP progress
ledger determines the age interval and once-only return; no timestep-sized
mass debit is invented by this adapter. A future directed impulse or recoil
model requires its own explicit accounting.

Birth occurs at the completed endpoint, with local synchronized gas velocity,
`tpp=texp`, `mp0=mp`, and zero feedback progress. Newborns receive no drift or
kick from the preceding interval. Native IDs are allocated above the maximum
live stellar/CDM ID without changing the stellar count statistic. This avoids
collisions between imported CDM IDs and the ordinary stellar count namespace.

The old roster remains immutable during KDK. After sources, merge to canonical
base lists and match all old bodies by identity. Require unchanged species,
positions and velocities, nonincreasing stellar mass and exactly unchanged
CDM mass. Only well-initialized stellar births may be added. Rebuild storage
and assignment lookup from the new roster; retain the chosen physical cloud
width and source role, including QUMOND's CDM exclusion.

Increment the physical-source generation and recompute both N/Q stages at
the same time. Replace current fields while preserving the previous
distinct-time potentials. There is no third gravity kick. Ordinary cache
reuse still rejects unexpected same-time source changes outside this API.

## Birth correction

The existing sampler allowed up to 90 percent of cell mass in stars, while
the gas commit retained at least 50 percent. Events above 50 percent could
therefore create mass. The correction retains the existing 50-percent gas
floor and applies that same cap to the accepted stellar count before credit,
including the legacy loading factor in the debit bound. The gas debit is
committed with its corresponding particle credit. The native path requires
zero loading because the extra loaded mass otherwise needs an explicit
destination. Newborn `indtab` is initialized even for reused particle slots.

Optional `RAR_STELLAR_BIRTH_LEDGER=1` records actual donor densities, accepted
particle mass, volume, synchronized velocity and gas specific internal energy
at each birth. Removal of gas internal energy under the existing isothermal
star-formation prescription is not a conserved stellar thermal reservoir.

## Validation scope

Compilation and runtime validation are pending at implementation commit.
Use an isolated Grammar CPU build and one bounded native source fixture with
actual births, SSP return, cooling and both refinement levels. Check return
mass/momentum/energy/elements, accepted birth debit/credit, age progression,
newborn absence from the preceding KDK step, roster/list integrity and the
post-source field against its actual source. A deliberately high-depletion
case must exercise the corrected cap. Stochastic birth histories cannot
support a deterministic pointwise half-timestep claim without matched events.
Keep the existing deterministic particle/thermal and legacy source regressions.
Preserve all outputs and the already measured nonzero AMR particle self force.

Sink admission stays closed. Sink clouds have repeated negative sink IDs and
are rebuilt by `create_sink`; the next sink implementation needs a separate
physical sink-mass adapter and explicit body/cloud motion and accretion budget.
