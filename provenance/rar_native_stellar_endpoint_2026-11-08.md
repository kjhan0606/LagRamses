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

Compilation and runtime validation were pending at implementation commit
`0248da63cd5b9e4a6a90c6ce4e1fc3ad5876e6ba`; measured follow-up is recorded below.
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

## Measured implementation follow-up

Grammar build `1279381` completed in 3m26s with PHASE0 stellar enrichment,
CPU FFTW and HDF5, without dust or RT. Source `0248da6` was pushed after
compilation and bounded execution; it was not described as fully validated
before the independent audits below.

Run `1279399` completed in 1m07s (53.1s simulation time). It starts with no
particles and evolves a smooth moving baryonic gas on 960 physical leaves,
levels 3 and 4, in a periodic box. Four intervals of 0.002 code time produce
599 actual newborns at the first endpoint, followed by 599 return events at
each later endpoint. Both levels contain births and feedback receivers.
The deliberately large `eps_star=1000` is a synthetic depletion stress, not
a proposed physical efficiency. The accepted star-mass quantum is 1e-5;
the maximum donor depletion fraction is 0.499901676, below the 50% bound.

Independent audit `1279402` passed. Audit `1279406` additionally verifies the
two gas gravity half-kicks and exactly one post-source generation per time
step, using the same preserved outputs without rerunning the simulation.

| Check | Measured result |
|---|---:|
| Birth mass | 0.383810 code mass |
| Accepted return events | 1,797 |
| Stellar mass returned | 0.000581735384183595 code mass |
| Maximum gas-plus-star mass difference over all nine fields | 6.67e-16 code mass |
| Maximum individual return mass mismatch | 1.09e-19 code mass |
| Maximum whole return gas-inventory mismatch | 8.54e-16 code units |
| Eleven-element return receipts | 19,767 |
| Maximum element request/accepted-deposit difference | 2.12e-22 code mass |
| Old-roster KDK residual | 0 at exported precision |
| Maximum gas gravity transaction residual | 5.21e-17 code units |
| Cooling/heating leaf records | 1,797 / 2,043 |
| Maximum cooling-call energy receipt mismatch | 5.38e-17 code energy |

The independent particle-cloud integration, N/Q subface equations, phantom
source and gathered forces pass at all nine field epochs. Both full dumps
(after steps 2 and 4) agree with the post-source gas, particle masses,
positions, velocities, birth masses/times, committed SSP ages and gravity
fields. Newborns are absent from the preceding KDK roster and remain at their
donor-cell centres with the donor's synchronized velocity. No third gas
gravity kick occurs at the stellar endpoint. The unitless small residuals
describe this numerical fixture, not the accuracy of a galaxy force model.

Regression run `1279401` and audit `1279404` pass the existing five-case
particle/cooling Newtonian and QUMOND checks, including matched half steps.
Legacy source run `1279403` and audit `1279405` pass SSP elemental receipts,
AGN physical-mass inventory and cooling reconciliation. Native sink/AGN
admission remains closed: the legacy regression is not evidence for native
sink coupling. No deterministic timestep-convergence claim is made for the
new stochastic birth history.

The effective namelist, result JSON, field checks, regression summaries,
input/output hashes and exact remote paths are preserved in
`provenance/results/rar_native_stellar_1279381/`. New run/audit scripts and
their helpers are under `provenance/tools/`. Both full snapshots and all
raw logs remain in `/gpfs/kjhan/RAR_QUMOND`; nothing was removed.
