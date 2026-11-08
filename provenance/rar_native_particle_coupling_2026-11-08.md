# Controlled native particle coupling

The current RAR implementation bundle connects gravitating particles to the
native leaf/subface QUMOND operator and the existing common-time gas update.
The approved Astra design review is continued under the user's instruction
that an Astra driver performs the review directly. This is a coupling change
within the fixed, serial, two-level, periodic CPU experiment, not a claim of
cosmological or galaxy-formation readiness.

## Source and time contract

Each particle has a fixed physical cubic cloud of width equal to the base
cell width. Exact periodic cloud/leaf overlap gives nonnegative weights that
sum to one. Deposit `mass * weight / leaf_volume` before mean subtraction;
gather cell acceleration with the same weights. This agrees with CIC on a
uniform grid of that width. It is a controlled smoothing convention, not a
calibrated stellar or dark-matter softening prescription. The octree query
visits intersected leaves without constructing a uniform force grid.

The initial live roster comes from owned linked lists after the root split.
Identity, mass and species remain fixed; particle creation, deletion, SIDM,
tracers, sinks and PBH evolution are rejected. Cache reuse requires unchanged
positions as well as the existing gas and mesh checks. The explicit source
role is baryon (stellar species), CDM (Newtonian only), or diagnostic. The last
role admits synthetic particles for numerical controls and is not evidence
for a CDM-free galaxy model.

One common timestep obeys the gas constraints and a prospective kicked-drift
bound of half the smallest leaf width. Existing particle kick/drift consumers
are bypassed in this opt-in path. After gas transport/cooling and root tree
merge, the old field supplies the first particle half kick; a full drift and
periodic relink precede the new-source solve. The new field supplies the
second particle half kick and the existing endpoint gas half source.
Particle velocities are synchronized at each completed common endpoint.
Canonical base lists are checked before and after relinking; merged child
aliases are not counted as additional particles.

## Bounded validation and limits

Use Grammar Slurm and immutable source/build identities. Preserve all runs.
The initial moving-particle cases pair Newtonian and QUMOND forces at the same
gas/particle mass, with a QUMOND half-timestep replay. Independently reconstruct
cloud overlaps, deposited mass and force gathers from exported geometry and
particle states; check saved snapshots, each KDK transaction, and crossings
of native refinement and periodic boundaries. Retain the existing gas-only
regression and source/field checks. Separately measure a small self/pair-force
sample: matched deposit/gather weights alone do not guarantee exact
self-force cancellation or total momentum conservation on an AMR mesh.

Results must distinguish implementation, successful execution, and numerical
validation. A paired Newtonian diagnostic with identical particle mass is
not the Lambda-CDM science comparison unless an explicit CDM halo is added.
Cooling remains the validated level-local Lie split. MPI, regridding,
subcycling, restart, particle mass exchange, stellar/AGN feedback and native
cosmological initial conditions remain outside this first particle bundle.

Status at implementation commit: integrated source; compilation and execution
results pending. The default-off path is retained.
