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

## Execution and independent checks

Numerical source `b86ee8b814f49518c4c70079a21476298cab2ed7` compiled in
Grammar job `1279259` (3 min 2 s). All following jobs use that immutable build.

| Execution | Audit | Result |
| --- | --- | --- |
| 1279263 | 1279264 | Execution completed; initial fixture missed the actual refinement interface. Preserved. |
| 1279266 | 1279269, 1279272 | Corrected interface positions. Initial audits exposed cancellation roundoff in reconstructed gas density. Preserved. |
| 1279266 | 1279275 | Three moving-particle cases passed independent checks. |
| 1279265 | 1279267 | Six gas-only cooling/heating and accuracy cases passed regression checks. |
| 1279277 | 1279280 | Five simultaneous particle/gas/cooling-heating cases passed. |
| 1279276 | in job | Six initial self/pair-force measurements completed. |

The actual refinement interfaces at the selected transverse coordinates are
`x=0.6875` and `0.3125`, rather than the nominal refinement-radius endpoints.
The corrected four-particle fixture starts at `x=0.6873, 0.9998, 0.3127, 0.0002`
with masses `0.002, 0.001, 0.002, 0.001`. Each N/Q/half-step case records two
refinement crossings and two periodic crossings. Deposited mass agrees within
`3.47e-18` code mass; independent force gathers agree within `1.78e-15` code
acceleration, including thermal cases. KDK transactions and independently
saved synchronized particle states pass. With cooling, halving the common
timestep changes particle position by `2.94e-8` of total travel and velocity
by `9.20e-10` relatively. These are short-interval replay differences, not
estimates of long-time or spatial convergence.

The independent gas reconstruction subtracts particle density from the
exported total. A mismatch of `1.03e-14` in one cell was roundoff from that
subtraction, not a different gas source. Its comparison now propagates an
absolute bound of `32 * eps * (abs(total_density) + abs(particle_density))`,
in addition to the existing `1e-14` gas-relative tolerance. Gas-only
tolerances, direct source/snapshot checks and momentum/work checks are
unchanged. No Fortran numerical change followed the implementation commit.

Run directories are under `/gpfs/kjhan/RAR_QUMOND/runs/rar_evolving_gas_JOB`;
audits are under `audits/rar_evolving_gas_JOB`. Each evolving case has only
two full snapshots. The self/pair runs have no full snapshots and explicitly
execute one step while analyzing the initial force only. Their first attempt
`1279271` exposed that `nstepmax=0` still performs one full step in this path;
the preserved failure led to the explicit `nstepmax=1` contract.

## Measured spatial limitation

A single synthetic cloud near the refinement interface has nonzero native
self acceleration. Its norm is `0.0477060 * G m / s^2` for Newtonian gravity
and `0.262895 * G m / s^2` for QUMOND, where `s=0.0625` is the fixed cloud
width. These normalizations are numerical cloud scales, not percentages of
an observed rotation-curve acceleration. Symmetric pairs cancel their total
force to relative `8.76e-16` (N) and `5.25e-14` (Q), but this symmetry does
not remove the single-cloud bias. Newtonian field superposition agrees to
`6.04e-10` relative maximum error; QUMOND fields are not superposed.

Therefore the source and time coupling is demonstrated within the controlled
domain, while the current particle force prescription is not ready for a
precision live-galaxy calculation. The next force work must reduce or bound
the assignment/AMR self-force error and test spatial convergence. Do not
silently subtract a single-particle QUMOND field: the nonlinear phantom
source does not permit that Newtonian correction. This limitation is
separate from the already tested prescribed, continuous baryon force path.

Compact results and reproduction tools accompany this note in
`provenance/results/rar_native_particle_bundle_1279259/` and
`provenance/tools/`. The broad RAR/BIMOND/FDM programme remains incomplete.
