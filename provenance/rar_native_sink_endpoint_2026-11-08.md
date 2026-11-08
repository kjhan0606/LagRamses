# Native physical sinks and synchronized AGN sources

Q-GOAL: connect physical black-hole gravity and motion, existing gas accretion,
and AGN deposition to the evolving native Lambda-CDM/QUMOND calculation.
The full goal also retains stellar return/birth, enrichment, cooling/heating,
cosmological initial conditions and subsequent BIMOND/FDM. Dust and ray
tracing remain excluded by the operator. Q-LEAN: reuse the existing Bondi,
sink creation/merging and AGN routines and the native fixed-cloud assignment;
do not introduce a new accretion law or an unrelated benchmark programme.
Driver review follows the operator's direct-decision instruction and the
existing Astra review's separate physical-sink adapter recommendation.

Physical `msink/xsink/vsink/idsink` constitute the gravitating sink roster.
Repeated negative-ID cloud particles remain gas-sampling proxies and are
excluded from the ordinary stellar/CDM gravity roster. Deposit physical mass
once and gather force using the same volume-integrated cloud weights. Evolve
physical bodies with the synchronized KDK and move the cloud coordinates
coherently; suppress their legacy gravitational motion and cloud-averaged
overwrite of physical body coordinates. The ordinary particle-list checks
must still detect duplicate, missing or foreign cloud slots.

At the completed common-time endpoint, retain the existing SSP return
transaction, then rebuild/merge/create sink sampling clouds, sample gas,
perform existing Bondi accretion once per active level and apply existing AGN
feedback once. Repartition lists before level-local operations. Cooling and
stellar formation follow. Reconcile both physical rosters and publish the
changed source at the same time without another gravity kick. Initial cloud
construction occurs once before the initial field. Ordinary non-opt-in
execution retains its scheduling.

Accretion retains the existing radiative-efficiency convention: gas removal
exceeds retained BH mass, with the remainder an explicitly recorded radiative
mass/energy channel. A gas+BH mass sum alone must not be asserted conserved.
AGN mode labels alone do not validate deposition. Record physical-body KDK,
gas/BH inventories around creation, accretion and feedback, and verify nonzero
accepted accretion and a feedback event against existing source receipts.
Existing stellar, cooling and particle regressions remain required. Use new
Grammar Slurm run directories and preserve all outputs.

Initial implementation is commit `e24e4a7`. Grammar build 1279487 completed;
seeded AGN-on/off executions 1279494 and 1279495 stopped at the first endpoint
because the legacy cloud rebuild found no exact central proxy. The scheduler
correctly marked the runners failed despite the executable's clean-stop exit
code zero. No successful native sink evolution or feedback result is claimed.

The bounded repair retains the known central particle slot across the serial
KDK step, requires its periodic displacement from the physical body to be at
most 32 machine epsilons of the box, records that displacement, and restores
exact equality before legacy cloud removal. This handles possible scalar/
array rounding without accepting a displaced body or selecting a nearest
cloud. The next execution must confirm the measured discrepancy and whether
this resolves the failure. Opt-in accretion receipts record donor states,
gross/retained mass, efficiency, radiative energy and retained momentum for
independent accounting. They do not change the accretion prescription.

Runtime verification is pending. Existing restrictions
on serial CPU, noncosmological periodic fixed two-level geometry, common
timesteps, restart and regridding remain explicit. This is not a completed
galaxy-formation or cosmological validation.
