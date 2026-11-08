# Native physical sinks and synchronized AGN sources

Status: the seeded physical-BH adapter, Bondi accretion, thermal/jet source
budgets, stellar sources and cooling passed the bounded paired audit
`1279553`. Tested Fortran revision is `a9375731598d3671f1d0c8d3ab55037b657c63ad`.
This is an implementation/coupling result, not a galaxy or cosmological result.

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

Build 1279509 / commit `a2a5c43`, runs 1279511 and 1279512, ruled out the
initial rounding hypothesis: central-proxy coordinate differences were
exactly zero and the same failure persisted. The unnecessary coordinate
repair was removed. The actual list precondition was violated: SSP leaves
particles partitioned by level, whereas `create_sink` begins by repartitioning
base lists and clearing all child headers. Merge to canonical base lists
before invoking `create_sink`; verify cloud mass before and after the source
transaction. This preserves the existing cloud construction and Bondi rules.
The no-sink stellar regression 1279513 completed on build 1279509; independent
audit 1279517 was submitted separately.

Opt-in accretion receipts record donor states,
gross/retained mass, efficiency, radiative energy and retained momentum for
independent accounting. They do not change the accretion prescription.

Existing restrictions
on serial CPU, noncosmological periodic fixed two-level geometry, common
timesteps, restart and regridding remain explicit. This is not a completed
galaxy-formation or cosmological validation.

Build 1279524 (`051b4e1`) completed. AGN-on/off runs 1279526/1279527
completed all four steps after the list repair. Independent audits 1279528
and 1279533 retained the completed physical outputs: the first exposed a
Fortran implicit-E exponent parsing issue, and the second exposed missing
AGN source-ledger rows. The native endpoint follows `update_time`, so the
legacy `nstep_coarse_old` restart guard suppresses new diagnostic rows after
step zero. Exempt the native path from that guard, retaining the per-step
duplicate guard; this changes diagnostic output, not AGN deposition.
Accretion receipts now use explicit three-digit exponents. The reader also
supports valid historical implicit-E notation. Rerun from the corrected
diagnostic revision before claiming complete physical-sink validation.

The no-sink stellar audit 1279517 passed. On revision `051b4e1`, particle/
cooling regression 1279529 and independent audit 1279530 passed; legacy
stellar/AGN/cooling regression 1279531 and audit 1279532 also passed. These
do not replace the physical-sink audit.

## Completed paired validation

Grammar build 1279535 completed in 3m30s. AGN-on/off runs 1279536/1279537
completed in 1m16s/1m13s, using 960 physical leaves on levels 3 and 4,
one seeded BH, 2109 sampling proxies and four common-time steps. Both form
599 stellar particles in the synthetic depletion stress. Each execution
writes exactly two snapshots; all original outputs remain on GPFS.

Audit 1279538 identified a limitation of the initial thermal-only verifier:
the actual AGN sequence is THERMAL, JET, JET, THERMAL. Audit 1279553 expands
the independent accounting to the actual existing jet model, without
changing or repeating the physical execution. It reconstructs the containing
donor from the endpoint gravity state, SSP returns and accretion receipts;
checks the prescribed loading cap, transported bulk momentum/kinetic energy
and material scalars; and compares receiver-minus-donor transactions with
whole-gas inventories. Existing loading leaves internal energy in the donor.
The net energy oracle uses the declared retained-mass convention,
`eAGN_(T or K) * 0.1 * retained_mass * (3e10/velocity_unit)^2`.
The radiated-rest-mass receipt separately uses `c=2.99792458e10 cm/s`.

Initial and first-drift field/particle/BH exports are byte-identical between
the controls. First gross/retained accretion masses match exactly. At
`t=0.002`, measured feedback energy is `0.0017568113658726` code units with
AGN enabled and exactly zero with it disabled. Later timesteps differ:
AGN-on uses `[0.002, 0.001561868348367044, 0.002, 0.002]`, ending at
`0.007561868348367045`; AGN-off ends at `0.008`. No unmatched-time endpoint
subtraction is used as a causal feedback estimate.

The two jet events load `1.2478212059814934e-10` and
`1.664217240023054e-9` code mass. Their maximum material-balance errors are
`6.391e-23` and `5.534e-25` code mass. All four measured feedback energies
match the corresponding declared source formula within the audit tolerance.
The gas+stars+physical-BH mass decreases by the radiative mass-equivalent
channel, as intended; after accounting for that channel, maximum mass errors
are `6.662e-16` (on) and `4.441e-16` (off). Physical-BH and stellar KDK
residuals are zero in the exported transactions. Gas gravity transaction
errors remain below `1.57e-16`. Each case verifies 3840 cooling/heating
receipts, 1797 stellar-return events and 19767 element receipts. Nine source
epochs per case independently satisfy the native N/Q operators below
`1e-9` relative residual, and both snapshots match the gas+star+physical-BH
source and force fields. Proxies are checked against the physical body and
are not counted a second time in gravity.

Evidence is retained in `results/rar_native_sink_1279535/`, including exact
inputs, executable/source hashes, paired audit, regression evidence, Slurm
accounting and earlier failure excerpts. Full run outputs remain under
`/gpfs/kjhan/RAR_QUMOND/runs/rar_native_sink_1279536` and `1279537`.
No sink formation/merger, deferred-energy replay, resolved jet morphology,
MPI/regrid/restart, cosmological IC, galaxy-RAR or calibrated accretion
result follows from these tests. Global gravitational energy conservation
and exact AMR self-force cancellation are not asserted. These limitations
do not erase the measured native coupling result; they define subsequent work.
