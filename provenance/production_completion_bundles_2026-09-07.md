# Remaining production implementation bundles

## Latest operator scope: comparison implementation closed

The operator approved a bounded closeout after explicitly objecting to the
repeated expansion of "final" work. The **selected single-rank/OpenMP comparison
execution implementation is closed**; unrestricted production/publication
qualification is not claimed. This decision supersedes the active NEXT-WORK
instructions and remaining-bundle counts below, which are retained as history.
The newly proposed spectrum-specific transport/mixture work is deferred, not
another mandatory completion bundle. No new simulation or test ladder is needed
for the unchanged, already exercised profile.

The fixed executable, inputs, local-only dependencies, launch recipe, evidence,
limitations and separate optional verification proposal are in the
[closeout handover](rt_feedback_dust_comparison_closeout_2026-09-07.md).
Do not automatically resume historical deferred work from this file.

## Historical implementation plan and amendments

User-directed sequence, updated 2026-09-07. Scope: RT, stellar/AGN feedback,
and dust in `/gpfs/kjhan/LRD_JWST` (kjhan0606/LagRamses).

1. **Real-source integration (in progress):** load actual stellar/BH particles,
   connect accepted accretion and stellar population evolution to feedback,
   primary RT and live dust; preserve source progress and unconsumed energy
   across restart without loss or duplication. Do not require spontaneous BH
   formation as a prerequisite for an existing-BH integration run.
2. **Physical input completion:** resolve the selected yield, stellar/AGN SED
   and dust material inputs and their runtime connections. Keep reference
   controls separate from scientifically approved production inputs.
3. **GPU/OpenMP runtime automatic allocation (user-added):** provide runtime
   backend selection and resource assignment for these high-level operators.
   Account for compiled capabilities, visible devices, MPI local-rank device
   ownership, memory headroom and work size. Preserve explicit CPU/GPU choices.
   Never silently fall back to a nonexistent or physically different CPU
   implementation; define state synchronization and safe fallback boundaries
   before enabling automatic switching. Update the namelist generators with
   any new namelist controls. Reuse existing backend infrastructure.
4. **Production execution closeout (formerly #3):** finish the selected profile's
   MPI/AMR/restart wiring, retain short operational checks, and fix the runnable
   configuration and build identity. This is not a promise to support every
   backend/model combination.

The user inserted #3 between the original #2 and #3; the remaining planned
count is therefore **four**, including the first bundle already in progress.
Do not create new bundles for each defect, test, or audit finding.

Only necessary build, short execution and focused regression checks accompany
implementation. **Propose a separate comprehensive verification plan after
all four implementation bundles are finished.** Do not add intermediate
audit gates. Production readiness remains subject to the eventual results;
implementation completion alone is not a scientific validation verdict.

## Execution status after the preapproval to work through bundle 4

The four bundles were worked on together where their dependencies overlap;
they have **not all reached production completion**. No new numbered bundle
or intermediate audit has been introduced.

| Bundle | Implemented / exercised | Remaining boundary |
| --- | --- | --- |
| 1 | Actual BH IC loading and accretion; actual STAR-particle photon integration; primary RT/live dust coupling; five AGN pending reservoirs and exact stellar table persisted in HDF5 | Full physical stellar mechanical feedback combined with the selected yield/population package is not certified by the reference run |
| 2 | Native stellar table adapter binds IMF mass limits, population, age/Z bounds and common spectral identity; existing physical admission checks remain enforced | Selected channel-resolved fate policy is still `review_only_unresolved`; no joint approved Chabrier stellar/AGN SED and live-dust material package is installed |
| 3 | Primary SNRT species/dust CUDA/OpenMP auto placement, MPI-local device assignment/sharer counts, memory/work-size selection and host synchronization | Feedback mechanical and thermal/IR operators remain their existing host implementations; toolkit-free CPU build and unimplemented CUDA counterparts are not claimed |
| 4 | NVAR=30 build; corrected blocked-layout and half-open source ownership; bounded primary/IR MPI packets; single-rank real-star/AGN/dust execution, GPU-to-OpenMP HDF5 restart, level-3/4 AMR control; two-rank restart completes after the bounded near-bath correction | Combined physical MPI+stellar+AMR qualification is not claimed; physical production remains contingent on bundle 2 |

Detailed run results, failures and final binary identity are recorded in
`real_source_integration_progress_2026-09-07.md`. Do not interpret a reference
control's `Run completed` as approval of a scientific source model. The existing
physical-package admission contract has `physical_node_inventory=[]` and
`physical_package_selected=false`; fabricating metadata is not completion.
The previously parked 40--120 Msun source choice stays explicitly unresolved,
not silently set to zero or removed from the approved IMF support.

## User-directed deferral and continued implementation

### High-mass model selection: operator approval, 2026-09-07

The operator accepted a consistent baseline plus limited, explicit namelist
subparameters instead of indefinitely blocking on competing source models.
This supersedes parking the **implementation** of the high-mass options;
it does not certify a source correction or permit broken conservation.

The first implementation adds `high_mass_preset` and
`high_mass_remnant_adjust_max_fraction` in `STELLAR_ENRICHMENT_PARAMS`, shared
with the namelist generator and mkrun CLI/GUI advanced editor. Omission means
`source_consistent` and zero correction allowance. The native endpoint
resolver in `stellar_yield_audit.f90` implements:

- `source_consistent`: same-source input required; keep the endpoint unchanged
  and reject a non-closing budget.
- `wind_only_collapse`: retain supplied wind, suppress terminal ejecta and
  terminal energy/momentum, put remaining baryonic mass in the remnant. This
  is an explicit comparison approximation, not a universal fate statement.
- `mixed_remnant`: keep both supplied ejecta channels and compositions; derive
  remaining baryonic mass and report its signed difference from the supplied
  remnant. The user must set a positive correction limit as a fraction of
  initial mass. No default correction tolerance has been invented.

All modes enforce finite/nonnegative gross material, element budget and mass
closure; errors publish no partial result. Raw data remain unchanged. The
endpoint helper consumes complete endpoints, NOT cumulative intermediate-age
rows. Its source-energy convention excludes directed bulk kinetic energy,
which the existing bridge adds separately per channel. Isotope-decay handling,
source identity and timing are not inferred by this helper. The source-match
input must ultimately come from the source adapter, not a user trust switch.

Native endpoint tests include an actual LC18 wind **total** at 60 Msun,
[Fe/H]=-3, 150 km/s: 18.06491167124404 Msun gives 41.93508832875596 Msun
remaining under the chosen wind-only collapse approximation. This is solely
a mass-budget check, not a physical elemental/energy/time realization.
Synthetic tests cover consistent, mixed, excessive correction, negative
remaining mass, and invalid elemental budgets. The existing namelist tests
cover selection propagation and atomic rejection. The GUI suite passes
21 tests with one display skip. Historical hydro output comparison now
exempts only the intentionally completed stellar namelist block; all other
outputs still match. The old mkrun block omitted required population fields.

Both gfortran and ifx bounds-checked endpoint executables pass in
`.high-mass-options.50u91r`; the expanded native configuration fixture passes.
The complete NVAR=30 SNRT/DUST_LIVE/HDF5/CUDA executable linked successfully:
`.production-native.YCek9O/ramses_high_mass_options3d`, SHA256
`621be64229c9bf6502b523c8636fd075da257024f1c1c52f3909c8d506c5b23f`.
Build log: `.production-native.YCek9O/build-high-mass-options.log`.
No simulation, external audit, commit or push was performed in this increment.
The existing unrelated generator deletions and staged `informat` addition
remain separate and preserved.

**Subsequent increment — native connection implemented:** the user requested
the source-node/time-history adapter and restart binding. These now operate
through the explicit `user_selected_model_v1` route, not by flipping the old
`snii_source_node_fate_consumer_available` approval flag. A validated per-table
history controls wind timing, discrete terminal events, common IMF integration
and existing particle-to-gas deposition. HDF5 persists the actual source values
and selection; MPI consensus and changed-input restart rejection are active.
Single-rank, MPI2 and same-input restart evolution completed on synthetic input;
all three presets passed the native cumulative/increment fixture. See
`real_source_integration_progress_2026-09-07.md` for paths and measured results.
The initial route was single-star/wind+SNII with optional AGB. The subsequent
operator-requested increment adds a matching effective binary SSP/SNIa/AGB-WD
route and optional linear-in-Z cumulative source mixtures within supplied Z
bounds. Real approved N100/Maoz SNIa now executes in the combined control;
DTD and interpolation semantics are bound to restart. PISN stays disabled.
This completes those engineering connections, not physical qualification of
LC18/BR26/AGB, Z extrapolation, or an all-channel physical production package.
No new audit bundle was created.

Latest amendment: the operator selected primary article/detailed yield-table
data over the discrepant AGB auxiliary array. The KL16 4 Msun, Z=.03 remnant
now selects .774 Msun (3.226 Msun expelled), retaining raw .744 as provenance;
this particular mass choice no longer blocks implementation. The printed
article contains only the 3.5 Msun example, so the .774 source is explicitly
the detailed yield header, not a falsely attributed article value. Other AGB
source-data issues and the 40--120 Msun issue remain deferred. See the updated
F-P1H-F plan for the precise decision and scope.

The user explicitly re-parked **both** the 40--120 Msun fate/yield discrepancy
and the KL16/CK22 AGB raw-data issues. They are follow-up scientific input work,
not immediate prerequisites for unrelated RT/feedback/dust implementation.
Do not reopen their source/literature audits, contact authors, fabricate missing
values, or broaden the approved physical domain while continuing this work.
Production claims involving these unapproved inputs remain restricted.

The next increment completed the native SNIa population-to-runtime binding
and a Draine optical/thermal table exporter to the existing live-dust input.
The already approved Kroupa/binary SNIa baseline is retained; it is not silently
applied to default Chabrier/single-star particles. The Draine optical input is
real, but its constant heat capacity remains an explicitly reference-only
control. A short native cooling run and HDF5 output succeeded; see the existing
progress record. No additional bundle or intermediate audit gate was added.

Continue the remaining stellar/AGN spectral-input and dust-material runtime
connections within bundle 2, using explicit model identities and the existing
native operators. Resolve scientific choices when evidence supports them;
do not relabel the reference runs as physical production or let the two parked
data problems halt independent implementation.

**Latest operator-approved LC18 increment:** actual table8 wind isotope ejecta
and table7 lifetimes now feed the native 40--120 Msun wind-only collapse path.
Missing time/energy/composition detail is represented by REQUIRED, explicitly
selected uniform-release/mean-composition/isotropic-thermalization assumptions
and a user-supplied wind speed, not by invented source data. Native endpoint,
SSP, live hydro and restart checks passed; changing speed rejects restart.
This completes the local comparison-model connection, not phase-resolved
wind/terminal physics or unrestricted publication approval. The former statement
that LC18 is wholly unconnected is superseded for this selected route only.
Actual AGB input and the corresponding WD supplier remain the next physical
input implementation. Evidence is appended to the existing progress record;
no additional bundle or intermediate external audit was created.

**Superseding AGB increment (operator approved items 1--3):** normalized KL16
gross ejecta now connect to exact M/Z/Y/overshoot lifetimes and native
simultaneous terminal-envelope/WD release. Explicitly selected AGB wind speed
drives native thermal energy. Actual LC18+AGB CPU hydro, mass conservation and
same-input restart pass. No extra audit bundle was introduced.
The SNIa-on combined route was also run, but correctly rejects: approved DTD
onset 40 Myr precedes the selected ordinary CO-AGB WD supplier's first event
at 68.89 Myr. Thus the connection is implemented, but combined physical
qualification is NOT complete. A consistent early binary/WD source or a
separately approved comparison DTD/population is the next scientific decision;
do not silently change the approved DTD or borrow future WDs. See the same
progress record for source hashes, tests, run paths and measured results.

The subsequent approved prompt-WD review found no ready, consistently
normalized early supplier for unchanged KL16 + 40 Myr DTD + N100. The
observational prompt bin is <420 Myr; the exact cutoff is a project model
choice, not independently resolved by that observation. Candidate binary
channels also require donor-mass transfer and avoidance of double-counted
KL16 returns. Prefer a bounded offline binary-history feasibility sample
before new native plumbing, and only after selecting its physical population
prescriptions. Do not launch a large BPS campaign or invent an early reservoir
to force a pass. Details and quantitative missing-fuel estimates are appended
to the same progress record; runtime parameters and approval IDs are unchanged.

**Offline feasibility executed:** operator approved the bounded comparison.
COSMIC 4.2.0 evolved 1620 binaries over four Z values plus four selected
follow-ups. Early CO WDs form at 21--28 Myr, and donor transfer is resolved.
However the selected fast explosions follow COSMIC's sub-Chandrasekhar
He-trigger prescription, not the approved N100 source. Thus the previous
"no executed early supplier" status is superseded, but there is still no
qualified prompt N100 supplier to splice into native feedback. Keep the
existing physics unchanged; select a compatible near-Chandrasekhar population
or separately authorize a sub-Chandrasekhar comparison before further source
plumbing. Small data, model choices and measured outcomes are recorded in
the existing progress log; no new audit bundle was opened.

**N100-preserving retention comparison:** selected/implemented an offline
Wang-2017 stable-burning + conditional Kato-Hachisu-2004 flash reference and
screened the saved fine histories. Stable, unsupported-flash and expansion
regimes are distinguished without inventing retention outside source support
or automatically turning off-centre ignition into N100. This is not coupled
binary re-evolution. A paired same-input COSMIC check also showed that the
representative early CO WD becomes a neutron star with the IFMR modification
disabled: the prior early-WD result is prescription-dependent. Before any
native splice, a declared remnant model and self-consistent donor/transfer/
retention/orbit history are still required. Do not substitute saved net mass
change rates for incident He transfer or continue old histories beyond their
original explosion. Native DTD/N100/source defaults remain unchanged.

## Current operator decision: effective SSP SNIa, not a BPS prerequisite

2026-09-07: the operator approved the recommendation with
"추천한 대로 정리합니다." This supersedes the preceding NEXT-WORK instructions
to complete a coupled near-Chandrasekhar binary population before resuming
high-level RT/stellar+AGN feedback/dust. It does not invalidate the measured
WD deficit or promote the COSMIC comparison into a physical N100 population.

**Status at the scope decision: model direction approved, implementation pending.**
Superseded by the completed implementation note below. At that decision the
native route still required causally supplied AGB WDs and still
rejects the selected 40--68.89 Myr deficit. No runtime default, source contract,
namelist, binary or existing output was changed by this documentation decision.

The next bounded implementation, within the existing source-integration work:

1. Add an explicitly selected effective-SSP accounting mode. Keep the approved
   empirical DTD and N100 event mass/yields/energy, including their existing
   IMF/population applicability and rate normalization. Debit ejecta from the
   particle's remaining stellar mass, not from a claimed KL16 WD inventory.
   Bind the new accounting identity to input validation, startup reporting and
   restart; update the namelist generators together if a main namelist field
   changes. Do not reuse the old physical approval as proof of the new mode.
2. Account for AGB, SNII/winds and SNIa returns together before committing a
   particle/gas update. Prevent duplicate returns and aggregate overdraw; do
   not add a fictitious binary reservoir on top of the full KL16 population,
   borrow future WDs, silently clip the DTD, or disguise a rejected transaction
   as successful feedback. Retain the strict WD-supply mode and its rejection
   tests as a separate comparison, not a bypassed or newly passed gate.
3. Use the existing native checks for combined mass/element/energy accounting,
   unchanged strict-mode rejection, and same-input restart/mode-change rejection.
   Then resume the existing RT/feedback/dust integration queue. Do not create
   another Python framework, a new gate hierarchy or a large BPS campaign.

This is a phenomenological SSP approximation, not demonstrated microscopic
progenitor consistency. Aggregate mass conservation cannot establish that the
same physical binary progenitors have been partitioned consistently across
single-star yields and SNIa. Document this limit in any production/publication
claim. Existing source-domain and wind/dust approximation limits also remain.

Medium-term, non-blocking work: self-consistent binary population normalization,
IFMR/remnant sensitivity, incident H/He transfer, retention/outflows and orbital
response, near-Chandrasekhar versus sub-Chandrasekhar fate/yields, and their
metallicity dependence. Preserve the offline scripts and measured histories;
do not resume this research without a separate task. No new external audit,
calculation, commit or push is part of this scope-consolidation decision.

**Subsequent implementation authorized and completed (2026-09-07):** the
operator's "진행 승인." authorized the bounded native connection above.
`mass_accounting='effective_ssp'` with its separate accounting approval now
selects the remaining-stellar-mass route in the external SNIa contract. The
default `strict_wd` path, DTD/N100 values, IMF applicability, main namelist and
historical strict restart identity are preserved. Effective identity rejects
mode changes on restart. No fictitious WD inventory is populated.

GNU bounds-checked and optimized Intel native tests pass, including actual
KL16+LC18 combined accounting over 0--13.7 Gyr within common source Z support.
Four-step effective-SNIa CPU hydro and checkpoint-2 restart complete with
conserved aggregate mass; strict-WD rejection and both restart-mode-change
rejections remain demonstrated. Details and limitations are appended to the
existing source integration record. This closes the effective SSP accounting
increment, not microscopic binary-population qualification or all-channel
physical production approval. Resume the existing RT/feedback/dust input and
integration queue; do not reopen the parked BPS work or add another gate ladder.

**Next integration increment:** real KL16/LC18 + effective SNIa now runs with
accepted BH accretion, reference AGN RT and real-opacity/reference-thermal dust.
A canonical-sink coordinate-sync traversal defect exposed on restart was fixed
locally; the continued run reproduces stellar/RT/dust/AGN-reservoir state.
This is a selected serial/OpenMP control, not joint physical source approval.
The oldest star is only 51.36 Myr, so nonzero AGB release remains separately
demonstrated, not combined here. The incompatible stellar SED stays disabled.
The existing integration record retains successful and failed cases and their
scope. Remaining physical input priority is a compatible stellar SED and dust
thermal material data; do not extend microscopic binary studies as a prerequisite.
