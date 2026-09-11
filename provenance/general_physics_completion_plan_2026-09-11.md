# General physical-model completion: scope and source assessment

Operator request: complete the general physical model after the named
comparison implementations. This is new physics, not permission to relabel
the comparisons or remove their admission checks. Workspace: GPFS LRD_JWST;
origin kjhan0606/LagRamses. Existing dirty implementation is preserved.

## Active operator decision: galaxy-scale effective model

Operator subsequently directs sequential execution of three steps: (1) mixed
v5 + DTD/N100, (2) compatible effective-model configuration, (3) coupled run
and closeout. [Galaxy calibration/resolution science](galaxy_calibration_science_plan.md)
is explicitly separate. No calibration campaign is launched by this approval.

Execution result: [three-step evidence and exact limits](effective_population_execution_2026-09-11.md).
Native v5/DTD connection and selected short coupled/restart profile pass.
The failed~69Myr CHIMES macrostep is now separately
[reproduced, repaired and verified](chimes_long_interval_repair_2026-09-11.md)
in the unchanged four-step coupled run, not deferred as microscopic physics.
No universal production
qualification or calibration follows from the bounded implementation result.

The operator has now accepted the recommended completion standard
(2026-09-11): literature-based approximations are admissible with explicit
domains, consistent mass/energy/element accounting, exclusive source/remnant
ownership and working RT/feedback/dust coupling. Direct microscopic
calculations move to the [long-term backlog](microscopic_physics_long_term_backlog.md).
This supersedes the proposed six-program scope and the clarification-wait
status below, retained as review history. No further approval of this same
scope is needed. It does not certify current code as complete or approve
every combination of existing comparison selectors.

Completion requires the selected effective model's actual source-to-consumer
connection, no double mass/energy/element release, and proportionate native
and coupled execution/restart evidence. Existing evaluated evidence is reused;
no per-exclusion audit or new monitoring framework. Domain limitations remain
explicit, including unsupported radiation, temperature, age, mass and Z.

Next native target remains mixed v5 + empirical DTD/N100, reconciling the
population identity and chosen mass-accounting contract together. Do not
merely delete v5 rejection checks or silently switch strict-WD/effective-SSP
semantics. Same-author stellar libraries and resolved binary evolution are
not prerequisites for a consistent effective population.

## Historical proposed scope (superseded as mandatory work)

Replace the six exclusions reported in
`remaining_physics_bundle_execution_2026-09-11.md` with source-supported
closures and native consumers where data permit. Preserve current defaults
and comparison selectors. Do not reopen already evaluated AMR/layout tests
or silently add nonideal MHD, dust--B or anisotropic CR--B projects.

“General” must mean an explicitly defined model domain, not exact physics
for arbitrary grain composition, stellar mass, metallicity and radiation.
Do not claim implementation complete while any requested consumer or its
necessary physical input is absent. No additional monitoring framework.

1. **PAH:** replace unit-yield atomization by size/charge/H-dependent daughter
   transitions with formation energies, optical/IR response and electron
   energy partition. Existing C24 neutral/cation tables cannot stand in for
   all daughters. Huo et al.2023 measured coronene-cation fragmentation near
   the carbon edge (279--300 eV), not a full 0--10 keV state network. Its data
   statement requires contacting the author for generated datasets. Use
   public supplementary data if sufficient; do not send mail automatically.
   <https://doi.org/10.1140/epjd/s10053-023-00763-w>
2. **Metallic Fe:** finite-time charge and hard-photon electron/fluorescence
   partition must share one conserved energy/charge transaction with gas and
   grain enthalpy. WDB2006 provides carbonaceous/silicate estimates, not an
   iron-specific keV calibration. Atomic cross sections alone do not specify
   solid-state escape probabilities. Retain UV and thermal-limit models until
   an explicit iron closure/data set is selected.
   <https://arxiv.org/abs/astro-ph/0601296>
3. **Stellar population:** obtain same-family evolution/postprocessing
   histories for a supported mass/Z subset before replacing mixed endpoints.
   NuGrid set1ext has useful published histories; it is not a continuous
   .08--600 Msun, Z=0-inclusive matched grid. Atmosphere matching must provide
   actual spectral energy integrals, not uniquely infer nine groups from
   five photon integrals. No synthesized missing evolutionary tracks.
   <https://nugrid.github.io/content/data>
   <https://download.nugridstars.org/set1ext/>
4. **SNIa:** evolve donor/WD transfer, retained mass, outflow and ignition
   consistently with binary state and event yields. Saved frozen-contact
   tracks are insufficient for general binary ownership. Preserve empirical
   DTD/N100 and explicitly frozen comparison until a concrete evolutionary
   implementation and channel-matched ejecta are available.
5. **Wind/pulses:** use actual time-resolved mass/composition/radiation and
   pulse event lists. Current PARSEC steady-wind histories are useful, but
   speed postprocessing cannot invent LBV mass loss. Woosley2017 contains
   pulse examples, not a matched full-grid per-pulse element history.
   <https://arxiv.org/abs/1608.08939>
6. **Radioactive energy:** nuclear branching/charged-particle/gamma spectra,
   source-age convolution and a spatial deposition/escape model are needed.
   Do not add total decay Q directly as gas heat or place MeV photons into
   the current <=10 keV SNRT groups. Local full deposition would be another
   comparison limit, not removal of the general-model limitation.
   <https://arxiv.org/abs/1309.0018>

## Execution decision under review

First determine which public data can actually remove an exclusion. Native
implementation follows only for that supported domain, with a coupled
source-to-consumer check and one relevant live/restart test. Missing physical
inputs are reported distinctly from coding work; no disconnected generic
table loader is counted as completion. A same-family NuGrid subset is a
candidate first substantive replacement, not yet an approved scientific
selection or an acquired/evaluated package.

This document records a proposal and checked limitations, NOT implemented
general physics. No new simulation has been launched for this expansion.

## Independent source checks and driver correction

Read-only local checks by Gibbs and Popper confirm useful source inputs but
not a universal replacement package. NuGrid HTTP directory access succeeds
(2026-09-11); `set1ext/set1.2/see_wind/` lists solar-Z models from 1 to25
Msun with gaps. No large archive was downloaded and no library was replaced.
Mixed published libraries are not intrinsically invalid: common IMF weights,
initial composition conventions, time support, return budgets and exclusive
remnant ownership are the requirements. A single-family replacement is an
optional scientific choice, not an automatic galaxy-production prerequisite.

The saved COSMIC seed169 remains detached at100Myr (R/RL=.280739).
Its forced contact therefore cannot be reinterpreted as simulated overflow.
An actual overflow follow-up ends with a .714100Msun WD rather than Ia.
Resolved donor+WD MESA input files have a published release:
<https://arxiv.org/html/1901.04512v2>, DOI10.5281/zenodo.2630887.
The paper distinguishes central/off-centre ignition but does not uniquely
identify the explosion subclass/yields. Adopting it would require new binary
evolution calculations and a declared mapping, not a native ledger patch.

Cleeves et al.2013 provides energy/range inputs and spatial escape equations,
but uses disk geometry and grey stopping approximations. Its old Fe60
half-life1.5Myr must not overwrite the pinned2.62Myr value. General galactic
charged-particle deposition is not established by adopting its disk formula.

## Fable review and driver disposition

Completed read-only review:
`general_physics_fable_2026-09-11.json` (successful result;
reported substantive model `claude-fable-5-1`). Verdict CONDITIONAL.
Q-GOAL/Q-LEAN: deriving six mandatory research programs from comparison
exclusions overexpands the galaxy-simulation objective. The driver accepts
that criticism: the list above is an assessment, not six new mandatory gates.

Candidate next native connection: mixed v5 + empirical DTD/N100. Inspected
`stellar_yield_audit.f90` and `snrt_parsec_source.f90`: both explicitly reject
SNIa with v5. However, the review's "code only" description is incomplete:
v5 binds single-population ID0/binary fraction0, while
`stellar_snia_population_contract.f90` requires binary SSP identity and
consistent fraction. Removing two checks alone would leave a false population
binding. Existing `strict_wd` and `effective_ssp` accounting paths both exist;
neither may be changed into the other silently. Strict-WD additionally needs
causal CO-WD supply at each event interval, not only nonnegative final mass.

Other review claims NOT adopted as established facts:
- A <1% Mg/Fe abundance effect does not bound radioactive ionization/heating.
- Declaring the Fe thermal upper limit does not make it a calibrated hard
  photoelectric model.
- The presence of other public PARSEC low-mass releases does not establish
  same-version tracks and radiation coverage for the selected11-Z package.
- Pulse-integrated shock energy does not prove resolved pulse timing is
  irrelevant at every intended resolution.

No physics restrictions, default model or population identities were changed.
No code implementation or production qualification is claimed for this new
request. A user clarification has been requested: galaxy-scale effective
closure completion versus explicit microscopic-model expansion. It materially
changes the required native consumers, data and scope. Recommended scope is
the former, retaining physically justified, declared subgrid approximations.
