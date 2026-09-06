# AGN accepted fuel and overlapping deposition — bundle plan

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Baseline: `8faa15821a2bf0526241b403d8d4204699662f40`
Status: Fable CONDITIONAL APPROVE on original draft; driver amendments below
recorded after source verification. Operator approved the amended package with
"승인" on 2026-09-06. Native implementation/tests complete; Opus CONDITIONAL PASS;
bounded repairs recorded in evidence, focused Opus closure PASS.
Evidence: `agn_accepted_fuel_overlap_bundle_evidence_2026-09-06.md`.
Review: `fable_agn_accepted_fuel_overlap_plan_audit_2026-09-05.md`.

## Objective and scope decision

Final goal: production/publication-ready native RT, stellar/AGN feedback and
dust. This package must change the actual Fortran source and deposition paths,
not produce another validation framework. It addresses two demonstrated AGN
energy errors together: photons funded by gas that was not actually accreted,
and mass-weighted heat whose normalization changes during overlapping jets.

The preceding discussion called this a common radiative/mechanical budget.
Source inspection shows that a fully closed simultaneous model additionally
requires a physical partition decision. This plan deliberately does NOT claim
that closure. It establishes trustworthy accretion receipts consumed by RT and
closes the existing mechanical deposition budget, while retaining exclusivity.
The operator must approve this bounded interpretation before implementation.

Do not silently subtract legacy heating from the existing escaped bolometric
SED fractions, renormalize omitted groups, cap a MAD efficiency, disable MAD,
or modify BH mass growth to invent a new wind efficiency. Those are physical
model changes, not bookkeeping corrections.

## Source evidence

- `sink_particle.kjhan.f90:accrete_bondi` clips requested gas withdrawal by
  the donor floor, multiplies the accepted amount by `(1-epsilon_r)` for BH
  retention, then removes the gross accepted mass from gas. Existing
  `dMBH_coarse`/`dMEd_coarse` sum *unclipped rate estimates*. Individual element
  densities are not removed with accreted gas, although total metal is.
- `grow_bondi` reduces actual retained mass and subsequently calls
  `kjhan_growspin`. The efficiency used to consume a later cumulative ledger
  need not equal that used at the original withdrawal.
- `snrt_ramses_driver` emits from the difference of cumulative
  `min(dMBH_coarse,dMEd_coarse)` and applies the current efficiency. A loose
  retained-mass upper bound does not make this an accepted-fuel receipt.
- `average_AGN` computes mass normalization before donor loading;
  `AGN_blast` multiplies it by receiver density modified by loading and earlier
  jets. The same mismatch affects thermal events and saved-energy replay.
- `AGN_blast` uses retained mass for its legacy thermal/jet recipes and
  approximate `c=3e10`; SNRT uses the exact cgs speed of light. Keep these
  conventions explicitly distinct until a physical common model is approved.
- `kjhan_growspin` DOES evolve MAD spin, including a low-state polynomial.
  This is not evidence of an explicit energy reservoir debited by each jet.
  Do not describe the code as lacking spin evolution or assert closed spin
  energetics from the presence of the spin update alone.

Scientific context: spin-powered jets and radiative accretion efficiencies are
distinct in the NewHorizon model; adding them or forcing a single cap is not a
neutral correction. See Beckmann et al., MNRAS 536 (2025), 1838,
https://doi.org/10.1093/mnras/stae2595 . This motivates keeping the physical
partition open, not claiming this bundle validates that literature model.

## One native implementation package

### A. Accepted accretion receipts, used by the real SNRT source

Capture the actual gross withdrawal, retained mass and efficiency inside the
existing successful accretion event. For each accepted event define

```
dm_retained = (1 - epsilon_event) * dm_gross
dE_radiative = epsilon_event * dm_gross * scale_m_event * c_exact**2
```

Use the contemporaneous efficiency already applied to retention (including
the current low-state suppression), not an efficiency recomputed after spin
evolution. Accumulate energy directly rather than subtracting nearly equal
large cumulative mass totals. Validate finite mass, units and efficiency before
mutation. Accept `0<=epsilon_event<1`; zero efficiency or zero accepted mass
creates zero energy. Publish receipts only for accepted events; threaded
reductions must retain actual event amounts. Do not call `clean_stop` inside
an OpenMP loop or cell lock: propagate an error through existing thread/reduction
boundaries and use collective consensus where applicable before stopping.

Use ONE per-sink pending radiative-energy reservoir, carried with the sink
arrays, not cumulative totals plus a second consumed-energy cursor. Initialize
it at zero in `init_sink` and creation; add accepted energy in `grow_bondi`
through the existing event/thread reduction. In `merge_sink`, sum old pending
energies into `igrp` alongside `dMsmbh_new` and `Esave_new`, then install the
new array with the existing sink-array commit. This makes reorder/merger
ownership explicit without a separate identity-remapping service. Preserve
existing valid unique `idsink` checks. No energy is created by seed formation
or BH-BH mass combination itself.

Keep the reservoir in physical energy units, converted using the accretion
event's units before accumulation, so later emission does not reinterpret old
fuel using a different unit scale. A single extra thread-reduction component
is sufficient for the accepted radiative energy. Accumulate it only when the
live SNRT path is requested; standalone legacy runs retain zero unused pending
energy and their existing mechanical recipe. The current serial/fresh-start
restriction remains; no new restart records or distributed qualification.

Give the native photon budget an energy-input path and wire the actual driver
to pending accepted radiative energy. Retain the old mass-based entry only as
a delegating compatibility interface where still used. Apply the current nine
escaped bolometric fractions exactly once, relative to that energy; do not add
a historical 0.5 factor. The unrepresented fraction remains unrepresented,
not invented photons or mechanical energy.

Clear the emitted pending energy only after successful all-group source
deposition. In the current source call there is no concurrent accretion; do
not turn this into a new global transaction. A source failure leaves it pending.
RT/chemistry transport failure
after successful source injection must not cause photon re-emission: source
commit and transport rollback are separate existing boundaries. Delete the old
`accounted_ids/accounted_inflow/retained_seen/retained_initialized` source
accounting and its supplied-regression/retained-bound skips. Do not retain
current-efficiency resolution as another gate on already accepted event energy.
Keep coarse supplied-rate ledgers for existing legacy mode selection and
diagnostics; they no longer fund SNRT photons. No generic multi-process
transaction service or persistent event database.

### B. Carry composition with actual swallowed gas

At the same accretion withdrawal, scale total metal and declared element-density
fields by the actual removed fraction, reusing the established scalar-layout
selection/validation from the entrainment bundle. Factor that layout selection
into one small shared routine: validate once in `grow_bondi` before accretion,
as well as in `AGN_feedback`. The latter is not called in live SNRT mode.
Scale selected densities by `1-dm_gross/(rho_before*cell_volume)` inside the
existing accretion lock. Preserve disjoint auxiliary fields under their existing
convention; reject overlapping maps before mutation.
This must support both legacy and channel-resolved storage without initializing
a yield source or treating `active_element` as a gas-transport mask.

Accretion removes proportional specific total energy under the existing model;
do NOT transplant the cold-entrainment rule that leaves donor internal-energy
density unchanged. Do not call the entrainment withdrawal helper (its 25% cap
is also inapplicable). No new NENER refusal in accretion: its existing auxiliary
energy convention is not being qualified here. No accretion/drag recalibration,
new species, or dust model.

### C. Freeze thermal and saved-energy weights across all overlapping events

Use the pre-loading density measure already used by `average_AGN` for both
denominator and numerator. Preserve existing spherical geometry, fallback,
temperature cap, and the separate volume-weighted jet-thermal convention.

Reuse the existing owner-local `ind_blast` and globally reduced actual `mAGN`;
do not add donor-payload storage or another normalization pass. At each receiver
cell, before entering its deposition loop, set `rho_ref=rho_current` and add
`mAGN/cell_volume` for every loaded donor entry whose `ind_blast` equals that
cell, using the existing neighboring-bin candidate list. Include all shared
donors and fallback-loaded events, not just the current thermal event. Freeze
that scalar throughout the cell's event loop. The numerator and old denominator
then refer to the same pre-loading gas despite preceding in-cell jet deposits.
Source inspection confirmed that between normalization and deposition the only
density mutation is the sequential loading; fallback deposition follows the
main cell traversal. No full-grid copy or cells-by-AGN matrix.

For each thermal/replay event, requested cell energies sum to its `EAGN` and
accepted gas energy plus newly deferred energy equals the requested amount.
Existing unresolved donor fallback returns the whole event exactly once.
Saved `Esave` is a remaining mechanical liability, never new fuel or photons.
Changing event order may alter capped placement/deferred amounts; do not claim
identical final thermodynamics. Each order must conserve its own event budget.

## Proportional verification and completion

Extend `snrt_agn_source_smoke`, `snrt_agn_efficiency_smoke` where applicable,
and `agn_feedback_deposition_smoke`; reuse `run_fp15_agn_efficiency.sh`.
A few combined cases cover:

1. Floor-limited/zero gas withdrawal and two events with different efficiencies:
   photons follow the sum of actual event energies, not supplied rates or the
   latest efficiency. Coarse reset, sink reorder/creation/merger, source retry and
   successful source followed by transport retry do not duplicate energy.
2. Accretion removes the selected gas constituent masses in proportion to the
   actual withdrawal, preserving donor fractions and existing auxiliary-field
   semantics; invalid layouts/receipts reject before mutation.
3. Multiple overlapping thermal, jet and replay events, including a shared
   donor, unequal cell volumes and an active cap: requested heat integrates to
   the source budget and gas-energy change plus remaining liability closes.
   Reverse source order and test closure, not identical capped gas states.

Reuse production helpers in these native tests. Inspect and compile the actual
callers in scratch with the supported SNRT/CPU configurations. Do not substitute
source-string assertions for behavioural checks. No new Python report pipeline,
AMR reader, large HDF5 asset, global infrastructure gate or long simulation.
Helper tests and caller compiles do not qualify live evolution or publication.

Completion is real driver consumption of accepted receipts, corrected gas
composition withdrawal and normalized overlapping native deposition, backed by
these bounded tests, one evidence record, and one Opus 5 end review (Fable only
if the primary fails). Fable reviews this plan for necessity, feasibility and
overinstrumentation. Findings outside this package go to the science roadmap,
not automatically to another blocking gate.

## Explicit follow-on, not secretly included

After this package, design and obtain approval for the physical simultaneous
RT/mechanical model: gross versus retained accretion convention, unresolved
wind/heat versus resolved radiation coupling, MAD spin energy, contemporaneous
mode history and the luminosity basis of the SED fractions. Then qualify the
joint mode and relax exclusivity only with evidence. A receipt shared in name
alone is not a closed physical energy partition.

No runtime activation, restart/migration certification, SED/obscuration/IR
approval, AGB-source reopening, new job, commit or push is authorized by this
planning record. Unrelated dirty and staged files are preserved.

## Driver disposition of Fable conditions

1. Adopt the simpler pending reservoir and explicit existing merger mapping.
   Independently confirmed `merge_sink` forms `igrp` and installs replacement
   arrays. The original plan's new independent energy ledger would not
   necessarily inherit the coarse-rate-array bug cited by Fable; the reason to
   prefer this amendment is smaller ownership logic and no cumulative-energy
   subtraction. Keep event-unit energy rather than Fable's suggested deferred
   code-mass conversion. Do not accumulate an unused legacy reservoir merely
   to reset it later. No other coarse-ledger repair is claimed.
2. Accept shared layout validation in BOTH callers and proportional selected
   scalar withdrawal; no transplant of cold loading or its cap. Keep existing
   auxiliary energy semantics without adding a new accretion NENER gate.
3. Accept reuse of `ind_blast/mAGN`; the needed donor information and ordering
   already exist. Drop the new sparse map and alternate snapshot-pass design.
4. Accept energy-input photon API, old-wrapper delegation, and removal of the
   post-event efficiency and cumulative supplied-mass source gates.
5. Keep native behavioural coverage of the new reservoir operations, using
   production helpers if factoring is useful, plus caller compiles/inspection.
   A compile alone cannot demonstrate retry or merger conservation. No new
   initialized-RAMSES harness or Python framework is required by this plan.
6. Do NOT accept the suggestion that source-side checks after gas withdrawal
   replace validation before mutation. A late photon rejection cannot repair
   invalid gas/BH changes. Keep bounded pre-mutation checks with error reduction,
   not `clean_stop` inside the lock. This is an implementation invariant, not
   an additional audit gate. Fable's approximate line-count estimates are not
   an implementation estimate; the actual native caller wiring remains required.

Fable judged the bounded precursor necessary and not a delay of the scientific
goal, conditional on simplification. These are driver dispositions of that
single planning review. The operator subsequently approved implementation;
completed work and the Opus closure are recorded in the evidence. The
simultaneous physical partition remains explicitly unapproved.
