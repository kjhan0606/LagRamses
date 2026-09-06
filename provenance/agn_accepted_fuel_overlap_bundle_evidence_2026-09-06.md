# AGN accepted fuel / overlapping deposition — implementation evidence

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`, baseline `8faa158`.
Authority: operator approved the amended 2026-09-05 plan with "승인".
Status: native implementation and bounded verification complete; Opus end review
CONDITIONAL PASS followed by bounded repairs and focused closure PASS.
Review: `opus5_agn_accepted_fuel_overlap_end_audit_2026-09-06.md`.
Closure: `opus5_agn_accepted_fuel_overlap_closure_audit_2026-09-06.md`.
Plan: `agn_accepted_fuel_overlap_bundle_plan_2026-09-05.md`.

## Actual implementation

- `agn_accretion_receipt` replaces the actual `accrete_bondi` inline gross/net
  computation. It clips requested gas by the original 0.75 initial-density
  floor BEFORE retention, then computes `(1-epsilon_event)*gross` and
  `epsilon_event*gross*event_mass_unit*c_exact^2`. Zero efficiency is valid;
  unity/nonfinite efficiency and invalid mass/units are rejected. The gas
  loses this same gross amount, without reconstructing it by division of net.
- `accrete_bondi` records accepted energy only when live SNRT was requested,
  inside the existing successful cloud-event record. The existing OpenMP
  per-sink reduction and MPI packed sum gain one energy component. Invalid
  event inputs report through per-thread errors; the lock is released and no
  invalid event's donor is modified. Collective failure occurs outside the
  lock/parallel loop. This is per-event safety, NOT batch rollback of earlier
  valid accretion events. Zero-rate events retain the original drag path and
  coarse supplied-rate diagnostics; no drag recalibration was made.
- `pm_commons` owns ONE `agn_pending_erg` array in physical erg. The active
  Makefile inherited `init_sink.f90` from cuRamses; a project-local override
  copies that initializer with zero initialization in both normal and
  allocation-only entry points. The inherited shared file is untouched.
  Creation zeroes only genuinely new sink slots. `agn_merge_pending` uses the
  actual `gsink` grouping before compaction, with map/finite/overflow checks,
  and the sink commit installs the merged reservoir. It replaces in-loop
  unchecked addition, not a separate identity/remapping service. Old unused
  tail entries are zeroed. No checkpoint format is changed.
- `snrt_ramses_driver` consumes pending energy via
  `snrt_agn_photon_budget_energy`, applying the existing escaped bolometric
  group fractions once. The old mass entry delegates to the energy entry.
  Invalid budget status is propagated, not silently treated as successful
  zero emission. Successful all-group deposition alone clears pending energy
  through `snrt_agn_source_commit`. Failed lookup/source deposition keeps it.
  The existing RT snapshot follows source commit, so transport retry neither
  restores fuel nor re-emits it. The old cumulative inflow and retained-mass
  cursors, regression checks and current-efficiency source gate are removed.
  Existing unique positive sink-ID checks remain. Coarse rate arrays remain
  legacy/diagnostic data, not the photon fuel source.
- Shared `agn_select_scalar_fields` selects legacy or channel-resolved storage
  without initializing yields. Both `grow_bondi` and `AGN_feedback` call it
  before their respective loops. `agn_accrete_scalars` scales total metal and
  every declared element density by the actual swallowed gas fraction, inside
  the cell lock. Existing auxiliary fields and hydro/MHD specific-energy
  convention remain with the caller. No cold-loading cap/energy rule is
  transplanted into accretion, and no new accretion NENER exclusion is added.
- `AGN_blast` restores the pre-loading receiver density using existing
  owner-local `ind_blast` and actual globally reduced `mAGN` entries in the
  existing neighbor bins. All shared/fallback-loaded donor masses are included.
  This scalar is private to the OpenMP cell iteration and frozen before its
  source loop. Thermal and saved-replay numerators therefore match the
  pre-loading mass denominator despite earlier in-cell jets. Jet thermal
  volume weighting, cap, geometry, nominal EAGN and whole-event fallback
  remain unchanged. No full-grid copy, donor map or extra normalization pass.

## Verification actually performed

The existing `simulation/snrt/tests/run_fp15_agn_efficiency.sh` passed at O0
with Intel runtime checks, compiling actual production helpers. Extended
`snrt_agn_source_smoke` and `agn_feedback_deposition_smoke` passed at O3 too.
No new Python test or runner was introduced.

Markers:

```
SNRT_AGN_EFFICIENCY_OK
SNRT_ACCEPTED_ENERGY_SOURCE_RETRY_OK
SNRT_AGN_SOURCE_OK
AGN_NATIVE_CELL_COUPLING_SMOKE_OK
SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK
```

New native cases:

- Requested mass 100, available accepted mass 4, retained mass 3.6; only the
  accepted mass funds photons. Next event epsilon 0.2 contributes its own
  energy, not a reinterpretation of the preceding epsilon 0.1 event.
- Exhausted donor floor and zero efficiency yield zero energy; invalid
  efficiency and scalar maps reject without helper mutation.
- Actual constituent withdrawal follows the gross fraction, leaving the
  selected donor abundances unchanged after proportional density removal;
  disjoint auxiliary fields are not overwritten.
- The actual merger helper sums two sources into one, preserves a zero-fuel
  new sink, reorders slots and rejects an invalid mapping.
- The energy-input API agrees with the mass compatibility entry. An actual
  failing multigroup deposit leaves intensity and fuel unchanged; successful
  retry clears fuel. A post-source snapshot/retry illustration then produces
  zero new photons. This last check is NOT a full RT-driver integration run;
  actual source-before-transaction ordering was inspected in the driver.
- Two thermal/replay events and two jet events overlap on three unequal-volume
  cells, with two loads sharing a donor. The fixture demonstrably has a
  pre/post-loading mass-normalization difference. Both source orders integrate
  each requested thermal budget and conserve gas energy plus nonzero deferred
  energy under an active cap; mass, momentum and transported species close.
  This exercises production depositors with caller-equivalent density recovery,
  not an initialized AMR/bin traversal. Cell-state equality across source order
  is not required; conservative capped placement can be order dependent.

Small synthetic absolute tolerances: 1e-12 for requested heat, 1e-11 for
overlap energy/mass/momentum/species, 1e-14 for receipt arithmetic; source API
comparison 1e-13 relative. Existing c_float deposition tolerances are retained.
These are focused arithmetic tolerances, not cosmological convergence limits.

Scratch: `build/agn_fuel_overlap.0ypNaI`. `pm_commons`, accepted-energy/source
helpers, project-local `init_sink` and `snrt_ramses_driver` compiled. Actual
`sink_particle.kjhan.f90` compiled for SNRT+PHASE0, CPU+PHASE0 and CPU legacy.
Compiler: `mpiifx -qopenmp -fpp -O0 -g`, NDIM=3, NVECTOR=500, NPRE=8,
NVAR=18, NENER=0, SOLVER=hydro, LONGINT, QUADHILBERT,
OUTPUT_PARTICLE_POTENTIAL, USE_FFTW. Scratch module interfaces precede `bin`.
No shared executable/module was overwritten. `git diff --check` passed.

## Limits / no new physical approval

The native path now has accepted radiative fuel and conservative overlapping
mechanical deposition, NOT a joint radiative/mechanical energy partition.
Legacy thermal/jet EAGN still uses its existing retained-mass/current-mode
recipe and approximate `3e10` c; radiation uses exact c at accretion-event
units. Simultaneous legacy AGN+SNRT remains prohibited. Live SNRT sinks remain
serial fresh-start only; pending energy is intentionally not restart-persistent.
MAD spin energy, SED physics, obscuration, winds versus radiative coupling,
auxiliary-energy advection, dust and larger convergence/live qualification are
not approved by these tests. No full binary link, MPI evolution, snapshot,
large job, AGB-source investigation, commit or push occurred. Existing unrelated
dirty changes and staged fixture moves were preserved.

## Opus disposition and bounded repairs

Opus verified all three actual Fortran deliverables and issued CONDITIONAL PASS.
The original review remains verbatim in its own record. Changes below address
that review, not a new physical partition or additional implementation bundle.

1. Finite negative constituents / total-metal density above gas density now
   return a distinct skip status from `agn_accrete_scalars`. The caller unlocks
   and skips the WHOLE event before gas, particle, sink or radiation receipts
   commit. Nothing is clamped: clamping would invent/remove species mass and
   was not approved. Such cells can therefore have reduced accretion until
   their advected composition is admissible; this is not a positivity fix.
   Structural maps, nonfinite inputs and invalid receipt units/efficiency still
   fail through the existing collective error boundary. Small integer notice
   counts reduce per call and print only the first warning of each type.
2. Did NOT mirror the AGN `ndim/=3` return into `grow_bondi`, which would silently
   turn off existing lower-dimensional accretion. Instead both scalar helpers
   accept the actual last hydro index; the shared selector and accretion pass
   `ndim+2`. Default five-field depositor behaviour stays unchanged. Native
   1D/2D index/withdrawal cases now pass. These are layout tests, not 1D/2D
   evolution qualification. The real `amr_step.jaehyun:770-793` call is already
   nested inside `if(hydro)` and `if(sink)`; no new hydro-disable gate is needed.
3. Opus identified a *possible*, not reproduced, zero initial-density reference
   in reception cells. For finite `rho_initial==0`, the receipt now uses current
   donor density for the existing 0.75 floor. The caller counts/reports this
   fallback. Negative/nonfinite references remain errors. This preserves a
   positive per-event residual, rather than using the zero reference to empty
   a donor or aborting. It does NOT reconstruct a missing original coarse/fine
   step snapshot: repeated fallback events can cumulatively exceed 25% of the
   original cell mass. Exact reception-history qualification remains a live
   follow-on, not a new AMR-copy framework in this package.
4. Keep the stated live-claim limit: the native fixture tests caller-equivalent
   density recovery with actual depositors; bin traversal was inspected and
   compiled, not evolved. A short live overlapping-source case is needed before
   claiming production/publication qualification. No run is launched here.
5. `bin_ledger` with its separate `patch/ledgertest` override is NOT qualified
   for this implementation; its initialization/build dependencies are separate
   and were already divergent. Use the active `bin/Makefile` / `patch/lagRamses`
   build only. Do not infer any other tree's readiness or silently migrate it.
6. Updated the stale inflow-cursor warning and packed-field comment. The energy
   API's bolometric luminosity output is used for validity checking only; it is
   not a new time-resolved source luminosity measurement.

After repairs, the full existing native runner passed again at O0/check-all,
and the amended depositor smoke passed at O3. All three sink caller compiles
and the SNRT driver compile passed again in the same scratch directory.
Added checks cover skip-without-clamping, over-unity total-metal skip, zero
reference/current-density floor, and 1D/2D scalar placement. No new Python
validator, large job or per-helper audit was added. One focused Opus closure
review addresses these conditions only, under the existing bundle cadence.

The focused Opus closure returned PASS, with no further local repairs. It
confirmed C1/C2 closed and the C3/C4/C5 limits correctly recorded. Its one
non-blocking note is retained: the zero-reference counter counts receipt
attempts, including those subsequently skipped for invalid composition; it
is not a count of committed accretion events. No accounting change was needed.
