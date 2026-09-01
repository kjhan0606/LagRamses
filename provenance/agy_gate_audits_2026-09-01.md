# Independent AGY gate audits

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: `/home/kjhan/local/bin/agy` (AGY/Antigravity CLI), version 1.1.21  
Model selected: `gemini-3.1-pro-high`

## Execution record

G0, G1, G2, G3, G4, G5, G6, and G7 were requested as separate AGY CLI
sessions. Each prompt restricted the scope to one gate and prohibited edits,
builds, job submission, simulation launch, deletion, copying, and commits.
The G3/G5/G6 retries used `--sandbox --dangerously-skip-permissions` because
headless mode otherwise auto-denied a read command; the prompts still limited
the agent to read-only inspection. No additional project artifact attributable
to these sessions was found in the project status after the audits.

The first all-gates umbrella prompt and the earlier Claude `--model agy`
attempt are not AGY audit evidence. The former was stopped when the requested
scope was clarified; the latter was an incorrect CLI/model interpretation and
was stopped.

## Gate scorecard

| Gate | AGY status | Independent local interpretation |
|---|---|---|
| G0 inventory/provenance | **BLOCKED** | Correct: unresolved hashes, dirty external checkout, incomplete environment/license metadata, and no production namespace refusal. |
| G1 native contract/algorithms | **PASS** | Closed within native-contract scope after six-point native/JAX differential re-audit; MPI/checkpoint qualification remains G5. |
| G2 physical stellar yields | **BLOCKED** | Correct: only test/legacy data; no approved full physical grid. |
| G3 SED/events/AGN ledgers | **BLOCKED** | Correct: candidate clamps, non-coeval ledgers, no SNIa DTD, and no approved obscuration. |
| G4 dust/RT/thermochemistry | **BLOCKED** | Correct: dust model incomplete, atlas edge clamp, and source-cell science convergence open. |
| G5 live coupling/HPC | **BLOCKED** | Correct status; the code has an opt-in live hook, but it is not production-qualified or active in the registered baseline. |
| G6 science qualification | **BLOCKED** | Correct: pilots and internal conservation passes cannot support science claims. |
| G7 release/publication | **BLOCKED** | Correct: no release package, manuscript, regeneration scripts, or final external audit record. |

No gate is promotable at this time. Passing unit/ledger tests below are
component results, not gate promotions.

## G0 — inventory, provenance, and identity

**AGY decision: BLOCKED.**

Evidence reported by AGY:

- `manifests/lrd_jwst_external_assets.json` has `p4_snapshot_hdf5` marked
  `not_migrated` with `sha256: null` and `snrt_cuda_smoke_executable` marked
  `missing` with `sha256: null`.
- The lagRamses checkout is `available_external_dirty`; the registered
  `output_00011` build identity is
  `65d080243d29f918490148352543269796462e30-dirty`.
- The manifest does not yet carry the required license, Python environment,
  and JAX-version fields for a complete reproducibility record.
- The `output_00011` baseline label is correctly transitional Phase-0 and not
  pure legacy.

Required closure:

- complete all source/output/executable hashes or explicitly mark the asset
  unavailable;
- record licenses, environment lockfile, JAX version, compiler/flags, and
  effective namelist;
- make the source/build identity clean and immutable for production;
- add a machine-readable production namespace check that rejects legacy and
  transitional assets and the embedded fallback.

## G1 — native contract and algorithmic correctness

**AGY decision: BLOCKED.**

Confirmed evidence:

- `stellar_source_increment.f90` implements cumulative
  `C(age+timestep)-C(age)` through `cumulative_difference`.
- `stellar_ramses_runtime.f90` updates `indtab` after deposition, but no
  explicit restart/idempotence test exists.
- `stellar_yield_interpolation.f90:157--158` clamps out-of-range queries to
  the channel minimum/maximum, contrary to the fail-closed production plan.
- Legacy feedback contains explicit `1d51` conversions, whereas the Phase-0
  path consumes `source%energy / scale_energy`; no one-event dimensional test
  proves the two conventions cannot be mixed.
- The native audit lacks a complete non-finite-value check and there is no
  native-Fortran/JAX differential grid test.

Required closure:

- verify the RAMSES `aexp**2` time conversion and enforce explicit `age_yr`/
  `dt_yr` boundaries;
- add interval-boundary, repeated-call, and restart-equivalence tests;
- make energy units explicit and test a single SN event in cgs and code units;
- replace endpoint clamping with a production error;
- add IEEE finite checks, field-map/NVAR/He/metal/delayed-cooling startup
  validation, and native/JAX differential tests.

## G2 — physical stellar yield model

**AGY decision: BLOCKED.**

Evidence:

- The legacy asset is H/O/Fe only, with SHA256
  `ba1099c5a4c3afe5e9ba28b3eb59d2e85fd3d40b7e7cb4ec30799eec00a5ac2e`.
- The canonical fixture is nine rows and has no production-complete grid,
  no approved source sidecar, and no SNIa/PISN data.
- The 11-element contract is a reduced-chemistry foundation, not a complete
  nucleosynthesis model.
- Wind/AGB/SNII sources, IMF/binary consistency, lifetimes, fallback/remnant
  model, and returned-mass/energy closure are not approved.

Required closure:

- select and cite physical source families for each stellar channel;
- generate complete mass--metallicity--age Cartesian coverage and an age-zero
  anchor;
- make IMF and binary/single-star assumptions consistent with the SED;
- close `initial = living + remnant + returned` and the elemental/energy/
  momentum ledgers;
- formally choose expanded nucleosynthesis or publish only a quantified
  11-element reduced-chemistry scope;
- attach an immutable provenance/approval sidecar and pass the existing
  fail-closed audit.

## G3 — SED, SNIa/PISN, and AGN ledgers

**AGY decision: BLOCKED.**

Evidence:

- The BPASS candidate records metallicity/age clamping, per-initial-mass
  normalization as an assumption, and `fesc=1`; it is explicitly marked
  candidate rather than science input.
- The current stellar candidate is at output 00011 while the AGN ledger is at
  output 00017; the strict merger correctly rejects this non-coeval pair.
- The retained AGN pilot also uses a different group count from the P0
  nine-group stellar ledger.
- SNIa DTD/event yields and PISN/PPISN eligibility are not implemented.
- AGN interval, inflow-versus-retained mass, obscuration, and complete
  energy/momentum closure are not yet approved.

Local cross-check found a provenance discrepancy requiring G0 resolution:
AGY reported SHA256 `aeb9b72fbc0716ec25836f431bf9a897bd52011c8fe5553ad185c293dc97e1c4`
for `feedback_transition_phase0_output_00011_bpass_stellar_photon_ledger.json`,
while the current `/gpfs` file hashes to
`939c7ab2305b1de6e44c7e021801a2a0bbf580acc4372374d544aff90c74954b`.
Neither value is accepted until the artifact is rehashed and the manifest is
updated.

Required closure:

- obtain strictly coeval stellar and AGN ledgers on the same group edges;
- remove or formally justify SED domain clamps and `fesc=1`;
- validate normalization per initial stellar mass;
- implement SNIa as a separate DTD convolution and PISN behind an explicit
  population/core-mass gate;
- deduplicate AGN coarse-step records and close accretion, radiative, jet,
  radiation-pressure, mass, and momentum ledgers.

## G4 — dust, RT, and thermochemistry

**AGY decision: BLOCKED.**

Evidence:

- P4/P5 dust tests pass technically, but the current Draine/WD01 sidecar is
  not an approved source-specific dust-to-metal/depletion model.
- Only absorption, dust heating, and an absorption-only momentum diagnostic
  are implemented; scattering, IR re-emission, and full radiation pressure
  are not production physics.
- The P4 thermal atlas high-edge clamp is documented near the output 00017
  epoch.
- Local source-cell convergence remains open: the limit-0.25 test reports
  maximum local `xHII` difference `0.0215405`, and the compact3 control reports
  `0.134382`; both are internal controls, not science promotion.

Required closure:

- approve the dust mixture, opacity source, dust-to-metal relation, depletion,
  and spectral weighting;
- implement or explicitly scope scattering, IR re-emission, and radiation
  pressure with separate ledgers;
- extend the thermal atlas to bracket all production epochs without clamping;
- establish physical source-size/birth-cloud and source-cell spatial
  convergence, plus angular/timestep/mesh/reduced-light-speed convergence.

## G5 — live RAMSES coupling and HPC qualification

**AGY decision: BLOCKED.**

AGY found no restart/AMR/MPI/performance qualification suite and identified:

- the embedded yield fallback is still reachable when
  `PHASE0_YIELD_TABLE` is empty;
- `snrt_ramses_driver.f90` maintains an `accounted_mass` array without a
  demonstrated cross-MPI atomic ownership protocol;
- radiation heating and delayed-cooling energy lack a production
  no-double-count proof;
- the live source-cell limiter is not recomputed from changing source/gas
  state;
- floating-point reduction choices may prevent bitwise reproducibility.

Important interpretation correction: `snrt_ramses_driver.f90` does contain an
opt-in in-process RT hook that updates `uold`, so the source tree has partial
live capability. However it is enabled only by `SNRT_RT_ENABLE=1`, requires a
CUDA device, and the registered transitional checkpoint log reports zero
`snrt_advance` and `snrt_diagnose` time. The current validated workflow is
therefore static/post-processed for science purposes; the existence of a live
hook does not constitute live production qualification.

Required closure:

- make approved table/sidecar refusal unconditional in production;
- implement and test AMR/MPI ownership, ghost/boundary deposition,
  subcycling, restart idempotence, and no-double-count ledger;
- qualify a clean binary at multi-rank production shape with memory/I/O/wall
  time measurements;
- document and test the live timestep/source-cell contract.

## G6 — science qualification and comparison

**AGY decision: BLOCKED.**

Evidence:

- The registered baseline is correctly identified as transitional Phase-0 and
  has executable SHA256
  `9e28e6445a8876290249230b616f5b91f11b345eac0d1bd94c96b721956d695c`.
- P4/P5 photon and thermal conservation passes exist, but source-cell spatial
  convergence is open and the inputs are zero-dust/AGN-pilot controls.
- The baseline lacks the RT/dust/source-ledger state needed for a matched
  production comparison.
- SNIa, PISN, AGB-delay, and approved yield/SED/dust inputs are missing.
- No uncertainty/sensitivity envelope supports paper observables.

Required closure:

- define observables and predeclared numerical/model acceptance thresholds;
- compare production and transitional baseline at identical state, epoch,
  resolution, and measurement definition;
- complete yield, IMF/binary, SNIa, AGN, dust, and RT sensitivity matrices;
- include source-cell convergence, uncertainty intervals, and independent
  replay; reject pilot PASS labels as science claims.

## G7 — publication package and release reproducibility

**AGY decision: BLOCKED.**

Evidence:

- `paper/README.md` states that no Paper-III manuscript source is migrated.
- No finalized release directory, figure/table regeneration scripts,
  effective production namelists, or complete release manifest exists.
- The external audit schedule is recorded as `not_registered`; a scheduled
  request is not evidence that a completed audit exists.
- Upstream G0–G6 are not closed.

Required closure:

- create a release directory containing source/executable/environment
  identities, approved sidecars and hashes, effective namelists, outputs or
  certified reductions, validation reports, and regeneration scripts;
- write the manuscript limitations section, including reduced chemistry,
  omitted channels, dust assumptions, source resolution, and static/live scope;
- record a completed final external algorithm/wiring/physics audit;
- run an independent clean-environment reproduction from the release package.

## G0 final re-audit after implementation

After adding the production-readiness overlay, CPU/JAX environment
fingerprint, fail-closed manifest auditor, and synthetic regression test, a
fresh AGY audit was commissioned for **G0 only**. AGY returned **BLOCKED**.
The detailed record is in
[`agy_g0_final_audit_2026-09-01.md`](agy_g0_final_audit_2026-09-01.md).

The local audit also returned `status=blocked` / exit code `2`; its report is
`simulation/snrt/data/g0_production_manifest_audit.json`. The implementation
is therefore recorded and tested, but G0 is not promoted. G1 was subsequently
handled as an independent native-contract gate and has now received an AGY
**PASS**; that does not waive the G0 blockers.

At that audit snapshot the blockers included the non-migrated 866 GB HDF5
payload, missing CUDA smoke executable, dirty source/repository identity,
incomplete required asset registrations and per-asset approval metadata, and
an unlocked environment fingerprint. Subsequent local work qualified the CUDA
smoke candidate and locked the CPU environment; the current local manifest
audit remains blocked on approval metadata, source identity, and the HDF5
reference/provenance policy.

## G1 final re-audit after differential-matrix expansion

The first G1 review was **CONDITIONAL PASS** because its JAX differential
covered one interior point only. The matrix was expanded to six queries:
exact low/high corners, mass/metallicity/age boundary faces, and an interior
point. The native Fortran runner and CPU JAX 0.11.1 differential both passed.

A fresh **G1-only** AGY re-audit (CLI session `4842`) returned **PASS**. The
full report is in
[`agy_g1_reaudit_2026-09-01.md`](agy_g1_reaudit_2026-09-01.md).

G1 is therefore closed within its scope. AGY explicitly kept true atomic
MPI/checkpoint behavior out of scope for G1 and assigned that qualification to
G5; physical yield-paper approval remains a G2 task.

## G2 audit after source-contract and closure scaffolding

A fresh **G2-only** AGY audit (CLI session `33756`) returned **BLOCKED**. The
project-local record is
[`agy_g2_audit_2026-09-01.md`](agy_g2_audit_2026-09-01.md).

The blocker is external and deliberate: the source matrix contains literature
candidates but no approved, checksummed physical full-grid asset for required
channels 1--3. The legacy table and nine-row fixture remain test/comparison
inputs only. The converter, configuration contract, fail-closed preflight, and
native population mass ledger are implemented and tested, but they cannot
manufacture scientific source data or an approval identifier. G2 therefore
cannot be promoted until the source package and its approval sidecar are
provided.

## G2 re-audit after population-ledger wiring

The G2-only re-audit was repeated after the population-ledger and SSP-driver
changes. AGY again returned **BLOCKED**, while explicitly confirming the
following resolved controls: one-time population living-mass derivation,
channel ejecta/returned closure, terminal-remnant ownership, fail-closed
SNIa/PISN activation, and passing G2 configuration/converter/preflight tests.
The detailed record is in
[`agy_g2_reaudit_2026-09-01.md`](agy_g2_reaudit_2026-09-01.md).

The sole G2 promotion blocker remains the absence of approved, checksummed
physical wind/AGB/SNII full-grid assets and their provenance/approval sidecar.
No synthetic or legacy input is being promoted in their place.

## G2 final re-audit after Python/native ownership symmetry

After the Python contract was brought into symmetry with the native ownership
check, AGY performed one final read-only G2-only audit. It again returned
**BLOCKED** and confirmed the ownership checks, population-ledger wiring,
deterministic converter, fail-closed SNIa/PISN paths, and G2 test results. The
final project-local record is
[`agy_g2_final_reaudit_2026-09-01.md`](agy_g2_final_reaudit_2026-09-01.md).

At that audit snapshot, a read-only search of the two relevant `/gpfs`
locations found only the synthetic/legacy yield files and their audit reports;
no physical wind, AGB, or SNII source package was present. G2 was therefore
waiting for the explicit source selection and provenance/approval input listed
in the final audit.

## G2 candidate-source re-audit after public-source staging

After the Limongi & Chieffi 2018 CDS and NuGrid Set1ext candidate package was
staged on `/gpfs`, a new read-only G2-only AGY re-audit inspected the source
files, acquisition manifest, candidate parser/report, and preflight. AGY
returned **BLOCKED**. The project-local record is
[`agy_g2_candidate_source_reaudit_2026-09-01.md`](agy_g2_candidate_source_reaudit_2026-09-01.md).

AGY independently verified the manifest hashes, G2 contract, native
population closure, terminal-remnant ownership, and generic converter guards.
It confirmed that neither candidate has the required age-resolved cumulative
history, canonical energy/momentum fields, complete runtime mass coverage,
disjoint channel semantics, or license/approval sidecar. This re-audit does
not promote G2 or authorize a simulation.

## Local checks run alongside the AGY audits

Passing component tests included `stellar_yield_asset`, `agn_photon_ledger`,
`p0_sed_closure`, `p0_photon_ledger`, `stellar_photon_ledger`,
`merge_photon_ledgers`, HDF5 staging/ingestion, native stellar catalogue,
Draine/dust opacity, P4/P5 dust runner, P0 smoke, P1 validation, P2/P3
validation with two virtual CPU devices, P8 sharding/angular controls, and
source deposition. The relevant negative/limiting tests also passed by
reporting open gates: `p5_source_cell_convergence`,
`p5_controlled_deposition`, and `p5_refined_mesh_convergence`.

These results certify component behavior only. They do not override the AGY
BLOCKED gate decisions because the missing scientific assets, approvals,
cross-domain ledgers, native restart/MPI tests, and production release package
remain unresolved.

## Critical path from the independent audits

1. G0: repair manifest/hash/environment identity and enforce namespace/fallback
   refusal.
2. G1: closed within native-contract scope; MPI/checkpoint qualification is
   retained for G5.
3. G2: approve sources and generate the physical complete yield grid.
4. G3: produce coeval approved SED/SNIa/AGN ledgers and optional PISN gate.
5. G4: approve dust/thermal physics and close spatial/temporal convergence.
6. G5: qualify live AMR/MPI/restart coupling and production performance.
7. G6/G7: perform science comparison, uncertainties, external audit, and
   reproducible release.
