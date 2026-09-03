# Stellar/AGN feedback and yield-table implementation plan

Status: active plan, 2026-09-01.  The stopped comparison run is recorded in
[`legacy_feedback_baseline.md`](legacy_feedback_baseline.md).  The plan is
deliberately fail-closed: the legacy three-species table is a comparison input,
not a canonical Phase-0 yield table, and no missing physical values will be
invented by an adapter.

The native/phase0 mirror contract is closed with an AGY PASS. The G1
implementation mirror, six-query Fortran/JAX differential matrix, and final audit are recorded in
[`agy_g1_reaudit_2026-09-01.md`](agy_g1_reaudit_2026-09-01.md). G2 code
scaffolding and review-only source staging are implemented. The current G2
AGY audit gives the engineering path PASS and the scientific path BLOCK
because no approved physical full-grid asset has been supplied. Its source-selection
and physical-contract records are
[`g2_source_selection_matrix_v1.json`](../simulation/snrt/config/g2_source_selection_matrix_v1.json)
and
[`g2_physics_contract_v1.json`](../simulation/snrt/config/g2_physics_contract_v1.json).
The re-audit is recorded in
[`agy_g2_reaudit_2026-09-01.md`](agy_g2_reaudit_2026-09-01.md); it confirms the
code-level controls but leaves G2 blocked on the physical source package.
The final ownership-symmetry re-audit is recorded in
[`agy_g2_final_reaudit_2026-09-01.md`](agy_g2_final_reaudit_2026-09-01.md) and
has the same physical-asset blocker. A subsequent candidate-source re-audit,
after public Limongi/NuGrid files were staged and hashed, is recorded in
[`agy_g2_candidate_source_reaudit_2026-09-01.md`](agy_g2_candidate_source_reaudit_2026-09-01.md);
it independently confirms that G2 remains blocked.
The subsequent lossless source adapters and internal-closure diagnostics are
recorded in
[`g2_source_adapter_review_2026-09-01.md`](g2_source_adapter_review_2026-09-01.md).
They emit zero canonical rows and do not change the gate decision.
The current CC BY 4.0 Huscher et al. 2025 AGB release is now staged and
audited as well; its license is resolved, but its population-table
normalization is not.
The CC BY 4.0 Boccioli & Roberti 2026 CCSN release is also staged and audited.
Its F23 branches close internally and improve candidate SNII coverage to
11--40 Msun, while its LC18 failed-model Wind inconsistency, missing
machine-readable energetics, and absent age history remain fail-closed.
The current G2 candidate-source and high-mass projection audits are recorded in
[`fp1_source_package_selection_plan_2026-09-02.md`](fp1_source_package_selection_plan_2026-09-02.md),
[`claude_opus5_g2_source_package_staging_final_audit_2026-09-02.md`](claude_opus5_g2_source_package_staging_final_audit_2026-09-02.md),
and
[`agy_g2_source_package_staging_audit_2026-09-02.md`](agy_g2_source_package_staging_audit_2026-09-02.md).
They preserve zero canonical rows and zero runtime deposition; the next step
requires an explicit project source/physics approval, not a synthetic default.
The manifest-scoped fingerprint and selection gates now record Sukhbold
W18/N20 as a review-only validation branch for the 40--120 Msun seam; they
keep `production_source_id=null` and runtime activation disabled. This is a
comparison choice, not a physical source approval, and the next conversion
step remains blocked until the population, age/decay, energy/momentum,
channel-ownership, and licensing decisions are explicitly closed.
The F-P2 DTD interval kernel and separate expected-event ledger are now also
compiled in both the native mirror and production source order, with
subdivision/restart and mass/energy/momentum closure tests plus a fail-closed
contract. Its `alpha=-1` shape is a mathematical fixture only; no physical
DTD normalization or SNIa event yield has been selected. DTD and event-yield
literature candidate matrices are staged and machine-audited, but neither
candidate set can authorize runtime activation.
The CC BY 4.0 Limongi et al. 2024 transition-fate evidence is now pinned and
audited as a reference rather than a yield source. It shows that the runtime
8 Msun SNII boundary is not a universal explosion threshold and leaves the
8--8.8 Msun edge as a non-interpolable fate-policy seam. The inherited
comparison catalogue has also been audited: all 42,342 stars have
`Z <= 1.18135e-9`, more than 4.43 dex below the lowest staged positive-Z
full-grid lower edge at the comparison maximum. No floor or solar
extrapolation is permitted.
Roberti, Limongi & Chieffi (2024) now provides a CC BY 4.0 sparse candidate at
`Z=0`, `3.236e-7`, and `3.236e-6`, but only for 15 and 25 Msun. Its four
missing zero-Z MRT columns, quarantined `025z600` mass-budget outlier,
unselected rotation population, and unresolved wind/terminal ownership keep
it from defining the production domain.
Heger & Woosley (2010) now supplies an official VizieR exact-`Z=0` candidate
with 120 masses from 10--100 Msun and 5,760 coupled energy--piston--mixing
coordinates. It improves primordial mass coverage but does not close the
8--10 Msun runtime edge or define positive-Z behavior, and none of its free
explosion/mixing dimensions is population-selected.

The 2026-09-02 Fable audit is independently reproduced in
[`fable_sn_agn_independent_reproduction_2026-09-02.md`](fable_sn_agn_independent_reproduction_2026-09-02.md).
It confirms the production **BLOCK**: the AGY/G1 PASS applies to the native
mirror, while `bin/Makefile` builds the unclosed `patch/lagRamses` runtime.
The independent reproduction is the authoritative disposition for F1--F17
and takes precedence over the earlier broad “G1 closed” wording below.

## Current evidence

- The legacy asset is `/gpfs/kjhan/Run_JWST/opt_run/yield_table.asc`.  It has
  3 species (H, O, Fe), no explicit channel axis, no remnant/energy/momentum
  columns, and no channel-resolved release history.
- The Phase-0 Fortran reader expects 32 canonical fields for 11 elements, but
  the current integration fixture has only one mass per channel and no SNIa or
  PISN rows.  The runtime integrates over fixed mass intervals, so that fixture
  cannot be a production asset.
- The first `/gpfs` audit is reproducible: the legacy file has 12,000 valid
  data rows in 60 blocks of 200, with SHA256
  `ba1099c5a4c3afe5e9ba28b3eb59d2e85fd3d40b7e7cb4ec30799eec00a5ac2e`, and is
  classified `legacy_only`.  The 9-row Phase-0 fixture has no row-level
  closure failure, but is rejected for its one-point channel mass grids,
  pending channel approval, and missing provenance sidecar.
- The machine-readable reports are
  [`legacy_yield_table_audit.json`](../simulation/snrt/data/legacy_yield_table_audit.json)
  and
  [`phase0_validation_yields_audit.json`](../simulation/snrt/data/phase0_validation_yields_audit.json).
- The current source path is a post-processed chain
  `lagRamses coarse-state/rate ledger -> AGN SED -> photon ledger -> static
  SNRT`.  It is not live RT--hydrodynamics feedback.  Dust absorption/heating
  is diagnostic; scattering, IR re-emission, radiation pressure, and a
  chemistry/RHD closure are not approved production physics.
- The new audit tool is
  [`audit_stellar_yield_asset.py`](../simulation/snrt/tools/audit_stellar_yield_asset.py).
  It classifies the legacy format, checks canonical rows, and gates the
  runtime coverage and provenance requirements in
  [`stellar_feedback_contract_v1.json`](../simulation/snrt/config/stellar_feedback_contract_v1.json).
- The batch preflight wrapper is
[`validate_stellar_yield_asset.sbatch`](../simulation/snrt/validate_stellar_yield_asset.sbatch);
  it requires `PHASE0_YIELD_TABLE` and stops before any production executable
  can consume an unapproved table.
- The G2 source matrix records Jost et al. and Boccioli & Roberti 2026 (CCSN),
  Karakas & Lugaro, Huscher
  et al. 2025, and Osborn
  et al. (AGB/binary AGB), and Keegans et al. (SNIa event yields) as candidates,
  not approvals. Limongi & Chieffi 2018, NuGrid Set1ext, Huscher, and
  Boccioli--Roberti candidate files are
  now staged under `/gpfs/kjhan/LRD_JWST/external/g2_candidates/` with the
  acquisition manifest and reproducible source audit. Huscher and
  Boccioli--Roberti are CC BY 4.0; license terms for the other selected
  sources, conversion, channel closure, and approval records remain incomplete.
  The OUP record confirms that the Jost supplementary data exist, but the
  dynamically hosted package has not been retrieved and fingerprinted; values
  are not reconstructed from the paper. The official Heger--Woosley VizieR
  package is staged for internal scientific review with public redistribution
  still unresolved.
- The production timestep driver executes the population ledger on the same
  cumulative SSP states used for each source increment.  It checks
  channel-level ejecta/returned closure, terminal-remnant ownership,
  nonnegative living mass, and consistency with the independently advanced
  RAMSES particle mass.  Production startup uses the same strict yield audit
  as the native oracle. SNIa/PISN activation remains fail-closed until their
  distinct event models exist.
- The Limongi and NuGrid source-specific adapters verify all manifest hashes
  and preserve actual source values and duplicates. Source-internal closure
  diagnostics pass: no Limongi source-yield sum exceeds initial mass; NuGrid
  total/pre-explosion rows close to the printed-table rounding level. These
  diagnostics do not supply missing age histories, channel semantics, energy,
  momentum, license, or approval.
- Primary-source review now resolves the source-level component meanings:
  Limongi set R separates wind and non-wind terminal ejecta over 13--25 Msun
  and is wind-only above 25 Msun; NuGrid supplies wind-only, wind+pre-SN, and
  wind+delayed-SN totals. Project runtime channel ownership, time distribution,
  energy, momentum, and approval remain separate blockers.
- The reduced-chemistry mass contract now preserves
  `untracked_ejecta = returned_mass - sum(ejecta_H..Fe)` without adding a
  thirty-third table field. Native validation rejects a negative residual,
  channel/population ledgers record it explicitly, and RAMSES adds it only to
  generic metallicity. Disabled individual metal fields no longer remove their
  mass from the total-metal source. G1 native/JAX regression and G2 preflight
  pass after this contract change.
- Limongi decay-horizon sensitivity is reproducible in
  `g2_limongi_decay_projection_audit.json`. The pinned decay matrix covers 307
  of 333 source nuclides directly; checksummed NUBASE2020 records resolve 22
  short-lived beta decays and retain four nuclides whose half-lives exceed the
  largest tested horizon. No source nuclide remains unresolved. A 1 Myr
  projection increases the non-IMF-weighted terminal Fe grid sum by about
  3.373 Msun, so no-decay parent-element mapping is unsuitable for Fe abundance
  claims. The projection horizon remains a project approval decision.
- The Huscher 2025 audit resolves a large portion of the AGB mass--metallicity
  source search: 120 single-star models span 0.8--7 Msun and ten metallicities.
  It also exposes a fail-closed normalization defect. Integrating the already
  IMF-weighted population Mdot table under its documented units gives
  1326--2895 Msun per claimed Msun formed. No `1e4` or other hidden factor is
  inferred; the age table remains validation-only until clarified. The
  single-star rows still lack 7--8 Msun, per-star age release, S/Ca/Fe,
  energy, and momentum.
- The Boccioli--Roberti 2026 audit resolves the strongest currently staged
  CCSN composition/remnant candidate. It verifies 206 models and full reduced
  element coverage. F23 single covers each integer mass from 11 to 45 Msun and
  its released post/wind components close, but it is solar-only. All 56 failed
  LC18 models have positive summary wind masses and zero Wind tables, so that
  branch is quarantined. No branch tabulates per-model explosion energy,
  canonical momentum, or age-resolved release.
- Limongi 2024 constrains the source-solar nonrotating terminal-fate seam but
  contributes no yield nodes. Its potential ECSN lower edge depends on the
  electron-capture core-mass criterion, and the Stockinger e8.8 event cannot
  be promoted into a population fate law.
- The inherited comparison population quantitatively requires an
  ultra-low/primordial-metallicity review: all 42,342 stars lie below every
  staged positive-Z full-grid node. The baseline does not select the future
  production domain, but it rules out silently flooring these stars into the
  current tables.
- The Roberti 2024 source places all baseline metallicities inside its
  zero-to-first-positive coordinate interval, but this is not full-grid
  coverage: the source has only two masses, no approved metallicity or mass
  interpolation, no selected rotation population, and no age-resolved or
  wind/terminal-partitioned release. The official MRT subset has 30 of 34
  source models and `025z600` remains quarantined.
- The Heger--Woosley 2010 source provides dense exact-`Z=0` terminal yields
  from 10--100 Msun and source kinetic energies at infinity. It does not
  authorize interpolation into the 8--10 Msun or positive-Z domains, and its
  piston, artificial mixing, fallback/remnant, omitted neutrino-wind, and
  redistribution policies remain blockers.

## Priority order and acceptance gates

### P0 — make the executed yield and source contracts physically executable

The first P0 work is compiled-tree closure, not physical-table import:

1. **P0.1 source identity/build parity (F2/F15):** use `patch/lagRamses` as
   the canonical production source and retain `simulation/snrt/native/phase0`
   only as a differential oracle. The production-linked harness must force a
   `bin/Makefile` build, verify binary linkage and startup smoke, and record
   source/config/tool/binary evidence; the native G1 runner cannot close this
   gate by itself. The gate also checks the source sidecar and records
   `NVAR/NENER/inener/imetal/idelay/ichem` in its JSON report. The build and
   smoke logs carry the binary SHA-256, evidence survives descendant commits
   through hash pinning plus an ancestor check, and the native shared profile
   is a bounded differential diagnostic rather than a fixed hash set. The
   implemented gate and its current PASS result are recorded in
   [`p0_source_parity_gate_2026-09-02.md`](p0_source_parity_gate_2026-09-02.md).
2. **P0.2 time and interval semantics (F3/F4):** fix the compiled age unit at
   the table boundary and make cumulative release telescope over variable
   timesteps, first interval, retry, and restart. **Status (2026-09-02):
   CONDITIONAL PASS.** The implementation and bounded native/production-linked
   tests pass; closure still requires the committed module set and an explicit
   cross-step physical-age/restart evidence record.
3. **P0.3 HDF5 stellar state (F5):** round-trip `ptypep`, `tpp`, `mp0`,
   `indtab`, and the stellar mass ledgers. **Status: CONDITIONAL PASS.** The
   linked writer/reader test preserves nonzero release state; final closure is
   limited to an uninterrupted-versus-restarted physical stellar continuation
   test. Generic hydro, gravity, AMR, ksection/CPU-box, and checkpoint-reader
   hardening are not P0.3 blockers and are tracked separately in
   [`long_term_hdf5_restart_validation_backlog.md`](long_term_hdf5_restart_validation_backlog.md).
4. **P0.4 fail-closed runtime (F7/F8):** remove embedded fallback, reject
   out-of-domain queries, and make IMF/channel windows and approvals explicit.
   **Status (2026-09-02): PASS.**  Implementation, production-linked build,
   production-binary negative execution, and Claude Opus 5 re-audit pass. See
   [`p04_fail_closed_runtime_gate_2026-09-02.md`](p04_fail_closed_runtime_gate_2026-09-02.md),
   [`claude_opus5_p04_fail_closed_runtime_audit_2026-09-02.md`](claude_opus5_p04_fail_closed_runtime_audit_2026-09-02.md),
   [`claude_opus5_p04_fail_closed_runtime_reaudit_2026-09-02.md`](claude_opus5_p04_fail_closed_runtime_reaudit_2026-09-02.md),
   and the mandatory
   [`population/DTD roadmap`](feedback_population_dtd_active_roadmap.md).
   The temporary
   rejection of `binary_ssp`, SNIa, and PISN is an admission control, not a
   decision to omit their physics.
5. **P0.6 executed field/species semantics (F14/F15):** validate the actual
   binary map, prove the `NENER=0` index relations, and close
   H/He/tracked/untracked/energy/delayed-cooling ledgers.
6. **P0.5 physical asset approval (F1):** only then promote a complete,
   checksummed wind/AGB/SNII grid with age-zero, remnant, energy, momentum,
   release, and provenance closure.

The native mirror tests remain prerequisites, but their PASS no longer closes
P0. They must be joined by a compiled-tree parity and runtime test.

#### P0 acceptance gates

1. Select and approve a physical yield source for wind, AGB, and SNII, with
   explicit IMF, progenitor-mass range, birth-metallicity grid, age grid, and
   release-time convention.  Add SNIa only after its binary normalization and
   delay-time distribution are specified; add PISN behind an explicit
   population gate.  The asset must contain a complete mass--metallicity--age
   grid covering the runtime integration ranges, an age-zero anchor, cumulative
   actual ejecta (not net yields), remnant mass, energy, and momentum.
2. Add immutable provenance: source citation/version, units, conversion
   procedure, assumptions, checksum, and an approval record.  The production
   gate must reject a missing or mismatched sidecar.
3. The native interpolation and conservation contracts are implemented and
   tested in the G1/G2 mirrors; retain those tests as prerequisites for every
   physical-table import and native runtime rebuild.
4. Complete conservation tests at source, cell, restart, AMR, and MPI boundaries:
   `returned_mass = sum(tracked_ejecta) + untracked_ejecta`,
   `initial = living + remnant + returned`, no gas update from `net_*`, and
   exact energy/momentum accounting. Keep the
   legacy table on a separate comparison path.

P0 is complete only when the audit tool passes a named, cited, approved asset
and a native runtime test demonstrates non-negative increments and conservation
for every enabled channel.

### P1 — approve channel physics and AGN bookkeeping

1. **P1.1 population and fate model (mandatory):** select the single/binary
   population basis, IMF normalization, binary fraction and parameter
   distribution, metallicity dependence, stellar lifetimes, and mutually
   exclusive channel/remnant ownership.  Prove that a population-integrated
   table is not convolved with the IMF a second time.  Replace the one-point
   synthetic channel fixture with age-dependent wind and AGB release histories;
   document whether each table is cumulative or an instantaneous rate and test
   the conversion.
   **Medium-term 40--120 M☉ seam:** do not close this interval with a universal
   ZAMS direct-collapse bin. Select a source-node or pre-supernova-structure
   fate resolver with explicit metallicity, rotation/binary state, engine,
   lifetime, mass-cell, and PPISN/PISN axes. Reject out-of-hull and all unsafe
   interpolation. Keep wind, terminal ejecta, failed-collapse envelope
   ejection, remnant, and pulse history as separate ledger components. This
   remains blocked until the source package, age/decay/energy/momentum closure,
   license, sidecar checksum, and approval id are complete. See
   [`fp1_mass40_120_literature_dossier_2026-09-02.md`](fp1_mass40_120_literature_dossier_2026-09-02.md)
   and the zero-node resolver contract.
2. **P1.2 core-collapse channel (mandatory):** specify SNII energy, ejecta,
   momentum, remnant, failed-explosion/fallback, and delayed-cooling semantics
   together.  Verify that delayed cooling receives only the intended channel
   contribution and that thermal/kinetic energy is not double counted.
3. **P1.3 SNIa DTD (mandatory):** implement SNIa as a distinct SSP event-rate
   convolution.  Select and pin the DTD shape, minimum delay, normalization per
   initial stellar mass, binary progenitor assumption, metallicity dependence,
   event yield/energy, and stochastic versus expectation-value realization.
   Test analytic DTD integrals, timestep telescoping, restart continuity,
   single-event closure, and population-integrated mass/energy closure.
4. **P1.4 PISN/PPISN population decision (mandatory gate):** select a
   stellar-evolution/fate source and implement eligibility in terms of the
   declared population, metallicity and core mass.  The reviewed scientific
   configuration may explicitly disable PISN, but only after this gate records
   why; no universal ZAMS interval or accidental interpolation may activate it.
5. **P1.5 AGN ledger (mandatory):** audit the AGN coarse-state ledger against
   the exact accumulator-reset event:
   deduplicate `(nstep_coarse, sink_id)`, resolve conflicts, and close
   accretion, radiative efficiency, `Lbol`, thermal/jet energy, mass, and
   momentum ledgers.  Record how the SED, escape fraction, obscuration, and
   source deposition consume that ledger.

Each P1 sub-gate receives a separate read-only Claude Opus 5 physics/code
audit before the next sub-gate starts.  P1 is complete only when the selected
population model, DTD, channel/fate rules, and AGN ledgers close over a restart
and the transitional baseline can be replayed without ambiguity.

### P2 — connect approved stellar/AGN spectra and dust

1. Approve the stellar SED/IMF/metallicity/age mapping (the BPASS ledger is
   currently a candidate with explicit clamps and `fesc=1`, not a science
   input) and define the AGN SED/escape/obscuration model.
2. Define dust-to-metal normalization and source-specific opacity units.  Add
   scattering, absorption, IR re-emission, and radiation-pressure tests before
   claiming a dusty LRD solution.  The staged Draine/WD01 table alone is not a
   complete dust mixture closure.
3. Add H2/metal chemistry and Compton/photo-heating coupling only after the
   elemental source ledger is approved.

### P3 — live coupling and production qualification

Implement and validate live RT--hydrodynamics/AGN accretion coupling, timestep
and subcycling policy, momentum feedback, AMR/MPI behavior, restart determinism,
and production-scale convergence.  No production rerun is authorized by this
plan until P0--P2 gates are green.

## Immediate implementation sequence

1. Obtain explicit source/IMF/chemistry decisions and complete the review of
   the staged physical candidates; the current G2 blocker cannot be resolved
   by synthetic data.
2. Resolve the remaining source decisions exposed by the implemented review
   adapters. Limongi's exact `[Fe/H]`--`Z` mapping and identical duplicate PSN
   rows, and NuGrid's identical duplicate coordinate and component meanings,
   are resolved at source level. Still required are Limongi unit/project
   approval, decay horizon, rotation weighting, and terminal-versus-wind
   ownership; the Huscher Mdot normalization and 7--8 Msun AGB gap; the
   Boccioli--Roberti LC18 failed-wind inconsistency, F23 metallicity/population
   weighting, the 8--8.8 Msun terminal-fate seam, and an
   ultra-low/primordial-metallicity full channel grid beyond Roberti's sparse
   two-mass candidate, including resolution of its MRT omissions and
   `025z600` quarantine; the Heger--Woosley 8--10 Msun and positive-Z gaps plus
   its explosion-energy/piston/mixing population selection and redistribution;
   and an
   approved cumulative age-release construction with separately sourced
   energy and radial momentum. Only an approved adapter may feed the generic
   canonical converter and immutable sidecar. The legacy table remains
   `legacy_only`.
3. Rebuild the `/gpfs` native runtime against the approved canonical table and
   run channel-by-channel closure, source-increment, and differential tests
   before any large job.
4. Commission the G2-only AGY re-audit after the physical asset and sidecar
   are present; do not advance to G3 while G2 is blocked.
5. Re-run the P0 comparison against
   `feedback_transition_phase0_10cMpc_h_z0_epsstar002`, labeling it
   `transitional_feedback_baseline`, not pure legacy.

The one-time Fable audit requested for 2026-09-02 09:00 KST is captured in
[`fable_sn_agn_feedback_audit_schedule.json`](fable_sn_agn_feedback_audit_schedule.json).
