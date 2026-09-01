# SNRT stellar/AGN feedback production and publication-readiness plan

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Status: planning baseline; production gate closed

This plan closes every item raised by the AGY-role, Fable, and
GPT-5.6-sol reviews. Starting with B2, Claude Opus 5 assumes the independent
gate-auditor role; earlier AGY reports remain historical records. It is
deliberately fail-closed: a passing Python pilot,
a synthetic yield fixture, or a successful RAMSES exit code is not sufficient
for a production or publication claim.

## 1. Definition of done

### Production-ready

The production path must satisfy all of the following:

- no embedded synthetic yield fallback, silent extrapolation, range clamp, or
  unrecorded environment-dependent source selection;
- every enabled channel has an approved, cited, checksummed source table with
  complete mass--metallicity--age coverage and a declared IMF/population model;
- ages, times, masses, energies, momenta, luminosities, and photon counts have
  one explicit unit contract at every language boundary;
- source, cell, AMR, MPI, restart, and stellar-population ledgers close;
- all enabled feedback channels are applied exactly once, with no overlap
  between stellar, delayed-cooling, AGN thermal, kinetic, radiation-pressure,
  and dust terms;
- the exact binary, source manifest, namelist, table hashes, compiler flags,
  and run metadata are recorded and reproducible on `/gpfs`;
- native RAMSES and JAX/SNRT tests pass at the intended production scale.

### Publication-ready

In addition to the production gate:

- the scientific scope of the chemistry is explicit: either expand the table
  for nucleosynthesis claims or call the 11-element model *reduced chemistry*
  and quantify the effect of omitted species;
- source-selection alternatives, IMF/binary assumptions, dust model,
  obscuration/escape prescription, resolution, timestep, angular quadrature,
  and reduced-light-speed sensitivities are reported;
- convergence and uncertainty intervals are shown for the observables used in
  the paper;
- the baseline comparison is reproducible and labeled correctly as the
  transitional Phase-0 feedback baseline, not as pure legacy feedback;
- the paper, data/metadata release, figures, tables, and exact launch scripts
  point to immutable source and executable identities.

## 2. Current baseline and non-negotiable labels

- `/gpfs/kjhan/Run_JWST/opt_run/yield_table.asc` remains a
  `legacy_comparison_only` asset: 12,000 rows, 60 blocks, H/O/Fe, no channel,
  remnant, energy, momentum, or release-history fields.
- The 32-field/11-element/five-channel contract and the 9-row canonical
  fixture remain integration tests until replaced by a complete physical grid.
- SNIa and PISN remain disabled until their independent population gates pass.
- The stopped `output_00011` run remains a transitional Phase-0 comparison
  checkpoint; it is not a pure legacy run, not a production run, and not a
  canonical SNRT input.
- Current P4/P5 dust and SED artifacts are pilots. A Draine/WD01 sidecar and a
  BPASS candidate do not by themselves constitute an approved LRD source
  model.

## 3. Gate structure

| Gate | Purpose | Required result | Promotion decision |
|---|---|---|---|
| G0 | inventory and identity | immutable manifest, source/build map, baseline labels | permits development only |
| G1 | native contract | units, field map, time interval, channel semantics pass | permits table generation |
| G2 | physical yields | approved full grid, provenance, IMF/remnant/closure pass | permits stellar source use |
| G3 | source spectra and event channels | stellar SED, SNIa DTD, AGN ledger, optional PISN pass | permits photon-ledger merge |
| G4 | dust and thermochemistry | opacity/depletion/IR/thermal closure and convergence pass | permits dusty static SNRT |
| G5 | live coupling | RAMSES--SNRT feedback, restart, AMR/MPI, no-double-count pass | permits hydro production |
| G6 | science qualification | convergence, sensitivity, uncertainty, baseline comparison | permits paper claims |
| G7 | release | reproducible package and final independent audit | publication release |

No gate may be marked passed because a downstream gate happened to produce
finite numbers.

## 4. Work packages

### G0 — freeze the inventory and provenance

1. Freeze the `/gpfs/kjhan/LRD_JWST` project manifest without overwriting the
   existing external-asset manifest. Add each source, table, executable,
   checkpoint, and generated ledger with absolute path, size, SHA256, license,
   and status (`candidate`, `test`, `approved`, `comparison_only`).
2. Record the exact lagRamses source commit, dirty state, compiler, compiler
   flags, linked modules, JAX version, Python environment, and effective
   namelist for every validation artifact.
3. Separate three namespaces in both code and reports: legacy comparison,
   transitional Phase-0, and production candidate. A production runner must
   refuse the first two namespaces.
4. Preserve all existing snapshots and pilot outputs; use new directories for
   promoted artifacts.

**Exit:** a manifest can recreate every audit result and no production command
can resolve to the legacy or embedded fallback implicitly.

### G1 — make the native contract unambiguous

#### G1.1 Verify and fix the four Fable findings

Before changing source, add minimal reproductions and inspect the exact binary
path:

1. **Age units:** establish whether `tpp`, `texp`, `indtab`, and table `age_yr`
   are code time, conformal time, physical seconds, yr, or Gyr. Use explicit
   names (`age_yr`, `dt_yr`, `age_gyr_internal`) and one tested conversion.
   Verify the `aexp**2` factor against the RAMSES time convention. Reject a
   mixed-unit call at the interface.
2. **Interval semantics:** define the source as
   `C(age + dt) - C(age)` for cumulative quantities, apply it once, and update
   the particle progress marker only after successful deposition. Test a
   restart at the interval boundary and a repeated call for idempotence.
3. **SN energy:** trace legacy and Phase-0 energy from source-table units to
   `unew`; explicitly encode whether `energy_erg_per_star` is cgs and whether
   any `1e51 erg` normalization is already included. Remove magic conversion
   factors and add a one-event dimensional test. Do not change the suspected
   `feedback.kjhan3.f90` paths until this test identifies the actual convention.
4. **Field metadata:** replace inferred `scalar_01..scalar_11` meanings with a
   versioned native map tied to `NVAR`, `imetal`, `ichem`, and the executable
   build. Prove He indexing, total-metal treatment, and delayed-cooling field
   ownership. Reject an unrecognized map at startup.

#### G1.2 Interpolation and table API

- Replace `find_bounds` endpoint clamping with explicit out-of-range errors in
  production mode; permit boundary equality only when it is an exact table
  node and record it.
- Preflight complete Cartesian corners for every enabled channel, detect
  duplicate coordinates, non-finite values, negative increments, and
  non-monotonic cumulative fields.
- Make interpolation policy explicit per axis and channel; no silent
  nearest-neighbor substitution or extrapolation.
- Keep `net_*` as diagnostics only. Gas deposition uses actual ejecta and the
  declared returned-mass field.
- Preserve reduced-chemistry mass as
  `untracked_ejecta = returned_mass - sum(ejecta_H..Fe)`. Deposit this residual
  into generic metallicity only, never a tracked element, and require the
  cumulative residual itself to remain non-negative and monotonic.
- Add native Fortran tests and Python/JAX differential tests on the same table,
  including boundaries and exact-grid queries.

**Exit:** G1 has a passing native contract test matrix and a production binary
that cannot start with an ambiguous field map, unit contract, or query outside
its declared domain.

### G2 — build and approve the physical yield model

#### G2.1 Source-selection review

Create a source matrix before selecting values. For every candidate record
mass range, metallicity grid, remnant/fallback prescription, mass-loss model,
explosion energy, release-time resolution, binary assumptions, IMF convention,
uncertainties, license, and citation.

The review must separately select and cross-calibrate:

- massive-star winds;
- AGB winds and thermal-pulse ejecta;
- core-collapse/SNII ejecta, remnant mass, fallback, energy, and momentum;
- SNIa event yields and event-rate model;
- PISN/PPISN only for a declared Pop-III or very-low-metallicity population.

Candidate literature families include modern CCSN grids such as [Jost et al.
2025](https://academic.oup.com/mnras/article/536/3/2135/7920781) and
[Boccioli & Roberti 2026](https://doi.org/10.1051/0004-6361/202557714), AGB grids
such as [Huscher et al. 2025](https://doi.org/10.3847/1538-4357/ae0199),
[Karakas & Lugaro 2016](https://doi.org/10.3847/0004-637X/825/1/26), or
Cristallo et al., and SNIa yields such as [Keegans et al.
2023](https://arxiv.org/abs/2306.12885). These are candidates for review, not
automatic approvals.

#### G2.2 Population and channel conventions

1. Select one IMF as the default and make it data/configuration driven. Remove
   the hard-coded Kroupa assumption from the runtime or prove that it matches
   the selected population-synthesis table.
   Population-weighted source tables must pass a physical-time integral and
   unit-mass-return test and must never receive a second IMF convolution.
2. Define non-overlapping mass/channel boundaries and lifetime mapping. Winds,
   AGB, and SNII must not release the same mass twice.
3. Define cumulative actual ejecta, living mass, remnant mass, energy, and
   momentum per initial stellar mass. Require
   `initial = living + remnant + returned` within a declared tolerance.
4. Generate a full mass--metallicity--age grid, including an age-zero anchor,
   all runtime mass ranges, sufficient corners for interpolation, and an
   explicit maximum age. One-point fixtures are retained only as tests.
5. Decide the chemistry scope. Either add the species needed for a full
   nucleosynthesis claim (at minimum assess Na, Al, Ar, Ni and any paper-critical
   species), or formally publish an 11-element reduced-chemistry model with
   omitted-species closure and sensitivity tests.

#### G2.3 Converter, closure, and provenance

- Write a deterministic converter from each approved source format into the
  canonical table; store the original file hashes and conversion code hash.
- Convert rates to cumulative quantities only once and test numerical
  integration against analytic release histories.
- Check row, grid, channel, population, mass, elemental, energy, and momentum
  closures in both source units and RAMSES units.
- Add a signed/approved sidecar containing citations, versions, units,
  assumptions, IMF, channel boundaries, metallicity definition, solar
  abundance set, remnant model, and checksums.
- Require the sidecar approval identifier in the batch gate.

**Exit:** `audit_stellar_yield_asset.py` and the native audit pass on the named
asset, all required channels are approved, and a channel-by-channel native
run reproduces the independent ledger within tolerance.

### G3 — event rates, SEDs, and source ledgers

#### G3.1 Stellar SED and photon ledger

1. Select a population-synthesis SED consistent with the approved IMF,
   metallicity grid, age convention, and binary/single-star assumption. The
   current BPASS candidate remains unapproved until its low-Z/young-age
   clamping, normalization per initial mass, and `fesc=1` assumptions are
   replaced or explicitly justified.
2. Extend or reject outside the SED domain; never silently floor metallicity
   or age. Integrate the intrinsic spectrum into the single pinned photon
   group table and store photon-number and energy closures.
3. Choose a physical birth-cloud/escape prescription, separate intrinsic
   emission from escaped emission, and report its uncertainty. Do not use
   `fesc=1` as a science default without an obscuration model.
4. Merge stellar and AGN ledgers only after coevality, source IDs, group edges,
   positions, normalization, and total emitted photons pass the strict merger.

#### G3.2 SNIa

- Implement SNIa as a separate event-rate convolution over the SSP age,
  metallicity, and DTD; it is not an ordinary 3--8 Msun yield channel.
- Select DTD shape, minimum delay, normalization, progenitor model, metallicity
  dependence, ejecta composition, event energy, and stochastic versus
  expectation-value treatment.
- Test single-event and SSP-integrated rates, restart continuity, and mass/
  energy closure. Keep the channel disabled until all inputs are approved.

#### G3.3 PISN/PPISN

- Keep PISN disabled for ordinary metal-enriched populations.
- If the science case needs it, implement an explicit population eligibility
  gate using metallicity, Pop-III status, core mass, PPISN outcome, and the
  selected stellar-evolution model. Do not encode a universal ZAMS 140--260
  Msun rule.
- Validate event rates and yields against the selected source family and test
  that PISN cannot be activated by an accidental mass interpolation.

#### G3.4 AGN ledger

- Use the pre-reset coarse-step diagnostic and record the exact interval,
  inflow versus retained-BH-mass convention, radiative efficiency, `Lbol`, SED,
  obscuration, escape fraction, thermal/jet energy, and momentum.
- Deduplicate `(coarse_step, sink_id)`, reject zero/multiple matches, and make
  restart replay bitwise or tolerance-reproducible.
- Prove that AGN mass, radiative energy, jet energy, radiation momentum, and
  deposited gas changes close without double counting.
- Declare whether the production product is static post-processing or live
  RAMSES coupling; a static ledger must not be described as live feedback.

**Exit:** approved STAR, SNIa/optional-PISN, and AGN ledgers use one group table,
one epoch convention, one source-ID namespace, and close independently and
after merging.

### G4 — dust, radiation transport, and thermochemistry

#### G4.1 Dust physics

1. Approve a dust mixture and opacity source/version. The current Draine/WD01
   sidecar is a technical candidate, not a complete astrophysical prescription.
2. Select a cell dust-to-metal relation, depletion pattern, metallicity floor,
   grain-size assumptions, and redshift/environment dependence. Store these in
   the snapshot and manifest, not as an implicit zero or constant.
3. Implement and separately ledger absorption, scattering, IR re-emission, and
   radiation pressure. If any process is excluded from the paper, state the
   limitation and test that its omission is not being counted as zero physical
   opacity by accident.
4. Ensure the source SED and opacity use the same group edges and spectral
   weighting. Dust photon absorption and heating must be excluded from H/He
   primary-absorption closure exactly once.

#### G4.2 RT and thermal closure

- Extend the thermal atlas to bracket every production epoch; high-edge
  clamping is a pilot-only behavior.
- Couple photoheating, dust heating, background heating/cooling, chemistry,
  and gas energy exactly once. Keep diagnostic heating separate until the live
  hydro contract is active.
- Establish source-cell spatial convergence, not only volume-mean convergence;
  replace point-source controls with an audited physical source-size or
  birth-cloud model.
- Run angular, spatial, timestep, reduced-light-speed, and thermal-subcycle
  convergence matrices with predeclared thresholds.
- Validate H/He/H2 and metal chemistry ranges, temperature bounds, positivity,
  photon conservation, and energy/momentum conservation in float64 reference
  runs and the target JAX precision.

**Exit:** the chosen dusty RT/thermal model has an approved opacity/chemistry
sidecar, no atlas or source-domain bound hits, converged source cells, and
closed photon/thermal/dust ledgers.

### G5 — live RAMSES coupling and HPC qualification

1. Build from a clean, named source commit with the selected compile-time
   modules. The executable must print its field-map, table hash, source mode,
   and feedback-channel status at startup.
2. Remove the embedded Phase-0 fallback from the production build, or make it
   compile-time unavailable. Require an explicit approved table and sidecar.
3. Implement atomic source-progress and ledger updates across AMR levels,
   MPI ranks, subcycling, and restart. Test ownership and ghost/boundary
   deposition so a source is applied once.
4. Define the live timestep/subcycling contract between RAMSES and SNRT. The
   source-cell limiter must recompute when luminosity, gas state, opacity, or
   refinement changes.
5. Verify that delayed cooling, thermal energy, kinetic energy, radiation
   pressure, dust force, and AGN jet terms are mutually exclusive or explicitly
   additive with a closed ledger.
6. Run a small native hydro gate, then a multi-rank AMR gate, then a restart
   equivalence test before any production volume. Record memory, I/O, wall
   time, and failure diagnostics.

**Exit:** live coupling passes conservation, restart, AMR/MPI ownership, and
performance gates on the exact production binary.

### G6 — science qualification and baseline comparison

- Re-run the registered transitional feedback baseline with its exact embedded
  build identity and compare it against the approved production model at the
  same initial state, output time, resolution, and measurement definition.
- Use controlled comparisons that change one ingredient at a time: yield set,
  IMF/binary model, SNIa DTD, AGN SED/obscuration, dust model, RT resolution,
  and live versus static coupling.
- Define primary observables before looking at results: ionization/thermal
  structure, dust heating, radiation force, metal distribution, stellar/AGN
  photon escape, and JWST-facing synthetic observables.
- Quantify numerical and model uncertainties; distinguish them from the
  legacy-to-production physical difference.
- Require independent replay from the manifest and an external final audit of
  algorithm, wiring, and physical justification.

**Exit:** all paper observables have a convergence statement, sensitivity
range, provenance, and reproducible baseline comparison.

### G7 — publication package

Prepare one release directory containing:

- source and executable commit IDs, compiler/JAX environment, and manifest;
- approved yield/SED/dust/AGN sidecars and all checksums/licenses;
- effective namelists and launch scripts;
- unit, closure, convergence, and restart reports;
- baseline and production outputs or certified reduction products;
- scripts that regenerate every paper table/figure;
- a limitations section covering reduced chemistry, omitted channels, dust
  assumptions, source resolution, and static/live scope.

The final paper wording must not call the 11-element table “complete
nucleosynthesis” unless the chemistry-expansion gate has passed.

## 5. Execution order and decision points

1. **Complete:** G0/G1 — identities and native/JAX contracts are frozen and
   audited.
2. **Now:** G2 — select sources and generate the first fully cited table; this
   is the only point at which a scientific yield asset may be promoted.
3. **Then:** G3 — implement approved SED, SNIa, AGN ledgers, and optional PISN
   gating; merge only coeval, same-grid ledgers.
4. **Then:** G4 — promote dust and thermal/RT physics after spatial and temporal
   convergence, not merely conservation.
5. **Finally:** G5/G6/G7 — live coupling, HPC qualification, uncertainty,
   baseline comparison, external audit, and release.

At each decision point, a failing gate produces a named `blocked` artifact and
the previous certified baseline remains untouched. No large RAMSES rerun is
authorized by this plan before G2--G5 are green.

## 6. Immediate deliverables

- `provenance/production_publication_readiness_plan.md` (this plan);
- a source-selection matrix and approval record for wind/AGB/SNII;
- native reproductions for age units, interval accounting, SN-energy units,
  and field metadata;
- a production-side table manifest with full-grid coverage and checksums;
- a native/JAX differential and conservation matrix;
- a paper-scope decision for 11-element reduced chemistry versus an expanded
  nucleosynthesis set.

The existing detailed implementation notes remain active in
[`feedback_implementation_plan.md`](feedback_implementation_plan.md), while
the three-model audit is preserved in
[`yield_table_multi_model_audit_2026-09-01.md`](yield_table_multi_model_audit_2026-09-01.md).
