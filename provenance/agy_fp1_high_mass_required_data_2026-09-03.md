# AGY opinion — F-P1 40–120 M☉ required data

Date: 2026-09-03
Model: `gemini-3.8-flash-high` via AGY
Mode: read-only analysis of `/gpfs/kjhan/LRD_JWST`

## Verdict

AGY gives a `CONDITIONAL PASS` for the fail-closed contract, resolver architecture, and runtime guards, but a `BLOCK` for production table emission and runtime activation. The scientific verdict is `BLOCK`.

The reason is that the staged W18/N20 evidence is useful comparison data, not an approved production fate law. Both engine branches contain positive- and non-positive-energy outcomes within 40–120 M☉, so a mass-only direct-collapse boundary is unsafe.

## Minimum per-node data

AGY requires the following for every physical source node:

- immutable source identity: source/version, article DOI, data DOI, archive URL, file SHA256, license and redistribution permission;
- coordinates: ZAMS mass and mass cell, explicit edge convention, metallicity definition and abundance set, rotation or declared marginalization, binary state or marginalization, and engine/branch ID;
- pre-supernova structure: compactness ξ2.5, μ4, C/O-core, He-core, and iron-core masses;
- fate: explicit outcome (successful CCSN with fallback, direct collapse variants, PPISN, PISN, or wind-only), explodability criterion, PPISN pulse count, and PISN disruption flag;
- timing: lifetime, lifetime source, and either an age-resolved cumulative wind history or an explicitly approved terminal-lumped approximation;
- wind: cumulative mass, all tracked-element yields, untracked residual, and release-age grid;
- terminal ejecta: total mass, tracked-element yields, and untracked residual;
- decay: decay projection, horizon, raw isotope count, duplicate-isotope resolution, and closure residual;
- remnant/ownership: baryonic and gravitational remnant masses, fallback, terminal owner channel, and a double-counting guard;
- energetics: final kinetic, diagnostic explosion, and thermal energies;
- momentum/deposition: source-frame vector, scalar launch momentum convention, runtime momentum behavior, energy mode, deposition contract, and untracked-residual policy.

## Coverage and interpolation

The approved domain must partition all of [40, 120] M☉ with source-node cells. AGY recommends at least four metallicity tracks spanning Z≈10⁻⁶–0.02 and at least three rotation values (0, 150, 300 km s⁻¹), or a formally convolved population distribution.

Mass interpolation across different outcomes is forbidden. The recommended rule is piecewise-constant source-node assignment. Continuous interpolation is allowed only along a single node’s monotone cumulative wind-age axis. Cross-source and cross-engine interpolation is forbidden. Queries outside the approved convex hull must fail closed; endpoint clamping is forbidden.

## Quantitative closure/reproducibility tests

AGY proposes these checks:

1. progenitor mass closure, relative tolerance 10⁻³ for rounded tables and 10⁻⁵ for native outputs;
2. tracked ejecta must not exceed returned mass, with residual no less than −10⁻¹² M☉;
3. untracked residual must be non-negative;
4. decay/baryon closure within 10⁻⁴ relative;
5. source-frame momentum vector must be zero within 10⁻¹² g cm s⁻¹;
6. cumulative mass/energy release must be monotone;
7. population IMF ledger closure within 10⁻⁷ relative to initial stellar mass.

It also calls for immutable input fingerprints, complete grid/hull tests, restart/time-subdivision reproducibility, and converter/resolver differential checks.

## Blocking versus later items

AGY classifies as audit-blocking physics/provenance: multi-Z and rotation coverage, per-node fate classifier, age-resolved winds or approved terminal lumping, complete decay and closure, event energy/momentum/deposition semantics, terminal remnant ownership, PPISN/PISN classification, redistribution rights, populated resolver nodes, deterministic conversion, and a hash-coupled approval sidecar.

Binary SSP distributions, reduced-chemistry sensitivity, and W18/N20 uncertainty quantification are classified as later sensitivity/publication work once a baseline is approved.

## Candidate assessment

Sukhbold (2016) has verified archive fingerprints, 9 high-mass nodes, W18/N20 outcome diversity, separated wind/explosive components, and useful exploding-model yield evidence. It remains blocked by solar-only, non-rotating coverage; incomplete/nonuniform high-mass wind records; incomplete decay inventory; no age-resolved winds; no canonical momentum; and redistribution restrictions.

LC18/Boccioli & Roberti has CC BY 4.0 licensing, four metallicities, three rotation values, modern GR1D+SkyNet physics, explicit component semantics, and good closure on successful explosions. It remains quarantined because failed-model summaries report nonzero winds while corresponding wind tables are zero, machine-readable energy and momentum are absent, high-mass sampling is sparse (40, 60, 80, 120 M☉), and the 8–13 M☉ transition is uncovered.

## Promotion prerequisites

AGY’s recommended promotion sequence is: acquire or curate a clean multi-Z/multi-rotation grid; ingest machine-readable asymptotic kinetic energy and an explicit momentum convention; formalize age-resolved winds or approve terminal lumping with a dynamical error budget; populate the resolver with physical source nodes; then regenerate the SHA256 sidecar with a non-null approval ID and change the fate map to approved.
