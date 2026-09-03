# F-P1 40–120 M☉ required-data comparison

Date: 2026-09-03
Scope: data and evidence needed to promote the massive-terminal-fate seam to production/publication readiness
Auditors: AGY (`gemini-3.8-flash-high`) and Claude Opus 5

Source records:

- [AGY opinion](agy_fp1_high_mass_required_data_2026-09-03.md)
- [Claude Opus 5 opinion](claude_opus5_fp1_high_mass_required_data_2026-09-03.md)
- [local high-mass review artifact](../simulation/snrt/data/fp1_high_mass_seam_review.json)

## Joint verdict

| Area | AGY | Opus | Adopted disposition |
|---|---|---|---|
| Resolver/fail-closed scaffolding | Conditional pass | Conditional pass | Keep active as review-only infrastructure |
| Scientific 40–120 M☉ fate law | Block | Block | No production fate map or physical rows |
| Production table emission/runtime activation | Block | Block | Keep disabled |
| Current W18/N20 evidence | Comparison evidence | Anti-monotone comparison evidence | Validation branch only |

Both auditors independently reject a mass-only direct-collapse boundary. The staged W18/N20 records contain both positive- and non-positive-energy outcomes in the same mass interval, and the engines disagree at some nominal masses. This evidence demonstrates why the current unresolved seam is necessary; it does not select a source or a fate prescription.

## Common minimum data package

The following is the intersection of the two audits and is the minimum production admission package.

1. **Identity, provenance, and rights**: source/version, article and data DOI, archive URL, immutable file bytes/SHA256, package fingerprint, retrieval date, license, research-use status, redistribution permission, citation/acknowledgement requirements, converter version/hash, and approval id.
2. **All resolver coordinates**: ZAMS mass and non-overlapping mass cell with the declared half-open edge rule; metallicity value and definition; solar abundance set; rotation or approved population marginalization; binary state/distribution or explicit single-star baseline; engine/branch id; and cell assignment rule.
3. **Pre-SN structure and fate**: pre-SN mass, He/CO/Fe core masses, compactness and definition, μ4/M4 if used, classifier/version/calibration/threshold provenance, explicit outcome, and pair-instability criterion/pulse/disruption fields.
4. **Time-dependent release**: lifetime and its definition/source, age-zero convention, ordered cumulative wind mass and 11-element/untracked composition versus age, or a formally approved terminal-lumping approximation with a quantified dynamical error budget.
5. **Terminal and remnant ownership**: terminal ejecta mass/composition/untracked residual; explicit zero-because-direct-collapse semantics; fallback; baryonic and gravitational remnant or a declared conversion convention; remnant type and owner; and wind/terminal/PISN double-counting rules.
6. **Decay**: complete or explicitly typed isotope inventory, decay projection/horizon, decay-data fingerprint, cross-segment duplicate resolution (including K-40), missing-nuclide policy, returned-mass closure, and rest-mass-loss treatment.
7. **Energy and momentum**: typed energy quantity (asymptotic kinetic vs diagnostic vs injected), binding-energy mapping where needed, approved injected-energy contract, source-frame vector convention, scalar launch-momentum convention, coupling/deposition mode, and advective policy. Momentum must not be inferred from energy.
8. **Reproducible admission evidence**: converter output, source-node fingerprints, closure results and tolerances, complete hull/cell coverage, resolver edge/out-of-hull behavior, Fortran/JSON mirror, sidecar hashes, and runtime digest binding.

Missing physical quantities must stay `null`; an absent record is not a zero record. In particular, a failed/direct-collapse node must be an explicit wind-only node with terminal ejecta zero and a remnant once source data support that claim.

## Coverage and interpolation policy

Both audits require source-node cells to tile [40,120] M☉ with no gaps or overlaps, and fail-closed behavior outside the full approved hull. Mass interpolation across different outcomes, cross-engine interpolation, and cross-source interpolation are prohibited. Only monotone cumulative wind interpolation along the age axis of one assigned node is admissible.

AGY’s minimum grid recommendation is at least four metallicity tracks spanning approximately Z=10⁻⁶–0.02 and at least three rotation values (0/150/300 km s⁻¹), or an explicitly documented weighted marginalization. Opus adds that the grid must reach the target population’s low-metallicity domain; endpoint clamping and solar extrapolation are forbidden. The exact production grid remains a design/physics decision, but the present solar-only, non-rotating Sukhbold branch cannot satisfy it.

## Closure and test set to adopt

AGY provides seven core identities; Opus expands them into a more operational test matrix. Adopt the union, with tolerances tied to source precision rather than one universal number:

- per-node progenitor mass closure, using a demonstrated rounded-table tolerance (Opus’s Sukhbold analysis supports about 7×10⁻³ relative) and a tighter full-precision tolerance (about 10⁻³ or the release-demonstrated bound);
- direct-collapse remnant identity and explicit derived-value uncertainty where applicable;
- non-negative tracked/untracked returned-mass residual and no tracked over-return;
- complete decay/baryon closure and duplicate-isotope rejection;
- wind monotonicity, age anchoring, telescoping, terminal-time constancy, and restart/time-subdivision reproducibility;
- exactly-once terminal energy/event addition and typed energy semantics;
- source-frame vector momentum closure; scalar momentum remains null without a deposition contract; no code path computes momentum from energy;
- population/IMF ledger closure;
- resolver query tests at node centres, both edges, and ±ε outside every axis, plus independent mutation of every resolver key;
- source, contract, converter, sidecar, Fortran/JSON mirror, and deployed-runtime digest re-hashing;
- same-progenitor wind consistency across engine branches and named handling of any source-confirmed discrepancy;
- cross-source wind consistency as an expected blocking test until the source semantics are reconciled.

The tolerances need a source-specific justification. Printed Sukhbold values must not be judged by the same absolute cap as full-precision Zenodo tables, and runtime numerical tolerances must remain separate from source-table rounding tolerances.

## Additional findings from Opus

Opus identified four current wiring risks that AGY did not call out explicitly:

- the Fortran approval token currently validates shape/non-emptiness, not the digest of the deployed map; a production implementation needs a compiled-in expected digest or startup hash comparison;
- the 32-field canonical row does not carry all resolver axes, so rotation/engine/criterion information could be silently dropped unless the schema is extended or explicit marginalization is recorded;
- the current SNII channel window `[8,40]` excludes the exploding evidence nodes above 40 M☉ and therefore needs an explicit ownership/domain decision;
- the canonical contract has no home for scalar radial launch momentum, so a deposition contract is required before it can be emitted.

Opus also flagged source-level checks for branch wind consistency and radioisotope reference-epoch consistency. These are to be added to the high-mass admission test set, not resolved by averaging or by inventing a decay epoch.

## Candidate disposition

- **Sukhbold W18/N20**: retain as a solar, non-rotating, review-only engine-comparison/validation branch. It is strong evidence against mass-only fate assignment but is not promotable as the production backbone.
- **Boccioli & Roberti 2026 / LC18-style release**: best long-term structural candidate because of CC BY 4.0, multi-Z, and multi-rotation coverage, but quarantined until the failed-model Wind-table anomaly is resolved with the authors or a corrected release is obtained. Energy, momentum, age-history, lifetime, and high-mass sampling still need to be made machine-readable/adequate.
- **Limongi & Chieffi 2018 CDS**: lifetime and wind-semantics comparison branch only; its wind-only failed-model scenario is not an approved fate law and its licensing/reconciliation remain open.

## Implementation priority derived from both audits

1. Obtain a redistributable, corrected source package and author clarification for failed-model winds.
2. Freeze and validate the per-node schema, including all resolver axes, explicit null/zero semantics, fate classifier, ownership, and energy/momentum types.
3. Build the complete source-hull/cell admission and deterministic converter checks, including the runtime digest binding and channel ownership decision.
4. Add age-resolved wind/lifetime data or approve terminal lumping with a quantitative error budget.
5. Add remnant/fallback, decay, PPISN/PISN, energy, momentum, and deposition closure evidence.
6. Populate physical nodes only after all blocking tests pass, regenerate the cryptographic sidecar, and request a bundled gate audit.

The local high-mass artifact and tests remain review-only. No candidate, physical value, fate law, or runtime activation was changed by this comparison.
