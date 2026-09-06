# Claude Opus 5 opinion — F-P1 40–120 M☉ required data

Date: 2026-09-03
Model: `claude-opus-5` via Claude CLI
Mode: read-only plan review of `/gpfs/kjhan/LRD_JWST`

## Verdict

Opus gives `CONDITIONAL PASS` for the existing fail-closed contract, resolver architecture, and runtime guards, and `BLOCK` for production table emission/runtime activation. The scientific verdict is an unqualified `BLOCK`.

The staged W18/N20 records are anti-monotone comparison evidence, not a project fate law. W18 has positive asymptotic explosion energy at 60 and 120 M☉ and non-positive outcomes at the other staged nodes; N20 is positive at 60, 80, 100, and 120 M☉ and non-positive at 40, 45, 50, 55, and 70 M☉. The same nominal mass can therefore have opposite outcomes under different engines.

## Engineering conditions before promotion

1. The runtime admission token is self-asserted: the Fortran path checks a policy string, 64 hex characters, and a non-empty approval id, but does not hash the deployed map or compare against an allow-list. Compile in the approved digest or hash the deployed map at startup and stop cleanly on mismatch.
2. The 32-field canonical row has only channel/mass/Z/age coordinates while the resolver requires twelve axes, including rotation, engine/branch, lifetime source, pair-instability criterion, and abundance set. Extend the schema or require an explicit, recorded marginalization/freeze rule; never silently drop an axis.
3. A failed/direct-collapse node must be an explicit wind-only row with zero terminal ejecta and a remnant, not a missing row that invites interpolation.
4. The current SNII channel window `[8,40]` excludes the exploding evidence nodes at 60–120 M☉. The channel domain/ownership contract must be extended or explicitly split before those nodes can enter the canonical table.
5. The canonical contract has no valid home for scalar radial launch momentum. Keep the vector firewall, but add an explicit deposition/momentum contract before emitting a non-null scalar.

## Minimum per-node data

Every node needs:

- identity and rights: source/version, article/data DOI, archive URL, per-file size and SHA256, package fingerprint, retrieval date, license, research-use and redistribution status, citation/acknowledgement requirements, converter version/hash, and approval id;
- all resolver coordinates: ZAMS mass and half-open mass cell, metallicity value and definition, abundance set, rotation or declared marginalization, binary state/distribution (with binary prescription if applicable), engine/branch id, and cell-assignment rule;
- pre-SN structure and fate classifier: total pre-SN mass, He/CO/Fe core masses, compactness ξ2.5 with definition, μ4/M4 if used, classifier/version/calibration/threshold provenance, and an explicit outcome;
- timing: lifetime, lifetime definition and source, age-zero anchor, and either ordered age-resolved cumulative wind data or an approved terminal-lumping approximation with a quantified dynamical error budget;
- wind: monotone cumulative mass, 11 tracked elements, untracked residual, and release ages;
- terminal: ejecta mass/elements/untracked residual plus component reference, ownership, and an explicit zero-because-direct-collapse flag;
- fallback/remnant: fallback, baryonic and gravitational remnant or a declared neutrino-loss convention, remnant type and owner;
- PPISN/PISN: core-based criterion, criterion input, pulse history with age/mass/composition/energy, final remnant, and explicit complete-disruption confirmation;
- decay: projection/horizon, raw isotope count, completeness status, duplicate-isotope resolution, decay-data fingerprint, missing-nuclide policy, closure residual, and rest-mass-loss treatment;
- energy: energy kind (asymptotic kinetic, diagnostic, injected, or central-engine deposited), final/diagnostic/thermal energy, binding energy if mapping is derived, and an approved injected-energy mapping id;
- momentum/deposition: source-frame vector, scalar launch convention, deposition contract, coupling/advective policy, and residual policy. Momentum must never be derived from energy in the converter.

Missing values must remain `null`; they cannot be inferred from integrated yields. A zero failed-model energy is an outcome and must be typed as such.

## Coverage and interpolation rules

Nine mass nodes at one solar, non-rotating point are not coverage of the twelve-dimensional resolver key. The approved cells must tile [40,120] with no gap/overlap and must fail closed outside the hull in mass, Z, rotation, and binary coordinates. Endpoint clamping is forbidden. At promotion, add a completeness check for the cell tiling.

Mass interpolation across different outcomes is forbidden because it manufactures a fractional explosion between collapse and explosion nodes. Interpolation is allowed only along one node’s own monotone age-resolved wind history. Cross-source and cross-engine interpolation is forbidden unless a named, approved population marginalization is introduced and independently tested.

AGY’s minimum of four metallicity tracks spanning roughly 10⁻⁶–0.02 and three rotation values (or a formally convolved distribution) is retained. Opus additionally notes that the staged positive-Z hull is above the target catalogue: the minimum demanded metallicity must be reached rather than clamped or solar-extrapolated.

## Closure and reproducibility tests

Opus recommends:

- per-node mass closure `|wind + terminal + remnant − ZAMS|/ZAMS ≤ 7×10⁻³` for printed Sukhbold values, or `≤10⁻³` for full precision; use the release-appropriate bound and report signed residuals;
- direct-collapse remnant identity `remnant = ZAMS − wind`, marked as derived with its uncertainty if no source remnant is supplied;
- exact wind agreement across engine branches for the same progenitor, or a named source erratum; do not average a discrepancy;
- radioactive wind epoch consistency before using isotope values;
- duplicate-isotope rejection (notably K-40), decay/baryon closure, and tracked/untracked non-negative residual closure;
- explicit energy typing and exactly-once terminal addition;
- zero source-frame vector momentum, scalar momentum null unless a deposition contract exists, and a source assertion that no code path computes momentum from energy;
- age anchor/monotonicity/telescoping, interval-boundary restart reproducibility, and population ledger closure;
- resolver queries at node centres, both edges, ±ε outside every hull axis, and independent mutation of all twelve axes;
- source-file, contract, converter, sidecar, Fortran/JSON mirror, and deployed-runtime digest re-hashing.

Opus notes that the current printed Sukhbold residual pattern supports a seam-specific relative tolerance around `7×10⁻³`, while full-precision releases should use their own demonstrated tolerance. Runtime and analytic IMF tolerances remain separate from source-table rounding tolerances.

## Blocking versus later

Blocking physics: approved fate law; multi-Z and rotation coverage; target-hull coverage; age-resolved winds or an error-budgeted lumping approximation; failed-model remnants/fallback; complete decay and isotope epoch semantics; injected-energy mapping; momentum/deposition; PPISN/PISN classification; and cross-source wind consistency.

Blocking engineering/provenance: runtime digest binding; schema/axis completeness; ragged failed-node semantics; channel-window ownership; scalar momentum home; node/cell completeness; seam-specific acceptance bound; N20 wind coverage; redistribution rights; decay-data rights; and the known LC18 failed-model wind inconsistency.

Later sensitivity/publication: W18/N20/GR1D engine uncertainty after one baseline exists, binary SSP distributions, reduced-chemistry sensitivity, terminal-lumping sensitivity after its error budget is approved, a second structure-based classifier, Pop-III/PISN science branch, and IMF-shape sensitivity.

## Candidate disposition

Sukhbold (2016) should remain a solar-metallicity, non-rotating review/validation branch. It has immutable archive fingerprints, nine high-mass nodes, parsed W18/N20 outcomes, component separation, six terminal-yield records, thirteen wind records, and useful closure evidence. It lacks Z/rotation axes, ages/lifetime, age-resolved wind, failed-node remnants, complete decay, radioisotope epoch consistency, injected-energy mapping, canonical momentum, N20 wind coverage at several failed nodes, redistribution permission, and an applied seam acceptance bound.

Boccioli & Roberti 2026/LC18 is the best long-term structural path because it is CC BY 4.0 and has four metallicities and three rotations. It is blocked by the release anomaly in which all 56 failed-model summaries report nonzero winds while their Wind tables are zero, by sparse 40/60/80/120 M☉ sampling, and by absent machine-readable energy, momentum, age history, and lifetime. Querying the authors about the failed-model Wind tables is the shortest credible next step.

The Limongi & Chieffi 2018 CDS branch is useful for lifetime and wind-semantic comparison, but its wind-only failed-model treatment is a scenario, not an approved fate law. Its redistribution terms and internal wind reconciliation also require resolution.

No physical value should be invented, no momentum inferred from energy, no lifetime inferred from integrated yields, and no engine comparison converted into a mass-only direct-collapse rule.
