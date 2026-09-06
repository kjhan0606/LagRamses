# F-P1: 40--120 M☉ stellar-fate literature dossier

Date: 2026-09-02  
Scope: production/publication-ready stellar population, RT source, and
feedback wiring for the high-mass interval.  This dossier is a review input;
it does not approve a yield source or emit canonical rows.

## Executive finding

The literature does not support a universal rule of the form
`40 <= M_ZAMS/M☉ <= 120 -> direct collapse`.  The outcome is model-dependent
and is controlled by at least metallicity, rotation, wind history, the
pre-supernova core structure, the explosion engine, fallback, and (at higher
core masses) pair instability.  The production representation must therefore
be a source-node/structure-based fate resolver, not a single ZAMS mass bin.

The current project map consequently remains fail-closed.  Its 40--120 M☉
entry means “not yet assigned a reviewed model”, not “all direct collapse”.

## Candidate strategies used in other work

| Strategy | What it resolves | What it would mean for SNRT | Current disposition |
|---|---|---|---|
| Sukhbold et al. 2016 | A solar-metallicity 9--120 M☉ progenitor grid with several calibrated neutrino-engine branches; per-model explosion, fallback, and remnant quantities | Look up fate and terminal quantities at a source progenitor node. Keep presupernova wind ownership separate from terminal ejecta. | Best direct evidence against a mass-only bin; candidate source, not approved. The staged archive still lacks an age-resolved wind history, complete decay inventory, canonical momentum, and a project redistribution approval. |
| Limongi & Chieffi 2018 | 13--120 M☉, four metallicities, three initial rotations, presupernova evolution and explosive yields | Use explicit `(M_ZAMS, [Fe/H], rotation, yield-set)` coordinates. Their set-R convention treats high-mass models above 25 M☉ as wind-only under a stated direct-BH scenario. | Useful comparison branch and wind/ejecta semantic reference, but set-R is a scenario, not a universal fate law; no continuous age-resolved release history or project event-momentum contract. |
| Fryer et al. 2012 | Rapid/delayed analytic remnant prescriptions with fallback and metallicity dependence | Fast population-level remnant sampler when detailed progenitor structure is unavailable; must be tagged as a population prescription, not a per-star yield source. | Possible uncertainty/sensitivity branch only. It has no PPISN treatment and explicitly leaves high-mass wind loss/metallicity as a major uncertainty. |
| Patton & Sukhbold 2020 / Ertl-style criterion | Explodability islands from pre-SN structural quantities such as the iron-core/mass-gradient diagnostics | Classify a progenitor using pre-SN structure supplied by the selected stellar grid; do not infer the class from ZAMS mass alone. | Recommended algorithmic shape for a future resolver, but threshold calibration is source/engine dependent and the current canonical rows do not contain the required structural coordinates. |
| Ugolini et al. 2025 | HYPERION explosions on the LC18 grid, with metallicity/rotation dependence, fallback, and an explicit direct-collapse criterion based on envelope ejection; PPISN/PISN handling through CO-core criteria | A modern grid-based branch can provide outcome and remnant mass at sampled `(Z, rotation, M_ZAMS)` nodes, with out-of-hull queries rejected. | Strong candidate for a future multi-metallicity fate branch; not approved until its source tables, wind ownership, event quantities, and table conversion are independently fingerprinted. |
| Woosley 2017; Heger & Woosley 2002 | PPISN/PISN as a function of helium/CO-core evolution, with pulses, partial ejection, complete disruption, or collapse | Add a separate core-based P-F3 resolver. A PPISN pulse is not an ordinary SNII row, and a PISN has no remnant owner. | Must be represented explicitly and separately; never inferred from a fixed 40--120 M☉ ZAMS interval. |

## Independent reproduction from the staged Sukhbold archive

The local archive was read without modifying or converting it.  The high-mass
tail of its published explosion-result files is non-monotonic:

* W18: 40, 45, 50, 55 M☉ fail in the listed diagnostic outcome; 60 M☉
  explodes at 0.65 foe; 70, 80, and 100 M☉ fail; 120 M☉ explodes at 0.67 foe.
* N20: 40, 45, 50, 55, and 70 M☉ fail; 60, 80, 100, and 120 M☉ explode at
  0.92, 0.54, 0.62, and 1.01 foe respectively.

These values are a local reproduction of the staged source files, not a claim
that either engine is the project truth.  They are sufficient to reject a
monotonic ZAMS-mass partition.  The local candidate audit still correctly
keeps the source out of production because it does not provide all required
SNRT fields and provenance permissions.

## Physical and implementation requirements extracted from the comparison

### Fate key and lookup semantics

The canonical resolver must accept an explicit source/model coordinate:

```text
(source_id, M_ZAMS, birth_metallicity, initial_rotation, optional_binary_state)
```

It must return a terminal outcome only when that coordinate is inside the
approved source hull.  Nearest-node substitution, cross-source interpolation,
metallicity/rotation extrapolation, and a mass-only fallback are rejected by
default.  If a population model intentionally marginalizes over rotation or
binary state, that marginalization must be declared and recorded as a model
choice rather than hidden in the lookup.

### Separate source components

The wind component is cumulative and age-resolved when used as a time source.
The terminal component is added once, after the terminal outcome is known.
Failed or direct-collapse models may have wind return but no terminal ejecta;
their remnant mass is still a required closure quantity.  A successful
fallback model must not count pre-SN wind material again as terminal ejecta.

### Pair-instability branch

PPISN/PISN classification is core-based and source-dependent.  PPISN requires
pulse timing, ejected mass/composition, and the final remnant semantics if it
is to drive time-dependent feedback.  PISN requires a complete-disruption
record and explicitly no remnant-owner channel.  If only a classification is
available, the model may be used for a review/sensitivity branch but not as a
canonical energy, momentum, or yield row.

### Required canonical fields

For each approved source node the converter must account independently for
initial mass, age/lifetime, returned and remnant mass, tracked and untracked
ejecta, energy, momentum/deposition semantics, isotope-to-element/decay
projection, source metallicity and rotation, and provenance/license.  Missing
fields remain missing; they must not be synthesized from a figure, a fixed
energy, or an arbitrary momentum direction.

## Recommended application after Fable review

1. Keep `40--120 M☉` unresolved in the production admission map, but label it
   `model_dependent` and attach the candidate strategy identifiers above.
2. Adopt `per_source_node_fate_lookup_v1` as the implementation target: a
   source-specific node table plus an out-of-hull fail-closed resolver.
3. Use Sukhbold as the first engine-specific solar-metallicity validation
   branch, without promoting its values until the missing wind, decay,
   momentum, and redistribution requirements are closed.
4. Add the LC18/Ugolini multi-metallicity branch only after its source package
   and semantics are independently approved.
5. Keep PPISN/PISN in the separate P-F3 gate; do not fold it into the SNII
   40--120 row.

Fable is asked to verify this choice and to identify which parts can be
implemented now without claiming scientific coverage.  A favorable audit can
authorize the resolver contract and candidate metadata; it cannot by itself
authorize third-party source redistribution or invent absent physical fields.

## Fable result and applied subset

Fable returned **CONDITIONAL PASS**.  The scientific gate remains blocked.  Its
recommendation was applied as a zero-node resolver contract, explicit
piecewise-constant source-cell semantics, fail-closed hull behavior, PISN
complete-disruption/no-remnant ownership, and stronger fate-map policy and
owner cross-checks.  The Sukhbold audit now parses W18 and N20 high-mass
outcomes and all 105 `implosions_W18` wind-only tables; its result is retained
as review evidence and emits no canonical rows.

One Fable note about the Limongi & Chieffi DOI was independently checked against
the article record and arXiv entry: `10.3847/1538-4365/aacb24` in this dossier
is correct.  The DOI correction is therefore not applied.

The remaining production blockers are source approval, metallicity/rotation
coverage, age-resolved winds, terminal ejecta/remnant/deposition semantics,
decay closure, PPISN/PISN policy, and a fresh production-linked build over the
new runtime.  The non-interpolating runtime, unresolved bucket, and
approval-id/map-checksum gate are now implemented, but they do not promote any
physical source values. No 40--120 M☉ physical values have been promoted.

## Primary references

* [Sukhbold et al. 2016](https://doi.org/10.3847/0004-637X/821/1/38) and the
  [MPA CCSN archive](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/).
* [Limongi & Chieffi 2018](https://doi.org/10.3847/1538-4365/aacb24).
* [Fryer et al. 2012](https://arxiv.org/abs/1110.1726).
* [Patton & Sukhbold 2020](https://arxiv.org/abs/2005.03066).
* [Ugolini et al. 2025](https://arxiv.org/abs/2501.18689).
* [Woosley 2017, Pulsational Pair-Instability Supernovae](https://arxiv.org/abs/1608.08939).
* [Heger & Woosley 2002](https://arxiv.org/abs/astro-ph/0107037).
