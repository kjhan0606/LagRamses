# F-P2 SNIa/HESMA bundle — independent AGY audit

- Model: `gemini-3.8-flash-high` through `/home/kjhan/local/bin/agy`
- Date: 2026-09-03
- Mode: read-only plan/sandbox audit; no files modified, no model selected, no runtime activated, no commit/push

## Verdict

**PASS for review-only contract/scaffolding; BLOCK for production runtime activation.**

The DTD integration and event-ledger kernels are judged mathematically sound, conservative, restart-independent, and identical between native and production mirrors. The checksum chain and fail-closed gates are judged complete. Production remains blocked because the physical SNIa event contract and source selection are unapproved.

## Core judgment

- The normalized power-law interval integral, exact log branch, support clipping, telescoping, and defensive input checks are correct.
- The event ledger scales per-event fields by expected event count, preserves untracked residual mass, and rejects tracked ejecta above returned mass.
- Native/production SNIa mirrors are byte-identical. Artifact, sidecar, manifest, tool, and audit hashes were reported consistent.
- The fail-closed chain rejects `use_snia=.true.` at support/configuration/driver layers, and selection fields remain null across comparison, estimator, packet, and contract audit.
- HESMA source-format evidence is kept separate from unresolved decay horizon, isotope aggregation, returned mass/remnant, energy, momentum, and population weighting.
- The 5% screen is explicitly diagnostic and cannot silently select a source. n300c has approximately 641% profile discrepancy and n1600c approximately 5.11%; estimator sensitivity is small.

## Severity-ranked findings

### Production blockers

1. Physical event fields and DTD normalization/bounds remain null by design and must be approved before promotion.
2. No progenitor family or population mixture weights are selected.
3. The 384-isotope to 11-element decay aggregation is not approved.
4. n300c and n1600c profile discrepancies must be resolved or explicitly excluded from any production mixture.

### Medium improvements

1. Replace the near-α=−1 subtractive power-law expression with an `expm1`/series-stable formulation for small nonzero exponents.
2. Reconsider the event-ledger tolerance scale `max(1.0, returned_mass, tracked_ejecta)` if inputs can be in units much smaller than one; use a contract-specific scale or relative floor.
3. Record the compositional incompatibility between the HESMA profile source and the Keegans/NuGrid candidate, which lacks H, He, C, and N.

## Required before production approval

- Select a source model or mixture and named approval.
- Approve DTD delay bounds and normalization.
- Approve decay horizon and isotope-to-project-element mapping.
- Set authoritative returned mass, terminal remnant, energy, and momentum conventions.
- Generate and hash-bind the canonical event asset.
- Activate the runtime gate only under the named approval.

## Review-only assessment

Correctly enforced. The current gate is safe for review-only evidence; the residual risk is scientific incompleteness of the physical event contract, not accidental activation.
