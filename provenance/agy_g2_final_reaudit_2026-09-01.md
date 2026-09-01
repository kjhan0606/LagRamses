# AGY final G2-only re-audit

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: `/home/kjhan/local/bin/agy` (AGY/Antigravity CLI)  
Model: `gemini-3.1-pro-high`  
Scope: G2 physical stellar-yield gate only; read-only audit.

## Verdict

**BLOCKED**

## Confirmed resolved controls

- The main Python contract and asset audit now carry the terminal-remnant
  ownership map and reject non-terminal remnant rows.
- The native audit enforces the same ownership rule.
- [`stellar_population_ledger.f90`](../simulation/snrt/native/phase0/stellar_population_ledger.f90)
  aggregates channel states, checks channel and total ejecta/returned closure,
  enforces remnant ownership, and derives living mass once at population
  finalization.
- [`stellar_enrichment_driver.f90`](../simulation/snrt/native/phase0/stellar_enrichment_driver.f90)
  wires the SSP channel integrations to that ledger; the SSP integrator carries
  channel remnant contributions.
- Enabled SNIa/PISN paths fail closed until separate event models/gates exist.
- `G2_CONFIGURATION_TEST_OK`, `G2_POPULATION_LEDGER_TEST_OK`,
  `YIELD_CONVERTER_TEST_OK`, and the expected `G2_PREFLIGHT_BLOCKED` result are
  reproducible. Runtime adapter syntax compilation also succeeds.

AGY confirmed these are code-level and synthetic controls only; they are not a
physical-data approval.

## Remaining blockers

The project contains no approved, checksummed physical full-grid asset for the
required wind, AGB, and SNII channels. The source matrix contains literature
candidates only. The legacy tables and synthetic fixtures remain correctly
refused by the production gate.

## Required input for G2 promotion

1. Select final wind, AGB, and CCSN/SNII source versions and assumptions.
2. Supply the original source files and convert them to the canonical 32-field
   complete mass--metallicity--age grid, including age-zero coverage.
3. Supply a matching provenance sidecar with source/conversion SHA256 values,
   license/provenance approval, units, IMF, population model, channel
   boundaries, metallicity/solar definitions, remnant model, and approval ID.
4. Mark channels approved only after both Python and native audits pass on the
   actual physical asset.

This audit does not promote G2 and does not audit other gates.
