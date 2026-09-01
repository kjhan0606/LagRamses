# AGY G2-only re-audit

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: `/home/kjhan/local/bin/agy` (AGY/Antigravity CLI)  
Model: `gemini-3.1-pro-high`  
Scope: G2 physical stellar-yield gate only; read-only audit.

## Verdict

**BLOCKED**

## Resolved code-level controls

- The native population ledger is wired through the SSP driver and computes
  living mass once after channel aggregation.
- Channel-level ejecta/returned closure and terminal-remnant ownership are
  enforced; non-terminal channels cannot contribute a remnant.
- SNIa/PISN activation is fail-closed until their separate event models and
  gates exist.
- The G2 configuration, population-ledger, converter, and preflight tests pass
  while retaining the expected blocked status for unapproved assets.

Evidence includes:

- [`stellar_population_ledger.f90`](../simulation/snrt/native/phase0/stellar_population_ledger.f90)
- [`stellar_enrichment_driver.f90`](../simulation/snrt/native/phase0/stellar_enrichment_driver.f90)
- [`g2_population_ledger_test.f90`](../simulation/snrt/native/phase0/g2_population_ledger_test.f90)
- [`run_g2_preflight.sh`](../simulation/snrt/tests/run_g2_preflight.sh)
- [`g2_preflight.json`](../simulation/snrt/data/g2_preflight.json)

## Remaining blockers

The source matrix contains only candidate literature references. No approved,
checksummed physical full-grid asset exists for required wind, AGB, and SNII
channels. The legacy `yield_table.asc` and nine-row validation fixture are
correctly refused and cannot promote G2.

## Exact input required for promotion

1. Select final physical source versions for wind, AGB, and CCSN/SNII.
2. Supply the source files and convert them into the canonical 32-field,
   complete mass--metallicity--age grid with age-zero coverage.
3. Supply a matching provenance sidecar with source/conversion SHA256 values,
   license/provenance approval, units, IMF, population model, boundaries,
   metallicity/solar definitions, remnant model, and approval ID.
4. Mark the channels approved only after the Python and native audits pass on
   those actual files.

This re-audit does not promote G2 and does not audit G0, G1, G3--G7, HDF5,
RT, or HPC.
