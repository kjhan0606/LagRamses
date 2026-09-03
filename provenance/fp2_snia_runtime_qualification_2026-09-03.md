# F-P2 runtime qualification — 2026-09-03

Status: **qualified review-only caller; runtime activation remains blocked**.

## Implemented in this bundle

- Added direct startup negative evidence against the linked
  `bin/ramses_final3d` for missing contract, mismatched source commit,
  mismatched approval id, missing thermal group, and valid-contract
  activation.  The first four return the contract failure path; the valid
  handoff reaches the independent production gate and returns before
  evolution.
- Unified the runtime contract environment variable as
  `PHASE0_SNIA_RUNTIME_CONTRACT` and included it in the contract audit.
- Added explicit SNIa linkage symbols and the Makefile dependency from the
  runtime accounting object to `stellar_enrichment_config.o`.
- Made the RAMSES bridge treat a zero element-field index as an explicit
  inactive-element sentinel in both `uold` and `unew` scatter paths.
- Corrected the production population gate so a binary SSP is structurally
  valid only with an explicitly enabled SNIa path, while a single-star SSP is
  valid only without it.  Both still require the reviewed F-P1 fate map, so
  this change does not open production activation.
- Split F-P2 audit semantics into `physical_baseline_ready` and actual
  `production_ready`.  The former is true for the approved Maoz/HESMA
  baseline; the latter is false while runtime activation is disabled.
- Added immutable source-binding diagnostics.  The declared `c6c8042...`
  staging revision exists and is ancestral to current HEAD, but the dirty
  worktree is explicitly reported as not production-code-bound.

## Evidence

The following all pass under `/gpfs/kjhan/LRD_JWST`:

- full forced production linked build and source parity;
- `run_fp2_snia_dtd_contract.sh`;
- `run_p04_production_negative.sh`;
- `run_stellar_feedback_policy_unit.sh`;
- `run_stellar_residual_deposition_unit.sh`;
- `run_g2_population_ledger.sh`;
- `run_g1_native_contract.sh`, including JAX 0.11.1 CPU differential;
- direct Python syntax checks for the new audit/negative harness.

## Remaining blockers

This is not an active SNIa evolution run.  F-P1's terminal-fate map remains
review-only; a hard-crash pending-event journal, distributed MPI/neighbour
deposition, and full metallicity-sensitive/net-yield population semantics
remain separate production/publication gates.  No commit or push is implied
by this evidence record.
