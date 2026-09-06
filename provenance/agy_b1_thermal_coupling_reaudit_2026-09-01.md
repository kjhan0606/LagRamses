# AGY B1 thermal-coupling re-audit

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: `/home/kjhan/local/bin/agy` (Gemini Antigravity CLI)  
Model: `gemini-3.1-pro-high`  
Scope: B1 thermal coupling only; read-only re-audit.

## Verdict

**PASS**

## Original-finding disposition

1. **Critical — fixed.** Format v3 removes the metallicity dimension. Host and
   JAX runtimes multiply the solar-metallicity table analytically by `Z/Zsun`.
   Host inputs reject negative/non-finite metallicity; JAX emits NaN for those
   invalid values rather than physical cooling/heating.
2. **High — fixed.** The B1 test exercises off-grid scalar and vector values,
   including `10^-2.35`, `0.37`, and `sqrt(10)`, on both host and JAX paths.
3. **Medium — fixed.** The generator continuously subtracts the coefficient at
   `T_CMB` without the Grackle two-dex cutoff, and HDF5 provenance explicitly
   records `grackle_2dex_cutoff_removed_for_continuity`.

## Additional findings

No new algorithm defect or B1 validation gap was found. AGY independently
confirmed the `(a, n_H, T)` artifact dimensions, runtime broadcast semantics,
legacy-format rejection, tests, and recorded hashes.

## Promotion disposition

The code, test, and artifact-provenance portion may clear
`independent_gate_audit_pending`. The separate upstream data-license status
remains `pending_explicit_data_license_confirmation` and independently blocks
production promotion. This PASS applies only to B1 and does not imply overall
production or publication readiness.
