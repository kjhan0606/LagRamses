# AGY G2 candidate-source re-audit — 2026-09-01

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: `/home/kjhan/local/bin/agy` (AGY/Antigravity CLI)  
Model: `gemini-3.1-pro-high`  
Scope: G2 physical stellar-yield gate and the newly staged candidate-source
package only; read-only audit.  
AGY original artifact:
`/home/kjhan/.gemini/antigravity-cli/brain/fe995b42-7635-42f6-bd0d-3347f7f67170/g2_audit_report.md`

## Verdict

**BLOCKED**

## Verified controls

- The Limongi & Chieffi 2018 CDS and NuGrid Set1ext candidate files are
  immutably recorded with hashes, and the acquisition-manifest check passes.
- The G2 physics contract requires a complete 32-field Cartesian grid over
  mass, birth metallicity, and age for channels 1 (wind), 2 (AGB), and 3
  (SNII), with no mass gaps.
- The native population ledger enforces
  `initial = living + remnant + returned` and single terminal-remnant
  ownership.
- The generic JSON-to-ASCII converter enforces canonical ejecta/net-yield,
  energy, and momentum fields and refuses undeclared rate-table integration.

## Remaining blockers

- Both staged candidates provide integrated/terminal yields and lifetime
  scalars rather than explicit age-resolved cumulative histories.
- Canonical per-star energy and momentum fields are absent.
- Limongi & Chieffi does not cover the complete runtime wind (0.8–120 Msun)
  and SNII (8–40 Msun) ranges; its isotope mapping, decay, and pre-SN
  wind/terminal-ejecta partition still require an explicit adapter.
- NuGrid has a duplicate `(5.0 Msun, Z=0.01)` coordinate and incomplete
  runtime mass coverage; its total/wind/pre-explosion partition also requires
  source-semantics review.
- No candidate has verified license status, a project approval identifier, or
  a compliant production provenance sidecar.

## Exact evidence paths

- `simulation/snrt/config/g2_physics_contract_v1.json`
- `simulation/snrt/config/g2_source_selection_matrix_v1.json`
- `external/g2_candidates/acquisition_manifest_v1.json`
- `simulation/snrt/data/g2_candidate_source_audit.json`
- `simulation/snrt/data/g2_preflight.json`
- `simulation/snrt/tools/convert_yield_rows_to_canonical.py`
- `simulation/snrt/native/phase0/stellar_population_ledger.f90`

## Minimum requirements for promotion

1. Supply complete physical Cartesian coverage for wind 0.8–120 Msun, AGB
   1–8 Msun, and SNII 8–40 Msun.
2. Supply explicit age-resolved cumulative histories; a rate-to-cumulative
   conversion must be source-specific and reviewed before the generic
   converter is used.
3. Supply valid per-star energy and momentum for every enabled channel.
4. Define isotope-to-11-element/decay mapping and disjointly partition
   pre-SN winds, AGB release, and terminal CCSN ejecta.
5. Verify the license, obtain an `approval_id`, and pass the Python/native
   audits with a compliant provenance sidecar.

This re-audit confirms the new candidate package and controls but does not
promote G2 or authorize a lagRamses simulation.
