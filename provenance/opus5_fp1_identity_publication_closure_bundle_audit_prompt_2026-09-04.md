# Claude Opus 5 F-P1 identity/publication closure bundle audit request

Act as the independent scientific and implementation auditor for the
completed F-P1 identity and publication closure bundle in the lagRamses/SNRT
project. Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, launch
jobs, select physical sources, contact authors, or redistribute CDS data.
Audit the current checkout at HEAD `5aeb6d3`; the implementation commit is
`25bd05f`.

The final project goal is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack for radiative transfer, stellar/AGN
feedback, and dust. F-P1 is an integrity boundary before physical yield/fate
activation. It intentionally leaves the real state review-only and blocked:
zero physical nodes, unresolved fate intervals `[0.8,1.0]` and
`[40,120] M_sun`, and false conversion/deposition/production/publication
flags. Do not penalize the bundle for not claiming to solve those later
physics/data gates.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` for an F-P1 defect or absent evidence that makes source identity,
mapping integrity, publication rights, or fail-closed behavior unsound. Use
`CONDITIONAL PASS` for material but non-blocking follow-up. Distinguish both
from later work that is correctly deferred.

This is an algorithm/architecture/wiring/scientific-legitimacy audit, not a
generic style or superficial bug sweep. Independently read the code, tests,
contracts, and generated evidence. Reproduce bounded fixture tests when they
are non-writing; if a runner regenerates tracked JSON, inspect the code or use
an isolated temporary copy. Do not trust report labels without tracing their
inputs.

## Acceptance contract

1. Passed source-identity validator reports require a valid package
   fingerprint. Positive admission consumes all nine identity-matched passed
   reports and binds the selected package SHA to executable source-identity
   evidence. The real selection remains null/blocked.
2. One shared canonical source-node mapping serializer owns exact hash bytes,
   deterministic numeric normalization, duplicate/invalid-coordinate rules,
   and the bindings for source-node contract, approvals, package, asset, row
   count, and node coverage. Conversion must verify exact admitted/generated
   mapping equality before writing. Proposal mode must perform no writes.
3. The derived-artifact publication gate must be code-owned and fail closed:
   hash-locked terms path/bytes, verified production license, source and
   derived redistribution permissions, citation/attribution, explicit
   approval and approval identity, candidate/artifact identity, and
   `review_use_only=false` are all required. A report label alone cannot
   promote an artifact.
4. LC18 successful-control and total endpoint counts must be symmetric and
   correctly interpreted as 48 positive/4 zero successful, 53 positive/3 zero
   failed, 101 positive/7 zero overall, while existing anomaly blockers stay
   active.
5. Evidence must be reproducible and scope-honest: tests/compilation pass,
   config/data hashes are invariant under deterministic regeneration, exact
   contract hashes are updated, and no runtime feedback/source activation is
   introduced.

## Required report format

Give a short verdict first, then findings ordered by severity. Every finding
must include exact file and line references, the impact on the final project
goal, F-P1 disposition (blocker versus later gate), and a concrete remedy.
State which acceptance items you independently reproduced and whether any
observed issue is a genuine contradiction or merely the intentional blocked
state. End by confirming that a next bundle requires a driver plan and Fable
approval before implementation.

## Primary files

`simulation/snrt/tools/fp1_source_node_mapping.py`,
`simulation/snrt/tools/fp1_gate_validator_registry.py`,
`simulation/snrt/tools/audit_fp1_physical_package_admission.py`,
`simulation/snrt/tools/convert_yield_rows_to_canonical.py`,
`simulation/snrt/tools/fp1_publication_rights.py`,
`simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`, their F-P1
tests, the F-P1 contract/config files, generated F-P1 reports, and the
provenance plan/completion records.
