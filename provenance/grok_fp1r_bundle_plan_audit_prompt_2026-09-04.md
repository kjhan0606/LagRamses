# Grok bundle-start plan audit: F-P1R evidence and semantics closure

Audit the driver's proposed F-P1R implementation bundle in
`/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`) as the active bundle-start plan
reviewer. Work read-only. Do not edit files, download data, launch HPC or
RAMSES jobs, select a physical source, contact authors, or create approval
artifacts. Inspect the plan and the current parent evidence at HEAD
`db1bb66`.

The final purpose is a production-ready and publication-ready lagRamses
high-level hydrodynamics stack focused on RT, stellar/AGN feedback, and dust.
F-P1R is only a pre-admission evidence/semantics/workflow hardening bundle.
It must preserve the intentional fail-closed state: zero physical source
nodes, unresolved fate seams `[0.8,1.0]` and `[40,120] M_sun`, no canonical
conversion or runtime deposition, and a blocked LC18 publication gate.

Return exactly one plan decision: `APPROVE`, `APPROVE WITH CHANGES`, or
`REJECT`. Evaluate the plan, not just the existing implementation. Check
whether R1–R4 are scientifically and technically justified, correctly
ordered, feasible as one bounded bundle, and aligned with the final project
purpose. Do not demand physical source selection, 40–120 M_sun fate
resolution, source licensing resolution, or runtime activation in this bundle.

## Proposed work to audit

1. **R1 positive converter evidence:** isolated, fully admitted synthetic
   conversion reaches the normal write path only in a temporary directory;
   matching output, mapping-content mutation, and mapping-hash mutation are
   tested; all mutable seams are restored; real repository evidence stays
   unchanged.
2. **R2 LC18 semantics:** distinguish parsed/quantized CDS Table 5 exact-zero
   endpoints from inferred physical zero; expose the `0.005 M_sun` half-bin
   and `physical_zero_inferred=false`; preserve the 4 successful and 3 failed
   zero split, 48/4, 53/3, 101/7, and the unresolved failed-release anomaly.
3. **R3 same-run freshness:** regenerate the high-mass seam report before
   dependent physical/fate admission, or run dependent audits after the final
   regeneration; prove that downstream evidence references the post-run hash.
4. **R4 publication API:** make the code-owned publication gate read and hash
   the locked terms bytes internally, parse the candidate record from those
   bytes, and refuse caller-supplied digest authority; keep synthetic locks
   test-only and current LC18 publication review-only.

## Required judgments

- Is R1 a sufficient and safe way to close the current positive-path evidence
  gap without weakening production contract checks?
- Are R2's naming and physical interpretation correct for rounded CDS data,
  and does the requested outcome split avoid conflating the failed BR26
  zero-wind release anomaly with a CDS endpoint quantization result?
- Does R3 establish a genuine same-invocation dependency/freshness guarantee,
  and is the proposed ordering compatible with the existing fixture runner?
- Does R4 improve provenance security without accidentally converting review
  terms into production permission or relying on mutable labels?
- Are the acceptance tests sufficient and non-destructive, and are all
  unrelated infrastructure/physics items properly deferred?

List mandatory changes, risks, and valid deferrals with exact file references
where possible. If the plan is approved with changes, identify the required
edits before implementation. End by stating that implementation may begin
only after the plan decision is acceptable, and that Claude Opus 5 audits each
completed implementation step.
