# Fable plan audit — next F-P1 admission-closure bundle

Review, in read-only mode, the proposed next implementation bundle in
`/gpfs/kjhan/LRD_JWST`:

`provenance/fp1_next_bundle_plan_admission_closure_2026-09-04.md`

Do not edit files, run jobs, download data, or contact authors. Read the plan,
the current roadmap, the current audit comparison, the AGY and Opus 5 audit
records, and the referenced source/config/tool files as needed.

## Final project purpose

lagRamses must become a production-ready and publication-ready high-level
hydrodynamics stack focused on radiative transfer, stellar/AGN feedback, and
dust. The present bundle is a prerequisite for physically trustworthy stellar
feedback source selection; it is not itself a source-selection or runtime
deposition bundle.

## Context

The preceding implementation bundle repaired the code-owned source-rights trust
root and independently cross-checked the Boccioli--Roberti LC18 failed-wind
anomaly against locally staged Limongi--Chieffi CDS data. AGY returned PASS;
Claude Opus 5 returned CONDITIONAL PASS. Independent reproduction confirmed the
Opus findings that fate admission does not currently couple its production/
publication result to the physical-package result, physical evidence hashes
are self-certified in an editable contract, the positive selection branch is
unreachable under the current review guard, direct LC18 tool invocation reports
but does not raise on phase violations, and CDS rights are not explicitly
labelled non-authoritative. The real state remains zero physical nodes and all
production/publication/conversion/deposition flags false.

## Questions to answer

1. Is the proposed bundle aligned with the final RT/stellar-AGN-feedback/dust
   purpose, or does it spend effort unrelated to high-level hydro physics?
2. Are the priorities and boundaries scientifically and technically justified?
   In particular, assess the order: commit/provenance boundary, admission
   coupling, non-circular physical-package selection, fail-closed LC18 evidence,
   then adversarial coverage.
3. Is the plan feasible as one implementation bundle without silently
   selecting a source or opening runtime feedback? Identify any part that should
   be split, deferred, or added.
4. Does the plan preserve the unresolved 40--120 M_sun and failed-wind issues as
   honest blockers, avoid inventing age-resolved winds/energy/momentum, and
   maintain publication provenance?
5. Are the acceptance tests sufficient to support a later physical-source
   implementation, and what exact acceptance condition is missing if not?

## Required response

Return a self-contained decision: `APPROVE`, `APPROVE WITH CHANGES`, or
`REJECT`. List mandatory changes with priority and rationale, optional changes,
out-of-scope items that should remain deferred, and a final statement whether
implementation may begin after the listed changes. Do not treat a schema-only
or synthetic test as physical source approval.
