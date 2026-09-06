# Grok replacement audit: F-P1 identity and publication closure plan

Act as the replacement for the Fable plan reviewer for the F-P1
identity/publication closure bundle in `/gpfs/kjhan/LRD_JWST` for
`kjhan0606/LagRamses`. Work strictly read-only: do not edit files, download
data, contact authors, launch HPC/production jobs, or create approval
artifacts. Inspect the plan, the current implementation, and its evidence.

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. The reviewed bundle is an integrity boundary
before physical stellar-yield/fate activation; it must not be judged as a
physical-source or runtime approval. The real state intentionally remains
fail-closed with zero physical source nodes, unresolved fate intervals
`[0.8,1.0]` and `[40,120] M_sun`, and false production/publication/conversion/
deposition flags.

Review the proposed/implemented B1--B4 work for scientific and technical
justification, ordering, feasibility, claim scope, and alignment with the
final purpose. Return exactly one decision: `APPROVE`, `APPROVE WITH CHANGES`,
or `REJECT`. If reviewing the already implemented tree, explicitly separate
plan approval from implementation findings.

Check in particular:

1. Selected package identity is bound to the passed executable
   `source_identity_and_rights` package fingerprint, and the validator
   registry cannot pass missing/malformed fingerprints or a mere list of gate
   IDs.
2. The canonical source-node mapping serializer is deterministic and shared
   safely by admission and conversion; mapping hash, package/asset hash,
   contract hash, approval identity, coordinates, and node coverage are all
   bound. Proposal mode must be non-writing.
3. The derived-artifact publication gate is a genuine code-owned rights
   boundary, not a mutable `review_only`/`publication_ready` label. The
   current Limongi CDS terms remain blocked unless explicit terms, license,
   attribution, and artifact approval evidence exist.
4. LC18 successful-control and total zero-wind accounting is symmetric:
   48 positive/4 zero successful controls, 53 positive/3 zero failed models,
   and 101 positive/7 zero overall, with the anomaly still unresolved.
5. Synthetic fixtures cannot modify real evidence; deterministic regeneration,
   contract hashes, and fail-closed state are preserved.

List mandatory changes, risks, and valid later-gate deferrals. Do not demand
that this bundle resolve the deferred physical 40--120 M_sun fate decision,
energy/momentum data, source licensing, author inquiry, or runtime activation.
End by stating that any next implementation bundle requires a driver plan and
formal review/approval under project governance.
