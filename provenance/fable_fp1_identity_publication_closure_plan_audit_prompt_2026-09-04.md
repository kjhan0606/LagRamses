# Fable plan audit: F-P1 identity and publication closure bundle

You are reviewing a proposed implementation bundle in read-only mode. Do not
edit files, download data, contact authors, launch HPC/production jobs, or
create approval artifacts. Read the plan and relevant committed source/data
only. The project is `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`).

The final purpose is a production-ready and publication-ready lagRamses
high-level hydrodynamics stack focused on radiative transfer, stellar/AGN
feedback, and dust. The proposed bundle follows the completed F-P1 admission
closure bundle and targets the audit findings F1 (selected package hash is not
bound to executable source identity), F2 (CDS-derived review artifacts lack a
technical publication gate), and F3 (successful-control zero-wind statistics
are incomplete).

Evaluate whether B1--B4 are scientifically and technically justified, ordered
correctly, feasible in one coherent bundle, and aligned with the final
purpose. Check especially:

1. whether the selected package hash is bound to the passed executable
   `source_identity_and_rights` fingerprint rather than a self-consistent
   editable value;
2. whether the proposed canonical source-node mapping serialization is
   deterministic and can be shared safely by admission and the converter;
3. whether the publication gate is a real executable rights boundary and
   cannot be bypassed by `review_only`/`publication_ready` labels, while the
   current Limongi CDS record remains blocked;
4. whether the synthetic positive fixtures are isolated and cannot change real
   evidence; and
5. whether the acceptance criteria preserve the current fail-closed state and
   do not silently promote missing wind, energy, momentum, or 40--120 M_sun
   physics.

Return a decision of APPROVE, APPROVE WITH CHANGES, or REJECT. List mandatory
changes, risks, and any items that should be deferred. Do not treat the plan
as physical-source or runtime approval.
