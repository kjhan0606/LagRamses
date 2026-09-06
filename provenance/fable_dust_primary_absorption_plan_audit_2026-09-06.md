# Fable audit: DUST-6 bundle plan

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Base: `578bb3b`
- Auditor/model: Fable / `claude-fable-5-1`
- Session: `488dc027-b326-4fec-9bac-0e5329b68ddd`
- Duration: 173 s
- Verdict: **CONDITIONAL APPROVE**

## Scope finding

The bundle is appropriately scoped as a native primary-photon dust
absorption/heating boundary. The existing CUDA path has a three-species H/He
finite-inventory ABI, while dust must be a fourth, non-depleting absorber.
The current bundle should therefore leave the CUDA ABI, live RAMSES driver,
AMR/MPI/restart behavior, and physical opacity-mixture approval unchanged.

## Required amendments

1. Do not transfer all H/He saturation residual to dust. Allocate the direct
   optical-depth fraction to dust, and transfer from a saturated H/He excess
   only `excess * (1-exp(-tau_dust))`. Return the remainder to the photon
   field so dust heating is not overestimated at an ionization front.
2. Define one group-wise closure including the returned remnant:
   `raw_removed = HHe_assigned + dust_assigned + returned + unassigned`.
   `returned` is a physical photon-field restoration; `unassigned` is only a
   numerical post-partition residual.
3. The zero-dust case must call the existing H/He partition contract and
   retain an exact regression. No separate Python validator, JSON schema, or
   HDF5 evidence asset is needed.
4. Keep the partition contract homogeneous in photon amount. Only the
   heating ledger converts physical photons/cm3 and mean eV over seconds to
   erg cm-3 s-1. Code-unit conversion belongs at a future RAMSES adapter.

## Acceptance conditions

The native smoke should cover proportional allocation, H/He saturation,
zero-dust regression, returned-photon conservation, zero/weak heating,
invalid input, and caller-state nonmutation, with GNU and Intel builds. The
result is a contract test with manufactured values, not approval of a
cosmological dust mixture or production readiness.
