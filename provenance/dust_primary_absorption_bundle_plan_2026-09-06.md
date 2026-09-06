# DUST-6: native primary-photon dust absorption contract

Project `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`.
Base commit: `578bb3b`; operator preapproved continued work on 2026-09-06.
Previous DUST-5 added and audited the native FP64 secondary IR operator, but
the live primary SNRT path still has only the three finite H/He absorbing
inventories. This bundle is the next high-level RT/dust task: define and
implement the fourth, non-depleting dust absorption channel and its physical
heating ledger. It is not live activation or a complete dust science model.

Plan review is Fable first; bundle-end review is Opus 5. One coherent native
coupling contract and one compact native smoke test are intended. No new
Python validator, JSON schema, HDF5 loader, global gate or unrelated AMR audit.

## Physical contract

1. Input dust absorption is an explicit caller-supplied cross section per H
   at the selected reference mixture, `sigma_dust(g)` [cm2/H], and explicit
   cell factor `dust_relative_abundance` [dimensionless]. No hidden
   dust-to-metal law, depletion law or redshift interpolation is introduced.
   The caller supplies `n_H` [cm-3] and path [cm], yielding
   `tau_dust(g)=sigma_dust(g)*dust_relative*n_H*path`.
2. Primary photon transport uses total optical depth
   `tau_total=tau_HI+tau_HeI+tau_HeII+tau_dust`. Dust is not a finite atom
   inventory and is never passed through the H/He reservoir cap. H/He
   assignment remains bounded by its three available neutral inventories.
   If that bound saturates, the residual of the already-removed photon
   absorption is assigned to dust only when dust opacity is present; if no
   finite absorber inventory is available, it is returned to the photon
   field. Physical returned photons and numerical `unassigned` residuals are
   separate ledger terms, so dust absorption is never silently counted as
   H/He ionization or made to disappear.
3. Implement `snrt_dust_coupling.f90` as a pure FP64, caller-owned-array
   contract: validate finite/nonnegative inputs, group shape, total-tau
   closure, and return dust tau, total tau, H/He assigned absorption,
   dust-assigned absorption, unassigned absorption and a per-cell heating
   ledger. Inputs/outputs have explicit physical units: absorbed photons
   [photons cm-3], mean group energy [eV], dt [s], output heating
   [erg cm-3 s-1]. Do not convert code `n_gamma/n_H` inside this generic
   contract; the RAMSES adapter will do that at its existing unit boundary.
4. Provide a native transport-facing adapter that accepts a raw removed
   photon amount per group plus component optical depths and reservoirs. Its
   output must satisfy, group by group,
   `raw_removed = H/He_assigned + dust_assigned + returned + unassigned`.
   `returned` is the physical remnant that the caller restores to the photon
   field; `unassigned` is only a numerical post-partition residual.
   This adapter is a reference/ABI boundary for the next CUDA integration;
   do not alter the existing three-species production ABI in this bundle.

## Implementation and evidence

- Add one native module and one extension to the existing native dust test.
  Test proportional dust/HHe allocation, H/He saturation transfer to dust,
  a sequential multi-group shared reservoir, zero-dust regression in the
  non-saturating regime, returned-photon and group-wise closure, zero/weak
  heating, invalid/negative/nonfinite/range inputs and exact caller-state
  nonmutation on failure. Compile with GNU and Intel.
- Add a small reference Python comparison only inside the existing
  `dust_ir_transport.py --native` fixture path if needed; do not add another
  test file or generate a new evidence dataset.
- Record status/units and one manufactured multi-group case. The case is a
  contract test, not a cosmological or astrophysical dust-mixture claim.
- Keep the DUST-5 IR module and live RAMSES driver unchanged. Do not use this
  contract to mutate `snrt_intensity`, H/He state, gas energy or momentum.

## Exit and handoff

Exit is an audited native dust absorption/heating boundary with local ledger
accounting to floating-point closure. The next bundle can add a new CUDA ABI that includes dust
as a non-depleting absorber, then wire it into the existing transaction. That
future step must prove FP32 allocation, total photon conservation and MPI/AMR
semantics. This bundle does not claim those.

Deferred: dust opacity source/mixture physical approval, depletion/evolution,
IR scattering, live persistent spectral state, gas/force coupling, native
CUDA integration, AMR/MPI/restart, reduced-c physical interpretation and
production/publication qualification.
