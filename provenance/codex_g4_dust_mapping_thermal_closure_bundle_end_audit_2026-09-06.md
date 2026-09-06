# Codex end audit: G4 dust mapping and thermal-ledger closure bundle

Date: 2026-09-06
Project: `/gpfs/kjhan/LRD_JWST`
Base: `cedf187`

## Verdict

**PASS within the bounded engineering scope; conditional for
production/publication.**

The array and yt staging boundaries no longer accept a partial
`metallicity_solar`/`dust_to_metal` mapping. The derived path uses the exact
cell product and preserves both source fields plus the origin in P5 output.
The existing source-bound dust and one-pass thermal ledger controls pass in
the same CPU/JAX bundle gate.

## Independent scope check

No physical DTM/depletion law, opacity-mixture approval, live RAMSES dust
activation, gas-energy or force injection, recursive IR transport, or
AMR/MPI/restart qualification is claimed. The analytic `Z=1e-6`, `DTM=1`
fixture is only a data-contract test. The existing IR status remains
`recorded_not_transport_reemitted`.

## Residual promotion work

G4 still requires source/mixture-specific physical review, an approved
depletion/DTM prescription, source-domain and thermal-atlas convergence, and
matched transport convergence before dusty static SNRT can be promoted. Live
RAMSES dust state and feedback remain G5 work.
