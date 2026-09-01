# SNRT algorithm and wiring audit

Status: recorded 2026-08-31  
Audit mode: scientific architecture review (algorithm, data flow, and physical
justification); this was not a defect-hunting or production-run audit.  
Reviewer role: Fable-requested independent review (the agent runtime displayed
the reviewer as `Nash`).

## Executive decision

The current implementation is internally coherent as a static, fixed-grid,
non-scattering, limited H/He S_N transport prototype.  The CPU P0/P1/P2-P3/P8
checks provide deterministic synthetic-regime evidence for that scope.

It is not yet scientifically sufficient for the full LRD problem described in
the P0 contract.  In particular, the repository must not claim full 9-group
RT, dust obscuration or IR re-emission, stellar plus live AGN radiative
feedback, radiation-pressure dynamics, H2/metal chemistry, or coupled RHD.

## What the audit found justified

- Carlson level-symmetric S4/S6/S8 quadrature, first-order upwind finite-volume
  transport, vacuum boundaries, and local exponential absorption form a
  consistent limited transport operator.
- The H/He conservative path assigns absorbed photons by species opacity and
  uses analytic H relaxation, backward-Euler helium chemistry, and explicit
  photon/recombination/photoelectron ledgers.
- Strömgren, shadow, crossing-beam, H/He residual, and sharding checks are
  useful evidence for the declared synthetic prototype regime.
- The source provenance boundary is correct: external sink diagnostic → AGN
  rate ledger → audited photon ledger → cell-centred source deposition.  SNRT
  does not infer luminosity from sink mass alone.

Relevant implementation anchors:

- [`quadrature.py`](../simulation/snrt/snrt_core/quadrature.py)
- [`transport.py`](../simulation/snrt/snrt_core/transport.py)
- [`conservative_primordial.py`](../simulation/snrt/snrt_core/conservative_primordial.py)
- [`sink_diagnostic.py`](../simulation/snrt/snrt_core/sink_diagnostic.py)
- [`source_ledger.py`](../simulation/snrt/snrt_core/source_ledger.py)

## Findings requiring implementation or an explicit scope gate

1. Group opacity must be SED-weighted.  A production group needs
   `∫N_E σ(E)dE / ∫N_E dE`, and the photoelectron excess energy must use the
   same absorber-weighted integral.  The new offline closure is implemented in
   [`primordial.py`](../simulation/snrt/snrt_core/primordial.py), serialized by
   [`p4_build_agn_photon_ledger.py`](../simulation/snrt/tools/p4_build_agn_photon_ledger.py),
   and consumed by P4-P6 runners.
2. The photon ledger must use the transport operator's exact absorbed loss,
   including source photons emitted within the step.  The ledger now rebuilds
   the same `pre_absorption → exp(-tau)` operation in
   [`ledger.py`](../simulation/snrt/snrt_core/ledger.py).
3. The RAMSES/HDF5 staging route needed a versioned conservative field map.
   The map, all-field leaf deposition, coverage check, density mass balance,
   and v2 canonical output are now implemented by
   [`p4_stage_hdf5_level15.py`](../simulation/snrt/tools/p4_stage_hdf5_level15.py).
   The current pilot map still uses explicit non-production constants for
   unresolved composition and chemistry fields, so the production gate remains
   closed until the real checkpoint is audited.
4. The existing `implicit.py` closure is local chemistry, not implicit
   transport.  Backward-Euler transport, matrix-free iteration, and DSA remain
   future work.
5. Dust now supplies absorption/heating plus an absorption-only momentum
   diagnostic.  Scattering, dust temperature, IR re-emission, and full
   radiation-pressure coupling are not wired into the live transport/hydro
   path.
6. The stellar SED input contract and offline CSV-to-photon-ledger converter
   are now implemented, but no production SED asset is staged and live RAMSES
   feedback is not implemented. Current AGN sources are post-processed,
   ledger-driven snapshots, not a live RT → accretion/hydro/feedback loop.

## External-data gate

The 866 GB HDF5 payload remains external/not migrated, while its small metadata
sidecars and both thermal-atlas files are now available under the active
`/gpfs` tree and have verified checksums in
[`lrd_jwst_external_assets.json`](../manifests/lrd_jwst_external_assets.json).
JAX CPU 0.11.1 and `h5py` are installed in the active `/gpfs` working copy.
The real HDF5/field-map preflight passes for populated levels 1--15; no
production calculation was launched.

## Go / no-go

- **Go:** fixed-source, static-grid, non-scattering H/He ionizing/X-ray
  transport prototype with explicit SED closure and accounting tests.
- **No-go:** LRD dust obscuration, IR observables, stellar/AGN radiative
  feedback, radiation-pressure dynamics, H2/metal chemistry, full P0 9-group
  claims, or quantitative RHD interpretation.

The remaining scientific gates are: stage and validate the adopted stellar SED
table (including its IMF and metallicity convention); replace the remaining
dust/chemistry pilot constants with audited HDF5 fields or an explicitly
approved subgrid prescription so the production contract passes; and then
choose and implement the dust/stellar/feedback scope explicitly.
