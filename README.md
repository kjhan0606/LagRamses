# LRD_JWST

Independent working repository for the HR5 high-redshift SMBH, dual-AGN, and
JWST little-red-dot (LRD) zoom-in programme. The active scientific scope is
recorded in `docs/lrd_zoom_simulation_design.md`.

## Repository policy

- `paper/` contains manuscript sources once drafting starts.
- `simulation/` contains reproducible input configurations, source ledgers,
  analysis tools, test fixtures, and Slurm submission scripts.
- `manifests/` records immutable locations and SHA256 hashes for ICs, outputs,
  catalogues, and other large external assets. Large data are never committed.
- `provenance/` pins external code and migration provenance.
- `docs/archive/legacy_dualagn_paper3/` preserves the superseded dual-AGN
  Paper-III planning record. It is not the active LRD/JWST work plan.

`lagRamses` remains an external dependency and is not vendored here. See
`provenance/lagramses.md`.
