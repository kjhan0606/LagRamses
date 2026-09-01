# Migration from HR5_dualAGN

Created 2026-08-31 from `/home/kjhan/BACKUP/HR5_dualAGN` without modifying the
source workspace.

Transferred:

- active LRD/JWST simulation design;
- SNRT source, configuration, input ledgers, tests, analysis tools, and Slurm
  submission scripts;
- pilot run configuration only;
- archived prior dual-AGN Paper-III planning documents.

Excluded deliberately:

- `lagRamses` and the RAMSES Bitbucket reference clone;
- virtual environments and third-party vendor trees;
- scheduler logs, transient status files, job directories, and `output_*`;
- raw sink catalogues and all other large simulation data.

Post-migration audit on 2026-08-31:

- Active working tree for continued implementation is
  `/gpfs/kjhan/LRD_JWST`; the source checkout is read-only provenance.
- Active scripts and derived metadata were relocated to the current
  `/gpfs/kjhan/LRD_JWST/simulation/` layout.
- Slurm and pilot submission scripts reference
  `manifests/lrd_jwst_external_assets.json` before starting.
- The 866 GB `output_00017/data_00017.h5` payload remains external by design;
  its source path, size, pending payload hash, and reference-only staging
  policy are recorded in the asset manifest. It is a raw RAMSES HDF5
  checkpoint, not a direct lagRamses runtime payload; the P4 adapter extracts
  the required AMR leaves into a smaller canonical SNRT input. A 320 GB
  interrupted transfer was quarantined and is not a valid input. Small output
  metadata sidecars, both thermal atlases, and the regenerated 32³ P4
  canonical inputs now live under `/gpfs` with verified checksums.
