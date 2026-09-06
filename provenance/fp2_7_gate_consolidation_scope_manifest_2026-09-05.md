# F-P2.7 checkpoint scope manifest

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Purpose: define the first checkpoint commit before the consolidated gate and
initialized-RAMSES smoke.

## Included in checkpoint scope

- Modified tracked source and build wiring under `patch/lagRamses/` and
  `bin/Makefile` for the F-P1/F-P2 stellar/AGN/RT/chemistry work.
- New native F-P2 Fortran/CUDA modules, focused native smoke sources, native
  test runners, SNRT contracts/configuration, and the associated
  `simulation/snrt` source/test/tool files.
- F-P1/F-P2 plans, implementation evidence, audit records, operational
  decisions, and this bundle's provenance index/scope manifest.
- The existing small JSON/CSV manifests and ledgers that are part of the
  source/asset contracts; no contents are rewritten by the checkpoint.

## Explicitly excluded

- Root-level compiler products matching `*__genmod.f90`, all `*.mod`/`*.o`,
  executables, and temporary build directories.
- `simulation/snrt/.venv/`, `build/`, `bin/` generated products, and cached
  Python/JAX outputs.
- `.quarantine_hdf5/`, including its 299 GB partial HDF5; no move or deletion.
- Existing large validation trees and external assets; manifests remain, data
  ownership is unchanged.
- Any unrelated source outside the RT/feedback/dust project scope.

## Commit and cleanup policy

The checkpoint is local only and exists to make source/binary evidence hashes
resolvable. It is not a GitHub push. No superseded JSON, prompt, native phase0
mirror, HDF5, virtual environment, or JAX tree is moved or deleted in this
bundle. Such changes require a separate explicit storage/archive decision.

After the checkpoint, the only implementation work in F-P2.7 is the single
native bundle gate, the small initialized-RAMSES smoke, and the production-
source transition for the four mirror-building runners. The bundle has one
implementation evidence record and one end audit.
