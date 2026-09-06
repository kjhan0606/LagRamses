# Primary/IR state transfer at grid creation and destruction

This bundle adds the state-transfer portion of AMR support. It does **not**
complete coarse-fine transport, time synchronization/reflux, load-balance
migration, or a live evolving-AMR qualification. The v3 dust startup restriction
`levelmin=levelmax` remains in place until those pieces are connected.

## Active source and VPATH

The user explicitly reminded us to verify VPATH precedence. The current build
has `PATCH=../patch/lagRamses` and
`VPATH=$(PATCH):../patch/cuda:../patch/oct_tree:../patch/cuRamses:../$(SOLVER):../aton:../hydro:../pm:../poisson:../amr`.
There is no earlier `refine_utils.f90` override: the active source is
`patch/cuRamses/refine_utils.f90`, **not** `amr/refine_utils.f90`.

Confirmed with
`make -C bin -n -W ../patch/cuRamses/refine_utils.f90 SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1 refine_utils.o`:
the compiler command explicitly uses `-c ../patch/cuRamses/refine_utils.f90`.
`nm bin/refine_utils.o` shows references to the new refine, coarsen and error
handling routines. No duplicate fallback copy was edited. The Makefile adds
the module dependency so parallel builds compile the hook module first.

## State semantics

`snrt_regrid.f90` transfers the actual persistent primary and IR payloads.
Grid creation copies parent **densities** to each child: no factor of eight
is applied to density, because each child's volume is one eighth of the
parent's. This is a first-order, positivity-preserving prolongation, not a
claim of second-order spatial reconstruction. The parent is retained as a
coarse representation; physical totals must sum current leaves, not parent
and children together.

Grid destruction averages photon/IR density over equal-volume children.
H II, He II and He III fractions use gas-density weights instead, preserving
their rho*x inventories. Zero total gas density produces neutral fractions.
Hydro remains responsible for its own conservative gas and dust mass/thermal
energy restriction; these radiation routines do not rewrite uold.

All children are read and validated before a coarsen update. Mixed missing/
present radiation children, invalid radiation or negative/nonfinite densities
are rejected. Retired child photon/IR/chemistry payloads are cleared while their
cell-ID/slot associations are retained for pool reuse and raw-checkpoint
compatibility. Refining an unrelated uninitialized parent cannot resurrect
old child radiation. `snrt_state_initialize` now extends its cell-slot lookup
when RAMSES grid capacity increases, preserving existing entries.

Calls are attached to the active `make_grid_coarse`, `make_grid_fine`, and
`kill_grid` routines. Restriction occurs before child topology is disconnected
or uold is zeroed. Only locally owned non-boundary grids are handled; balance
operations are excluded because migration is a separate unresolved path.
Early initialization with no SNRT slots is a no-op. A transfer failure uses
MPI_Abort rather than rank-local MPI_Finalize, which could strand other ranks.

## Verification

Retained root: `/gpfs/kjhan/LRD_JWST/.snrt-regrid.Gwcrnl`.

`snrt_regrid_native` links the same production module, AMR/hydro and CUDA
objects as the simulation, with a small Fortran main instead of RAMSES's main.
It directly invokes the native transfer APIs on one parent and eight children;
it is not a full live-grid refinement calculation and produces no snapshots.

Build and run:

```
make -s -C bin -j4 SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1 \
  EXEC=/gpfs/kjhan/LRD_JWST/.snrt-regrid.Gwcrnl/ramses snrt_regrid_native
SNRT_ALLOW_REFERENCE_CONTROL=1 \
SNRT_DUST_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/dust_native_reference_control_v3.nml \
OMP_NUM_THREADS=1 .snrt-regrid.Gwcrnl/ramses_regrid
```

Result: `SNRT_NATIVE_REGRID_PASS`. Tested: exact uniform prolongation;
nonuniform child photon/IR restriction; density-weighted H and He fractions;
retired payload clearing; reuse under an unrelated parent; negative photon
and gas-density rejection without parent mutation; preservation of existing
cell mappings during capacity growth. Primary restriction allows its FP32
rounding; IR and ion checks use double precision.

Test executable SHA256:
`25b28f704cb8e03db99c0669755625b6c639206f0dd6f3ec610ef9d4c7fb16bc`.

`bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh` passes with
GNU, including raw checkpoint roundtrip and negative/NaN/Inf rejection. This
checks that the state-map change does not break that persistence path.

The full live-profile RAMSES binary also builds, SHA256:
`17049b76bfe90a10b552155b01034593641bc77de2a519a8ae63aff7c7cb2d6a`.
`pattern2` repeats the preceding bundle's two-rank nonuniform primary/IR
restart from a copied checkpoint. All 452608 primary/IR entries and all 512
dust mass, dust energy and gas-energy entries are bitwise identical to the
previous binary's outputs. This is a fixed-mesh regression, **not** runtime
evidence that nonzero regridding has been exercised by RAMSES itself.

The effective `pattern2/effective.nml` retains a fixed level-3 hydro mesh,
`nrestart=1,nstepmax=4,nremap=0`; no stars/AGN, gravity or gas cooling. Output
policy: `noutput=1,aout=2,tout=1e30,foutput=4,fbackup=1000000`, HDF5 input/output.
The schedule times are not reached; two additional steps write one 4.2 MiB
checkpoint, with 173 TiB free at preflight. Original checkpoints are preserved.

No simulation namelist schema changed. Next work remains coarse-fine
transport/reflux and radiation restriction at level synchronization, followed
by a live populated-AMR refinement/coarsening test. Physical dust source
approval and the other previously recorded production limitations are not
waived by this state-transfer result.
