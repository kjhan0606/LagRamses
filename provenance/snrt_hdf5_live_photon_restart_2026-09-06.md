# Live radiation persistence and dust absorption restart

Project `/gpfs/kjhan/LRD_JWST`, remote `kjhan0606/LagRamses`.
Base commit `510f170`; date 2026-09-06.

## Implementation

`snrt_hdf5.f90` connects to the active lagRamses HDF5 dump and hydro-restore
callers. Each `/snrt/level_N` dataset stores grid-major cell records in the
same file ordering as hydro. A record contains presence, HII/HeII/HeIII and
80-by-9 directional photon densities (724 values). Photon floats are
represented exactly as doubles on disk. Cell-to-slot registration occurs
using the restored hydro grid layout rather than serialized RAMSES memory
cell identifiers. Native neutral H is reconstructed from HII.

Same-rank restoration uses per-rank file offsets. Variable-rank restoration
uses the existing AMR `varcpu_grid_file_idx` mapping and 64-grid streaming
chunks; that branch is implemented but not runtime-qualified by this record.
Its replicated chunk reads and the double-precision on-disk photon layout
also need scaling measurements before large production deployment.

The `/snrt` header binds format version, record width, spectral source hash,
commit, approval/status, edge hash, fraction semantics and secondary manifest.
DUST_LIVE additionally binds the numerical opacity/thermal contract, including
heat capacity. Missing radiation data, identity mismatch, malformed extent
and invalid cell payloads trigger a collective failure. Checkpoint payload
restoration is fatal on invalid data; it is not a recover-and-continue API.

This adds no simulation namelist field. `mkrun.py` and
`ramses_nml_generator.py` therefore require no schema change for this bundle.

## Actual executable evidence

Run root `/gpfs/kjhan/LRD_JWST/.snrt-hdf5-live.bSYICX`.
Final binary `ramses_bound3d`, SHA-256:
`228078f765181e058ebebcbf12a135d8aac994b3a97fa965347fe017172df12b`.
Built with Intel MPI/ifx, `DUST_LIVE=1 SNRT=1 USE_CUDA=1 HDF5=1`.
Build logs, effective namelists, run logs and checkpoints are retained there.
The `run_case.sh` runner records the exact reference contract environment.
The default `bin/ramses_final3d` was not overwritten; shared build objects
belong to the live profile, requiring `make -B` when switching profiles.

All positive cases used a noncosmological 8-cubed hydro box, one MPI process,
one A10 GPU (device 1) and one OpenMP thread. There are no active sinks or
stellar sources. The existing analytic reference dust contract remains
NONPRODUCTION.

- `bound`: two source-free steps; one 3.5 MiB HDF5 checkpoint containing
  radiation and nonzero dust. Output policy `noutput=1, aout=2, tout=1e30,
  foutput=2, fbackup=1000000`. Storage estimate under 4 MiB; GPFS had 173 TiB
  free before the run.
- `bound_seeded`: an independent copy of that checkpoint was intentionally
  modified for a numerical test: all 80 directional bins of group 1 received
  `1e-4` code photon density in all 512 leaf cells. The group representative
  energy is 0.1 eV, below the H/He thresholds. This is a test initial radiation
  state, not a physical AGN source prescription. One resumed step uses
  `nrestart=1,nstepmax=3,foutput=3` and writes one new checkpoint.
- `bound_continuous`: starts from the identical modified checkpoint and
  advances two steps with `nstepmax=4,foutput=4`.
- `bound_resumed`: restores the one-step result using
  `nrestart=2,nstepmax=4,foutput=4` and advances the final step.
- `rejected`: a copy of the one-step checkpoint has its stored dust heat
  capacity doubled. The final binary rejected it with exit 10 and
  `SNRT HDF5 radiation checkpoint rejected`; no new dump was created.

Every positive run reported radiation restore (where applicable), coupled
RT transaction/closure PASS and `Run completed`.

## Numerical comparison

HDF5 comparisons used the final executable's outputs, not a Python transport
surrogate. After one seeded step, directional group-1 photon density decreases
from approximately `1e-4` to `6.79966e-5`; mean dust energy grows from
`1.25e-10` to `1.97669e-8` in code energy-density units. Dust mass and gas total
energy remain unchanged. Angular entries are photon counts per directional
bin, so conservation sums bins without applying quadrature weights again.

Using `scale_nH=0.76`, `scale_d=1.66e-24`,
`scale_v=3.0857e21/3.1556926e13`, the maximum cellwise relative difference
between lost photon energy and dust energy gain is `6.6949988e-7`, below
`2e-6`. The input checkpoint photons are converted to FP32 before taking the
loss, matching the live storage representation.

At the end of two steps, the uninterrupted and restarted paths are bitwise
identical for all 370688 leaf radiation/chemistry record values, all 512 dust
mass values, all 512 dust energy values and all 512 gas total energies.
The existing native spectral/checkpoint runner also passed
`SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK` after the state API addition.

## Remaining scope

This establishes actual serial HDF5 radiation persistence and nonzero primary
dust absorption heating. It does not qualify varying MPI count, AMR evolution
or load balancing, AGN pending-fuel restart, binary-format radiation restart,
physical dust/SED approval, IR re-emission or full production dusty feedback.
The existing live AGN serial/fresh-start restriction remains in place. MPI
file-map tests and populated-radiation AMR lifecycle handling remain necessary
before broader simulation-ready promotion.
