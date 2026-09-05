# F-P2.5 native H/He thermochemistry — implementation evidence

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Base checkout: `bd0411bf9c2a` (`main`, worktree intentionally contains the
operator's in-progress bundle changes)
Status: implementation plus first-audit remediation complete; bundled Claude
Opus 5 follow-up audit returned `CONDITIONAL PASS`, with its record-only
conditions closed in the contract and remediation records.

## Native implementation delivered

- `patch/lagRamses/snrt_thermochemistry.f90` now loads the pinned 14-file,
  258-energy FS2010/21cmFAST table set through a strict Fortran namelist
  contract, performs the reviewed bilinear interpolation and boundary
  semantics, and normalizes the five deposition channels.
- The same module advances H II, He II, and He III. It partitions each
  transport-returned group absorption against H I/He I/He II inventories,
  charges unavailable secondary-ionization energy back to heat, and solves
  the bounded electron-density case-B recombination closure. He III uses
  `2 alpha_H,B(T/4)`; He II includes the reviewed radiative plus dielectronic
  case-B coefficient.
- `snrt_nlte_coupling.f90` computes species-resolved H I/He I/He II optical
  depths and their total. `snrt_ramses_driver.f90` passes group/species optical
  depth and the start-of-step H I/He I/He II inventory to the species-aware
  CUDA cap, partitions its returned absorption, updates the persistent state,
  and adds only the native gas-heating ledger to RAMSES total energy.
- `snrt_state.f90` carries the new fractions in checkpoint version 6 as
  double-precision authoritative state, binds the FS2010 source/commit/manifest
  identity, and rejects invalid or H I/H II-inconsistent payloads. The
  checkpoint smoke tests both spectral and secondary-identity rejection and
  round-trips all intensity and H/He state fields.
- `snrt_transport_step.f90` preserves the species inventory across groups and
  transport substeps; CUDA redistribution is masked by positive group opacity.
- `bin/Makefile` includes the thermochemistry module and its dependency graph
  in the SNRT production link.

## Native evidence

The following are native Fortran/shell checks; no Python or JAX is required.

```text
bash simulation/snrt/tests/run_snrt_native_thermochemistry.sh
  SNRT_NATIVE_THERMOCHEMISTRY_OK
  SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK

FC=mpiifx bash simulation/snrt/tests/run_snrt_native_thermochemistry.sh
  SNRT_NATIVE_THERMOCHEMISTRY_OK
  SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK

bash simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
  SNRT_CUDA_MULTIGROUP_OK relative_budget_error=5.596245E-09
  SNRT_NATIVE_CUDA_MULTIGROUP_ALL_OK

bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh
  SNRT_SPECTRAL_CONTRACT_OK
  SNRT_CHECKPOINT_OK
  SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK

FC=mpiifx bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh
  SNRT_SPECTRAL_CONTRACT_OK
  SNRT_CHECKPOINT_OK
  SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK

make -C bin SNRT=1 USE_CUDA=1 ramses
  exit status 0; linked bin/ramses_final3d
```

The FS2010 runner also reports `OK` for all fourteen vendored table files
against `TABLE_MANIFEST.sha256`. The thermochemistry smoke covers the pinned
200 eV, `x_HII=0.1` reference, the 99.9/100.1 eV continuity check, the table
floor behavior, raw table-channel closure at load, all three species
partition/inventory caps including opaque-zero masking, unavailable-target
heating, non-negative recombination ledgers, the H/He simplex, the exact He
III coefficient identity, the hot-temperature He II dielectronic term, and
the photoelectron energy residual. The CUDA smoke verifies group-order
inventory consumption and the absence of low-energy absorption by an
opaque-zero species.

The linked binary exports the native entry points
`snrt_thermochemistry_advance_cell`,
`snrt_nlte_primordial_optical_depth_groups`,
`snrt_state_checkpoint_read/write`, and
`snrt_ramses_advance_level`.

## Source hashes at evidence capture

```text
621f7d816a85286947da720a00f4203f34bcc9e28598439865b6a1e2f0b0abf2  patch/lagRamses/snrt_thermochemistry.f90
90c96f0caa579d1db814afac2daad7c65882e9cce290f9fbc756b9717274f1ae  patch/lagRamses/snrt_thermochemistry_smoke.f90
58e1760e46c4151656606f1fe6306ed67b4306973dca5f74ec4beca54d20f6dd  patch/lagRamses/snrt_thermochemistry_loader_smoke.f90
4c06658fae1caba76f6ad8e150874e00eeff867e0b42964d99014be32253feb5  patch/lagRamses/snrt_nlte_coupling.f90
fa89f25655f9496d8bb05afb201219da24a409212440fa150eab5da6d6dbd102  patch/lagRamses/snrt_ramses_driver.f90
fa53625d4ffd1950aeb65b55bd77abf035b292e1eb36eb75501f0fceb523bb94  patch/lagRamses/snrt_transport_step.f90
6deec7b53397f74c8d7a65fa3042c845cd0dd544dbdd0742b826b303a4f5024d  patch/lagRamses/snrt_state.f90
c16e54d46c663c0d369e8630eaa9b7a79ad888009fb1a1c75027e7daa50e9a56  patch/lagRamses/snrt_checkpoint_smoke.f90
9534e026876955d790c42ce9b5323538e56fa830d99836e6a4fb6de9ec38b9af  patch/lagRamses/snrt_cuda_kernels.cu
86b52d8b9c4bac2bfa85a6d9ea894833a7f7a5d0945ec1c25b43391e1d4621cb  patch/lagRamses/snrt_cuda_multigroup_interface.f90
0121964d552278beffbc820385e3a77263ab4813caf11bb31d984cec85b688a5  patch/lagRamses/snrt_cuda_multigroup_smoke.f90
12ea9422d786394277332b53e1ccab3a264d3a5aaede21dc29271183b160345f  bin/Makefile
df3c3d15da2e9af68def935dcabe5bd270318082d458c53c848538b911e4e9a8  simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
558cbc6114ba31d45dd068fc19f131499fecf366e09cba80d6a7e648aa43854c  simulation/snrt/config/snrt_secondary_table_contract_v1.nml
10153794818a40eb38468cbd208652b441eaa16712904a8e0ba3aa914358135f  simulation/snrt/data/furlanetto_stoever_2010/TABLE_MANIFEST.sha256
3c916d00c06fd8bce8655c2f74f3522b2edc75396657f11a7f92e2d2cf221cfc  bin/ramses_final3d
```

## Scope boundary

No RAMSES evolution or production hydro/AGN run was launched. The checked-in
reference-control SED remains a wiring control, not a physical SED approval.
This bundle does not yet provide the global implicit opacity/chemistry fixed
point, collisional/background/Compton/metal cooling, dust/radiation pressure,
recombination-line transport, stellar photon coupling, HDF5 restart wiring, or
publication-scale RT/feedback convergence.
