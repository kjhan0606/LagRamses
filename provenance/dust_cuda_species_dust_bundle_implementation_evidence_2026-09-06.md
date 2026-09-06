# DUST-7 implementation evidence: fourth-species CUDA dust boundary

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Base commit: `ad935dc` (DUST-6 native primary-photon coupling)
- Scope: native CUDA transport boundary only; no live RAMSES activation
- Plan: [DUST-7 plan](dust_cuda_species_dust_bundle_plan_2026-09-06.md)
- Plan audit: [Fable conditional approve](fable_dust_cuda_species_dust_bundle_plan_audit_2026-09-06.md)

## Implemented boundary

`snrt_cuda_multigroup_rt_step_species_dust_c` is a new `bind(C)` entry
point. The existing `snrt_cuda_multigroup_rt_step_species_c` symbol and its
three-species kernel body are unchanged. The new caller contract is:

- Fortran arrays are contiguous `(cell,group,species)` and are consumed in
  CUDA as `(species,group,cell)` for the H I/He I/He II component arrays.
- `tau_total`, the three H/He component optical depths, and `tau_dust` are
  FP32 group/cell arrays. The wrapper rejects a component mismatch beyond
  `8*FLT_EPSILON*max(tau_total,component_sum,FLT_MIN)`.
- H/He inventories are finite and depleted in group order. Dust has no
  inventory and is non-depleting in this boundary.
- Outputs are independent `raw_group`, `absorbed_hhe_species`,
  `absorbed_dust_group`, `returned_group`, `absorbed_group` (assigned total),
  and `absorbed` (assigned cell total), plus updated H/He inventory. CUDA
  has no unassigned output; the contract sets it to zero by omission.
- The legacy FP32 `0.99995` inventory guard band is returned to the photon
  field, while only true finite-inventory H/He excess is eligible for the
  bounded dust transfer `HHe_excess*(1-exp(-tau_dust))`. The weak-dust
  fraction uses `expm1f`.
- Invalid state/direction/neighbor/tau/inventory input is rejected by a
  device validation pass before any host state, inventory, or output is
  copied back. Failure leaves those caller arrays unchanged.

The transport state is updated by the existing upwind and absorption kernels;
the new kernel only performs the fourth-component cap and ledger. The
zero-dust branch retains the old cap arithmetic and reproduces the legacy
scalar/group reduction order, allowing a bitwise regression of state,
absorption, and H/He inventories. No dust heating, gas energy, momentum,
AMR/MPI transaction, restart, opacity-mixture choice, or live driver wiring
is included; those remain DUST-8 or longer-term work.

## Native CUDA evidence

The existing native smoke was extended in place, preserving its two prior
markers and adding dust-only groups 1--4, mixed H/He+dust groups 5--9,
finite-inventory saturation, returned-photon closure, zero-dust bitwise
equivalence, host-side checks of the direct-dust plus finite-excess split,
a non-saturating mixed group, and invalid-input transaction checks. The
invalid tests separately cover negative dust tau and an all-nonnegative
component-total mismatch. The runner now reports
the GPU model and accepts an overridable `CUDA_ARCH` (default `86`); a
missing runtime device prints `SNRT_NATIVE_CUDA_UNAVAILABLE` and exits
nonzero, so a compile-only result cannot be reported as a pass.

Commands were run from the project root:

```text
simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
FC=/opt/ohpc/pub/compiler/gcc/13.2.0/bin/gfortran simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
CUDA_VISIBLE_DEVICES=-1 simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
```

The Intel mpiifx/ifx run and GNU gfortran run both completed on an NVIDIA
A10 (`sm_86`) with these markers:

```text
SNRT_CUDA_MULTIGROUP_OK relative_budget_error=  5.365561E-09
SNRT_CUDA_MULTIGROUP_SPECIES_MIX_OK relative_budget_error=  2.483527E-09
SNRT_CUDA_MULTIGROUP_SPECIES_DUST_OK relative_budget_error=  1.192093E-08
SNRT_CUDA_MULTIGROUP_SPECIES_DUST_NON_SATURATING_OK
SNRT_CUDA_MULTIGROUP_SPECIES_DUST_ZERO_DUST_BITWISE_OK
SNRT_CUDA_MULTIGROUP_SPECIES_DUST_INVALID_INPUT_OK
SNRT_NATIVE_CUDA_MULTIGROUP_ALL_OK
```

The deliberate no-device run emitted:

```text
SNRT_NATIVE_CUDA_UNAVAILABLE
```

and exited with status `90`, without printing the all-OK marker. `git
diff --check` passes. The device line is informational: `nvidia-smi` is
queried independently of CUDA's visibility mask, while runtime availability
is determined by `snrt_cuda_available_c` and the smoke itself.

The bitwise legacy comparison and the `fmaxf` removal clamp assume the
transport step supplies non-negative photon amounts, as required by a
CFL-respecting physical state. DUST-7 records caller photon amounts and
FP32 closure tolerances; it does not claim FP32/FP64 identity. The
`0.99995` guard and FP32 residual reconciliation remain explicit DUST-8
comparison items.

## Source hashes at evidence capture

```text
5376df67d586c4c777c5d8601779f1430c887ed959ce0ef770235ce22e2ec3c1  patch/lagRamses/snrt_cuda_kernels.cu
4e2fb01bc875a57cd92aa8cbdf699aca92284be8c846d6fa194bf15e42aa1732  patch/lagRamses/snrt_cuda_multigroup_interface.f90
d8c099e36a80b23263668a5952286a8ffa526dfdffe94c27d78163fa27a3c165  patch/lagRamses/snrt_cuda_multigroup_smoke.f90
0d110a61199546f210514c1b472487a380f7582c1b00b3406d13c2b93c5a09f1  simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
```

## Disposition

DUST-7 is a native, auditable fourth-species CUDA reference boundary with
legacy ABI preservation and actual GPU evidence. It is not production or
publication approval. DUST-8 may wire these ledgers into the FP64 RAMSES
transaction, where it must prove FP32/FP64 conservation, primary
photon/H/He/dust heating closure, and rollback across the live AMR/MPI path.
