# DUST-7: fourth-species CUDA dust transport boundary

Project `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`.
Base commit: `ad935dc` (DUST-6 native primary dust coupling).
The operator has pre-approved continued bundle work. The project objective is
production/publication-ready high-level RAMSES RT with stellar/AGN feedback
and dust physics; this bundle is one implementation boundary toward that
objective, not a production qualification.

## Purpose and scope

DUST-6 defines the FP64 reference partition and heating ledger, but the live
CUDA transport still has a three-species H/He ABI. DUST-7 will add a separate
fourth-species CUDA ABI for primary dust absorption while preserving the old
three-species entry point byte-for-byte. It will establish a cross-language
transport ledger suitable for a later RAMSES transaction adapter.

This bundle does not select or approve a dust opacity mixture, change the
existing driver, activate live RAMSES dust heating, or add AMR/MPI/restart
qualification. Those are DUST-8/live-integration and longer-term science
work. No new Python validator, JSON/HDF5 asset, global gate or diagnostic
framework is planned; extend the existing native CUDA smoke only.

Fable pre-bundle review: **CONDITIONAL APPROVE**. The review confirmed that
deferring live wiring to DUST-8 keeps FP32 kernel validation attributable. Its
conditions below are part of this amended plan.

## ABI and physical contract

Add a new `bind(C)` entry point alongside the existing ABI. Its inputs are
explicit FP32 caller arrays:

- total group optical depth `tau_total`;
- the three H I/He I/He II component optical depths;
- dust optical depth as a separate group array; and
- the three finite H/He inventories. Dust has no inventory array and is never
  depleted.

The host wrapper rejects null/mismatched pointers, invalid dimensions,
negative timestep CFL, non-finite/negative input arrays, and a total optical
depth inconsistent with the sum of the four supplied components beyond
`8*FLT_EPSILON*max(tau,component_sum,FLT_MIN)`. A device validation flag is
checked before any output/state update. The new kernel retains the existing
upwind/absorption layout and consumes each cell's H/He inventory in group
order. The legacy kernel and its shared implementation body are not
refactored.

For each cell/group, after raw transport removal:

1. Set `HHe_target = raw_removed*tau_HHe/tau_total` and allocate the H/He
   target with the existing bounded active-set rule.
2. Set direct dust allocation to `raw_removed-HHe_target`.
3. Preserve the existing FP32 `0.99995` inventory guard band: the guard-band
   slice is returned, not treated as physical excess. If H/He truly
   saturates, transfer only `HHe_excess*(1-exp(-tau_dust))` to dust and return
   the rest to the photon field, matching the DUST-6 reference contract.
   Use `expm1f` in the weak-dust branch.
4. Scale directional removal by the assigned total, add the returned remnant
   to the surviving photon state, and update only H/He inventories.

The new ABI returns raw removal, assigned H/He species by group, assigned dust
by group, returned photons by group, assigned total by group/cell, and the
updated H/He inventories. CUDA `unassigned` is identically zero; the
independently testable closure is
`raw_removed = HHe_assigned + dust_assigned + returned` within FP32 tolerance.
The assigned-total output is explicitly not to be fed to the host H/He
partition as raw removal.
The CUDA ledger is transport-only; conversion to physical heating remains at
the DUST-6/Fortran unit boundary. No dust energy is written to gas or
momentum in this bundle.

## Implementation and evidence

- Add the new Fortran `bind(C)` declaration and CUDA implementation without
  changing the old C symbol or its three-species memory contract.
- Reuse the existing `snrt_cuda_multigroup_smoke.f90` and runner. Add a
  manufactured dust case with proportional allocation, H/He saturation,
  returned-photon closure, shared nine-group inventory monotonicity, dust-only
  groups below 13.6 eV, bitwise zero-dust equivalence to the existing ABI,
  component-tau mismatch rejection and non-negative outputs. Keep the two
  existing legacy cases and their recorded markers unchanged.
- Compile with the existing GNU/Intel Fortran paths and `nvcc`. If a CUDA
  device is unavailable, print a distinct unavailable marker and exit
  nonzero; compile-only status is not the DUST-7 exit condition and DUST-8
  cannot start from it. Make the tested CUDA architecture overridable and
  record the device model in the evidence.
- Record ABI layout, units, tolerance and device availability in one bundle
  evidence file. Do not add a second audit or a new dataset for this bundle.

## Exit and handoff

Exit is a compiled, auditable fourth-species CUDA reference boundary with
FP32 closure evidence, dust-only-group behavior and the legacy ABI preserved.
DUST-8 may then add the
optional transport adapter and live transaction wiring, where it must prove
FP32-to-FP64 agreement, primary photon/H/He/dust heating closure and
transaction rollback across AMR/MPI boundaries.

Deferred: opacity-mixture provenance and evolution, photoelectric heating,
IR scattering/recursive state, live RAMSES activation, AMR/MPI/restart,
reduced-c interpretation, production cost and publication qualification.
