# F-P2.5 first-audit remediation record — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Related audit: `claude_opus5_fp2_5_native_hhe_thermochemistry_bundle_end_audit_2026-09-05.md`
Status: implementation remediation complete; Opus follow-up was `CONDITIONAL PASS`.
The three record-only closure conditions from that follow-up are applied below.

## Dispositions

- HIGH-1: closed in the native path by explicit absolute/relative inventory
  tolerances, tolerance-sized clipping in `snrt_partition_absorption`, a
  named optional `unassigned_absorption_code`, and a `0.99995` FP32 guard band
  in the species-aware CUDA cap.
- HIGH-2: closed in the production adapter by passing a
  `(leaf,group,species)` optical-depth tensor and `(leaf,species)` inventory to
  CUDA. One cell thread consumes the shared H I/He I/He II inventory in
  deterministic group order across all groups and all transport substeps.
  Both CUDA and Fortran redistribution paths require positive group opacity.
  The old scalar-budget entry point remains only for compatibility benchmarks;
  `snrt_ramses_driver` uses the species-aware ABI.
- MODERATE-3: closed by checking raw `fion+fheat+fexc` against unity for all
  258×14 FS2010 table entries before setting the table-loaded flag.
- MODERATE-4: closed by exposing radiative and dielectronic He II terms and
  checking their sum and non-negligible dielectronic contribution at `10^5 K`.
- MODERATE-5: closed by checking that secondary H ionization in a saturated H
  II cell is zero and its unavailable energy is routed to heat.
- MODERATE-6: closed by storing authoritative H II/He II/He III and H I mirror
  fractions in `real(dp)` and rejecting a non-simplex or inconsistent payload.
- MODERATE-7: closed by checkpoint version 6 metadata for the FS2010 source id,
  upstream commit, and manifest identity, plus a pre-mutation mismatch smoke.
- LOW-10: closed by adding hydro module prerequisites to the driver target.
  LOW-9/11/12/13 remain non-blocking optimization/diagnostic follow-ups.

## Native verification

All checks below ran from `/gpfs/kjhan/LRD_JWST`; no RAMSES evolution or large
job was launched.

```text
bash simulation/snrt/tests/run_snrt_native_thermochemistry.sh
  SNRT_NATIVE_THERMOCHEMISTRY_OK
  SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK

FC=mpiifx bash simulation/snrt/tests/run_snrt_native_thermochemistry.sh
  SNRT_NATIVE_THERMOCHEMISTRY_OK
  SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK

bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh
  SNRT_SPECTRAL_CONTRACT_OK
  SNRT_CHECKPOINT_OK
  SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK

FC=mpiifx bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh
  SNRT_SPECTRAL_CONTRACT_OK
  SNRT_CHECKPOINT_OK
  SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK

bash simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
  SNRT_CUDA_MULTIGROUP_OK relative_budget_error=5.596245E-09
  SNRT_NATIVE_CUDA_MULTIGROUP_ALL_OK

make -C bin SNRT=1 USE_CUDA=1 ramses
  exit status 0; linked bin/ramses_final3d
```

The CUDA smoke ran on an NVIDIA A10 (compute capability 8.6). It verifies
global photon conservation, group-total versus directional reduction,
non-negative state/inventory, group-order inventory consumption, and zero
absorption in groups whose only nonzero opacity belongs to an unavailable
species.

## Remaining boundary

The Opus follow-up confirmed the species-aware CUDA/Fortran ABI and the
checkpoint/table hardening. It found no further source remediation required;
the three record-only conditions were:

1. document the diagnostic-only `unassigned` receiver limitation;
2. document that the ascending-group greedy depends on the nested eligibility
   invariant enforced by `validate_species_table`, including its priority bias;
3. reconcile the spectral-contract transcript (the full transcript includes
   `SNRT_SPECTRAL_CONTRACT_OK` before the checkpoint markers).

All three records are now present in this remediation record and
`SNRT_NATIVE_GROUP_CONTRACT.md`. F-P2.5 may close as a conditional pass, but it
still does not approve a physical AGN or stellar SED, a live RT+feedback
evolution, global implicit opacity/chemistry convergence, cooling receiver
integration, dust/radiation pressure, HDF5 restart integration, or the
40--120 M_sun yield seam.
