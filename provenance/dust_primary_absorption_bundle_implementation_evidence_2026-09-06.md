# DUST-6 implementation evidence: native primary dust coupling

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Base commit: `578bb3b`
- Scope: native reference contract only; no live RAMSES/CUDA activation
- Plan: [DUST-6 plan](dust_primary_absorption_bundle_plan_2026-09-06.md)
- Plan review: [Fable conditional approve](fable_dust_primary_absorption_plan_audit_2026-09-06.md)

## Implemented boundary

`patch/lagRamses/snrt_dust_coupling.f90` provides four FP64 caller-owned
operations:

1. Explicit `tau_dust = sigma_dust[cm2/H] * dust_relative * n_H[cm-3] *
   path[cm]`, with finite/non-negative admission checks.
2. Total optical depth from the three H/He components plus dust.
3. Group-wise partition of a raw removed photon amount. Dust is a
   non-depleting absorber. H/He is sent through the existing finite-inventory
   partition; direct dust allocation is opacity-weighted, and only
   `HHe_excess*(1-exp(-tau_dust))` is transferred from a saturated H/He
   target. The remaining physical remnant is returned to the photon field.
   This transfer is an explicitly bounded sub-step approximation; Opus notes
   a mild high bias in the resulting dust heating for that residual path and
   it is not a dust-mixture claim.
4. A separate heating ledger converting
   `(photons cm-3)*(mean eV)*1.602176634e-12/dt[s]` to `erg cm-3 s-1`.

The group closure is
`raw_removed = HHe_assigned + dust_assigned + returned + unassigned`.
`returned` is a physical photon-field restoration; `unassigned` is only a
non-negative floating-point residual. The existing three-species CUDA ABI,
primary transport driver, `snrt_intensity`, H/He state, gas energy and
momentum are unchanged.

## Native evidence

The existing `snrt_dust_ir_smoke.f90` was extended rather than adding another
test framework. `dust_ir_transport.py --native` now compiles, in dependency
order, the thermochemistry, DUST-5 IR operator, DUST-6 coupling module and
the smoke with both available compilers.

Command (from the project root):

```text
simulation/snrt/.venv/bin/python simulation/snrt/tests/dust_ir_transport.py --native
```

Result: PASS for GNU gfortran 13.2 and Intel ifx 2025.3. The original DUST-5
differential remains at machine precision:

| compiler | stored/escaped/reprocessed relative errors | temperature error | balance |
|---|---:|---:|---:|
| gfortran | 1.11e-15 / 0 / 3.33e-16 | 4.44e-16 | 7.407e-10 |
| ifx | 6.66e-16 / 4.44e-16 / 2.22e-16 | 0 | 7.407e-10 |

The smoke also reports `NATIVE_DUST_COUPLING_OK` and exercises proportional
allocation, saturation transfer, a sequential three-group shared reservoir,
exact zero-dust regression against `snrt_partition_absorption` in the
non-saturating regime, returned-photon conservation, zero/weak heating,
negative/overflow/shape/no-optical-depth input rejection, and caller reservoir
nonmutation on failure. In the zero-dust saturation regime, the existing
H/He routine rejects an over-inventory request while this adapter deliberately
returns the physically unassignable photon remnant; that semantic difference
is explicit and is not an exact regression claim. The default non-native regression also passes:
`DUST_IR_TRANSPORT_TEST_OK`, `DUST_SPECTRAL_IR_TEST_OK`.

The inner H/He routine retains its existing inventory-tolerance convention,
now anchored to the pre-partition H/He inventory. The outer DUST-6 closure
uses a separate 256-host-ulp check without an absolute-one floor; this keeps
allocation tolerance and ledger closure tolerance distinct.

`git diff --check` passes. Source hashes at evidence capture:

```text
608a3dc9c49cf0ed3ded86f7c8eb13f9691f9a9ec5ec1257a690a99b041275ae  patch/lagRamses/snrt_dust_coupling.f90
d372c2c793458228bfaddef14f86065ae0e427ff03f0eff49f4acbe021e66b44  patch/lagRamses/snrt_dust_ir_smoke.f90
d6fc98bb121bd1f9e86a822f5cf7306731df603fd7678a30a57117f04e9d82cb  simulation/snrt/tests/dust_ir_transport.py
```

## Disposition

This is an auditable native coupling boundary, not a production or
publication approval. The next bundle may define a fourth-species CUDA ABI
and wire it into the transaction, subject to FP32-vs-FP64 conservation and
AMR/MPI evidence. Dust opacity-mixture provenance, depletion/evolution,
photoelectric gas heating, IR scattering, persistent spectral state, live
driver integration, restart and cosmological qualification remain deferred.
