# DUST-9 implementation evidence: source-bound mapping and thermal receiver

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Base workspace commit before this bundle: `b30f004ed5179e46940a7281cbe7e541547c21e4`

## Implemented boundary

- `patch/lagRamses/snrt_dust_receiver.f90` adds a native FP64 boundary for:
  - numerical validation of source-bound group edges, per-H cross sections,
    representative energies, binding status, source identity, and SHA-256
    tokens;
  - explicit `metallicity_solar * dust_to_metal` cell mapping;
  - `tau_dust = n_H * path * dust_relative_abundance * sigma_per_H` per cell
    and group;
  - staged persistent dust thermal energy/temperature using caller-supplied
    volumetric heat capacity;
  - commit only after complete state validation.
- The native boundary deliberately receives identity/hash tokens from an
  upstream validated sidecar loader.  It does not parse JSON or treat a token
  check as a replacement for file/hash verification.
- `bin/Makefile` now includes `snrt_dust_receiver.o` in the SNRT production
  module graph.  The live driver remains explicitly `ZERO_SCAFFOLD`; no
  dedicated RAMSES dust state exists in `hydro_commons`.

## Native smoke

Command:

```text
simulation/snrt/tests/run_snrt_native_dust_mapping_receiver.sh
```

Result with the available compilers:

```text
SNRT_NATIVE_DUST_MAPPING_RECEIVER_OK binding=1 mapping=1 opacity=1 thermal=1 closure=1 rollback=1
SNRT_NATIVE_DUST_MAPPING_RECEIVER_IFX_PASS
SNRT_NATIVE_DUST_MAPPING_RECEIVER_OK binding=1 mapping=1 opacity=1 thermal=1 closure=1 rollback=1
SNRT_NATIVE_DUST_MAPPING_RECEIVER_GNU_PASS
SNRT_NATIVE_DUST_MAPPING_RECEIVER_RUN_PASS
```

The fixture covers accepted and rejected source binding, explicit abundance
mapping, optical-depth arithmetic, positive absorbed-energy staging, exact
energy closure, commit, positive absorption with zero abundance, zero-dust
no-op, and rejection of an invalid candidate without mutating persistent
state.

## Production gate

The existing consolidated native gate was rerun after adding the module to the
production graph:

```text
STAGE production_build status=PASS elapsed_s=196.833
STAGE agn_partition_reference status=PASS elapsed_s=1.773
NATIVE_SYMBOLS_CHECK count=5 status=PASS
STAGE dust_ledger_receiver status=PASS elapsed_s=0.397
STAGE thermochemistry status=PASS elapsed_s=0.923
STAGE spectral_contract status=PASS elapsed_s=1.520
STAGE transaction_mpi status=PASS elapsed_s=4.811
STAGE cuda_multigroup status=PASS elapsed_s=5.104
STAGE production_negative status=PASS elapsed_s=6.176
STAGE diff_check status=PASS elapsed_s=0.097
SNRT_BUNDLE_GATE_PASS
```

The full gate's commit line identifies the clean base commit because the
DUST-9 source was tested from the working tree before commit.  The production
build nevertheless consumed the new object through the modified Makefile;
the focused DUST-9 runner is the source-level evidence for the new API.

Source SHA-256 before commit:

```text
snrt_dust_receiver.f90 454cd34b9059ad4049f1a62a586b2962738bf6078db3aafcc629c248d1756505
snrt_dust_receiver_smoke.f90 19522c7d33c102e39108e36592ca9a71ae9a716b2d6f9a70ca0b9c788c0af618
run_snrt_native_dust_mapping_receiver.sh 26a43120f5bcdbbefa75248215034a23252ccb8eeb48a4d16485eb4bb60e6337
```

## Explicit non-claims

This bundle does not activate nonzero dust in the live RAMSES driver.  It does
not add a dust field to `uold`, a dust momentum receiver, IR re-emission,
scattering, grain growth/destruction, restart/migration payload, or a physical
temperature-dependent grain heat-capacity table.  Those remain required
gates for production/publication dust feedback.
