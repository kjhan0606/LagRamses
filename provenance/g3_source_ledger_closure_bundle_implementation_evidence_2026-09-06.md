# G3 source-ledger closure bundle implementation evidence

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Implementation

`simulation/snrt/tools/merge_photon_source_ledgers.py` now emits the explicit
mixed source kind, group count, source-ID policy, and source epoch disposition.
It rejects source-ID collisions before writing output, rejects component
ledgers with different scale factors, and marks missing scale-factor metadata
as `missing_source_scale_factor_not_production_eligible`. A coeval scale factor
is copied into the merged metadata only after all components declare the same
value.

The aggregate remains a candidate ledger. Its dust contract continues to
require a source-weighted aggregate dust sidecar; a STAR-only or AGN-only
sidecar cannot silently be reused.

## Gate

The consolidated gate is:

```text
simulation/snrt/tests/run_g3_source_ledger_closure.sh
```

Executed from `/gpfs/kjhan/LRD_JWST` with the project CPU JAX environment
(`jax==0.11.1`, `jaxlib==0.11.1`):

```text
MERGE_PHOTON_LEDGER_TEST_OK sources=2 groups=9 closure=aggregate
MERGE_PHOTON_SOURCE_LEDGERS_TEST_OK components=2 sources=3 mixed_dust_gate=1
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
G3_SOURCE_LEDGER_CLOSURE_PASS tests=3 backend=cpu
```

The gate covers normal aggregate closure, exact group-edge agreement, CSV /
metadata total agreement, duplicate source-ID rejection, mismatched epoch
rejection, source-SED identity binding, and dust source-binding integrity.
No RAMSES run, physical source activation, or production data modification was
performed.

## Disposition

G3 source-ledger closure is **PASS within scope**. The physical SED/yield,
obscuration, DTD, dust mixture, and live feedback gates remain unresolved and
are not promoted by this result.
