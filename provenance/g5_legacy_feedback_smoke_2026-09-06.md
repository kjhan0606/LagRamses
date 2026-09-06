# G5 legacy-feedback comparison smoke

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)

## Execution

The bounded comparison used two MPI ranks and the CPU production-linked
binary:

```text
binary=/gpfs/kjhan/LRD_JWST/bin/ramses_final3d
binary_sha256=5eaef4bdc321abbe2494a481d229b7170221c0be623d572f722b3f326253bff2
effective_namelist=/gpfs/kjhan/LRD_JWST/.legacy-feedback-smoke.SSbLrO/effective.nml
effective_namelist_sha256=5d09686cd0482d3300188596339c552a0a82fb1e3949bc36555dcc9cac1bcb34
mpi_ranks=2
nstepmax=2
SNRT_RT_ENABLE=unset
```

Result: `LEGACY_FEEDBACK_SMOKE_PASS`. The log contains
`Stellar feedback mode: legacy`, repeated finite-state `NaN_CHK` diagnostics,
and `Run completed`. The intentional future output policy produced no
`output_*` directory. The case directory and log are retained under
`/gpfs` for local inspection.

## Interpretation

This certifies the small legacy comparison path as an executable RAMSES
hydro/feedback smoke. It does not certify physical channel-resolved yields,
stellar/AGN source spectra, SNIa runtime activation, or dust; those remain
blocked by the production manifest and their individual contracts.
