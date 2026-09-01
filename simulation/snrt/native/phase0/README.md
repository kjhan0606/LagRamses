# Phase-0 native contract mirror

This directory is the `/gpfs` working mirror used for G1 contract tests. The
initial files were copied byte-for-byte from
`/home/kjhan/BACKUP/lagRamses/patch/lagRamses` on 2026-09-01; the external
checkout was at commit `c8e17039a02a7822cc2abf312203bc9ec78b9a43` with local
modifications. The external checkout is not modified by G1 work.

The mirror is not yet a production RAMSES source tree. Any source changes
here must be diffed against the recorded external files, compiled in an
isolated test, and explicitly integrated into a clean production checkout
before promotion. `stellar_ramses_runtime.f90` is a patched integration
candidate; its syntax was checked against the existing RAMSES module files,
but it has not yet been linked into a clean production executable.

The G1 contract uses the following canonical units:

- table age: `age_yr` in years;
- runtime age and timestep: `age_gyr` and `dt_gyr` in Gyr;
- RAMSES code time conversion: one explicit `code_time_to_age_gyr` factor;
- returned mass/ejecta: solar masses per initial star;
- energy: erg per initial star;
- momentum: g cm/s per initial star.

The source increment is defined as the cumulative difference
`C(age_gyr + dt_gyr) - C(age_gyr)`, and progress is committed only after a
successful deposition. No net yield is a gas-mass source.

The exact RAMSES layout used by this candidate is recorded in
`simulation/snrt/config/stellar_ramses_field_map_v1.json`. A different
`NVAR`, ATON setting, or delayed-cooling layout requires a new map version and
must not reuse this binary.
