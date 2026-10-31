# Combined stellar–AGN M5 CPU smoke — 2026-11-01

## Result

Grammar Slurm job `1252166` completed with exit code `0` on `grammar036`
(`normal`): 1 MPI rank, 8 OpenMP threads, CPU-only CHIMES, 32 GiB requested,
5:33 Slurm elapsed. The runner completed four coarse steps with M5 order 5,
live dust, and CHIMES. It reported
`M5_COMBINED_SOURCE_ATTRIBUTION_PASS overlap_commits=1`, followed by
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS`. One M5 commit received positive
stellar and thermal-AGN photon receipts together:

- Stellar photons: `7.71699705e66`
- Thermal-AGN photons: `1.56436650e70`
- AGN packets: `1`
- Momentum assertion cells: `256`
- Raw RAMSES outputs: `0`

The test used the `effective_population` profile and the `thermal_seed_reference`
AGN IC. The effective stellar package identifies PARSEC mixed low-mass evolution
with the Kroupa IMF and is marked `reference_control_not_production_approval`.
The AGN IC starts with a sink mass of about `2.45e11 Msun`; this deliberately
strong source seed is a wiring control, not a physically representative seed.
This run establishes source-path overlap and CPU coupling only. It does not
validate physical source normalization, production calibration, or GPU parity.

## Build and retained evidence

The CPU executable was built from LagRamses commit
`ee15548d6cf18064afdda507a95b15126879271e` with `CHIMES_CUDA=0`, the reference
CHIMES library, and no CUDA runtime dependency. The executable SHA256 is
`c0dd0a197db704ab696fe8936c9c3af01d347ffb491a2eb2ca88ca9d84e58877`; the
CHIMES library SHA256 is
`6f7abbd129f5382a7191b9352f0d1a4992b9faf0190ac42bba870461220cab37`.

The effective namelist is
`.m5-effective-dust-chimes-grammar-cpu-1252166/effective.nml` (SHA256
`01e623647752e38b00ea6dde9fb2be95925c923a81ac5453e5dc4f93c67a4104`). It uses
one future scheduled output (`noutput=1`, `aout=2.0d0`, `tout=1.0d30`) and
`foutput=fbackup=1000000`; the runner confirmed no raw output directories.
Inputs and data hashes are retained in the run directory's `inputs.sha256`,
with the compact result, source attribution, AGN state, and run log.

Slurm recorded 26,391,636 KiB peak RSS for the batch (build included) and
729,520 KiB for the application step. RAMSES reported 90.965 seconds for the
four-step run, of which 89.585 seconds (98.5%) fell in the aggregate `cooling`
timer. This is a small CPU reference smoke, not a controlled performance
benchmark; `cooling` includes the coupled M5/material path and is not a
standalone CHIMES measurement.

## Outstanding GPU verification

This CPU run was submitted while the A40 queue was blocked by the per-user GPU
QoS cap. Source review then found that test command `rhs` exercises direct
cell-level RHS parity, while the production `level_queue` comparison is in
test command `tile`. The GPU script now runs both modes; the corrected script
was committed as `9803e3a`. The original pending jobs `444877` and `444882`
were canceled before execution and replaced with `444910` (after-any
`444740`) and `444911` (after-success `444910`). Parent job `444740` remains
pending on `QOSMaxGRESPerUser`. This CPU result does not replace either GPU
verification.
