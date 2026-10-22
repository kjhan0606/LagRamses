# Evaluated raw snapshot cleanup

The following 37 small test dump directories were evaluated and selected
for removal. They contain only this agent's Grackle validation outputs.
Effective inputs, logs, executable identity and results.json are retained.
No active restart depends on these dumps; regenerated runs are possible,
but removed raw snapshots themselves have no retained backup here.

- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r1/omp1/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r1/omp1/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r1/omp1/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r1/omp1/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp1/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp1/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp1/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp1/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp2/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp2/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp2/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/omp2/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/mpi2/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/mpi2/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/mpi2/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/mpi2/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr1/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr1/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr1/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr1/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr2/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr2/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr2/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr2/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed1/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed1/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed1/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed1/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed2/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed2/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed2/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/amr_fixed2/output_00004`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/checkpoint/output_00001`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/checkpoint/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/restart/output_00002`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/restart/output_00003`
- `/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-validation-20261006-r2/restart/output_00004`

## Follow-on validation outputs evaluated and removed

After the 2026-10-06 extended validation bundle was analyzed, only the HDF5
raw-output directories below were removed. Effective namelists, logs, binary
checksums, test sources and scalar results remain in their run directories and
in `grackle_neq_validation_2026-10-06.md`.

- `.grackle-neq-dynamic-amr-20261006/omp1/output_00001`–`output_00008`
- `.grackle-neq-dynamic-amr-v2-20261006/omp1` and `/omp2`, each
  `output_00001`–`output_00008`
- `.grackle-neq-dynamic-amr-v3-20261006/omp1` and `/omp2`, each
  `output_00001`–`output_00008`
- `.grackle-neq-dynamic-amr-v4-20261006/omp1` and `/omp2`, each
  `output_00001`–`output_00008`
- `.grackle-neq-dynamic-amr-v5-20261006/omp1/output_00001`–`output_00010`
- `.grackle-neq-mpi-amr-20261006/omp1` and `/mpi2`, each
  `output_00001`–`output_00004`
- `.grackle-neq-delayed-v2-20261006/output_00001`–`output_00002`
- `.grackle-neq-delayed-control-v2-20261006/no_reservoir/output_00001`–
  `output_00002`
- `.grackle-neq-cosmo-ic-snapshot-20261006/output_00001`
- `.grackle-neq-postfeedback-contract-20261006/base/output_00001`–
  `output_00008` (step 1 was used; the others were unused duplicate raw dumps)
- `.grackle-neq-postfeedback-contract-20261006/restart/output_00001`–
  `output_00008` from the unintended fresh-start attempt (its log is retained)
- `.grackle-neq-postfeedback-contract-20261006/restart/output_00001` (corrected
  synthetic restart input) and `output_00002` (evaluated continuation)

The first delayed-control setup failed before producing a dump and was kept
for its log/input record. These removals are limited to the enumerated test
outputs; no ICs, executable, logs, or user simulation data were removed.

The following additional UV-switch dumps were evaluated with AMR leaf-cell
filtering and removed; effective namelists, binary/table hashes, logs, and the
scalar findings are retained in `.grackle-neq-uv-live-pair-20261006` and the
validation report:

- `.grackle-neq-uv-live-pair-20261006/uv_on/output_00001`
- `.grackle-neq-uv-live-pair-20261006/uv_off/output_00001`

After the 128^3 matched cosmological trajectories completed and the scalar
state checks were captured, their only final HDF5 dumps were also removed:

- `.grackle-neq-cosmo-extended-414735/real_ic/output_00001`
- `.grackle-neq-cosmo-extended-414735/uv_off/output_00001`
