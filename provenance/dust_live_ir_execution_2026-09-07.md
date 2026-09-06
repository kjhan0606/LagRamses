# Native live dust IR: fixed-mesh execution and restart

This bundle connects the v3 sidecar to actual RAMSES material and radiation
updates. It supersedes the startup prohibition in
`dust_ir_input_connection_2026-09-07.md`, but does not claim general production
readiness or approved physical dust data.

## Implemented path

`snrt_ramses_driver` completes its primary RT/chemistry iteration, then stages
IR transport, absorption and finite-capacity emission. IR starts from the
**pre-primary** dust energy and receives the accepted primary absorption once.
It must not start from the receiver's already-heated trial energy. IR opacity
is temperature-independent in this contract, so this stage does not change
the converged primary opacity/chemistry solve.

`snrt_dust_live` retains FP64 IR energy (erg/cm3 per normalized direction) on
the primary cell-slot map, independently of primary photon groups. Trials use
local arrays; persistent IR and RAMSES dust energy are written only after the
primary transaction commits. IR uses CFL substeps when required, with a fixed
primary heating rate across those substeps. Dust mass and gas energy are not
modified by this IR receiver.

The angular handoff converts the primary weights (sum 4*pi) to weights summing
one. Primary photon entries remain directional-bin counts; their conservation
ledger uses a sum, not another angular weighting. IR entries use normalized
angular integration. An initial live test rejected the unnormalized weights;
that error was fixed at the handoff, not by weakening the IR operator check.

The HDF5 record uses format 2 for v3 dust: original 724 primary/chemistry
entries, then `nir * 80` IR entries in Fortran `(group,direction)` order. The
two-node fixture therefore has width 884. Record ordering follows the same
RAMSES grid/child ordering as before. Actual IR quadrature values are already
bound in the dust attribute, and `ir_energy_units` describes the normalization.
Version-2 dust retains format 1 and width 724. IR payloads are checked for
negative/nonfinite values before primary restore; radiation with absent primary
cell state is rejected.

Current startup support is explicitly one rank, noncosmological fixed mesh
(`levelmin=levelmax`, `nremap=0`), HDF5 build/output and HDF5 restart input.
The adapter admits only local or classified physical faces, never an MPI,
coarse-fine or unmapped face disguised as vacuum. The tested domain is periodic.
Reference inputs print NONPRODUCTION. MPI, AMR migration and cosmology remain
implementation work, not silently approximated features.

## Real execution evidence

Retained root: `/gpfs/kjhan/LRD_JWST/.dust-live-ir.4XJVco`.
Final executable: `ramses3d`, SHA256
`8f1a62d5f73ac1dad1b56d42864550b922c64ddba0f081c42ba3ee2801f1b75e`.
Build: `make -s -C bin -j4 SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1
EXEC=/gpfs/kjhan/LRD_JWST/.dust-live-ir.4XJVco/ramses`.
Objects were already in this same live profile; the default executable was
not overwritten. Rebuild all objects when changing profiles.

The tracked fresh-run fixture is `simulation/snrt/config/dust_live_ir_smoke.nml`.
It uses 512 periodic cells, uniform stationary gas, synthetic dust mass density
`1e-26 g/cm3`, capacity density `1e-24 erg/cm3/K`, initial dust temperature 20 K,
and the v3 reference sidecar's 10 K background. No stars, sinks/AGN, gravity or
gas cooling are active. `run_case.sh` in the retained root records the exact
environment: canonical reference primary group and secondary table contracts,
v3 dust, reference opt-in, SNRT level 3, GPU 1, one OpenMP thread.

All effective namelists retain `noutput=1`, `aout=2`, `tout=1e30`, with
`fbackup=1000000`; scheduled times are not reached. Two-step runs use
`foutput=2`; four-step endpoints use `foutput=4`. Each successful case writes
one new 4.2 MiB checkpoint. Available storage before launch was 173 TiB.

- `cooling`: 2 steps; nonzero IR, material cools from 20 K to the background.
- `continuous_fixed`: fresh 4 steps, final `output_00001`.
- `resumed_fixed`: restart from `cooling/output_00001`, advance steps 3 and 4,
  final `output_00002`. All 452608 primary/IR record entries and all 512 dust
  mass, dust energy and gas total-energy entries are bitwise equal to the
  continuous run. Integrated material+IR energy drift is `1.5681e-10` relative
  to the initial material energy. Periodic escape is zero.
- `pattern`: a **copy** of the cooling checkpoint seeds group-1 primary bins
  with FP32 `1e-11 * (1 + 0.25*x_grid + 0.01*child_index)` and scales initial IR
  by the same spatial factor (1.03125 to 1.28875). Two further steps test
  nonuniform transport, primary absorption and re-emission. Summed lost primary
  energy is `3.110559397334588e-20 erg/cm3`; gained IR+material energy is
  `3.110559378320851e-20`, relative discrepancy `6.11265e-9` (<2e-6 FP32-primary
  allowance). Gas energy and dust mass remain bitwise unchanged. Radiation
  remains finite/nonnegative and IR varies spatially.
- `final_pattern`: final executable repeats the pattern case and produces
  bitwise-identical primary/IR and material arrays, including the new units
  attribute. This is the final-binary execution evidence.

The first four-step trial exposed roundoff when converting a background
`C*T` to/from code units: 10 K became 9.999999999999998 K and was rejected.
The native range checks now allow 64 machine epsilons at the table bounds;
any energy correction remains counted in the closure check. A native regression
accepts a one-ULP-below-background state, checks its energy change <1e-13, and
still rejects 9.99 K without material mutation. The first restart trial also
omitted `informat='hdf5'`; the fixture and explicit IR startup validation now
require it. Failed trials are retained, not reported as passes.

`JAX_PLATFORMS=cpu simulation/snrt/.venv/bin/python
simulation/snrt/tests/dust_ir_transport.py --native` passes with GNU and Intel,
including the new background-roundoff regression and the existing weak/stiff
emission and rollback tests. Python inspects HDF5 outputs and exercises the
existing reference test; the live evolution is the native RAMSES executable.

No new RAMSES namelist fields were introduced: generator schemas are unchanged.
Remaining production blockers include approved physical opacity/heat-capacity/
SED inputs, MPI/AMR IR exchange and persistence, cosmological radiation evolution,
and the previously recorded high-level source/migration limitations. These are
not waived by the fixed-mesh reference result.
