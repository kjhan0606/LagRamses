# Uniform CPU BIMOND verification

The two-sector checks exercise the `alpha=-1, beta=1` nonrelativistic
equations using tagged ordinary and twin particles. The implementation
currently requires a uniform periodic CPU FFTW grid with CIC assignment.
HDF5 restart requires the same MPI size and ordering. AMR, GPU, and coupled
hydro are not yet admitted.

Build the lagRamses patch with `NDIM=3`, `USE_FFTW=1`, `HDF5=1`,
`USE_CUDA=0`, and `SNRT=0`. Python analysis requires NumPy and h5py.
Use a compute-node allocation for the calculations. The RAR project routes
CPU calculations through Grammar Slurm.

Inside an allocation, create a new run directory and generate inputs:

```sh
mkdir /absolute/new/run
python tests/bimond/verify_bimond_two_sector.py prepare /absolute/new/run \
  --template tests/bimond/source_template.nml
```

The generator creates `equal32`, `pair32`, `pair64`, and `shift32`.
Run the compiled executable with `effective.nml` from each case directory,
using one MPI rank and one OpenMP thread. Each case requests one static
output. Inspect `output_policy.json` and available storage before launching.
The expected four snapshots and ASCII ICs fit within the reserved 1 GiB.

After all runs finish:

```sh
python tests/bimond/verify_bimond_two_sector.py assess /absolute/new/run
```

The verifier independently deposits each particle species and solves both
discrete field equations. Static RAMSES stores interpolated accelerations
in particle velocities, enabling a check of force routing. Identical sector
densities must cancel the differential and phantom sources exactly.

`verify_bimond_evolution.py` assesses four-step `equal16` and `pair16`
outputs. The shared `prepare` function accepts explicit case tuples for
the evolution fixture. Set `static=.false.`, `nstepmax=4`, and
`courant_factor=0.05` in **MOND_PARAMS**. HYDRO_PARAMS is not read when
hydro is disabled. The verifier reconstructs both half-kicks and drift from
saved timesteps, particle identities, and sector fields. A second argument
enables the supplemental `shift64` and `mpi2_pair32` comparisons against
the original static run.

Measured results are retained in
[`provenance/bimond_cpu_2026-11-07`](../../provenance/bimond_cpu_2026-11-07).
The source revisions, cluster jobs, limitations, and remaining acceptance
conditions are recorded in the
[implementation document](../../provenance/bimond_quasistatic_implementation_design_2026-11-06.md).

## Restart check

`verify_restart.py prepare RUN --template tests/bimond/source_template.nml`
generates a four-step displaced-particle run in `RUN/full`. After running
the executable in that directory, `verify_restart.py resume RUN` copies
the second checkpoint into a new `RUN/resume` directory and prepares
three invalid restart cases. Run the continuation in `RUN/resume`, then
use `verify_restart.py assess RUN` to compare final particle states and
both sector fields by identity and physical cell coordinates.

The `wrong_a0` and `missing_identity` cases must report a BIMOND model
identity error. The `wrapped_tag` case must reject particle type 286
before conversion to the one-byte representation. Each invalid case
must stop before writing a new snapshot. The source checkpoints are
copied, never modified in place.
