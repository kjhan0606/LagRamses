# Generating an isolated CPU Makefile

`configure` is an additional Python 3 script. It never edits the existing
`bin/Makefile` or any Fortran/C/C++ source. It generates a new Makefile in a
separately chosen build tree, links source directories into that tree, and
keeps objects/modules and the new executable in its own `bin` directory.
Simulation inputs, checkpoints and outputs belong in a different run directory
chosen by the user; `configure` never starts a simulation.

```sh
./configure --build-dir "$HOME/.astro/build/lagRamses"
make -C "$HOME/.astro/build/lagRamses/bin" -B
# The executable is .../bin/ramses_configured3d.
```

The default/auto profile detects Intel MPI `mpiifx` or legacy `mpiifort`, matching
this fork's existing supported compiler. It compiles/links a temporary MPI +
OpenMP Fortran probe before writing the generated Makefile. It does not launch
MPI processes. The initial profile retains the fork's source/patch order, 3D
NVAR=18 CPU MPI/OpenMP build and uses its iterative Poisson (`USE_FFTW=0`), without
activating optional CUDA, HDF5, RT or CHIMES. Physics is not newly calibrated or
validated by generating a Makefile.

`--fc /path/to/mpiifx` selects a wrapper explicitly. An experimental GNU flag
profile is available only with `--compiler gnu`; a successful probe does not
establish that this fork's full solver compiles or that its science is validated
with GNU. The configuration JSON records this distinction. `--no-probe` is
intended for offline Makefile review and records that no compiler check occurred.

```sh
./configure --compiler intel --no-probe --build-dir /tmp/lagramses-review
# Or explicitly probe a GNU MPI toolchain; full solver validation is separate:
./configure --compiler gnu --build-dir /tmp/lagramses-gnu-build
```

A nonempty build directory is reused only when its configuration marker belongs
to this source tree; unrelated files/directories and conflicting source links
are refused. Rerunning configure updates only its generated files. Use `make -B`
to rebuild objects after compiler/profile changes. Builds are serialized to
respect the existing Fortran module order. `lagramses-configure.json` records
source identity, original/generated Makefile checksums, compiler version and
probe status. No solver or simulation is run by configuration alone.

Astrafek invokes this configure file after the user's explicit GitHub download
request, then builds using the generated Makefile. Astrafek ships only its own
executor/configuration generators, not lagRamses source or binaries.
