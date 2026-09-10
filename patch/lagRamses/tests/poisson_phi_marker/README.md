# Poisson checkpoint marker smoke test

Run `fresh.nml` for two coarse steps.  The initial output must omit
`POISSON_PHI_VALID`, while `output_00002` must contain it because the static
mesh has completed a Poisson solve without a subsequent topology change.

Run `restart.nml` in the same directory. The legacy true setting must emit
the deprecation warning, and even a valid marker must select the cold
predictor. No `Poisson warm start from restored phi` message may appear.
Removing the marker must also select the predictor. Checkpoint phi and its
validity marker are still written/read; their presence no longer enables
warm initialization. The same policy applies to binary, variable-rank, and
HDF5 readers.

The C0 checkpoint-69 regression reproduced a timestep collapse with saved
phi used in the fixed AMR boundary RHS. Tightening epsilon from 1e-4 to
1e-10 did not fix it. Rebuilding boundaries while retaining the saved active
guess restored the cold timestep to 8.32e-11 relative agreement. Until a
separately validated boundary-safe warm path exists, all restarts use cold
initialization, including old namelists requesting warm start.

`topology_change.nml` exercises a live level-7 mesh.  Any output written
between a topology change and that level's next solve must omit the marker.
