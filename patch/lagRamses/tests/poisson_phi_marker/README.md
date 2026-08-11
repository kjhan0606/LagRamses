# Poisson checkpoint marker smoke test

Run `fresh.nml` for two coarse steps.  The initial output must omit
`POISSON_PHI_VALID`, while `output_00002` must contain it because the static
mesh has completed a Poisson solve without a subsequent topology change.

Run `restart.nml` in the same directory.  The log must report that the marker
was accepted and that restored phi was used as the first multigrid guess.
Removing the marker before the restart must instead select the predictor.

`topology_change.nml` exercises a live level-7 mesh.  Any output written
between a topology change and that level's next solve must omit the marker.
