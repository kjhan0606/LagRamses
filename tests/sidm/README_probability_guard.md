# SIDM probability guard

For three-dimensional PIC runs with a constant, time-independent cross section
and elastic, non-dissipative SIDM, enable the optional collision-only guard in
`&SIDM_PARAMS`:

```fortran
sidm_probability_guard=.true.
sidm_type='constant'
sidm_a_type='none'
sidm_inelastic=.false.
sidm_fdiss=0.0d0
sidm_courant=0.1d0
```

The guard defaults to false for compatibility with other SIDM modes. It checks
that the cross section is finite and positive, `0 < sidm_courant <= 1`, and
`sidm_npart_min >= 2`; unsupported modes are rejected. Before each collision
sweep it computes an all-pair rate bound from the owned leaf-cell particle
masses and velocity ranges, uses one common MPI substep, and rechecks the
bound after velocities change. The hard probability limit is the smaller of
0.1 and `sidm_courant`. If necessary, the collision interval is split, with
a maximum of 4096 substeps. The gravity/drift step is not repeated. Visit
counts in the estimator accumulate across substeps; they are not a unique
particle census. The next-step constraint uses the maximum physical rate,
not a substep probability interpreted as a whole-step probability.

The guard aborts on invalid inputs or a sampled probability above the guarded
limit. With the guard off, invalid sampled probabilities outside [0,1] still
abort, but the old next-step predictor does not guarantee a 0.1 ceiling for
the current sweep. Enabling the guard is a numerical safeguard, not a
validation of the collision operator, thermal transport, or any physical
simulation setup. Check convergence and conservation for the intended model.

Run the local deterministic regression with
`bash tests/sidm/run_sidm_subcycle_test.sh`; MPI and full RAMSES runs require
separate integration tests.
