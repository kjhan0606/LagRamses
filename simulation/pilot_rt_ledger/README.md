# RT source-ledger pilot

This is an I/O-contract pilot, not a scientific production run. It restarts
the unchanged CDM calculation at `output_00016` (`z=3.799`) and writes one
HDF5 output at the first subsequent coarse step.

The pilot changes only output controls:

1. `sinkprops=.true.` writes per-coarse-step `sink_*.dat` diagnostics with
   instantaneous Bondi and Eddington rates, retained accreted mass, and the
   active radiative efficiency.
2. `outformat='hdf5'` writes `/particles` and `/sinks`, including stellar
   age/metallicity fields and sink accretion state.
3. The next scheduled expansion factor is set to the restart expansion factor
   so the first step produces `output_00017.h5`; subsequent full outputs are
   suppressed during the short pilot allocation.

`prepare_pilot.sh` creates only local symlinks to the immutable restart and
the already-built executable. `run_pilot.sh` has no automatic resubmission.

Success requires `/particles`, `/sinks/dMBH_coarse`, `/sinks/dMEd_coarse`,
`/sinks/dMsmbh`, and `/sinks/eps_sink` in the output HDF5 file, as well as a
`sink_*.dat` instantaneous-rate diagnostic. SNRT still performs the declared
stellar and AGN SED-to-photon-group conversion.
