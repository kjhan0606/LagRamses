# Dust IR input connection (not live IR completion)

The version-3 native dust sidecar now carries an independent IR quadrature:
node energies and integration weights in eV, absorption cross section per H
in cm2, and background temperature in K. It does not reuse the nine primary
photon groups or their mean energies. The background must be a thermal table
node, matching the native IR initializer. Arrays are bounded to 256 nodes,
positive and finite, with strictly increasing energies.

Version 1/2 parsing and reference opt-in behavior are retained. IR inputs with
an older version are rejected rather than silently ignored. Version 3 requires
an IR status consistent with the opacity/thermal approval class for runtime
admission. Candidate IR remains inspection-only. Reference data require the
existing explicit opt-in and are not scientific approvals.

The new synthetic v3 fixture has nine primary groups but two IR nodes. The
native contract test loads this file, initializes `snrt_dust_ir` with its
arrays, and advances a finite-capacity, source-free cooling cell. Positive
emitted photons, decreasing dust temperature and material+radiation+escape
energy closure within 1e-10 are required. The bolometric inspection entries
are calculated from the same synthetic Planck quadrature; they are not a
separate physical dust model.

Validation:

- `bash simulation/snrt/tests/run_snrt_native_dust_contract.sh`: Intel ifx and
  GNU pass, including both v2/v3 with opt-in absent, 0, 1, invalid, and overlong.
- The active HDF5 dust binding appends IR numerical inputs only for v3;
  v2 retains its existing attribute layout. This new branch is compile-tested,
  not yet exercised in a live v3 restart.
- Build affected active objects with
  `make SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1 snrt_dust_contract.o snrt_hdf5.o read_params.jaehyun.o`.

RAMSES explicitly rejects v3 live execution until its IR transport and
persistent IR state are connected. It must not accept the spectrum and run
absorption-only while silently dropping re-emission. The isolated operator
test is not evidence of live AMR/MPI IR transport. This startup restriction is
to be removed as part of that connection, not treated as a completed feature.

Remaining work: live material/IR transaction and cell mapping, IR restart and
AMR/MPI state handling, physical opacity/heat-capacity/SED approval. Production
status is not claimed. No RAMSES simulation namelist fields changed, so neither
`mkrun.py` nor `ramses_nml_generator.py` needs a schema update for this sidecar.
