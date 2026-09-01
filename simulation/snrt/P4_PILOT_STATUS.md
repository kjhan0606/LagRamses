# P4 pilot status: output 00017

The RAMSES job reached wall time only after `output_00017` was written. The
job log explicitly reports `HDF5 output done.` The `levelmax = 17` entry in
`info_00017.txt` is an allowed refinement ceiling; the actual mesh hierarchy
ends at level 15 because it contains no level-16 or level-17 grids.

`p4_attach_pilot_sources.py` now rebinds an audited source ledger to a static
gas input without changing gas fields. It requires the photon metadata, checks
exact agreement with the configured P0 edges and CSV totals, and records hashes
of every input. The stage-4 `p4_coeval_static_rt_input_agn9.h5` rebinds the
coeval output-00017 gas cube to the nine-group AGN ledger; its gas/source scale
factors agree to `2.5e-16`. Scientific eligibility still follows the gas-field
contract rather than source coevality alone.

A coeval science input now requires an adaptive-leaf resampler for the complete
output 00017 HDF5 hierarchy through its populated level 15. The implementation
is complete and passes a synthetic AMR conservation test in
`tests/p4_hdf5_staging.py`: it streams all populated levels, uses only
`son_flag = 0` leaves, and requires unit coverage and a density mass-balance
check. The field map also records whether each quantity is a snapshot dataset,
an explicitly derived quantity, or a pilot constant.

The real output 00017 HDF5 payload remains external because it is
866,729,878,508 bytes; its `info_00017.txt`, hydro descriptor, and build
sidecars are now copied under `/gpfs/kjhan/LRD_JWST`, and both thermal atlases
are present locally with verified source SHA256 values. A real-value preflight
passes without loading the HDF5 payload: levels 1--15, six mapped datasets,
HDF5 `gamma=1.6666667`, and ten ledger sources are present.

The field map now decodes the raw conservative variables: momentum is divided
by density for velocity, total energy is reduced to thermal pressure, and
`uold_6` is converted to `Z/Z_sun` using the producing 0.02 solar mass
fraction. Dust, non-equilibrium ionization, and H₂ remain explicit
non-production placeholders, so the production-contract gate remains closed.
