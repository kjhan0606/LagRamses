# P4 pilot status: output 00017

The RAMSES job reached wall time only after `output_00017` was written. The
job log explicitly reports `HDF5 output done.` The `levelmax = 17` entry in
`info_00017.txt` is an allowed refinement ceiling; the actual mesh hierarchy
ends at level 15 because it contains no level-16 or level-17 grids.

`p4_attach_pilot_sources.py` creates a narrowly scoped transport/chemistry
pilot by preserving the completed `output_00016` 32-cubed gas input and adding
the audited instantaneous `output_00017` AGN photon ledger. Its output
metadata records both scale factors and their difference. The artifact is for
S_N interface and numerical tests only, never for physical conclusions.

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
