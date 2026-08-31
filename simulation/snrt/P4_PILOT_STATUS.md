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
output 00017 HDF5 hierarchy through its populated level 15.

The adaptive-leaf resampler is now complete. It streamed all populated levels,
used only `son_flag = 0` leaves, and covered each 32-cubed analysis cell with
volume weight one to floating-point precision. The resulting coeval static RT
input uses output 00017 gas and the output 00017 instantaneous sink photon
ledger at `aexp = 0.208497764676753`.
