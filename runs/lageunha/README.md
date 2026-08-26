# Lageunha SIDM run metadata

`snapshot-2026-08-24` is a lightweight copy of the configuration and status
metadata below `/gpfs/kjhan/Hydro/Sidm` on Lageunha. The snapshot contains 301
files and occupies about 2.7 MB.

Included files are namelists, shell and Python scripts, JSON diagnostics,
launcher return codes and process identifiers, job-control files, and RAMSES
`info_*.txt` files. Compiled executables, initial conditions, particle data,
mesh data, full logs, and other large outputs are intentionally excluded.

The snapshot is evidence of the 2026-08-24 state. Never use the snapshot as a
restart directory and never assume that editing a copied namelist changes a
live run. The canonical data remain on Lageunha.
