# Supporting-infrastructure long-term backlog

This project remains focused on RT, stellar/AGN feedback, dust, and their hydro
coupling.  The items below are not active physics gates unless a high-level
coupling test demonstrates that they block or corrupt those models.

- Normalize generic RAMSES parameter-error exits so all `clean_stop` failures
  return a nonzero scheduler-visible status.  P0.4 already blocks a missing
  stellar namelist; only the legacy global exit status remains.
- Add a generic MPI/OpenMP fault-injection framework.  P0.4 implements the
  required feedback-local parent-thread abort without redesigning unrelated
  error paths.
- Broaden generic HDF5, AMR, rank-layout, and low-level hydro validation only
  through their dedicated backlogs, not through feedback gate closure.
