# P4-C: AGN instantaneous-rate ledger

The RAMSES pilot enables `sinkprops=.true.`. Its `sink_*.dat` file is a
30-record Fortran-unformatted diagnostic, distinct from `sink_*.out` restart
state. It stores `dMBHoverdt`, `dMEdoverdt`, `dMsmbh`, and `eps_sink` with the
code-unit conversion factors in the same file.

```bash
python tools/p4_build_agn_rate_ledger.py \
  --diagnostic sink_000123.dat \
  --output data/agn_rate_000123.csv
```

The converter writes the instantaneous Bondi and Eddington rates, their
non-negative minimum, the active radiative efficiency, and
`L_bol = epsilon_r * mdot_inflow * c^2`. The calculation treats the selected
rate as inflowing mass rate, which is the convention used by the active
Bondi/Eddington limiter.

This is an AGN rate ledger, not a photon ledger. Conversion to
`q_group_N_s` requires a declared intrinsic AGN SED and photon-group energy
bounds, then feeds `P4_SOURCE_LEDGER.md`. It must not be replaced by a
conversion from sink mass or `dMsmbh` alone.

## Preferred coarse-step JSON source

The active `patch/lagRamses` writer appends `agn_coarse_state_v1.jsonl` before
feedback resets its coarse mass accumulators. It records the instantaneous
Bondi and Eddington rates, effective radiative efficiency, code-computed
bolometric luminosity, and sink position. Select exactly one matching coarse
step by the hydro snapshot expansion factor:

```bash
python tools/p4_build_agn_rate_ledger.py \
  --agn-coarse-json agn_coarse_state_v1.jsonl \
  --aexp 0.208497764676753 \
  --output data/agn_rate_000017.csv
```

This path preserves the active jet-mode effective efficiency. It rejects a
selection that matches zero or multiple coarse steps rather than guessing a
rate interval.
