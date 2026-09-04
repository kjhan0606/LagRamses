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
non-negative minimum, both the raw and effective radiative efficiencies, and
`L_bol = epsilon_eff * mdot_inflow * c^2`. The calculation treats the selected
rate as inflowing mass rate, which is the convention used by the active
Bondi/Eddington limiter.

The CSV columns are deliberately named `raw_radiative_efficiency` and
`effective_radiative_efficiency`; the latter is the luminosity convention.
For the coarse JSON source, the intermediate `radiative_efficiency` column is
the resolved base coefficient and the raw sink value is retained separately.
The coarse JSON path also emits `efficiency_status` and
`efficiency_contract_ok`.  Any row with `efficiency_contract_ok=false` is
rejected before the CSV is created, so configuration or initialization
divergence cannot silently enter a physical or publication artifact.  The
canonical reader keeps such rows readable for ledger diagnostics but does not
promote them into `AgnCoarseState`; the converter rejects a requested coarse
step containing one.  The
legacy `sinkprops` path has no mode-resolved contract; its corresponding CSV
fields are blank and `efficiency_contract_source` says
`legacy_sinkprops_mode_unresolved`.
Native `sinkprops` input has no mode-resolved effective field, so its
review-only rate ledger records raw == effective and does not approve an AGN
SED or obscuration model.  All rate conversions use a 365.25-day Julian year.

This is an AGN rate ledger, not a photon ledger. Conversion to
`q_group_N_s` requires a declared intrinsic AGN SED and photon-group energy
bounds, then feeds `P4_SOURCE_LEDGER.md`. It must not be replaced by a
conversion from sink mass or `dMsmbh` alone.

The current parameterized pilot converter
`tools/p4_build_agn_photon_ledger.py` defaults to the pinned P0 nine-group
table in `config/p0_photon_group_edges_ev.txt`. Its Sazonov-style pilot SED is
defined from 10 eV upward; the two lower P0 groups therefore carry explicit
zero-photon/zero-gas-opacity closure entries rather than an invented low-energy
extrapolation. The configured `[5.6,11.2] eV` group is explicitly marked
partially supported and integrates only `[10,11.2] eV`; all higher groups are
fully supported. The retained five-group outputs can be reproduced only with
`--legacy-five-groups`; they must not be merged with a nine-group ledger.
The canonical nine-group regeneration, exact edge gate, threshold-boundary
repair, and short production-runner integration are recorded in
[`AGN_NINE_GROUP_VALIDATION.md`](AGN_NINE_GROUP_VALIDATION.md).

## Preferred coarse-step JSON source

The active `patch/lagRamses` writer appends `agn_coarse_state_v1.jsonl` before
feedback resets its coarse mass accumulators. It records the instantaneous
Bondi and Eddington rates, raw/effective radiative efficiencies,
code-computed bolometric luminosity, and sink position. Its pre-reset and
instantaneous markers are required by the canonical reader. Select exactly one matching coarse
step by the hydro snapshot expansion factor:

```bash
python tools/p4_build_agn_rate_ledger.py \
  --agn-coarse-json agn_coarse_state_v1.jsonl \
  --aexp 0.208497764676753 \
  --output data/agn_rate_000017.csv
```

This path preserves the active jet-mode effective efficiency. It rejects a
selection that matches zero or multiple coarse steps rather than guessing a
rate interval. Identical restart duplicates are collapsed by key
`(nstep_coarse, sink_id)`; conflicting same-key payloads fail closed.

The bounded source/ordering audit can be run with a ledger from the same
production run:

```bash
python tools/audit_agn_coarse_ledger.py \
  --input agn_coarse_state_v1.jsonl \
  --helper ../patch/lagRamses/snrt_agn_efficiency.f90 \
  --output data/agn_coarse_ledger_audit.json
```
