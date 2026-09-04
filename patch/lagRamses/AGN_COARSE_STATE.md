# AGN coarse-step state ledger

`AGN_feedback` calls `dump_agn_coarse_state` before it deposits feedback or
resets the coarse-step accretion accumulators.  Rank 1 appends one JSON object
per sink and coarse step to `agn_coarse_state_v1.jsonl` by default.

The `&physics_params` controls are:

```fortran
agn_coarse_dump = .true.
agn_coarse_dump_file = 'agn_coarse_state_v1.jsonl'
```

Each row is keyed by `(nstep_coarse, sink_id)` and records the expansion
factor, redshift, sink mass, gas angular momentum, BH spin vector and
magnitude, physical BH angular-momentum vector and norm in cgs, gas--spin
angle, Bondi/Eddington/actual inflow rates, raw/resolved/effective
radiative efficiencies, instantaneous bolometric luminosity, saved feedback
energy and the last Bondi gas context.  It also records the accumulated
coarse-step Bondi-supply, Eddington-limit, and accreted BH masses.  The
`feedback_mode` classification uses their coarse-step Eddington ratio, as does
`AGN_feedback`; `feedback_energy_deferred` reports whether energy remains in
the saved reservoir.  Rates and luminosity are reported in code and cgs or
solar units as appropriate.

The machine-readable boundary is explicitly marked with
`ledger_phase="pre_feedback_pre_reset"` and
`source_interval_kind="instantaneous_pre_reset_state"`.  The
`raw_radiative_efficiency` field is the sink-array value.  The
`radiative_efficiency` field is the shared helper's resolved base coefficient,
and the distinct `effective_radiative_efficiency` field is the value used in
`L_bol` and the SNRT photon budget.  `efficiency_status_name`,
`efficiency_mode`, and `efficiency_contract_ok` expose initialization or
configuration divergence; a spin-enabled zero `eps_sink` is therefore not
silently promoted to physical parity.  Non-promotable rows remain readable in
the append-only diagnostic ledger so one initialization/configuration issue
does not make later coarse records inaccessible, but the canonical state
reader excludes them from `AgnCoarseState` promotion and the P4 converter
refuses the requested artifact before opening its output.  An idle/MAD-
quenched sink may have effective efficiency zero and zero luminosity.  All rate
conversions use the declared Julian year of 365.25 days.

The disk fields evaluate the coherent accretion-episode prescription used by
`kjhan_growspin` at the dumped current state: episode mass, viscous time, warp
radius, and (when enabled) self-gravity radius.  They are derived diagnostics,
not a persistent sub-grid disk reservoir or a history of every internal spin
substep.  `disk_model_valid=false` and JSON `null` are emitted when zero spin,
zero accretion, or another degenerate state makes that model undefined.  The
spin--gas angle is likewise `null` when either angular-momentum vector has zero
norm.

Independent restarts can append a repeated key.  The canonical Python reader
and `audit_agn_coarse_ledger.py` collapse an identical semantic duplicate
(including JSON formatting/key-order changes) and treat a conflicting
duplicate as a provenance error.  Because the record has no run UUID or dump
counter, a rewind that produces a different payload also fails closed.  A
durable crash journal and cross-coarse-step deferred-energy re-emission policy
are separate follow-ups; this ledger is not an interval-integrated photon
ledger.  Open, write, flush, or close failures stop the simulation rather than
silently losing requested diagnostics.
