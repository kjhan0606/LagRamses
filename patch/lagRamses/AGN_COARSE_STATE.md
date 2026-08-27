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
angle, Bondi/Eddington/actual inflow rates, effective
radiative efficiency, instantaneous bolometric luminosity, saved feedback
energy and the last Bondi gas context.  It also records the accumulated
coarse-step Bondi-supply, Eddington-limit, and accreted BH masses.  The
`feedback_mode` classification uses their coarse-step Eddington ratio, as does
`AGN_feedback`; `feedback_energy_deferred` reports whether energy remains in
the saved reservoir.  Rates and luminosity are reported in code and cgs or
solar units as appropriate.

The disk fields evaluate the coherent accretion-episode prescription used by
`kjhan_growspin` at the dumped current state: episode mass, viscous time, warp
radius, and (when enabled) self-gravity radius.  They are derived diagnostics,
not a persistent sub-grid disk reservoir or a history of every internal spin
substep.  `disk_model_valid=false` and JSON `null` are emitted when zero spin,
zero accretion, or another degenerate state makes that model undefined.  The
spin--gas angle is likewise `null` when either angular-momentum vector has zero
norm.

Independent restarts can append a repeated key.  Downstream consumers should
deduplicate identical `(nstep_coarse, sink_id)` rows and treat conflicting
duplicates as a provenance error.  Open, write, flush, or close failures stop
the simulation rather than silently losing requested diagnostics.
