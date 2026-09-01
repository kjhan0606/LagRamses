# P4 pilot AGN source selection

`p4_select_pilot_agn_sources.py` selects rows from an audited instantaneous
AGN rate ledger using the fixed, non-wrapping P4 high-density cube in
`p4_high_density_manifest.json`.

The resulting candidate ledger is deliberately not a photon-source ledger.
It preserves source position, sink mass, instantaneous Bondi/Eddington-limited
inflow rate, radiative efficiency, and bolometric luminosity.  A scientifically
audited AGN SED, photon-group integration, and escape/obscuration prescription
must add `q_group_N_s` later.  `snrt_core.source_ledger` will reject this
candidate file until those photon-number columns exist.

The selection is valid only for a sink diagnostic and a hydro cube at the same
simulation epoch.  It does not infer luminosity from a sink mass or from a
restart checkpoint accumulator.
