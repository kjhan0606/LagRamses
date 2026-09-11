# Bundle 3: explicitly mixed lower-mass comparison

Operator approved `근사모형, 명시적 비교로 진행`. This is an opt-in
comparison, not a full SSP or same-evolution yield/spectrum claim. No new
routine review gate: Main includes this design in its aggregate 3/4/Fe5 review.

Model `parsec_mixed_lowmass_truncated_v1`. New converter consumes the existing
upper PARSEC package, without editing Pauli's converter. Eleven target Z
branches are retained as comparison coordinates, not newly computed models.
Terminal-source coordinates and each projection are recorded separately.

* 2--7: nearest actual KL16/Fishlock terminal source in M and Z, scaled by
  target/source initial mass. Source lifetimes are retained (Fishlock's
  existing declared Raiteri fit). No implied PARSEC terminal yields.
* 7<M<9: scaled 7-Msun envelope proxy with that source lifetime. Its retained
  core is explicitly kind 3 (unclassified proxy), never a strict CO-Ia donor.
* 9--12: nearest Sukhbold solar stable-yield node, fractional mass scaling,
  existing residual-remnant convention and Raiteri lifetime. Solar material
  is explicitly a proxy at other target metallicities.
* 13: actual Sukhbold terminal source; time-dependent PARSEC 12-Msun radiation
  scaled by 13/12, explicitly not an original 13-Msun track.
* >=14: existing upper physical payloads retained (gross ejecta, remnant,
  kinetic energy, radiation and timing). Net-yield diagnostics are explicitly
  unavailable/zero throughout v5: no shared birth-composition convention is
  asserted for the mixed terminal sources. Gross ejecta are not zeroed.

PARSEC lower radiation uses actual Q5 constraints and existing within-band
Planck reconstruction, actual L/Teff only for gaps within the available track.
Stop at min(track end, selected terminal age); emit exactly zero thereafter.
Missing late luminosity is an omitted source, not gas heat, stored radiation
or recovered nuclear-fuel closure. No time-stretched or constant late plateau.

Material total lifetime budgets come from the selected terminal package.
Candidate pre-terminal PARSEC wind histories use nonnegative loss and surface
composition. Reduce the entire candidate history by one recorded scalar if
needed to fit each gross element AND the untracked reservoir. The terminal
release is the remaining total ejecta. Thus wind+terminal+remnant equals
initial mass, with no duplicated pre-AGB/stellar winds. Explicit wind speeds
set kinetic energy from the **scalar-reduced** cumulative wind mass:
`Ewind(t) = 0.5 * Mwind_reduced(t) * vwind**2` in cgs. Never retain the
unreduced wind kinetic energy. The residual AGB envelope uses the same
explicit speed. Missing luminosity never supplies mechanical work.

## Native contract / ownership

History namelist v5 extends `stellar_high_mass_history` with node arrays
`terminal_channel`, `remnant_kind`, `radiation_stop_age_yr`.
Terminal channel 2=AGB/envelope, 3=CCSN/failed/direct, 5=pair.
Remnant kind 0=CO,1=CO(Ne),2=ONe,3=unclassified retained-core proxy,
4=compact massive-star remnant. Fate 0 only for AGB; massive fate 1..5 retains
v4 meanings. Channel 2 owns AGB remnants; channel 3 owns massive remnants
(including PPISN remnants), channel 5 has ejecta/energy only.

All enabled generic source windows (1,2,3,5) are 2--600, with non-owner nodes
returning zero. Shared nearest-node cells and fractional budgets use full
0.08--600 IMF normalization. No source below2, but that mass stays living /
unresolved in the population ledger. Preserve the existing40-Msun seam for
upper-package compatibility in BOTH radiation and material cell selection.

Version5 is mutually exclusive with **every SNIa model** (`enable_snia=false`,
single-star population, zero binary fraction required by native admission).
Bundles3 and4 cannot be combined in one population. The strict CO inventory
exclusion smoke is a table-level classification check only, not a live Ia
pathway. In particular the 7<M<9 retained-core proxy never supplies Ia WDs.

SED wrapper version5: population_binding='common_imf_mixed_evolution',
radiation_population='PARSEC_mixed_lowmass_comparison',
spectral_tail_policy='zero_after_track_or_terminal_v1'; other Q/E metadata and
reference-control semantics unchanged. Native node payload magic
`SNRT_POPULATION_QE_INTERVAL_V2`, followed by nn,nr, model_id, then per-node
`M Z terminal_age_Myr radiation_stop_Myr terminal_channel remnant_kind fate nknot`
and existing21 positive interval moments. Last radiation knot equals stop,
NOT necessarily terminal age. Consumers retain their existing public calls.

Loader API is `parsec_sed_load(path,ierr,expected_version)` with optional
integer last argument. `expected_version=4` requires V1 nodes;
`expected_version=5` requires V2 nodes. Legacy two-argument native calls remain
valid. In-memory material names are `hm_terminal_channel`, `hm_remnant_kind`,
`hm_radiation_stop` (Gyr); history ages on disk remain years and SED ages Myr.
Node scalar row ranges are `hm_first(channel,node)`, `hm_count(channel,node)`.
`mixed_source_weights(table,population,n_mass_bins,weights,ierr)` is shared
by the material integrator and spectral binder, not a second population model.

Canonical v5 material rows are ordered Z, M, channel, age. Every node has
channels1,2,3,5, each with a zero origin and a terminal-age endpoint; empty
channels have zero payloads. The channel column is an integer token. There
is no channel4 source. Admission and per-release condensation use the indexed
ranges rather than global quadratic searches.

979 nodes for eleven branches: v5 cap1024, v4 cap512 unchanged. Row ranges
are indexed per node/channel for native evaluation and condensation; no
new per-particle state or MPI layout. Bind all new node fields in radiation
identity and material restart identity (Main-owned runtime serialization).

Main owns: v5 dispatch/reporting in snrt_stellar_source; v5 admission in
stellar_enrichment_driver pair guards and resolved bucket; runtime identity;
frontends, full build and integrated run. This worker owns converter, native
table/audit/interpolation/SSP/PARSEC consumer, physical fixture and evidence.

## Focused verification

Actual-source native calls, not disconnected Python-only validation:
early/late radiation and cutoff; wind/AGB/SNII/pair increments; AGB remnant
ownership and strict CO inventory exclusion; full-IMF mass accounting;
same-age Z mixture and interval telescoping; malformed ownership/duplicate/
negative-budget rejection. Regress v4 with its existing actual package.
Private native builds only here; no MPI, live launch or commits.
