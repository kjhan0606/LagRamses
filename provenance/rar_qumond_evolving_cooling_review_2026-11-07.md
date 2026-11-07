# Astra review: evolving native gravity with cooling/heating

Read-only Astra follow-up at1386db5 under the substantive coupling review
requirement. Q-GOAL approves direct integration of existing cooling/heating
into actual evolving gas toward the RAR Lambda-CDM/QUMOND comparison and later
BIMOND/FDM. Q-LEAN retains the existing full-dt level-local Lie split; no new
thermal solver or Strang scheduler is needed for this bounded integration.
This is not a Fable verdict or full feedback/enrichment approval.

Place original RAMSES cooling after each level's transport/old-density half
source, before final hierarchy restriction, new-density native solve and new
force half source. Cooling leaves rho/p unchanged and changes internal energy;
the final gravity kick preserves internal energy. They commute at fixed density
apart from roundoff. Keep old covered coarse uold until its Godunov commits;
do not upload fine cooling into old parent stencils prematurely.

Allow cooling only for qh_evolve; keep its rejection in frozen publication.
Require original/non-Grackle, isothermal and delayed cooling off, Jeans/eEOS
floors off, and existing PIC/star/sink/neq/RT/SGS exclusions. Declare abundance,
shielding and UV settings. Noncosmological init_time fixes aexp1, so HM tables
are fixed z0, not a separately selected aexp_ini epoch. Metal cooling uses the
advected rho_metal/rho/.02, not enrichment. Disabled T2_star must set its floor
exactly zero rather than evaluate zero times an undefined normalization.

Three existing diagnostic assumptions must change: bracket EVERY level cooling
call with step/time/level identity and actual owned leaf volumes; avoid unsafe
particle-list allocation/traversal in gas-only global inventory; and group by
physical call rather than constant aexp or hard-coded32768 cells. Old/new-half
energy joins must reconcile the intervening accepted cooling increment.
The selected coolfine1 does not call enforce_eeos_after_sink: after the energy
assignment it only updates an optional delayed tracer. Initial temperature
clipping may still distinguish the effective table RHS and actual energy
change, so record that separately and do not call it pure emissivity.

One evolving interface bundle suffices: active cooling and heating with declared
units/background/metallicity, matched cooling-off and same-endpoint half-dt
replay. Require active exchange on both levels, full per-call leaf coverage,
unchanged cooling-call mass/momentum/passives, accepted/clipping/whole-call
energy reconciliation, positive finite endpoints and existing moving-source
field checks. Keep default/frozen regressions; no disc survey rerun. Lie split
has not established second-order full physics evolution. No total gas-plus-
QUMOND energy conservation, source calibration or galaxy readiness is claimed.
