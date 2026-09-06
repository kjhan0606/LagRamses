# Fable audit: F-P1 40--120 M☉ fate resolution

Date: 2026-09-02  
Scope: scientific/algorithmic validity and RT/feedback/population-ledger
wiring.  This was a read-only audit; no project files were edited by Fable.

## Verdict

**CONDITIONAL PASS** for the resolution strategy and contract metadata.
The scientific F-P1 gate remains **BLOCKED**, which is the correct state.

The rejection of a universal 40--120 M☉ direct-collapse bin is justified.
The proposed source-node/structure-based resolver is the correct shape, but
the source, table conversion, and runtime admission are not yet production
approved.

## Scientific validity

Explodability is non-monotonic in ZAMS mass and depends on metallicity,
rotation, wind history, pre-supernova structure, fallback, and the explosion
engine.  Supporting references include Heger et al. 2003, O'Connor & Ott
2011, Ugliano et al. 2012, Ertl et al. 2016, Sukhbold et al. 2016, Sukhbold,
Woosley & Heger 2018, Müller et al. 2016, and Burrows & Vartanyan 2021.

The solar-metallicity Sukhbold outcome pattern must not be transferred to the
low-metallicity JWST population without a metallicity axis.  At low Z, retained
helium cores can enter the PPISN regime near this interval (Woosley 2017;
Heger & Woosley 2002; Farmer et al. 2019; Woosley 2019; Renzo et al. 2020).
PPISN/PISN classification should therefore be evaluated inside the same fate
call, before core-collapse classification, even though its production gate is
kept separate as F-P3.

Failed supernovae are not necessarily ejecta-free: a failed collapse can eject
the hydrogen envelope with a weak transient.  The outcome enum must include
“direct collapse with envelope ejection”, rather than treating all failed
models as zero terminal ejecta.

For this project, multi-metallicity grids (Limongi & Chieffi 2018; Ugolini et
al. 2025; Roberti et al. 2024; Heger & Woosley 2010) are the production target.
Sukhbold is a solar-only engine-validation branch.  Patton & Sukhbold 2020 is
an algorithmic structure-based shape requiring source/engine calibration.
Fryer et al. 2012 is a sensitivity prescription, not a complete per-star
yield source.

## Algorithm and wiring

Mandatory resolver axes are source id/version/checksum, ZAMS mass, birth
metallicity plus its definition and solar abundance set, initial rotation as a
value or declared marginalization, engine/branch id, mass-cell assignment and
edge convention, lifetime source, and pair-instability criterion id.  Binary
state becomes mandatory when `binary_ssp` is selected.  Mass-loss and fallback
tags may be optional source-specific metadata.

The resolver must reject out-of-hull queries.  It must not use nearest-node
substitution, cross-source interpolation, metallicity/rotation extrapolation,
or a mass-only fallback.  A mass-only direct-collapse rule is explicitly
invalidated by the non-monotonic source outcomes.

The largest implementation gap is the existing runtime interpolation layer:
linear blending of a failed node with an exploding node creates an unphysical
fractional explosion and remnant.  The resolver needs a declared
piecewise-constant source-node/mass-cell mode, enforced in the runtime and
covered by an audit.

Wind return must remain cumulative and age-resolved when used as a time source;
terminal ejecta is added once after the terminal outcome is known.  Failed or
direct-collapse models may retain wind return and a remnant but no terminal
ejecta.  PPISN pulse histories and full disruption require distinct outcome
semantics.  A PISN has no remnant owner; PPISN and PISN cannot be represented by
one boolean channel.

The current fail-closed policy string is useful but not sufficient as a future
admission mechanism: it must be coupled to the approved map checksum, sidecar,
approval id, zero unresolved intervals, owner mappings, and source contract.
Legacy mode currently bypasses the fate check and must be banned for production.
The physics contract's PISN owner flag also needs to agree with the feedback
contract and both Fortran mirrors.

## Safe to apply now

* Correct the PISN owner semantics and keep PPISN/PISN pending F-P3.
* Add a fixed metallicity axis, mass-cell assignment rule, edge inclusivity,
  and pair-instability criterion metadata.
* Harden the fate-map audit against every policy-flag mutation and cross-check
  the feedback and physics owner tables.
* Extend the Sukhbold audit to W18, N20, and implosion result files; do not
  hand-enter outcome vectors.
* Create a resolver contract with zero physical nodes, explicit outcome enums,
  and source-hull rejection.
* Add a ledger unresolved-mass diagnostic and couple future gate opening to a
  sidecar checksum and approval id.

The Fable session could not open the staged Sukhbold tarballs, so the dossier's
W18/N20 numbers remain local evidence requiring tool-based reproduction.  No
canonical 40--120 M☉ row, remnant, energy, momentum, lifetime, decay value, or
source coverage is approved by this audit.

## Blocking conditions

1. An approved source with a metallicity axis covering the runtime domain.
2. Piecewise-constant node-cell mode implemented and audited.
3. A terminal owner for 40--120 M☉ defined consistently in all contracts and
   Fortran defaults.
4. Explicit remnant and envelope-ejection semantics per node.
5. Age grid resolving each terminal lifetime step within a declared tolerance.
6. Age-resolved wind history, or an approved quantified terminal-lumped model.
7. Decay horizon and inventory closure.
8. Terminal momentum and deposition contract.
9. F-P3 PPISN/PISN outcome and ownership decision.
10. Gate coupling to approval id, map checksum, sidecar, and a production ban
    on legacy mode.
11. MPA redistribution permission or a verified alternative license.
12. Tool-based reproduction of W18/N20 outcomes.

## Required next implementation bundle

1. Contract and metadata fixes, including a zero-node resolver contract.
2. Fate-map audit hardening and Sukhbold result-file parsing.
3. A Python reference resolver with hull rejection, non-monotonic synthetic
   tests, and edge-inclusivity tests.
4. Table-format and Fortran piecewise-constant path, native mirror parity,
   and JAX differential.
5. Ledger unresolved bucket and production-gate coupling.
6. Source conversion only under a recorded physics approval id.

