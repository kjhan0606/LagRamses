# F-P1 40--120 M☉ literature bundle: application record

Date: 2026-09-02  
Project: `/gpfs/kjhan/LRD_JWST`

## Decision

Fable returned **CONDITIONAL PASS** for the literature-based strategy and
confirmed that the scientific F-P1 gate must remain blocked.  We therefore
applied only the contract and audit subset that cannot manufacture physical
coverage:

* added the `per_source_node_fate_lookup_v1` resolver contract with zero
  physical nodes;
* required piecewise-constant source-node mass-cell assignment, explicit edge
  semantics, source-hull rejection, and no nearest-node/cross-source/
  metallicity/rotation interpolation;
* added outcome enums for failed collapse with envelope ejection, PPISN, and
  complete-disruption PISN with no remnant owner;
* aligned PISN ownership in the feedback and physics contracts and both
  Fortran configuration mirrors;
* hardened the F-P1 map audit to reject unsafe policy mutations and owner/range
  disagreement with the G2 physics contract;
* added a runtime diagnostic-only Kroupa unresolved-mass bucket to the
  population ledger.  On the current 0.08--120 M☉ support it reports
  approximately 5.04% for 0.8--1 M☉ and 6.75% for 40--120 M☉; the bucket is
  excluded from returned/living/remnant closure and is never deposited or
  classified as feedback;
* added a production runtime selection for piecewise-constant source-cell
  mass assignment, with explicit half-open edge semantics;
* added a four-artifact SHA256 admission sidecar and mutation-tested audit
  coupling its map intervals and approval id to the source and physics
  contracts;
* extended the Sukhbold candidate audit to parse W18/N20 high-mass result
  files and all 105 `implosions_W18` wind-only tables.

## Verification evidence

The Sukhbold tool-based reproduction confirms the non-monotonic high-mass
pattern: W18 has positive explosion energy at 60 and 120 M☉, while N20 has
positive energy at 60, 80, 100, and 120 M☉ among the audited 40--120 M☉
nodes.  These are candidate-source evidence only; no canonical rows were
emitted.

Passed checks:

* `FP1_POPULATION_FATE_TEST_OK`
* `FP1_FATE_RESOLVER_TEST_OK`
* `G2_SUKHBOLD2016_CANDIDATE_TEST_OK`
* `G2_POPULATION_LEDGER_RUN_OK`
* `STELLAR_POPULATION_CONTRACT_UNIT_OK`
* `stellar feedback policy: PASS`
* `stellar yield audit contract: PASS`
* `stellar yield fail-closed policy: PASS`
* `stellar residual deposition: PASS`
* `FP1_IMF_JAX_DIFFERENTIAL_OK rows=8` (CPU fallback; no CUDA jaxlib)
* JSON validation and `git diff --check`

The source-parity checker remains blocked by the pre-existing
`production_linked_build_evidence` requirement; this is not treated as a
scientific approval.  A full production build and source conversion remain
pending.

## Still prohibited

No 40--120 M☉ canonical yield row, remnant mass, energy, momentum, lifetime,
decay projection, PPISN/PISN event history, or cross-metallicity coverage was
invented or promoted.  The existing linear yield interpolation path is not
allowed to resolve this seam; the runtime now selects the piecewise
source-node mass-cell contract, but physical source-node records are still
absent.

## Next bundle

Select a licensed metallicity/rotation source package, stage its immutable
source-node records and provenance, and independently reproduce the
age/wind/terminal/remnant/energy/momentum semantics before any physical node
is converted into the canonical table. Keep F-P1 blocked until that package
and its approval id replace the current zero-node sidecar.
