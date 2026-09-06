# Claude Opus 5 end-of-bundle audit — F-P2.1 source-closure verification — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Auditor: Claude Opus 5 (`claude-opus-5`), read-only restricted session
`c3b5158b-a278-45f1-b5e5-20a3c2fdbf19`. The auditor reported no repository
modification, job/build/simulation launch, or shell hash recomputation. The
local `/gpfs` tests and recorded hashes were checked independently in the
implementation evidence.

## Overall verdict

**CONDITIONAL PASS.** The three F-P2.1 repairs are genuine: explicit and
parameterized-pilot metadata use disjoint paths; the source-weighted Draine
closure has an independent offset-grid quadrature and closed-form test; and
source-SED, group-edge, Draine-table, and builder hashes are recorded,
re-hashed, and enforced fail-closed by both runners before output creation.
The null-identity v1 reference control remains usable and is rejected against
a bound source ledger.

There is no blocker. The condition is that the provenance chain and canonical
coverage need one further closure-integrity bundle before this path is treated
as a fully closed engineering gate or promoted.

## Findings

### MAJOR-1 — builder provenance omits the integration module

`simulation/snrt/tools/build_draine_dust_opacity.py:368-371` records only its
own `builder.sha256`, and `simulation/snrt/snrt_core/dust.py:164-176`
re-hashes only that file. The actual source integration also depends on
`simulation/snrt/snrt_core/sed.py` (`read_lbol_photon_sed` and
`integrate_photon_sed_groups`), which is not in the recorded/enforced closure
manifest. The AGN ledger similarly hashes its wrapper but not the complete
SED/primordial dependency set.

Disposition: does not block F-P2.1's narrower “builder code hash” repair, but
blocks a claim that the complete closure algorithm is provenance-bound. Put
in the next closure-integrity bundle.

### MAJOR-2 — v2 sidecar payload is not integrity-bound

`simulation/snrt/snrt_core/dust.py:72-235` re-hashes four input files and
structurally validates the arrays, but does not verify that
`absorption_cross_section_per_h_cm2` and
`absorption_weighted_energy_ev` are reproducible from those authenticated
inputs, and the sidecar has no self-hash. The current negatives mutate hash
strings, not the arrays (`simulation/snrt/tests/source_sed_dust_closure.py`),
so arbitrary opacity values with valid provenance strings could reach the
kernel.

Disposition: outside the explicitly scoped R3 hash-enforcement claim, but it
blocks treating v2 as proof that the transported opacity came from the
declared table. Add a payload self-hash or loader-side recompute gate and an
array-tampering negative test next.

### MAJOR-3 — dust status is unconstrained free text

`simulation/snrt/snrt_core/dust.py:221-223` accepts any non-empty `status`,
and `p4_run_transport_pilot.py` / `p5_run_thermochemical_pilot.py` stamp it
into HDF5 attributes. A sidecar could therefore label itself
`production_approved` without a vocabulary gate.

Disposition: carried over from the prior audit and outside F-P2.1's repair
scope. Replace it with a fixed enum before any promotion claim.

### MINOR findings

- The validated dust metadata/table/builder hashes are not all propagated to
  P4/P5 output attributes; add `dust_opacity_metadata_sha256`,
  `dust_source_table_sha256`, and `dust_builder_sha256`.
- The canonical nine-group validator is structurally pilot-only: it hard-codes
  the parameterized 10-eV status vocabulary and intervals, has no dust
  criterion, and validates a zero-dust control. The repaired explicit path is
  covered only in the temporary-directory test.
- Escaped photon totals are folded into the CSV `q_group_*` and
  `total_photon_rate_s` names rather than carrying an explicit escaped name;
  the sibling per-`L_bol` fields are correctly separated.
- Explicit support intervals are derived from configured edges, not the local
  non-zero support of the SED. This is truthful for the enforced full-edge
  table coverage, but an empty group can still advertise full table support.
- The AGN builder writes its CSV before late `aexp` validation, so a failure
  can leave an orphan CSV without metadata.
- The canonical validator records `git_head` without a clean-tree attestation;
  the current artifact's file hashes are internally consistent but the commit
  label is the pre-change HEAD.
- Group-edge parsing is duplicated in three modules and the ledger/sidecar
  duplicated energy-fraction fields are not cross-checked.

These are not blockers for F-P2.1 as scoped, but the first four should be
included in the next closure-integrity bundle where practical.

## Audit answers

1. The explicit AGN path is now truthful and distinct from the Sazonov pilot;
   the audit confirmed the absence of pilot citation/normalization/support
   claims and the intrinsic/escaped bookkeeping test.
2. The Draine v2 mathematics is correct and the independent numerical test is
   genuinely separate from the production integration routine. Its coverage
   is narrower than a full dependency audit because the chosen power-law case
   is especially interpolation-friendly.
3. The four recorded provenance legs are re-hashed and fail closed in P4/P5;
   v1 null-identity reference controls remain intentionally available. The
   complete closure-code dependency set and sidecar payload are not yet
   bound.
4. Group edges, source support, units, empty groups, and identity are
   consistent across the repaired single-source path. Remaining issues are
   vocabulary, dependency breadth, and canonical explicit-path coverage.
5. The canonical pilot artifact closes its reference-control gate at 27/27,
   but that validator does not yet exercise the repaired explicit-SED path.
6. Mixed STAR+AGN admission, physical SED/escape/obscuration selection, dust
   scattering/IR/grain physics, live RHD, production convergence, and
   publication claims remain correctly deferred.

## Genuinely closed

Truthful explicit/pilot metadata separation; intrinsic versus escaped moments;
mathematically correct source-photon-weighted Draine v2 closure independently
checked on offset grids and in closed form; four recorded and loader-rehashed
provenance legs with fail-closed P4/P5 wiring; the preserved null-identity v1
reference control; and canonical pilot nine-group artifacts with an internally
consistent 27/27 validation graph.

## Remains deferred or conditional

Complete closure-code dependency/payload integrity, fixed status vocabulary,
full output provenance propagation, and canonical explicit-path coverage are
the next engineering closure items. Mixed STAR+AGN aggregate dust admission,
physical stellar/AGN SED and obscuration approval, the `[40,120] M_sun` yield
seam, unimplemented dust physics, live RT–RAMSES coupling, production
convergence, and publication claims remain later science/production gates.

## Recommended next coherent bundle

**F-P2.2 — closure-record integrity and explicit-path canonical coverage:**
record and enforce a complete closure-code dependency manifest; bind or
recompute-check the sidecar payload; constrain dust status to an enum; expose
validated dust provenance hashes in P4/P5 outputs; add a canonical explicit-SED
nine-group artifact/validator branch; and record working-tree cleanliness.
This recommendation is not started pending driver approval.
