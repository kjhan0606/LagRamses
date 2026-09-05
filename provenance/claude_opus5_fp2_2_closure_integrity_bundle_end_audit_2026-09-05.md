# Claude Opus 5 end-of-bundle audit — F-P2.2 closure-record integrity and explicit canonical coverage

Date: 2026-09-05
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5
Audit session: `f36174ed-6d4c-4072-a5c1-6baf16efabdd` (launcher PTY `96379`)
Mode: read-only; file reads/search only; no repository edits, jobs, builds, or
shell/hash recomputation were performed by the auditor.

## Verdict

`CONDITIONAL PASS`

The core F-P2.2 implementation is substantially sound. I1 (complete closure
dependency manifests) and I2 (payload, status, and P4/P5 output provenance)
are implemented at the runtime boundary and fail closed. I3 (explicit
canonical coverage) is not yet closed because the canonical asset manifest,
explicit artifact test, and working-tree attestation do not provide a
reliable synchronization gate. F-P2.2 must therefore remain conditional;
the following repairs belong to the next approved bundle.

## Findings requiring closure before completion

### MAJOR-1 — explicit canonical asset manifest drift

`data/agn_nine_group_external_assets.json` contains a stale SHA-256 for
`data/p4_coeval_static_rt_input_agn9_explicit.json`. The live digest recorded
in `data/agn_nine_group_explicit_validation.json` differs, while the sidecar's
internal ledger/metadata bindings and the other explicit asset digests are
consistent. The validator checks only four expected external assets although
the manifest declares ten, so this drift was not detected.

Required disposition: validate every declared active-mode asset against the
manifest, regenerate the stale entry, and add explicit reproduction coverage.

### MAJOR-2 — explicit artifact freshness comparison is vacuous

`tests/agn_nine_group_explicit_artifact.py` compares fresh and canonical
criteria, but a passing validator makes both criteria vectors all-true. Its
working-tree assertions are likewise only type/length checks. Consequently it
does not perform the current-file hash sweep implemented by the pilot artifact
test, and the implementation evidence overstates what the explicit test
checks.

Required disposition: give the explicit test a complete provenance-hash sweep
and compare each canonical file digest to the live file.

### MAJOR-3 — working-tree attestation enforcement is asymmetric

The attestation values are truthful and internally consistent: both artifacts
record the same HEAD, dirty state, and status digest. However, the pilot test
compares a whole-repository porcelain digest and is invalidated by unrelated
provenance/scratch changes, while the explicit test enforces no meaningful
freshness condition. The attestation should either enforce HEAD plus a
truthful scoped clean/dirty state or treat the status digest as disclosure and
scope it to the relevant tree.

Required disposition: make the attestation criterion contentful and apply the
same meaningful rule to both pilot and explicit artifact tests.

## Additional findings and disposition

### MAJOR-4 — source-closure quadrature control (defer to science/closure gate)

The explicit SED closure integrates only table nodes and group boundaries,
without refinement or a convergence estimate. The canonical constant control
has only two nodes per group and uses linear interpolation of `q(E)` in energy,
which is not the exact implied `1/E` behavior; the canonical energy-fraction
check is tautological when all edges are table nodes. Add refinement,
convergence evidence, and interpolation-convention identity binding in a later
source-closure bundle. This is not a reason to promote the synthetic fixture
to a physical AGN SED.

Minor items to address with the F-P2.3 engineering repair include: making the
working-tree criterion non-tautological; deriving explicit-mode expectations
from SED support/integrals rather than the ledger under test; embedding a
synthetic/non-physical marker in the explicit ledger; enumerating and checking
`source_sed_contract.status`; cross-checking duplicated group fractions;
adding an explicit section and reproduction commands to the AGN validation
documentation; and refreshing the asset-manifest creation metadata.

The following remain correctly deferred: mixed STAR+AGN dust admission,
physical SED/escape/obscuration selection, dust-to-metal normalization,
`[40,120] M_sun` yield resolution, scattering/IR/grain physics, live RHD,
production convergence, and publication claims. The AGN CSV late-write issue,
unified edge parser, and root-relative manifests are later engineering items.

## Questions answered by the audit

1. The AGN and source-weighted dust manifests cover the reachable closure
   dependencies and enforce exact role/path/hash sets. Absolute resolved paths
   remain a portability limitation.
2. The canonical payload hash is deterministic, order-independent, rejects
   NaN, binds the serialized values/arrays, and has a meaningful transported-
   opacity tamper test. It proves payload integrity since build time, not
   independent reproducibility from authenticated inputs.
3. Dust status vocabulary is schema-consistent, and P4/P5 propagate validated
   metadata, payload, source-table, and builder hashes. v1 reference controls
   retain empty source-bound fields.
4. The explicit validator genuinely routes through the repaired source-bound
   closure and distinct explicit artifacts, but does not yet provide numerical
   quadrature validation, canonical source-bound dust coverage, or independent
   status derivation.
5. Attestation is useful and honest as recorded, but synchronization is not
   closed because of MAJOR-1 and MAJOR-2 and the asymmetric enforcement in
   MAJOR-3.

## Closed in this bundle

Complete exact-role live-rehashed manifests; canonical payload self-hashes with
negative tamper tests; fixed schema-consistent dust statuses; symmetric P4/P5
provenance propagation; a separately routed explicit nine-group branch; and a
truthful (though not yet sufficiently enforced) working-tree attestation.

## Recommended next bundle

`F-P2.3 — canonical asset synchronization and source-closure quadrature
control`:

1. Check every declared active-mode asset, repair/regenerate the stale digest,
   and add explicit static/transport reproduction commands.
2. Give both artifact tests the same complete provenance-hash sweep and replace
   whole-repository porcelain equality with a stable, meaningful attestation.
3. Make the attestation and explicit-mode criteria independent of the artifact
   under test.
4. Add quadrature refinement/convergence and bind the interpolation convention
   into the source identity.
5. Add non-physical labeling, status enum validation, duplicate-fraction
   checks, and explicit validation documentation.

No next bundle was started by this audit.
