# Approved next F-P1H bundle: trust-root repair plus failed-wind cross-check

Date proposed: 2026-09-03
Fable approval date: 2026-09-04
Status: historical approved plan; Part A and Part B' implementation complete;
the pre-retirement bundle-end AGY/Claude Opus 5 audit records are retained;
no AGY audit is pending or authorized
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Purpose fit

The project objective is a production-ready and publication-ready lagRamses
high-level hydrodynamics physics stack centered on radiative transfer, stellar
and AGN feedback, and dust. This bundle repairs the publication provenance of
the stellar-feedback source and advances physical lifetime, wind, and pre-SN
structure evidence. It does not reopen completed RAMSES topology, header,
HDF5, restart, or generic distributed-runtime work.

## Part A -- mandatory source trust-root repair

1. Add a code-owned lock profile that pins the admission/source IDs, release
   root, article/data DOI, Zenodo record ID and filename, CC-BY-4.0 identity,
   exact non-empty five-file inventory, strict integer byte counts, per-file
   SHA256 and published MD5 values, version-record SHA256, and expected
   composite fingerprint.
2. Require both all requirement booleans and no blockers for a pass in the
   validator and registry. Convert every malformed-input/runner exception into
   a controlled fail-closed result or admission error.
3. Reject all symlinks and non-regular files, unsafe path components,
   duplicate/extra/empty inventories, wrong scalar types, null identities, and
   invalid calendar dates.
4. Rename the requirement to `hash_locked_local_source_mirror`; do not claim
   operating-system immutability.
5. Derive redistribution from the pinned Zenodo license evidence and pinned
   attribution. The editable local terms file is reported but not
   verdict-bearing.
6. Require publication false while the F-P1 sidecar is production-blocked and
   accept only the exact pinned artifact set at repository-relative paths
   confined beneath `simulation/snrt`.
7. Add a temporary-fixture adversarial matrix covering every AGY/Codex bypass
   and prove the genuine staged fingerprint is unchanged before/after tests.

Part A accepts only when the genuine staged package has no blockers and
reproduces
`3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`,
all adversarial cases fail closed by name without uncaught exceptions, and all
higher production/publication/runtime gates remain false.

## Part B' -- source-derived failed-wind cross-check and inquiry packet

1. Read only manifest-pinned local bytes from Boccioli--Roberti `LC18.zip` and
   Limongi--Chieffi CDS `table5.dat`/`table7.dat`. Reuse existing source parsers
   rather than introducing a second interpretation of the table formats.
2. Join the 108 Boccioli--Roberti summary rows one-to-one with the 108 CDS
   `(initial mass, metallicity label, rotation)` rows. Record archive member,
   summary wind mass, Wind-table element sum, CDS cumulative wind
   `M_initial-M_total(PSN)`, PSN lifetime, signed differences, exploded state,
   and table7 structure fields or explicit nulls.
3. Use the 52 successful rows as a release-internal control: their Wind tables
   must be nonzero and their summary-to-table residual distribution must be
   reported without inventing a tolerance. The independent pre-implementation
   check showed that all 52 summary-to-CDS differences exceed the nominal
   0.005 M_sun table5 half-bin (maximum 0.5842 M_sun), so CDS agreement is a
   measured cross-source question, not an acceptance assumption. Record all
   56 failed rows as unresolved when the release Wind table is zero while its
   summary wind is positive; retain zero/quantized CDS values exactly.
4. Check that cumulative phase ages are strictly increasing and total masses
   non-increasing across the 3--8 unique phase points supplied per model after
   exact duplicate collapse. List duplicates and violations rather than
   modifying or filling values.
5. Require 96 table7 matches and 12 explicit null structure records. Do not
   infer missing compactness, binding energy, or envelope fields.
6. Generate an internal JSON review artifact and an unsent markdown inquiry
   packet with exact checksums, affected rows, and questions about failed Wind
   tables, summary wind provenance, machine-readable explosion energies, mass
   cuts, and the BHSN exception list.
7. Declare the CDS catalogue rights status accurately, emit zero canonical
   rows, leave all physical-package blockers untouched, and add the cross-check
   test to the F-P1 runner.

## Deferred coordinate validator

Do not register `fp1.coordinate_hull_and_population.v1` yet. It becomes useful
only after an approved target metallicity domain and rotation/binary decision
exist. At that point it must quantify uncovered source cells relative to that
approved domain and test every branch/axis endpoint and ±epsilon while
rejecting cross-source/engine interpolation, nearest-node substitution,
clamping, and flattened branch unions.

## Verification and stopping conditions

- Focused rights/registry/fate-sidecar adversarial tests.
- Exact 108/108, 52/56, and 96/12 cross-check counts with phase monotonicity.
- `run_fp1_population_fate_contract.sh` passes.
- `run_g2_preflight.sh` reaches its expected fail-closed
  `G2_PREFLIGHT_BLOCKED` state.
- JSON syntax, Python compilation, and `git diff --check` pass.
- No writes under `external/g2_candidates`, no download, and no author contact.
- Physical-node inventory remains empty and no canonical conversion, runtime
  deposition, production, or publication approval is enabled.
- The historical bundle-end AGY/Claude Opus 5 review was handled before AGY's
  retirement. Do not commission or retry AGY. Future completed steps use
  Claude Opus 5; a conditional or negative Opus result may use the existing
  GPT-5.6 Sol adjudication path, and the resulting next-bundle plan is reviewed
  by Grok under the current governance.
