# F-P1 source-package selection and staging plan — 2026-09-02

Project: `/gpfs/kjhan/LRD_JWST`  
Status: **review-only; no production source selected**  
Current G2 bundle: **candidate audit and high-mass review projection complete;
Opus engineering verdict conditional, scientific verdict blocked**

## Recommended staging order

1. **Sukhbold et al. 2016** is the first engine-specific validation branch
   for the 40–120 M☉ seam. The staged W18/N20 result files demonstrate that
   high-mass outcomes are non-monotonic in ZAMS mass. This branch is solar-
   metallicity validation only; its MPA archive terms permit internal
   non-commercial research with citation but do not authorize third-party
   redistribution.
2. **Boccioli–Roberti 2026 LC18 release** is the first multi-metallicity/
   rotation comparison branch because its Zenodo package is CC BY 4.0 and its
   LC18 coordinates cover 13–120 M☉, four source [Fe/H] values, and three
   rotations. It remains quarantined: failed-model wind summaries and Wind
   tables disagree, and per-model injected energy, canonical momentum, and
   age-resolved cumulative release are absent.
3. **Limongi–Chieffi 2018 CDS** remains an independent structure/wind semantic
   comparison. Its source-author metallicity mapping and recommended set-R
   semantics must be retained exactly; it does not by itself cover the
   runtime 0.8–120 M☉ source ranges or supply the complete injection ledger.
4. A modern multi-metallicity terminal-fate package such as Ugolini et al.
   2025 may be evaluated later, but it is not locally staged and no values are
   to be reconstructed from the paper.

This ordering separates a useful validation branch from a production choice.
It does not imply that Sukhbold, LC18, or any other candidate is the project
truth.

## Required staging record for each candidate

Before a candidate can feed a physical source-node table, retain the exact
source package and an immutable record of:

- citation, version/data DOI, license and redistribution terms;
- per-file SHA256 and a composite package fingerprint;
- source metallicity definition and solar abundance set;
- ZAMS mass, rotation, binary-state, engine, and source-hull coordinates;
- lifetime and age-release convention;
- separate pre-terminal wind, terminal ejecta, fallback/envelope ejection,
  remnant, PPISN/PISN, decay, energy, momentum, and deposition semantics;
- the conversion code revision and a channel-specific closure report.

Missing fields remain missing. In particular, a source's diagnostic explosion
energy, binding energy, or integrated yield is not silently converted into an
injected energy, momentum, or age-resolved source history.

## Promotion gate

The current `fp1_fate_admission_sidecar_v1.json` remains
`blocked_review_only`. Physical nodes may be added only after the candidate's
source-specific audit proves exact source-node keys, half-open mass-cell
coverage, mutually exclusive wind/terminal/remnant ownership, and all required
source quantities. The sidecar must then be regenerated with the package and
conversion fingerprints and a named approval id. Until that event, the
runtime's piecewise source-cell mode and unresolved bucket are safety guards,
not a physical fate model.

## Completed review-only execution in this bundle

The aggregate candidate audit now propagates acquisition-manifest and inline
parser integrity failures to a fatal top-level result. The Sukhbold projection
also couples the declared component to its actual source column, retains the
four source-package fingerprints, and parses the available high-mass W18/N20
yield tables plus the W18 implosion-wind records. The current projection
contains 26 low-mass review records and 19 high-mass review records; it emits
zero canonical rows and allows zero runtime deposition. Missing mass nodes
remain explicit and no interpolation is performed.

Evidence:

- `simulation/snrt/data/g2_candidate_source_audit.json`
- `simulation/snrt/data/g2_source_package_fingerprint_audit.json` — all 11
  staged manifest candidates and 65 files pass the manifest-scoped composite
  fingerprint check; these hashes identify the staged file sets only.
- `simulation/snrt/data/g2_sukhbold2016_candidate_audit.json`
- `simulation/snrt/data/g2_sukhbold_channel_projection_review.json`
- `provenance/claude_opus5_g2_source_package_staging_audit_2026-09-02.md`
- `provenance/claude_opus5_g2_source_package_staging_final_audit_2026-09-02.md`
- `provenance/agy_g2_source_package_staging_audit_2026-09-02.md`

## Next executable work

Resolve the remaining source-selection decisions and then, only with an
explicit project physics approval, design the source-specific conversion
contract. The immediate scientific prerequisites are an age-resolved wind
history (or quantified terminal lumping), decay horizon/projection, transition
seam ownership, explosion-engine uncertainty, canonical momentum/deposition
semantics, and MPA redistribution permission. The 40--120 M☉ review projection
is evidence for model dependence, not a production fate law.
