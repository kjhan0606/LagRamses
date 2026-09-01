# Independent stage-4 re-audit — AGN nine-group source closure

Model: `claude-opus-5`
Mode: read-only
Repository HEAD during audit: `f792791069810c92c3e23c2339cc16370cef15b3`

## Verdict

**PASS**

Stage 4 may close. M9 and Mo2 are closed on the production path; F1–F5 are
genuinely remediated and were verified independently rather than by trusting
the remediation's own claims. B3 field convergence is honestly deferred and
correctly not a prerequisite here.

## F1–F5

| # | Finding | Status | Independent evidence |
| --- | --- | --- | --- |
| F1 | Stale `primordial.py` hash in He artifact | **CLOSED** | The JSON diff changes only the stale source hash. All physical values are byte-identical, the report SHA is current, and the full tree contains the old source hash only in the deliberately historical five-group pin. Errors remain `heii=0.00168485`, `heiii=0.00932733`. |
| F2 | `[5.6,11.2] eV` mislabelled; in-band and bolometric power conflated | **CLOSED** | Group 2 is `agn_sed_partially_supported_10ev_to_upper` with `[10.0,11.2]`; groups 0–1 are below support and 3–8 are fully supported. Independent SciPy quadrature recovered `4.3784571829%` of represented SED power and `0.7946053108%` of candidate bolometric luminosity, agreeing within `8.6e-9` and `2.6e-8` relative. |
| F3 | Hard-coded group indices and threshold slices | **CLOSED** | Both X-ray groups are found by unique exact edge matching and species masks derive from production Verner thresholds. Four alternative edge tables retained a complete disjoint partition; moving the 13.6 eV edge changed the H I forced-zero count from four to three as required. |
| F4 | No external asset manifest or deposit status | **CLOSED** | Four ignored HDF5 dependencies have path, kind, status, size, SHA256, and role. Six tampered-manifest variants all failed closed, including one falsely claiming a completed deposit. The baseline reproduced the tracked artifact byte-for-byte. |
| F5 | Candidate identity, luminous count, and gas/zoom hashes ungated | **CLOSED** | Ordered IDs, kinds, positions, uniqueness, and source count are gated. The ledger has 10 rows and exactly 5 luminous rows. All seven static-sidecar hashes close. Every one of 90 CSV photon entries equals `Lbol * fesc * rate` to zero ULP. |

## M9 and Mo2 closure

- **M9 group span — CLOSED:** nine groups through 10 keV are present in the
  config, metadata, CSV, static HDF5 `(10,9)`, and runner control `(9,)`.
- **M9 edge agreement — CLOSED:** enforced at four layers and confirmed with a
  mutation test.
- **M9 fields — honestly deferred:** timestep and spatial convergence remain
  the next B3 stage.
- **Mo2 — CLOSED:** cross section is exactly zero in H I groups 0–3, He I
  groups 0–4, and He II groups 0–5. None of the three old leaked values remains.

The `0.1045 Myr` full-CFL failure is preserved with
`validation_passed=False` and H-ledger L1 `7.5731e-3`, matching the report.

## New low-severity observations

- **N1 — disclosure:** the nine-group table represents 18.148% of candidate
  bolometric luminosity; 81.85% lies outside `[10 eV,10 keV]`. Keep this
  explicit in final publication-facing limits.
- **N2 — parameterization:** `nu Lnu(13.6 eV)/Lbol = 0.1` implies a 2–10 keV
  bolometric correction near 126, versus an observational range roughly
  20–50 at these luminosities. The normalization is already a tunable input,
  but its scientific choice needs justification or sensitivity analysis.
- Additional observations were non-blocking: inert geometric means in zero
  groups, incomplete transport-HDF5 regeneration commands, duplicated SED
  minimum constants that do fail closed, literal slices retained only in a
  test, and gas-payload semantics not encoded in the validator. The auditor
  independently found all 14 non-source datasets and attributes bit-identical
  and reproduced `sources/cell_index` from candidates and the zoom cube.

## Validation class

Photon/photoelectron ledgers, group totals, and fixed-point residuals are
arithmetic closure and solver self-consistency, not physical field validation.
No historical five-group field result is presented as current. The auditor ran
23 tests, including all nine required tests: **23 green, 0 red**.
