# Final stage-3 audit — RSLA and refinement

Model: `claude-opus-5`
Mode: read-only
Repository HEAD during audit: `ca90a391296e4fbd99d183df3850de10c537cef4`

All four requested tests passed. The auditor independently reproduced every
v3 extrapolation value from the canonical matrix:

- inverse-ĉ intercepts `1.0073980090245733`, `1.0070156249072968`, and
  `1.0069731377831546`; padded upper bound `1.007822880265992`;
- photon-storage intercepts `1.0081343557466564`, `1.007069553822794`, and
  `1.0069459259966635`; padded upper bound `1.0093227854966493`;
- selected coordinate `photon_storage_fraction`;
- production RSLA term `0.017269837730221192` and total envelope
  `0.019764600042460632`, 98.823% of the 2% gate.

It independently confirmed that the four envelope terms sum exactly, photon
storage tracks the free-streaming estimate across the matrix, excluding the
non-asymptotic `0.001c` point is conservative, the P1 helium statement and
history are physically accurate, B2 is the same fixture, and all source and
cross-artifact hashes close.

## Prior-finding closure

| Finding | Status |
| --- | --- |
| F1, false helium attribution | **CLOSED** |
| F2, finite `0.03c` reference | **CLOSED** |
| F3, hidden mesh allowance | **CLOSED** |
| F4, stale B2 HEAD | **CLOSED** |
| N1, coordinate sensitivity omitted | **CLOSED**, subject to G1/G2 documentation below |
| N2, wrong old-P1 history | **CLOSED** |
| N3, undeclared escape tolerance | **CLOSED** |
| N4, insufficiently independent artifact test | **CLOSED** |

## New findings

- **G1 — LOW, conditional:** the predeclared acceptance-contract bullet still
  called the envelope “inverse-ĉ extrapolated”, although photon storage is the
  binding coordinate. The report body and implementation were correct.
- **G2 — LOW–MEDIUM, disclosed caveat:** with only 1.2% gate headroom, the
  coordinate set must be declared closed. The two physically motivated
  coordinates pass; arbitrary extension to a coordinate such as `(1/ĉ)^2`
  need not pass and is not the expected leading RSLA scaling.
- **G3 — LOW:** the intermediate v2 JSON itself was overwritten rather than
  preserved, although its complete numerical outcome and audit are retained in
  tracked provenance and v3 is strictly more conservative.

## Verdict

**CONDITIONAL PASS**

All numerical, source-binding, and prior-finding checks passed. Final PASS was
withheld solely until the acceptance bullet names the actual two-coordinate
gate and the coordinate set is explicitly closed for this version.
