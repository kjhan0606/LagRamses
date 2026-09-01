# AGN nine-group source-ledger validation

Date: 2026-09-02
Stage status: internal gates **PASS**; final Claude Opus 5 re-audit **PASS**

## Scope and acceptance contract

This stage closes the source-side M9 and Mo2 findings from the 2026-09-01 RT
audit. It promotes the canonical AGN source ledger from the retained five-group
pilot layout to the configured P0 nine-group layout, including the `2–10 keV`
tail. It does not promote a transported temperature or ionization field: that
requires the next timestep/spatial-resolution gate.

The hard gates require:

- exactly ten configured edges and nine CSV photon groups;
- bit-exact equality among `config/p0_photon_group_edges_ev.txt`, metadata
  `group_edges_ev`, and every serialized group interval;
- the declared interval convention `[lower, upper)`, except the final interval
  is closed on the right;
- a positive `[2000,10000] eV` photon rate and positive H I/He I/He II
  microphysics closure;
- exact zero cross section for every group wholly below each absorber
  threshold;
- agreement of CSV and metadata group totals and exact candidate-to-ledger
  identity, source-kind, and position rows;
- current generator, candidate-ledger, edge-table, chemistry, static-input,
  source-rebind, and short transport-control hashes;
- current gas-input, gas-sidecar, zoom-manifest, and external-HDF5 manifest
  hashes;
- a nine-group production-runner control passing photon, H/He, fixed-point,
  photoelectron-energy, root, and finite-value gates.

The legacy five-group converter remains available only through the explicit
`--legacy-five-groups` control. It is not accepted by this canonical gate.

## Boundary integration repair

The old SED quadrature selected both ends of every group. At an ionization
threshold located exactly on a group upper edge, the non-zero right-hand
Verner value entered the below-threshold trapezoid and produced a small but
formally false opacity.

The repaired integration keeps the continuous photon-number integral intact:
each group quadrature includes its two mathematical endpoints, but each
absorber numerator begins at `max(group lower edge, absorber threshold)` and
ends at `min(group upper edge, fit maximum)`. A group ending at a threshold
therefore has an empty absorber interval and exactly zero opacity, while the
next group owns the above-threshold integral. This avoids the missing-last-cell
error that a naive `E < upper` sample mask would introduce.

The canonical closure now has exactly zero:

- H I cross section in groups 0–3, through `[11.2,13.6] eV`;
- He I cross section in groups 0–4, through `[13.6,24.59] eV`;
- He II cross section in groups 0–5, through `[24.59,54.42] eV`.

The parameterized AGN SED begins at `10 eV`. Groups 0 and 1 are therefore
explicit zero-photon controls. Group 2 is serialized on the configured
`[5.6,11.2] eV` interval but is explicitly marked partially supported and is
integrated only over `[10,11.2] eV`; no spectrum is extrapolated below
`10 eV`. Groups 3--8 are fully supported.

## Restored hard-X-ray group

For the ten selected output-00017 sinks and the declared unobscured Sazonov-
style pilot SED, the restored `[2,10] keV` group has:

| Quantity | Value |
| --- | ---: |
| Photon rate | `3.5299649e51 s^-1` |
| Photon-weighted mean energy | `4023.5946 eV` |
| Photon-number ratio to `0.5–2 keV` | `0.219425` |
| Energy-power ratio to `0.5–2 keV` | `1.01523` |
| Fraction of total emitted photon number | `2.70656e-4` |
| Fraction of supported/in-band SED energy power | `0.0437846` |
| Fraction of candidate bolometric luminosity | `0.00794605` |

The group is rare by photon number but carries 4.38% of the represented
`10 eV--10 keV` SED power, 0.795% of the ten candidates' total bolometric
luminosity, and approximately the same energy power as the `0.5–2 keV` group.
It therefore cannot be silently omitted from diffuse pre-heating and
penetration studies.

Its H I, He I, and He II averaged cross sections are respectively
`3.14794e-25`, `9.58031e-24`, and `8.24702e-24 cm^2`; the corresponding mean
photoelectron excess energies are `2564.90`, `2558.20`, and `2537.36 eV`.

## Static-input and runner wiring

The canonical ledger was rebound to the existing coeval gas cube without
overwriting the historical five-group HDF5. The new external staged input is
`data/p4_coeval_static_rt_input_agn9.h5`, SHA256
`c25720164947df6a7f64b01d898401cf00e1420e74281f233323d92d2316a8b2`, with
source shape `(10,9)`. The ledger retains all ten candidate rows, of which five
have positive bolometric luminosity and five are explicit zero-luminosity
controls. Its tracked sidecar binds the original gas, gas metadata,
photon CSV, photon metadata, edge table, and source-rebind tool.

The ignored HDF5 dependencies are enumerated with exact size, role, status,
and SHA256 in `data/agn_nine_group_external_assets.json`. The gate fails closed
if they are missing or altered. A publication archive identifier remains an
explicit final-release requirement; the local assets are reproducible but are
not represented as already deposited.

A `0.001 Myr`, S4, float64 production-runner control with 32 opacity iterations
passes:

| Gate | Result |
| --- | ---: |
| Photon arithmetic ledger | `1.42037e-16` |
| H ledger L1 | `1.76661e-8` |
| Maximum fixed-point residual | `2.33031e-10` |
| Photoelectron-energy ledger | `8.49515e-22` |
| Electron-root failures | `0` |

The control HDF5 is
`data/p4_validation/p4_agn9_stage4_0p001myr.h5`, SHA256
`8db00cbab5699523e9652e3c3c05bf033f69b4445d5af7bffc1c3ea7f11feeef`.
A separate full-CFL one-step probe (`0.1045 Myr`) failed the H ledger at
`7.5731e-3`; this is preserved locally as evidence of the already declared
source-cell timestep problem and is not used as a passing artifact. The next
B3 gate must quantify timestep and spatial convergence with the nine-group
input.

## Reproduction and provenance

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
.venv/bin/python tests/agn_photon_ledger.py
.venv/bin/python tests/p0_sed_closure.py
.venv/bin/python tools/validate_agn_nine_group_ledger.py \
  --output data/agn_nine_group_validation.json
.venv/bin/python tests/agn_nine_group_artifact.py
```

- Canonical validation JSON SHA256:
  `2a66b4f3893c2ccf52b67cc2f4cbe10ac9102cf717cbbefec2e0fe6e16ecde75`.
- Validator SHA256:
  `749d26595a02470b4cf4c9b1733b45a59235c9777eb73cb4120b3f72ce922fd1`.
- Canonical photon CSV SHA256:
  `bfdab21cffc3fb9a7d02e7d6a3f2e892fd6e7d87076022ecb6b2fda902e19a4e`.
- Canonical photon metadata SHA256:
  `506bb47abd2481d7caab1f5ac3fa0451af10e0ec5907ec61ce03738639e9a53a`.
- External-asset manifest SHA256:
  `0dea881608b7fdbfc6045db08446e74f674ecac23e60efe9528b514fa61095a4`.

The initial independent audit is recorded in
[`provenance/claude_opus5_agn_nine_group_audit_2026-09-02.md`](../../provenance/claude_opus5_agn_nine_group_audit_2026-09-02.md).
It returned `CONDITIONAL PASS`; its F1--F5 findings are the remediations
documented above. The focused
[`claude-opus-5` re-audit](../../provenance/claude_opus5_agn_nine_group_reaudit_2026-09-02.md)
independently closed all five findings and M9/Mo2 with final `PASS`. Its two
new low-severity observations—full bolometric coverage and SED-normalization
sensitivity—remain explicit obligations for the final publication gate.

The SED and unit conversion are source-ledger closures, not an LRD obscuration
model. The intrinsic normalization and escape fraction remain explicit
systematic parameters.
