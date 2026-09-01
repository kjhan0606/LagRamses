# P5 thermochemistry validation

> Historical status (2026-09-01): the full-duration tables in this document
> predate B2 and used the retired atom-inventory attenuation cap with four
> opacity iterations. They remain convergence diagnostics, not current-solver
> science products. New P5 runs require at least 20 opacity iterations and
> record cap activity (which must remain zero) plus the maximum fixed-point
> residual. See
> [`B2_PRODUCTION_SOLVER_VALIDATION.md`](B2_PRODUCTION_SOLVER_VALIDATION.md).

Status: conservation gate passed; the tighter source-cell timestep treatment
substantially improves global and local convergence, but the local source-cell
convergence gate remains open.

## Current-solver FS2010 effect control (2026-09-02)

A matched 0.1 Myr, 32 cubed P5 pair remeasures the effect of secondary
ionization with the current uncapped solver. Both runs use S4, reduced light
speed `0.01c`, Courant 0.1, point deposition, a source-cell target of 0.25,
58 source subcycles, float64, and 32 opacity iterations. The OFF and FS2010 ON
runs both pass their internal P5 gates: maximum fixed-point residual is
`2.33e-10`, the worst H/He L1 ledger error is `1.55e-10`, and thermal closure
is `3.54e-7`.

FS2010 changes volume-mean xHII by `+2.24653e-8` and volume-mean temperature
by `-7.52274e-3 K`; the maximum local changes are `2.43194e-4` and `25.1009 K`.
All H I, He I, and He II secondary channels are non-zero. The small mean effect
is expected because the control ends at `<xHII>=0.978` and `<T>=5.58e6 K`:
FS2010's ionization-energy fraction near 200 eV falls from about 0.376 at
`xHII=1e-4` to 0.0022 at `xHII=0.9` and 0.0008 at `xHII=0.99`. It is a measured
property of this ionized, hot, UV-dominated AGN control, not evidence that
secondary ionization is generally negligible. The canonical regression bands
are tightened to `1e-8 < delta<xHII> < 5e-8` and
`-0.02 K < delta<T> < -0.002 K`. The source-bound report is
`data/p5_secondary_ionization_validation.json`, checked by
`tests/p5_secondary_ionization_artifact.py`. This pair is an effect and wiring
control, not a spatial-resolution or full-duration science promotion.

Both runs also require an internal multi-group photoelectron-energy ledger
below `1e-12`, zero electron-root bracket failures, and the explicit excitation
policy `radiative_line_escape_not_returned_to_gas`. The measured ON ledger is
`1.31e-17`; OFF is exactly zero.
The full-duration runs below were executed on `/gpfs/kjhan/LRD_JWST` with the
staged 32^3 P4 input, the audited AGN photon ledger, the P4 thermal atlas, S8
angular transport, zero dust, float64 arithmetic, and 6.3697581 Myr of
evolution. Short controls are explicitly marked.

## Changes validated

- The P5 runner now accumulates absorbed photons, unallocated photons, H/He
  transitions, recombinations, photoheating, background energy, and thermal
  residuals over the complete run rather than only the final outer step.
- The backward-Euler thermal solve brackets the root nearest the previous
  temperature. This is required because the tabulated cooling curve is not
  globally monotone; a full-range bisection can jump between thermal branches.
- Time-averaged opacity was retained in these historical controls; B2 has
  since replaced their capped transfer with an absorbed-rate fixed point.
- Thermal-atlas bound hits are recorded and are a hard failure for a promoted
  P5 result. None of the accepted conservation runs below hit a bound.
- The P5 runner now accepts the P4 dust-opacity sidecar and records dust
  absorption, dust heating, and absorption-only momentum separately; a
  non-zero-dust integration test passes with the staged Draine closure.
- Local source emission and absorption are now integrated with the exact
  constant-source response `phi(tau)=(1-exp(-tau))/tau`, rather than
  attenuating all newly emitted photons for the full step.
- The P5 runner supports a static source-cell limiter through
  `--source-cell-photons-per-neutral`; the source-limited runs below used a
  target of 1.0 and 0.25 while retaining static JAX control flow.
- The P5 runner supports an opt-in `compact3` controlled source deposition. It
  uses the normalized `[1/4, 1/2, 1/4]^3` kernel, renormalizes at boundaries,
  conserves each source group in host-side float64, and derives the limiter
  from the deposited source field. It is a numerical control, not a physical
  source-size model.

## Conservation results

| Run | outer CFL | thermal subcycles | time-average iterations | mean xHII | mean T [K] | absorbed photons | unallocated fraction | max H/He ledger | thermal closure | gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `C0.4_N4` | 0.40 | 4 | 4 | 0.106673 | 1.03835e7 | 1.81020e69 | 1.08e-17 | 2.18e-19 | 4.28e-8 | PASS |
| `C0.1_N1` | 0.10 | 1 | 4 | 0.106431 | 1.03835e7 | 1.80623e69 | 1.09e-17 | 6.34e-19 | 4.18e-8 | PASS |
| `C0.05_N1` | 0.05 | 1 | 4 | 0.065732 | 1.06261e7 | 1.89373e69 | 1.15e-17 | 5.78e-18 | 5.46e-8 | PASS |

The three artifacts have finite fields, valid H/He fractions, zero thermal
bound hits, photon-primary closure below `1.3e-16`, and
`validation_passed=True`.

## Convergence decision

The `C0.4_N4` and `C0.1_N1` states are close because they have nearly the same
effective transport timestep. Halving the effective timestep again to
`C0.05_N1`, however, changes the mean xHII from `0.106431` to `0.065732`.
The pair has mean absolute xHII difference `0.04160` and maximum local
difference `0.94210`. This is a timestep-sensitive point-source/cap coupling,
not a conservation failure: all ledgers still pass.

The P5 state is therefore not promoted to a science result, and stellar
sources and live feedback remain deferred. The dust wiring is now complete as
an opt-in sidecar path, but the Draine reference closure does not by itself
select a source-specific stellar/AGN spectrum or cell dust-to-metal model. The
next numerical implementation gate is a source-cell treatment that also
converges the unresolved source cell itself.
The current limiter makes the volume-averaged field much less timestep
sensitive, but the source cell and its immediate neighbours still show a
maximum difference of 0.1344 between the C0.1 and C0.05 runs. The next gate is
therefore a tighter source-cell limiter or controlled source deposition on a
refined mesh, followed by the same conservation matrix.

## Source-cell timestep results

The new runs use the same 32^3 mesh and 6.3697581 Myr duration as the controls,
with `--source-cell-photons-per-neutral 1.0` and four time-averaged-opacity
iterations. The limiter selected 13 source subcycles for C0.1 and 7 for C0.05;
the maximum initial source emission per actual substep was 0.9564 and 0.8881,
respectively.

| Run | outer CFL | source subcycles | mean xHII | mean T [K] | max H/He ledger | thermal closure | gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `source_limit1_C0.1_N1` | 0.10 | 13 | 0.027458 | 1.08653e7 | 2.27e-18 | 4.45e-7 | PASS |
| `source_limit1_C0.05_N1` | 0.05 | 7 | 0.027078 | 1.08676e7 | 1.16e-17 | 5.06e-7 | PASS |

Their field differences are mean absolute xHII `4.1464e-4`, maximum xHII
`0.1344`, and mean absolute temperature `2.283e3 K` (relative mean-temperature
difference `2.09e-4`). The maximum occurs at the dominant point-source cell
`(14,17,19)`, where xHII is 0.8556 versus 0.9900. The volume-weighted field is
therefore much better behaved, but the unresolved source cell is not yet a
converged science observable.

## Tighter source-cell limiter

The same full-duration pair was repeated with
`--source-cell-photons-per-neutral 0.25`. This selected 50 source subcycles at
C0.1 and 25 at C0.05; both runs passed the photon, H/He, thermal, and
finite-field validation gates.

| Run | outer CFL | source subcycles | mean xHII | source-cell xHII | mean T [K] | max H/He ledger | thermal closure | gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `source_limit0.25_C0.1_N1` | 0.10 | 50 | 0.0223196 | 0.9755265 | 1.08929e7 | 5.28e-18 | 5.36e-6 | PASS |
| `source_limit0.25_C0.05_N1` | 0.05 | 25 | 0.0223218 | 0.9970670 | 1.08929e7 | 5.91e-18 | 5.36e-6 | PASS |

The pair has mean absolute xHII difference `2.3907e-6` and maximum difference
`0.0215405`; both occur at the dominant source cell `(14,17,19)`. Relative to
the limit-1 pair, the maximum local difference falls from `0.1344` to `0.02154`
and the volume-mean difference falls from `4.1464e-4` to `2.3907e-6`. This
confirms that the remaining discrepancy is primarily source-cell temporal
stiffness. It does not yet establish spatial convergence of a point source
deposited into one unresolved cell, so the science gate remains open pending a
declared acceptance threshold or refined/controlled source deposition test.
The artifact pair is checked by
[`tests/p5_source_cell_convergence.py`](tests/p5_source_cell_convergence.py).

## Controlled source deposition

The unresolved point-source control was repeated with
`--source-deposition-mode compact3` and the same source target `0.25`. This
spreads each source over the local 3-by-3-by-3 kernel while preserving its
total group luminosity. The limiter then selected 8 source subcycles at C0.1
and 4 at C0.05; both outputs passed the photon, H/He, thermal, and finite-field
validation gates.

| Run | outer CFL | source subcycles | mean xHII | source-cell xHII | mean T [K] | max H/He ledger | thermal closure | gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `compact3_C0.1_N1` | 0.10 | 8 | 0.00882267 | 0.8321923 | 1.09628e7 | 1.40e-17 | 1.58e-7 | PASS* |
| `compact3_C0.05_N1` | 0.05 | 4 | 0.00885106 | 0.9665739 | 1.09628e7 | 1.40e-17 | 1.56e-7 | PASS* |

*`PASS` denotes the internal conservation/stability gate only. The pair has
mean absolute xHII difference `2.8383e-5` and maximum difference `0.13438`,
again at `(14,17,19)`. The maximum is larger than the `0.0215405` obtained by
the point-deposition limit-0.25 pair, while the mean ionization state also
changes materially (`0.00884` versus `0.02232`). Thus `compact3` is useful as
a controlled wiring and luminosity-conservation test, but it does not close
the source-cell science gate and must not be interpreted as a physical source
profile. The result is checked by
[`tests/p5_controlled_deposition.py`](tests/p5_controlled_deposition.py).

## Factor-2 refined-mesh control

To test whether the unresolved source cell is primarily a coarse-mesh effect,
the 32³ static input was prolonged to 64³ with piecewise-constant gas fields.
Each coarse source was split equally among its eight fine children, conserving
the full group luminosity and the source centre of luminosity. This is a
synthetic resolution control, not a reconstructed high-resolution hydro
snapshot. The pair below covers 0.5 Myr and uses point deposition on the fine
grid with source target `0.25`.

| Run | grid | outer CFL | source subcycles | mean xHII | fine source-block xHII | mean T [K] | max H/He ledger | thermal closure | gate |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `refined2_C0.1_N1` | 64³ | 0.10 | 25 | 0.000373783 | 0.9906207 | 1.15388e7 | 2.85e-17 | 1.88e-6 | PASS* |
| `refined2_C0.05_N1` | 64³ | 0.05 | 13 | 0.000372319 | 0.9822387 | 1.15388e7 | 4.42e-18 | 1.95e-6 | PASS* |

*`PASS` denotes the internal conservation/stability gate only. The fine
2³ block corresponding to coarse source cell `(14,17,19)` has mean absolute
xHII difference `1.4640e-6` over the full field and maximum difference
`0.008382`; the same maximum appears after coarse 2³ block averaging. This is
an improvement over the 32³ point-deposition limit-0.25 pair (`0.0215405`),
so the refinement control supports the interpretation that mesh resolution is
part of the source-cell discrepancy. It does not close the production science
gate because the gas prolongation is synthetic and the duration is only
0.5 Myr. The pair is checked by
[`tests/p5_refined_mesh_convergence.py`](tests/p5_refined_mesh_convergence.py).

The 0.5 Myr C0.4/N4 grouping check also passed the conservation gate with
`source_cell_subcycles=13` and `effective_subcycles=52`; it is not used as a
replacement for the full C0.1/C0.05 convergence pair because different outer
step remainders give a different local final-step partition.

## Artifacts

- [`FS2010 OFF 0.1 Myr`](data/p5_validation/p5_fs2010_off_s4_c0p1_0p1myr_source_limit0p25_remediated_n32_f64.h5), SHA256 `02ebe242f7fafea46edbb3bcc52a25d61cc26bb3e3f3d17794c28fcb6a822eea`
- [`FS2010 ON 0.1 Myr`](data/p5_validation/p5_fs2010_on_s4_c0p1_0p1myr_source_limit0p25_remediated_n32_f64.h5), SHA256 `a03e5ec5906fe2659b82e1fa6e9d162d805fd403aa3e09e9263e48fcc0b45d83`
- [`FS2010 matched-pair report`](data/p5_secondary_ionization_validation.json), SHA256 `0e2edc468e13bc1bd89dffcef0ce8d09fb9d919595a353e058fd270d081d5fb1`
- [`C0.4_N4`](data/p5_validation/p5_coeval_s8_sub4_6p37myr_conservative_timeavg4_f64.h5)
- [`C0.1_N1`](data/p5_validation/p5_coeval_s8_courant0p1_sub1_6p37myr_conservative_timeavg4_f64.h5)
- [`C0.05_N1`](data/p5_validation/p5_coeval_s8_courant0p05_sub1_6p37myr_conservative_timeavg4_f64.h5)
- [`source_limit1_C0.1_N1`](data/p5_validation/p5_coeval_s8_courant0p1_sub1_6p37myr_source_limit1_f64.h5), SHA256 `fe1e7e9485525e5b87fcb8e02a407af686c8ecc948e54df0431147d27ffaedb9`
- [`source_limit1_C0.05_N1`](data/p5_validation/p5_coeval_s8_courant0p05_sub1_6p37myr_source_limit1_f64.h5), SHA256 `55d90d7cd2c80770c3eee50434ae2342d10b0e5fe57b617d782690fdde3e9809`
- [`source_limit1_C0.1_N1` short 0.5 Myr](data/p5_validation/p5_source_limit1_courant0p1_0p5myr_f64.h5), SHA256 `38ea35911f84e864845ff900b22c9f460e38478b6b90170771a74023f0586`
- [`source_limit1_C0.05_N1` short 0.5 Myr](data/p5_validation/p5_source_limit1_courant0p05_0p5myr_f64.h5), SHA256 `5a8536b38a96dae63a5e404fe5260374f33cd15d7cfa30b2ac5023ffaee1e7e0`
- [`source_limit1_C0.4_N4` short 0.5 Myr](data/p5_validation/p5_source_limit1_courant0p4_sub4_0p5myr_f64.h5), SHA256 `533f249b36cfac29001f11bc004673c663a1266e68b18754b00f13593e61c73c`
- [`source_limit0.25_C0.1_N1`](data/p5_validation/p5_coeval_s8_courant0p1_sub1_6p37myr_source_limit0p25_f64.h5), SHA256 `23d0892d4734c9b851ba793dd62f08c18490ba71a068f6f6d507f08c82f92e32`
- [`source_limit0.25_C0.05_N1`](data/p5_validation/p5_coeval_s8_courant0p05_sub1_6p37myr_source_limit0p25_f64.h5), SHA256 `a014a3349c9f2cc40af35d56a665dd596afc0e7aa090cf0ced384d8fa980f7c4`
- [`compact3_C0.1_N1`](data/p5_validation/p5_compact3_courant0p1_sub1_6p37myr_source_limit0p25_f64.h5), SHA256 `e99e1c4bbf5f2591b60e13c1f5cc76b0b0cb2d30763e73e93ff5ac8f2abc8138`
- [`compact3_C0.05_N1`](data/p5_validation/p5_compact3_courant0p05_sub1_6p37myr_source_limit0p25_f64.h5), SHA256 `2ecc0354cd12925e77446b4f7973968817e49907e9caadd345b748aca43a6b75`
- [`refined2_static_input`](data/p5_validation/p4_coeval_static_rt_input_refined2.h5), SHA256 `b817d5e2b3352d6da175f8c604221ec6409312eeb4a4fd3e640272abaf71b4de`
- [`refined2_C0.1_N1` short 0.5 Myr](data/p5_validation/p5_refined2_point_courant0p1_sub1_0p5myr_source_limit0p25_f64.h5), SHA256 `6ababe261f134454e7ef884bf6e342b482f1ce10ee4669e2899b056c50aae287`
- [`refined2_C0.05_N1` short 0.5 Myr](data/p5_validation/p5_refined2_point_courant0p05_sub1_0p5myr_source_limit0p25_f64.h5), SHA256 `89e3b70485c3fa7a2896db755449d8e699b0b2abd3d5f3c505a37605c6a1b5a7`
- Earlier pre-fix and float32 outputs remain in the same directory as retained
  controls and must not be used as promoted P5 results.
