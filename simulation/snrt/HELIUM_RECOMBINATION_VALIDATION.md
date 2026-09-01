# Case-B helium recombination validation

Date: 2026-09-01

## Physics contract

All H/He ionization solvers now use one case-B contract. The radiative rates
follow Hui & Gnedin (1997):

- `alpha_HeII,B = 1.26e-14 lambda_HeI^0.75 cm^3 s^-1`, with
  `lambda_HeI = 2(285335 K)/T`;
- `alpha_HeIII,B(T) = 2 alpha_HII,B(T/4)`.

The standard He II dielectronic contribution is added separately to the
radiative case-B coefficient. The former mixture of a case-A He II radiative
rate and `4 alpha_HII,B(T)` for He III has been removed from every caller.
The primary reference is [Hui & Gnedin (1997)](https://doi.org/10.1093/mnras/292.1.27),
Appendix A.

## Temperature-resolved coefficients

| T [K] | total alpha_HeII,B [cm3/s] | alpha_HeIII,B [cm3/s] |
| ---: | ---: | ---: |
| 10000 | 2.616130035e-13 | 1.544760722e-12 |
| 20000 | 1.555560645e-13 | 9.085606647e-13 |
| 40000 | 9.442148745e-14 | 5.183631355e-13 |
| 100000 | 6.570365334e-13 | 2.337400265e-13 |

At each temperature, the case-B He II radiative rate is below its case-A
control. The non-monotonic total He II rate at 100000 K is the separately
retained dielectronic contribution, not a return to case A.

## One-zone gate

[`tests/helium_recombination.py`](tests/helium_recombination.py) evolves pure
He II and pure He III at fixed `n_e = 1 cm^-3`, zero photoionization, and each
of the four temperatures above. The tabulated reference values come from the
published Hui--Gnedin fits and are independently re-evaluated with NumPy rather
than by calling the production coefficient functions.

The first test integrates three recombination times using 512 backward-Euler
substeps and compares with `exp(-3)`. Because `alpha*dt = 3/512` is fixed by
construction, its identical temperature outputs are a cross-module consistency
check, not four independent stiffness measurements. The maximum relative error
is `0.00879332` for both ions, below the predeclared `0.02` threshold.

A second test uses one fixed physical interval, `2.0e12 s`, at all four
temperatures. Its final fractions therefore differ with the temperature-
dependent rates. Maximum relative errors are `0.00168485` for He II and
`0.00932733` for He III, again below `0.02`.

The thermal counterpart is also gated in
[`tests/b1_thermal_coupling.py`](tests/b1_thermal_coupling.py): radiative H II,
He II, and He III recombination cooling per event must remain below
`1.5 k_B T`, while the separate He II dielectronic cooling-to-rate ratio must
equal its matched Hui--Gnedin coefficient ratio. He III cooling is implemented
as `beta_HeIII,B(T) = 8 beta_HII,B(T/4)`; this closes the 4x prefactor defect
found in the first Opus 5 audit.

The machine-readable result is
[`data/helium_case_b_recombination_validation.json`](data/helium_case_b_recombination_validation.json),
SHA256 `c51cb89dce2311b1f68f8301efa9ddbcfa9fbece91cbb453d809f0e56c0ba305`.
The artifact records the JAX version and hashes of the one-zone test, B1 test,
coefficient source, cooling source, and implicit solver. The fail-closed
[`tests/helium_recombination_artifact.py`](tests/helium_recombination_artifact.py)
contract test rejects stale hashes, non-finite values, failed error thresholds,
or a temperature-degenerate fixed-time result.

## Regression closure

The B1 thermal, B2 zero-density, P2/P3 implicit/sharding, and canonical B2
transport gates pass after the change. The H-only B2 physical payload is
unchanged; its regenerated artifact differs only in runtime/provenance fields.
Full transport-coupled He-front convergence remains reserved for the planned
Coupled H+He gate.

The initial independent audit is recorded in
[`provenance/claude_opus5_he_recombination_audit_2026-09-01.md`](../../provenance/claude_opus5_he_recombination_audit_2026-09-01.md).
Its conditional findings were corrected. The
[`claude-opus-5` final re-audit](../../provenance/claude_opus5_he_recombination_reaudit_2026-09-02.md)
closed all four conditions with a `PASS`; this is the stage-1 close gate.
