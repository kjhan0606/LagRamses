# P4 transport conservation validation

> Historical status (2026-09-01): the HDF5 results below were generated with
> the retired atom-inventory attenuation cap. They remain diagnostic controls
> and must not be presented as outputs of the current B2 solver. The current
> transport-coupled validation is
> [`B2_PRODUCTION_SOLVER_VALIDATION.md`](B2_PRODUCTION_SOLVER_VALIDATION.md).

## Scope

This is a fixed-temperature, zero-dust, AGN-only static S8 pilot on the
canonical output-00017-derived input. It is an interface and numerical
validation artifact, not a production RHD or LRD-observable result. Stellar
sources, live AGN feedback, radiation pressure, dust opacity/IR re-emission,
and non-equilibrium metal chemistry are not active.

Dust activation is now wired through the explicit
[`P4_DUST_OPACITY.md`](P4_DUST_OPACITY.md) sidecar contract, but this retained
validation artifact intentionally remains a zero-dust control.

The native RAMSES implementation now uses the same canonical nine-group
contract, including the `[2000,10000] eV` hard-X-ray group. This is a native
Fortran wiring correction documented in
[`SNRT_NATIVE_GROUP_CONTRACT.md`](SNRT_NATIVE_GROUP_CONTRACT.md); it does not
retroactively promote the historical HDF5 result below to a nine-group live
RAMSES result. The checked-in native numbers are a reference-control closure,
and `SNRT_RT_ENABLE` still requires an explicit contract path, CUDA, and the
separate physical-source approvals.
The native contract records whether group fractions are intrinsic or escaped;
the resolved-domain injection gate accepts only `escaped`, and the checked-in
reference control requires `SNRT_ALLOW_REFERENCE_CONTROL=1`.

The native nine-group wiring does not yet close the emission-side
photon-number-weighted mean against the H absorber-weighted heating excess,
and it still lacks secondary ionization and recombination. For the new
`[2000,10000] eV` group these are explicit science-gate limitations, not
properties of the historical HDF5 result; see
[`SNRT_NATIVE_GROUP_CONTRACT.md`](SNRT_NATIVE_GROUP_CONTRACT.md) for the
quantified gap and follow-up gate.

The gas staging audit is internally consistent: the 32^3 mesh has cell width
`3.707491434780762e21 cm` (1.2015 kpc), side length 38.4485 kpc, AMR coverage
within `2.3e-16` of unity, and density mass-balance error
`2.10e-16`. Its mean `n_H` is `0.145638 cm^-3`, corresponding to
`2.4320e68` H nuclei or `2.0457e11 Msun` of hydrogen. This is a selected dense
galaxy-scale subvolume, not a cosmological mean-volume cell.

## Historical algorithmic correction

The initial direct-attenuation run removed photons from the radiation field
faster than the finite H/He inventory could absorb them: its unallocated
primary-photon fractions were 0.936 (fesc=1) and 0.842 (fesc=0.1). Those
files are retained as invalid diagnostic controls.

The retained P4 artifacts applied a local photon-conservative attenuation cap. It
limits gas absorption by available H/He atoms, including the requested
secondary-ionization demand; excess photons remain in the radiation field.
The dust component is left uncapped. The local source-plus-absorption operator
also now uses the exact constant-source response
`phi(tau)=(1-exp(-tau))/tau`: source photons are not all attenuated as if they
were injected at the beginning of the step. The validation runner uses
float64 for absolute cgs photon inventories and records the actual kernel
absorption in the finite-volume ledger. When dust is activated, its absorbed
photon number is recorded separately and excluded from the H/He primary-
absorption closure.

## Retained pre-B2 full-run result

Both historical source-exact reruns use 71 S8 steps over 6.3697581 Myr, reduced light
speed `0.01 c`, transport CFL `0.4`, 24 implicit recombination iterations, and
`photon_conservative_absorption=true`.

| Case | Emitted photons | Absorbed | Escaped fraction | Final in-domain fraction | Mean x_HII | x_HII >= 0.5 | Unallocated primary | Photon ledger error | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| fesc=1 | `2.16165e69` | `3.94924e68` | `2.97659e-2` | `7.87539e-1` | `0.538974` | `0.534363` | `3.18e-17` | `1.74e-16` | PASS |
| fesc=0.1 | `2.16165e68` | `1.30596e68` | `6.88015e-3` | `3.88973e-1` | `0.0969925` | `0.0897217` | `2.91e-17` | `3.89e-16` | PASS |

H/He ledger relative errors are below `1.6e-17` in both cases; all recorded
datasets are finite, ion fractions remain in bounds, and H/He recombination
densities remain non-negative. The mean absolute H II-fraction difference
between the two escape cases is `0.441976`.

## Artifacts

- `data/p4_transport_validation/p4_coeval_fesc1_s8_6p37myr_pcap_f64_sourceexact.h5`, SHA256 `a0ca30caeaca56bc79d9e8c160497e026ff3a7b0635b6072880602d7841d7516`
- `data/p4_transport_validation/p4_coeval_fesc0p1_s8_6p37myr_pcap_f64_sourceexact.h5`, SHA256 `c3943b1fd44fe7cb1a2c962ac48831e2cd4eb0cc5a568501084deacd7e1312d2`

These cap-era files and the previous v3/v2 files are retained controls. None
is a current-solver reproducibility artifact; P4/P5 science runs must be
regenerated with the B2 C2-Ray-style closure and its recorded limiter/fixed-
point diagnostics.

The HDF5 files include per-step source/absorption/boundary-escape inventories,
aggregate photon residuals, H/He transition ledgers, and numerical stability
attributes. The source SED and escape fractions remain explicitly labeled
pilot assumptions.
