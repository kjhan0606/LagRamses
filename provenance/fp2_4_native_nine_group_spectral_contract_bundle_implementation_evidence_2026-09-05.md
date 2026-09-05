# F-P2.4 native nine-group spectral contract bundle — implementation evidence — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Bundle: F-P2.4, native G3 spectral/source boundary

Parent: F-P2.3 canonical asset synchronization and source-closure quadrature

Work location: `/gpfs`

Status: implementation, native remediation, direct evidence, and Claude Opus 5
closure audit complete. F-P2.4 is PASS and closed.

## Delivered implementation

The new native module
[`patch/lagRamses/snrt_spectral_contract.f90`](../patch/lagRamses/snrt_spectral_contract.f90)
owns the canonical `snrt_ngroups=9` and the ten exact group edges. It loads a
strict namelist from `SNRT_GROUP_CONTRACT`, with no JSON parser and no implicit
default. The loader checks source identity, source/edge SHA256 fields, commit
binding, interval convention, group means, energy-fraction closure, H/He
threshold support, and absorber-weighted excess-energy arrays. Missing or
malformed input leaves the contract unloaded.

The checked-in reference-control input is
[`config/snrt_group_contract_reference_control_v1.nml`](../simulation/snrt/config/snrt_group_contract_reference_control_v1.nml).
Its source metadata is the current P4 pilot JSON (`d2326e4...62909d`), and its
edge digest is the canonical file digest
`d28f78f1703730c6c0b9a7d183edfe0c5e6337979e737ce002a572b66fc53ff1`. The input
is explicitly `reference_control`, not an approved physical AGN SED. The
contract records `fraction_semantics`; resolved-domain injection accepts only
`escaped` fractions. A `reference_control` contract additionally requires
`SNRT_ALLOW_REFERENCE_CONTROL=1`, while candidate and intrinsic-fraction
contracts remain non-runtime-admissible.

`snrt_state` now allocates nine groups through the native contract dimension,
and checkpoint version 4 stores the source identity, edge digest, interval
convention, fraction semantics, and status. A checkpoint read rejects a different or
non-runtime-admissible spectral identity. Existing RAMSES HDF5 backup/restore
call sites still do not invoke these optional SNRT checkpoint routines; that
integration is retained as a later G5 task.

`snrt_ramses_driver` now loads the contract once when the latched
`SNRT_RT_ENABLE` path is first entered. It reports the source identity and
represented/unrepresented energy fractions and returns before changing the RT
state if the contract is absent, invalid, or candidate-only. The old four-group
source table is gone, the `[2000,10000] eV` group is included in the transaction,
and the unexplained `0.5` group-energy multiplier was removed. Group photon
energy is now the resolved radiated energy times the declared group fraction.

The H chemistry boundary accepts the contract's absorber-weighted H I excess
energy. Its positional mean-energy behavior remains for monochromatic
benchmark callers. The native live chemistry scope is still H-only; the He
tables are loaded and threshold-validated but not silently activated.

## Native evidence

The reproducible Fortran smoke is
[`tests/run_snrt_native_spectral_contract.sh`](../simulation/snrt/tests/run_snrt_native_spectral_contract.sh).
It was run with both GNU Fortran and the production `mpiifx` compiler. The
reference-control run reported:

```text
SNRT_SPECTRAL_CONTRACT_OK
```

The smoke covers:

- namelist loading through `SNRT_GROUP_CONTRACT`;
- unset environment, missing file, malformed namelist, unsupported version,
  malformed identity, edge-digest mismatch, unknown fraction semantics,
  candidate status, intrinsic fraction, and reference-control opt-in paths;
- nine groups and ten boundaries;
- exact canonical edge values and source fraction closure;
- H/He threshold-safe opacity/excess support and upper-bound rejection;
- checkpoint identity match and edge-digest mismatch rejection;
- represented-Lbol source closure without hidden rescaling;
- canonical edge mismatch rejection;
- sub-threshold H opacity rejection; and
- out-of-band representative-energy rejection.

The same runner builds a separate real-state checkpoint smoke. It allocates a
small eight-leaf, nine-group payload, writes the version-4
header/identity/cell/intensity/neutral records, rejects the file under a
candidate contract before mutating the state, then restores the reference
contract and round-trips every payload record. The final marker is:

```text
SNRT_CHECKPOINT_OK
```

The complete spectral runner passed with both GNU Fortran and `mpiifx`, ending
with:

```text
SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK
```

The extended NLTE smoke, including the absorber-weighted excess-energy path,
reported:

```text
SNRT_NLTE_COUPLING_OK tau=  1.498962E+01 heating_rate=  1.401905E-27
```

The changed CUDA multigroup smoke was compiled and run on an NVIDIA A10 with
`ngroup=9`:

```text
SNRT_CUDA_MULTIGROUP_OK relative_budget_error=  5.582516E-09
```

The production module graph was compiled with:

```text
make -C bin SNRT=1 USE_CUDA=1 ramses
```

It linked `bin/ramses_final3d` successfully at `2026-09-05 10:17:38 +0900`
with size `12161608` bytes. The executable symbol table contains the native
spectral loader and the linked `snrt_ramses_advance_level` entry. No RAMSES
evolution or production feedback run was launched.

`git diff --check` is clean for the bundle files. The source tree remains a
shared dirty development tree containing pre-existing user changes; no
unrelated files were reset or overwritten.

## Opus conditional-pass remediation

The bundle-end audit returned **CONDITIONAL PASS**. Its native-path findings
and the remediation are recorded in
[`claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_end_audit_2026-09-05.md`](claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_end_audit_2026-09-05.md).

- C1/F3 is closed by conservative native upper bounds: cross sections are no
  greater than `1e-17 cm^2`, and absorber-weighted excess energy cannot exceed
  the group upper edge minus the species threshold.
- C2/F4 is closed by the direct loader rejection executable and the nonzero
  payload checkpoint round-trip described above.
- C3 is closed as documentation: the native mean/excess energy residual and
  the missing secondary-ionization/recombination channels are explicitly
  quantified in `SNRT_NATIVE_GROUP_CONTRACT.md` and carried to later G3/G4
  science gates.
- F6/F7 were hardened in the same remediation by recording fraction semantics,
  blocking intrinsic fractions at resolved-domain runtime, and requiring a
  reference-control opt-in.

No production RAMSES evolution or physical-source approval was added by this
remediation.

Key bundle hashes:

```text
patch/lagRamses/snrt_spectral_contract.f90  dea72a4128d58fdcd1d3a590d4e9442c49b1c45e6c1e130687b3ecf2603514b9
patch/lagRamses/snrt_state.f90              591e7efd72c5a0595a7224331fa25d53c2f0ecf51af1ccc21938e47a7403cf88
patch/lagRamses/snrt_ramses_driver.f90      9482ee63938792a7755ec3f4269f65bcf4daeb7240d43d75b2ec23253bd5faf4
simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
                                            5825dae4d55d5f4448880c0ebd5b9b727ecfa97cbcf573b609b6e9bfa6e6ca94
simulation/snrt/tests/run_snrt_native_spectral_contract.sh
                                            986fc7fb3c72b891fa9ec5fab4245ba12d7a4044cd2d6161252dc6ef5cb6e07a
patch/lagRamses/snrt_spectral_contract_smoke.f90
                                            e009d518da7936c1e1f2e982caa2cc26571e845a4b3429732b96d1ead4e23f62
patch/lagRamses/snrt_spectral_contract_loader_smoke.f90
                                            6016b3efd79ef845697fd25bace496055b2412201d7e30488bf20c1ab7ff1d1b
patch/lagRamses/snrt_checkpoint_smoke.f90
                                            51b800f7d66459fdf3a7e8aa9b84e5b648f7c96b5c5c1176eb2d2e7e2ac6be97
simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md
                                            19cba11e8ffa12f3caf3ba3e7e2d89fee91ff9048b049e9922a779f093a8c8ad
```

## Boundaries and remaining gates

This evidence closes the native nine-group/source-contract wiring scope only.
It does not promote the parameterized AGN pilot, approve a stellar SED, add
stellar photon emission, activate He chemistry, dust opacity, radiation
pressure, IR re-emission, or connect SNRT state to RAMSES HDF5 restart. It also
does not resolve the 40--120 M☉ physical-yield seam. Any production run still
requires an approved production contract, CUDA, and the separate physical
source, thermochemistry, coupling, and science qualification gates.

The end-of-bundle audit prompt is the paired
`claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_end_audit_prompt_2026-09-05.md`.
The recorded initial audit is
`claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_end_audit_2026-09-05.md`;
the intermediate follow-up is recorded by Opus in its read-only plan, and the
final closure audit is
[`claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_closure_audit_2026-09-05.md`](claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_closure_audit_2026-09-05.md).
The final verdict is **PASS**; C1, C2, and C3 are closed for this bundle.
