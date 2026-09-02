# P0.1 source-of-truth and build-parity gate

This gate prevents the native/phase0 contract PASS from being presented as a
PASS for the compiled `lagRamses` executable.

The production source is the `patch/lagRamses` tree selected by
`bin/Makefile`. The G1 runner intentionally builds a separate RAMSES-
independent mirror under `simulation/snrt/native/phase0`, which is retained as
a differential oracle. The production-linked harness targets the actual
`bin/Makefile`, which consumes the source-order sidecar, and records a forced
rebuild, binary linkage symbols, startup smoke, and source/config/tool hashes.
The source-order sidecar is checked against the effective make object order.

Run the read-only gate with:

```text
python3 simulation/snrt/tools/validate_stellar_source_parity.py
```

Use `--require-pass` in a production qualification path. It currently returns
`STELLAR_SOURCE_PARITY_BLOCKED blocked=production_linked_build_evidence` because
the production-linked harness has not yet produced fresh build evidence. The
14 differing native/production hashes are differential diagnostics, not a
P0.1 failure: the native mirror is explicitly not the production source of
truth. Existing `build/g1_native` objects are not evidence and are never used
by this validator.

For transparency, the 14 differing entries are 11 byte-different modules and
3 modules absent from the native mirror (`stellar_cell_deposition.f90`,
`stellar_ramses_bridge.f90`, and `stellar_ramses_mapped_bridge.f90`);
`stellar_yield_backend.f90` is identical. These counts are diagnostics only.

The build evidence is accepted only when the log has exactly one compile line
for each of the 16 required objects, the observed `-D` values match the
contract, the required active-branch flags are present, embedded yields are
absent, the expected binary exports the configured stellar runtime/bridge
symbols, and the no-argument startup smoke returns the expected usage line
with RAMSES's intentional `clean_stop` exit code 0.
The validator and evidence recorder hashes are also pinned, so changing the
gate tooling invalidates prior evidence. The build and smoke logs both carry
the exact production-binary SHA-256; the recorder derives the forced-rebuild
claim from the logged `-B` command instead of asserting it unconditionally.
Evidence remains valid across a later descendant commit: source, build-input,
tool, binary, and log hashes remain exact, while the recorded build `HEAD` is
required only to be an ancestor of the current `HEAD`.

The gate is invoked before the native G1 test by
`tests/run_g1_native_contract.sh`, so a future native contract PASS cannot
bypass compiled-tree parity. The production-linked path is
`tests/run_p0_production_linked_contract.sh`; it is opt-in and requires
`P0_BUILD=1`. It runs `make -C "$ROOT/bin" -B ramses`, records the binary,
Makefile, source, and compile-contract evidence, and then reruns the gate.
For differential-only work, set `P0_DIAGNOSTIC=1` on the G1 runner; it emits
`G1_NATIVE_DIAGNOSTIC_ONLY` and never claims P0 closure.

The native/production shared-file profile is a bounded differential diagnostic:
the configured shared list must partition into byte-identical, byte-different,
or native-absent entries, and every production entry must resolve. The exact
baseline counts are reported for drift diagnostics but are not a gate, so a
native mirror may legitimately converge toward production. The build writes
its generated Makefile source to the ignored
`bin/write_makefile.generated.f90` and does not delete the tracked
`bin/write_makefile.f90`.

The mirror-only modules `stellar_native_units.f90`,
`stellar_progress_contract.f90`, and `stellar_population_ledger.f90` are
explicitly declared as native differential-test scaffolding. They are not
production-linked sources and are not silently treated as part of the shared
production contract.

The current implementation record is
`provenance/p0_source_parity_gate_2026-09-02.md`; its expected status is
`STELLAR_SOURCE_PARITY_BLOCKED` until the opt-in production-linked build has
completed successfully.
