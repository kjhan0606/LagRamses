# Fable plan review: driver record

This is a driver summary, not a verbatim transcript. Claude CLI model `fable`
was invoked read-only with Read/Grep/Glob, permission-mode plan, without
`--bare`, against `snrt_secondary_spectrum_plan_2026-09-10.md`.

Verdict: approve with bounded corrections. Q-GOAL: necessary runtime physics,
closing the nonlinear FS2010 mean-energy approximation. Q-LEAN: eight scratch
doubles, no new persisted state/tables, and existing tests are lean.

Adopted findings:

- Latch the actual model name; distinguish native binary versions as well
  as HDF5 identity. Cross-model rejection is deliberate reproducibility
  policy, not an assertion that the old state has an incompatible shape.
- MPI agreement must compare model indices, not a band-enabled Boolean.
- Integrate over direction-summed accepted nodes using each species' cap;
  do not call interpolation once per direction.
- Excitation is in the local chemistry result but has no transported driver
  carrier. Here it represents escaping line energy, not stored radiation.
- Sum secondary demands, reserve primaries first and cap against remaining
  targets. The implementation extends the existing routine with explicit
  optional integrated moments rather than adding the suggested sibling:
  validation/recombination/state-publication code is shared without copying.
- Double accepted counts match FP64 absorption energy. FP32 public counts
  remain tolerance-checked diagnostics. Step-start xi is consistent with the
  existing nonlinear trial scheme. The callback only reads loaded tables.
- Test the mean-energy counterexample above the low-electron-energy band;
  a band where every electron is below10eV would be identically all heat.

The reviewer also confirmed cost and transactional failure propagation are
acceptable. No new completion criteria for unrelated groups were adopted.
Driver end evaluation will use the actual native/live evidence, not this
plan verdict as proof of implementation.
