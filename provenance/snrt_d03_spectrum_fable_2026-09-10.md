# Fable plan disposition — D03 spectral runtime

Driver summary, not a verbatim transcript. Read-only CLI fable plan review
returned REVISE (small specification items), then proceed. Q-GOAL yes;
Q-LEAN lean. No files/jobs/subagents/messages were produced by the reviewer.

Adopt: old two spectral modes retain no-live-dust/no-CHIMES restrictions;
new mode requires DUST_LIVE and rejects SNRT_CHIMES, non-D03 optics, Fe,
PAH, sublimation and relative motion. New nodal artifact/manifest does not
rewrite the v1 monochromatic/IR include or identity. Bind exact grid energies
and content SHA; dust endpoints use actual energies, not H/He nextafter.
Pass four columns and normalized angular weights; uncapped dust count/energy
uses existing outputs. Actual captured energy goes through the receiver's
deposited-spectrum argument. Kind3/native10/HDF offset and nodal SHA bind
restart. Skip the old primary scattering operator and supply zero legacy
group dust tau to the new spectral operator.

Trim adopted: no generator changes because no namelist semantics/key changes
and those scripts do not expose SNRT_SPECTRAL_MODEL. Update NATIVE_RUNTIME.

Driver decisions: retain the plan's sequential node absorption then elastic
delta-isotropic scattering (not the reviewer's suggested combined extinction
partition); this is the declared existing operator-splitting family and
conserves actual N/E, but does not model multiple scatter/absorb histories
inside one step. Rejected H/He packets retain E/N and are NOT transferred to
dust again, unlike the old fixed operator. This difference is explicit.

The measured source-resolution amendment selects128 instead of64 nodes for
the new mode (the reviewer referenced the original64-node text). No second
review is needed for this bounded numerical-resolution repair. The old64
models stay unchanged; zero-dust parity is with the same128-node kernel,
not a false bitwise128/64 equivalence. Monochromatic endpoint packets and
native node sums supply the extinction reference. Offline dense Mie checks
are evidence, not another native model/table or additional gate.

All necessary specification changes are adopted. Proceed under the existing
operator preapproval; driver performs the bounded end evaluation.
