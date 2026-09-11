# Selective plan review and driver disposition — 2026-09-11

Read-only Fable review of four approved comparison plans. Q-GOAL and Q-LEAN
were requested first. The transcript below is preserved, not adopted as an
independent source of physical facts.

Driver disposition: retain the named approximations and unchanged defaults.
Apply the Fe surface ionization-energy accounting clarification and the
mixed-grid wind-energy scaling/explicit Ia exclusion. The main-owned runtime
SNIa integration, v5 wrapper, source identity and pair guards have since been
edited; build and execution evidence are recorded separately when completed.
In particular, the review's claim that every UV Fe photoelectron is below
9 eV must **not** be used as an all-heat shortcut: negative charge and
photodetachment can invalidate that bound. Preserve the distributed-electron
secondary treatment and validate against actual charge-dependent energies.
The review is not a completion certificate for unfinished implementations.

## Original review

Write is disabled in this session, so the review is delivered here in full.

## Verdicts

| Plan | Verdict |
|---|---|
| Fe photon comparisons | CONDITIONAL |
| PARSEC mixed low-mass | CONDITIONAL |
| SNIa frozen He hybrid | PASS |
| PARSEC bistability | PASS |

## 1. Fe photon comparisons — CONDITIONAL

**Q-GOAL.** Aligned. The UV cycle adds real live coupling: Fe photoelectric gas heating through escaping-electron energy, H+/C+ surface neutralization, a single photon debit, and Fe sensible heating. The thermal limit adds a bounded hard-photon Fe heating comparison. Neither claims keV/Auger physics.

**Q-LEAN.** Clean. Charge populations are worker-local scratch, no hydro carriers, one kernel plus the existing IR transaction. The model gates and the 12-entry identity already exist in `dust_mass_physics.f90:394-396` and `:521-538`, with the driver guard at `snrt_ramses_driver.f90:334`.

**Verified, not fabricated.**
- The fixed groups already reach 10000 eV on the pinned Werner09/Henke plus size-Drude dielectric. The keV-group absorption is tabulated atomic data, not a Drude extrapolation.
- HD2017 is Hensley & Draine 2017 (arXiv:1611.08607). The execution record confirms it supplies metallic-Fe UV yields, photodetachment, and collisional charging with a 13.6 eV cutoff.
- The FS2010 grid is pinned at 10 eV and above, but the CHIMES photo adapter routes sample-less terms to all-heat (`snrt_chimes_photo.cpp:91-101`). Photoelectrons from ≤13.6 eV photons on Fe cannot exceed about 9 eV, so feeding distributions is harmless and consistent.

**One indispensable change, documentation only.** State where the ionization potential released at the surface goes. "Full accommodation, no recombination photon" implies the ion's ionization energy minus the effective work function is credited to the Fe sensible reservoir and debited from the CHIMES binding ledger. The closed-cycle identity the smoke must check is then explicit: photon energy plus ion kinetic energy plus ionization potential equals escaping-electron kinetic energy plus the 2 kT_dust outgoing atom plus Fe sensible heat. As written, "accounted separately" leaves the destination undefined, which is the one place a silent loss or double booking can enter.

**Authorized limitations, no change.** Fe-only photoelectric heating with no C/silicate or PAH counterpart makes the gas-heating delta Fe-specific. The ≤300 K Fe restriction rejects rather than clips cells where full retention pushes Fe above 300 K, which is existing behaviour in `snrt_dust_live.f90:519`. UV-only removal of Fe from transport opacity is a first-order splitting, conservative because there is a single debit.

## 2. PARSEC mixed low-mass — CONDITIONAL

**Q-GOAL.** Aligned. It extends live radiation and material sources to 2–14 Msun on the same IMF cells: intermediate-mass FUV/optical output over 10–1000 Myr and lifetime-timed AGB ejecta feeding dust condensation. Proxies are explicit and no matched-SSP claim is made.

**Q-LEAN.** Clean. Node arrays only, no per-particle state. The v5 cap of 1024 and v4 cap of 512 are already enforced in the loader and audit.

**Verified consistent.** The v5 admission and ownership logic in `stellar_yield_audit.f90:366-473` enforces wind plus terminal plus remnant equals initial mass per node, remnant only on the owner channel, kind 3 for 7<M<9, and fate 0 only for AGB. The V2 loader requires the last knot to equal the radiation stop and the stop not to exceed death. Nearest-node weights use full 0.08–600 normalization with a mass-closure check, and the 40 Msun seam is retained in the node selector. Non-owner channels return zero, the Ia channel returns zero, and the CO inventory excludes kind≠0. The v5 dust path condenses channels 1–3 only, which matches the legacy path, so there is no PISN-dust regression.

**Two indispensable changes, documentation only.**
1. State that v5 is mutually exclusive with any SNIa model. Native admission already rejects an enabled Ia channel for v5 in both `stellar_yield_audit.f90:109` and `snrt_parsec_source.f90:153`. Bundles 3 and 4 therefore cannot be combined, and the listed "strict CO inventory exclusion" check is a table-level check, not a live Ia pathway.
2. State that wind kinetic energy is computed from the scalar-reduced wind history, not the unscaled PARSEC history. Reduced mass with unreduced energy would exceed half v-squared per unit mass. The audit only checks energy monotonicity and cannot catch this.

**Already in the plan, still unwired, required.** Runtime identity binds fate only for v4, the pair-channel guards in the enrichment driver reject version≠4, and the SED wrapper reader stops at version 4. All are Main-owned and listed.

## 3. SNIa frozen He hybrid — PASS

**Q-GOAL.** Aligned as a named comparison: one delta-function prompt channel at 86.14 Myr with unchanged N100 ejecta and energy, replacing empirical counts under existing effective accounting. Narrow but coherent.

**Q-LEAN.** Clean. One event row in the ordered sidecar, no new namelist or environment variable, and an empty identity helper for the empirical model.

**Verified.** The contract module implements the plan exactly: effective-only, expectation-only, baked binary fraction, unity IMF factor, weight equal to events per initial Msun, half-open CDF with no age tolerance, W17 stable/central check at both the seed mass and 1.35 Msun, and N100 identity binding with zero remnant. The arithmetic checks out: 85.9619 Myr plus (1.40046 − 1.03448)/2e-6 yr gives 86.1449 Myr, matching the fixture. The converter checks a detached seed, an event within the retained history, and the donor budget, and adds no COSMIC WD accretion. The runtime currently does not pass birth Z, so the frozen model fails closed with error 78 until Main wires it. That is safe.

**Required, already listed.** The three Main-owned steps are indispensable: pass birth metallicity in both interval calls, append the event identity to the phase-0 SNIa identity (event age and Z are otherwise unbound on restart), and call the event-source validator in the contract loader.

**Authorized limitations.** Only birth Z=0.01 is admitted at roundoff, so any SSP with another birth Z aborts the run. The comparison is usable only in fixed-Z idealized setups. The weight carries COSMIC's illustrative IMF and binary fraction, not the runtime IMF. Both are stated, not fabricated.

## 4. PARSEC bistability — PASS

**Q-GOAL.** Aligned. Velocity-only sensitivity of wind kinetic energy for H-rich 10–27.5 kK states, with mass, yields, remnants, fates, and terminal ages unchanged.

**Q-LEAN.** Clean. One helper reusing the existing phase law, one model ID added to the native accepted list, no API, field, namelist, or Makefile change.

**Verified.** The helper scales the H-rich branch by 1.3/2.6 below 22.5 kK and ramps linearly to 27.5 kK, with the Z domain and precision-grid rejection as stated. Vink, de Koter & Lamers 2001 does give 2.6 hot-side and 1.3 cool-side with a 22.5–27.5 kK bracket, and the ramp is correctly labelled as a comparison closure. The energy quadrature in the plan is exactly the code. The 51-character model ID is accepted in `snrt_parsec_source.f90:50-53` and fits the 64-character field. The test file and the private evidence directory with the named logs exist.

**No conditions.** The note honestly states that native SED-loader admission of the new ID is not yet tested by the retained executable.

## Not done

The plan file could not be written because the Write tool is disabled in this session. The review above is the complete deliverable. No files were edited and no jobs were launched.
