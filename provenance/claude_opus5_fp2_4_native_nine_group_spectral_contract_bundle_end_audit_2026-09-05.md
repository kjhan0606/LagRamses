# Claude Opus 5 end audit — F-P2.4 native nine-group SNRT contract

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Auditor: Claude Opus 5

Date: 2026-09-05

Mode: read-only; no repository edits, jobs, or network access.

## Verdict

**CONDITIONAL PASS**

The native nine-group wiring is real and correct. The audit found no array
ordering, dimension, unit, CUDA ABI, or double-counting defect in the delivered
Fortran/CUDA path, and no stale four-group remnant in `patch/lagRamses`. The
conditions below are mechanical closure items and do not require a redesign.

## Conditions for closing F-P2.4

1. **C1 / F3 — value upper bounds.** The spectral validator required an upper
   bound on cross sections and absorber-weighted excess energies. Positivity
   alone would admit a decimal/transcription error. The follow-up remediation
   bounds each excess energy by the top of its group interval measured from the
   species threshold, with numerical tolerance. The initial remediation used
   `1e-16 cm^2`; the follow-up audit found that too permissive and required the
   final `1e-17 cm^2` ceiling now present in the source.
2. **C2 / F4 — direct rejection and restart evidence.** The original smoke
   only executed the successful loader and a pure identity comparison. The
   follow-up adds executable tests for unset/missing/malformed/version,
   malformed identity, edge digest, unknown fraction semantics, candidate
   status, and reference opt-in paths. It also writes and reads a real
   version-4 checkpoint containing a small nine-group state payload and tests
   identity rejection before state mutation.
3. **C3 — native science limitations.** The native documentation must name and
   quantify the emission-mean versus absorber-weighted heating gap and the
   missing secondary-ionization/recombination channels. These are subsequent
   science gates, not defects in the nine-group wiring.

The above conditions were implemented in the post-audit remediation and are
recorded in the paired implementation evidence.

## High-severity science findings carried forward

### F1 — emission/heating energy residual

Emission converts group energy to photon number using the photon-number-
weighted mean `Ebar_g`. Native H heating uses the Verner-absorber-weighted
photoelectron excess. Because the state carries photon number only, the
per-absorbed-photon energy does not close exactly. For the reference-control
contract, the approximate emission energy versus `13.6 + H excess` is:

| group | interval (eV) | emission mean (eV) | heating energy (eV) | relative gap |
| --- | ---: | ---: | ---: | ---: |
| 5 | 13.6--24.59 | 17.66 | 16.45 | 7% |
| 6 | 24.59--54.42 | 34.38 | 30.42 | 12% |
| 7 | 54.42--500 | 106.64 | 68.87 | 35% |
| 8 | 500--2000 | 869.63 | 625.81 | 28% |
| 9 | 2000--10000 | 4023.59 | 2578.50 | 36% |

Both averages are individually defensible; restoring `mean - 13.6` is not an
acceptable correction. A later G3/G4 gate must add an absorbed-energy residual
or a second energy state/mean with an explicit conservation convention.

### F2 — missing native secondary ionization and recombination

Native H chemistry currently counts one primary H ionization per absorbed
photon. The newly added group 9 carries a reference-control H excess of
`2564.90 eV`; in mostly neutral gas a photoelectron at that energy can produce
many secondary ionizations. The native path currently deposits that channel as
heat and one primary ionization only. No recombination update exists in the
native SNRT modules. Secondary ionization, recombination, and thermochemical
closure remain later G3/G4 implementation gates.

## Other findings and dispositions

- **F5:** source digests are declared identity fields and are not hashes of the
  loaded namelist numeric block. The current evidence must not claim stronger
  content-addressing than that. A numeric-block digest can be considered in a
  later provenance hardening task.
- **F6:** fraction semantics was added to the native contract and checkpoint
  identity. `escaped` is required for resolved-domain runtime; `intrinsic`
  remains inspectable but is blocked until an upstream escape conversion is
  explicitly approved.
- **F7:** reference-control execution now requires
  `SNRT_ALLOW_REFERENCE_CONTROL=1`; production status does not use this
  opt-in.
- **F8:** group-loop ordering is currently inert because the all-group CUDA
  neutral cap and the Fortran available-neutral budget share the same bound,
  but the invariant is documented for a future transport/chemistry audit.
- **F9:** ghost rows copied from the transport device buffer are currently
  harmless because the caller resets the full work buffer and reads owned rows
  only. A future hardening task may memset or document this ABI assumption.
- **F10:** source photons use step-start attenuation and step-start optical
  depth across substeps in the native path. This is a later transport accuracy
  and convergence gate, not a nine-group wiring defect.

## Gate disposition

- Canonical ten boundaries/nine groups: **PASS**.
- `[2000,10000] eV` allocation, transport, source transaction, and CUDA ABI:
  **PASS**.
- Fail-closed status/identity/edge/value validation: **PASS after C1/C2
  remediation**.
- Source energy closure and removal of the old unexplained `0.5` factor:
  **PASS**; independently checked against the ledger photon rates.
- Checkpoint identity binding: **PASS after C2 round-trip evidence**.
- Native versus Python/reference physics equivalence: **not claimed**; H-only,
  no-secondary/no-recombination status and F1 are explicitly carried forward.

## Not approved by this audit

This audit does not approve the parameterized pilot AGN SED, a physical stellar
SED, a production hydro run, live RT+feedback coupling, dust/radiation
pressure, or publication science results. No RAMSES evolution was launched.
