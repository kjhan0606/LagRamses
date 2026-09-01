# Yield-table scientific and implementation audit

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`

## Auditor provenance

- **AGY**: the exact `agy` model alias was unavailable (`Unknown model`). The
  delivered report was produced by a fallback subagent instructed to act in
  the AGY role. It must not be represented as an authentic AGY model call.
- **Fable**: completed through the Claude CLI with `claude --print
  --model fable`; read-only audit, exit code 0, no edits or commits.
- **GPT-5.6-sol**: completed as an actual `gpt-5.6-sol` subagent review.

## Evidence reviewed

- Legacy table: `/gpfs/kjhan/Run_JWST/opt_run/yield_table.asc`.
  It contains 12,000 data rows in 60 blocks of 200 and only H/O/Fe. SHA256:
  `ba1099c5a4c3afe5e9ba28b3eb59d2e85fd3d40b7e7cb4ec30799eec00a5ac2e`.
  It is suitable as a legacy comparison asset, not as the new production
  yield table.
- New contract: 32 fields, 11 tracked elements, five release channels, and
  fail-closed validation gates.
- Current canonical fixture: 9 rows, with one mass point for each of wind,
  AGB, and SNII. Row-level closure tests pass, but the mass/metallicity/age
  grid is not production-complete; SNIa and PISN are absent and approvals and
  provenance sidecar are missing.

## AGY-role report

**Verdict: RED — publication use is blocked.**

- The schema and channel separation are a reasonable implementation
  foundation, but the actual scientific yield content is not yet selected,
  cited, or approved.
- The 11-element set is adequate as a reduced chemistry set for the intended
  LRD RT/dust use, but is not complete nucleosynthesis. Na, Al, Ar, Ni, Li, F,
  isotopes, and detailed s/r-process products are not represented.
- Wind, AGB, and SNII source choices, remnant/fallback treatment, explosion
  energy, and AGB release-time behavior remain unapproved.
- SNIa is not implemented. It requires a separate SSP convolution with a
  delay-time distribution and metallicity/progenitor-dependent event yields;
  it must not be treated as an ordinary 3--8 Msun IMF-integrated channel.
- PISN is correctly disabled for now, but a future implementation needs
  explicit Pop-III/core-mass/PPISN eligibility rather than a universal ZAMS
  mass-window switch.
- The IMF, remnant ledger, and `initial = living + remnant + returned`
  conservation proof must be made consistent across the table and runtime.

## Fable report

**Verdict: the architecture is strong, but the production gate must remain
closed because scientific yield content is missing.**

Fable additionally identified the following items for source-level
verification before any fix is claimed:

1. Possible age-unit mismatch: table documentation uses years while the
   runtime query appears to use Gyr (`stellar_ramses_runtime.f90:168`).
2. Possible interval convention mismatch in the source increment path,
   apparently using `[age, age+dt]` (`stellar_source_increment.f90:57--58`).
3. Possible missing `1d51` conversion in legacy kinetic-SN energy paths
   (`feedback.kjhan3.f90:1319,1987`).
4. Possible metadata inconsistency between the three legacy species and the
   11-species deposit path, including He indexing and the NVAR-17/delayed-
   cooling interaction.

These are audit findings, not confirmed fixes or accepted defects yet.

Fable also judged the current BPASS-based path to be a candidate rather than
   an approved source: it still has metallicity/age clamping and `fesc=1`
   assumptions, and lacks a publication-ready provenance/approval record.

## GPT-5.6-sol report

**Verdict: RED — the framework is testable, but the yield model is not yet
scientifically complete.**

- The legacy H/O/Fe table must remain comparison-only; the 9-row canonical
  fixture is a plumbing test, not a physical yield grid.
- A real production table needs a complete mass--metallicity--age grid,
  literature provenance, strict interpolation/extrapolation policy, and
  closure of returned mass, remnant mass, energy, and momentum.
- The runtime's hard-coded Kroupa IMF and the candidate BPASS binary IMF/range
  are not yet demonstrably normalized to the same population. This is
  especially important for very massive stars and any future PISN treatment.
- SNIa DTD normalization, metallicity dependence, event energy, and ejecta
  composition are missing. Binary evolution can materially change the
  massive-star and AGB yields, so a single-star table cannot silently stand in
  for a binary population.
- The 11 tracked elements should be explicitly described as reduced chemistry,
  not as a complete nucleosynthetic inventory.
- SED, dust, scattering/IR re-emission, and live RT/hydro coupling remain
  separate downstream tasks; the current yield table does not establish them.

## Consensus and implementation order

All three reviews converge on the same operational decision:

> Keep the legacy run/table as a comparison baseline. Do not call the new
> yield implementation publication-ready until the scientific source data,
> approvals, and conservation/units ledger are complete.

Priority order:

1. **P0** — choose and approve wind/AGB/SNII sources; generate the complete
   grid; attach citations and checksums; resolve units, IMF normalization,
   remnant accounting, and mass/energy closure; verify the four Fable code
   findings.
2. **P1** — implement and validate AGB delay handling and SNIa DTD/event
   yields; make stellar and AGN energy/momentum ledgers unambiguous.
3. **P2** — add approved SED/dust physics, including the intended absorption,
   scattering, and IR-re-emission assumptions.
4. **P3** — validate live coupling to RT and RAMSES hydro at production scale.

Relevant recent literature checks include [Jost et al. 2024](https://academic.oup.com/mnras/article/536/3/2135/7920781),
[Keegans et al. 2023](https://arxiv.org/abs/2306.12885),
[Pepe et al. 2025](https://arxiv.org/abs/2412.07845),
[Osborn et al. 2025](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/using-binary-population-synthesis-to-calculate-the-yields-of-low-and-intermediatemass-binary-populations-at-low-metallicity/9B60FA78A037F92C51F5ECE9869B5592),
and [Gabrielli et al. 2024](https://academic.oup.com/mnras/article/534/1/151/7746769).
