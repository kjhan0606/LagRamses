# AGY F-P1 high-mass internal-gate re-audit, round 1

Date: 2026-09-03
Model: AGY / `gemini-3.8-flash-high`
Mode: read-only inspection of `/gpfs/kjhan/LRD_JWST`

## Verdict

- Top level: **CONDITIONAL PASS**.
- Internal F-P1H-A--E controls: **PASS**.
- Physical 40--120 M☉ gap: **BLOCK**.
- F-P1H-F, production, and publication: **BLOCK**.

AGY marked every first-audit remediation verified fixed. It confirmed computed
high-mass evidence and mutations, the 84-field node validator, actual-file
converter/mapping hashes, derived package gates, aggregate runner wiring,
transactional namelist identity, both driver guards, and corrected wording.

It independently reproduced 18 W18/N20 outcomes, 12 failed outcomes, six
terminal records, maximum relative residual `0.003301399145327014`, four common
wind nodes, the `1.0000010999533515e-4 M☉` 100 M☉ Mg difference, eight of eight
K-40 duplicate records, and 12 radioactive warnings.

The only new finding was low-severity environment coupling: the SNRT virtual
environment resolves Python through `/home/kjhan/miniconda3`, so a sandbox that
unmounts `/home` cannot run the shell runner. Normal project execution passed.

AGY explicitly confirmed the narrow statement that internal controls fail
closed while physical nodes, runtime consumer, linked build, physical source
admission, production, and publication remain blocked.
