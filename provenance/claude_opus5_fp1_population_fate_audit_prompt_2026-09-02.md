# Claude Opus 5 F-P1 audit request

Work read-only in `/gpfs/kjhan/LRD_JWST`; do not edit files, commit, launch a
simulation, or expand into generic HDF5/AMR/ksection/CPU-box work.  Use the
`claude-opus-5` model's own scientific and code judgment.

Audit only F-P1 stellar-population and fate semantics.  Start from
`provenance/fp1_population_fate_contract_2026-09-02.md`, then independently
inspect the production files under `patch/lagRamses`, their native mirrors,
the tests named in that report, `simulation/snrt/config/g2_physics_contract_v1.json`,
`simulation/snrt/config/stellar_source_identity_v1.json`, and the relevant
candidate-source/fate records.  You may run read-only tests and builds, but do
not launch RAMSES time integration.

Judge these criteria separately:

1. Source-basis units are unambiguous and an SSP-normalized table cannot be
   convolved with the IMF a second time.
2. IMF mass support and normalization are mathematically correct, explicit,
   configuration-driven, and independently tested.
3. Population-model and binary-fraction parsing is transactional and
   inconsistent or unsupported binary physics fails closed.
4. Enabled channel windows lie inside IMF support, phase ownership is
   non-overlapping by contract, and non-terminal channels cannot contribute a
   second remnant.
5. The population ledger correctly closes initial, living, remnant, returned,
   tracked, and untracked mass and is linked into the production source tree.
6. SNIa and PISN cannot accidentally use the ordinary IMF channel path.
7. The production build/source identity and test claims in the F-P1 report are
   reproducible.
8. The report does not overclaim physical source, IMF, lifetime, binary, or
   fate-model approval where the staged evidence is incomplete.

Return a concise but detailed Markdown audit with an overall verdict of PASS,
CONDITIONAL PASS, or BLOCK.  Distinguish code-contract defects that should be
fixed now from missing physical model/source decisions that keep the scientific
F-P1 gate open.  Put RT, dust, AGN, SNIa DTD implementation, PISN implementation,
and generic infrastructure suggestions explicitly out of scope unless they
directly invalidate an F-P1 claim.
