# LC18 prompt elemental decay: bounded implementation and driver evaluation

Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
This closes the immediate LC18 Ni/Co-to-Fe source defect as an explicit
comparison, NOT all medium group4 or live radioactive transport.

## One plan review and driver disposition

Fable reviewed [the original live-isotope proposal](stellar_radioactive_inventory_plan_2026-09-10.md)
read-only using `claude -p --model fable --permission-mode plan`, with
Read/Grep/Glob, no --bare, jobs, edits, subagents or external contacts. Exit0.
This section is a driver summary, not a verbatim transcript.

Q-GOAL: correcting missing Ni56/Co56 daughter Fe directly improves native
elemental feedback. Q-LEAN: five new advected carriers, Bateman source-time
convolution and new MPI/HDF5 plumbing are excessive for that correction.
Verdict PROCEED with corrections: source-side projection, unchanged default,
reuse the LC18 native fixture and combined AGB/Ia source integration.
Driver adopted this reduced implementation, without another approval gate.

The review suggested an explicit1Myr projection. Driver instead uses the
prompt cutoff below: the pinned ICRP matrix has Fe60 half-life1.5Myr whereas
the pinned NUBASE2020 evaluation has2.62Myr. A finite-horizon matrix call
would use the old rate. This correction needs no full decay-matrix rebuild.
The claim that every galaxy timestep exceeds a year is NOT assumed.

## Explicit model

`prompt_t12_le_100yr_baryonic_v1` fully cascades states with half-life<=100yr
at stellar release, stopping at stable/longer-lived daughters. This is a
coarse-grained source abundance prescription, NOT a fixed-time solution,
resolved young-SN composition or live ISM isotope evolution. Long-lived
Al26 and Fe60 remain in their parent elements; there are no isotope tracers.
Ni56/Co56 become Fe, Co60 becomes Ni, H3 becomes He. Already-decayed KL16,
Fishlock and PARSEC inputs are not processed again.

Ground-state half-lives use [NUBASE2020](https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt),
the existing pinned file SHA256
`1585a5eea86c5e17e90307c7e6e786d060049c4039e392a261ff6db977df9859`.
Branches/metastables use existing pinned radioactivedecay0.6.1 ICRP107 data;
22 missing fast ground states use the existing NUBASE100%-beta classification.
Four missing ultra-long-lived parents stay retained. Laboratory EC rates
are not an ionization-dependent plasma calculation. No Al26 isomer split
is fabricated from LC18's single yield label.

All333 source labels close in342 traversed states. Nine pinned branch sums
deviate from unity by at most5e-6; normalize their relative probabilities
and record every original sum, rejecting deviations>1e-5. Use A*m_u baryonic
mass, not atomic rest-mass, and account for emitted alpha baryons as He4
if a short alpha branch is encountered. In this actual source/cutoff,
long-lived alpha parents are retained; the native test does not claim an
active alpha-yield signal. No gamma/lepton heat, CR energy or SNRT photons
are added. No second gas mass or explosion energy is created.

Wind is projected separately from `(table8-table9)` terminal ejecta at
13--25Msun; table8 is wind-only at30Msun and above, as before. The total
source mass, remnant, release age/shape and energy are unchanged. Only the
eleven-element partition changes; generic metals are still derived from
the total minus H/He, not from an assumption that Fe is the entire metal.
Net LC18 columns remain unavailable diagnostics. The historical review-only
decay contract and CLI remain unchanged; new export requires explicit choice.

## Usage and native connection

Add `--decay prompt_t12_le_100yr_baryonic_v1` to
`simulation/snrt/tools/build_lc18_native_wind.py`, or
`--massive-decay prompt_t12_le_100yr_baryonic_v1` to
`build_kl16_lc18_native.py`. Default `as_tabulated_no_decay` is unchanged.
Select the produced table/history through existing `PHASE0_YIELD_TABLE`
and `high_mass_history_path`. The history/source identities change with the
physical projection. No runtime Python, RAMSES NML key, carrier, VPATH,
namelist-generator or checkpoint-layout change is required. Existing native
source loading, IMF/Z interpolation, deposition and identity checks apply.

## Measured evidence: PASS for this bounded correction

Retained work directory `.lc18-prompt.smztor/`:

- `physical-test.log`:333 positive baryon-closing isotope maps; Ni/Co chain,
  H3 conversion and long-lived retention; unknown model rejection. All36M/Z
  nodes at each rotation0/150/300 preserve mass/age/energy/net columns;
  same423/464/451 rows, no added age sampling. Pre-edit no-decay combined
  yield and history files remain byte-identical.
- `run_native.sh`, `gnu/{build,baseline,prompt}.log`: existing combined
  Fortran fixture rebuilt from current sources with GNU bounds checks and
  invalid/zero/overflow traps. All36 high-mass nodes,16 actual CCSN endpoints,
  Table5 wind timing,73 AGB nodes, native age/Z/IMF split-additivity and six
  metallicity full effective-Ia DTD mass closures pass. Raw material/energy/
  remnant/AGB equality to baseline is checked in Fortran.
- Solar-Z20Msun terminal Fe changes from0.012714845500128406 to
  0.084632412883759375Msun. This is a single actual node, not a claimed
  universal factor across the mass/Z grid.
- The existing strict40Myr DTD/selected WD-supply incompatibility remains
  correctly rejected; prompt elemental decay does not repair microscopic Ia.
- `git diff --check` passed. No additional MPI/restart campaign: this is an
  offline input change consumed by the already-tested native source path.
  No new simulation raw output/checkpoint was generated; nothing to delete.

Final combined artifacts in `prompt-final/`:
`yields.dat` SHA256 `f2ddbce660c2ec4d40272fec6032b6e93ebeac16411aff2a37c8c2f6ccac4efe`;
`history.nml` SHA256 `2e2d484032b6af552ff9254becc273009551db1bfd1f75fcf3c55c2214533bc1`.
Projection SHA256 `1af7b9ae18348f06166e53acc8a128879db3dd0318bf6b40662b6399b978dc7e`.
The earlier `prompt/` candidate includes local library installation paths in
its identity; use `prompt-final/`, whose nuclear identity binds bytes only.

Live long-lived decay and dust-lattice transmutation remain unimplemented;
they were not silently replaced by a claimed full isotope solution.
