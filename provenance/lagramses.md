# LRD/JWST lagRamses integration provenance

- Integrated repository: `git@github.com:kjhan0606/LagRamses.git`
- Integration baseline: `db9116018677e9ff37d85cdadc46d64ad528c70a`
- Historical migration pin: `d2271e87d837dd94b88a3609618521cbad2cfaab`
- Local development checkout at migration: `/home/kjhan/BACKUP/lagRamses`
- Default stellar-feedback policy: `channel_resolved`
- Historical reproduction policy: `legacy`

The LRD/JWST workflow and lagRamses source now share this repository. Each
submitted run must record the resolved commit, compiler, build switches, and
external yield-table manifest alongside its effective namelist.

## Sink-creation flag status

Commit `dcf5fb327bb380ff982743a6564586ce4a3e1bde` completes the active-patch
`create_sinks` wiring. Its parent `1e1d572` gates only the new-sink density
scan and `kjhan_make_sink` calls, while retaining merge, cloud, and Bondi
maintenance for existing sinks. Commit `dcf5fb3` reads the standard
`&SINK_PARAMS / create_sinks=...` field into the active `pm_parameters`
module and reports the resolved value at startup. Commit `db91160` enables
sink creation by default.

This source result does not certify older executables. The registered Phase 0
comparison binary was built from `65d0802-dirty`, before both fixes. A future
production run must build at or after `dcf5fb3`, include an explicit
`create_sinks` value in its effective namelist despite the later default
change, and record the startup report.

## Registered comparison run

The stopped Phase 0 RAMSES run on `lageunha` is registered as the external
transitional/development-feedback comparison baseline (its exact runtime mode
is not serialized in the run):
[`provenance/legacy_feedback_baseline.md`](legacy_feedback_baseline.md).
