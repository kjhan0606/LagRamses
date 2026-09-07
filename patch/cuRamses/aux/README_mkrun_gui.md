# Run configuration wizard

From the repository root, run `python3 mkrun.py --mode gui`. The older
`--gui` spelling remains a compatibility alias. With no arguments,
`mkrun.py` retains its terminal wizard. `--help` works without Tkinter or a display.
GUI mode requires Python's Tkinter module and a graphical desktop or forwarded
display; an unavailable dependency/display produces an actionable error and exit
status 2. No packages are downloaded automatically.

## Fixed RT / feedback / dust comparison

In either `python3 mkrun.py` or `python3 mkrun.py --mode gui`, select
**RT/feedback/dust comparison** in the **Run mode** stage and explicitly accept
the reference-only model. The default remains cosmological DMO; ordinary hydro
generation is unchanged. The comparison skips cosmology/IC/advanced editing and
uses the existing non-cosmological four-step profile, not a new physical model.

Choose a **new, nonexistent output directory** (type a new child directory if
the directory chooser selects an existing parent). The preview contains:

- `<name>.nml`: full existing profile including dust IC, HDF5 I/O,
  `create_sinks=.false.`, Kroupa/binary feedback and effective SSP SNIa.
- `ic_sink`, `<name>.history.nml`, `yields.dat`: copied local inputs; the
  namelist and environment point at the destination copies.
- `<name>.env.sh`: explicit BPASS, DL01, AGN/secondary/SNIa contracts and
  single-rank/OpenMP settings. Sourcing it exports settings but launches nothing.
- `README.txt`: build/dependency requirements, limitations, output budget and
  a manual launch command requiring a separate run/storage review.

The current local yield exports and verified executable must exist under the
repository's preserved `.agb-physical.4LAOTJ` and `.bpass-native.v0ZwR6`
directories. Missing files produce an error before any output is written; a
GitHub clone alone is insufficient. No source is synthesized/downloaded, and
no GPU, cosmological or production qualification is inferred. See the
[closeout handover](../../../provenance/rt_feedback_dust_comparison_closeout_2026-09-07.md).

No main RAMSES namelist field was added: SNRT and SNIa auxiliary contracts are
connected by the environment file, not invented entries in `RUN_PARAMS`.

## Shared wizard and save behavior

Next/Back navigates whole configuration stages: run files, DMO/hydro, dark
matter, gravity, cosmology, AMR, zoom, IC pipeline, output epochs, hydro and
advanced settings. Each stage shows its fields together in two columns.
Model choices appear as side-by-side radio buttons and reveal related fields
on the same page. Yes/No buttons also appear together. Long pages scroll
while navigation stays visible.

Use Update fields after changing a value that controls other fields.
Next validates the complete stage; newly exposed fields are shown for review
before proceeding. Back restores the whole previous form. Changing a branch
discards fields that no longer apply. Generation still uses the terminal
wizard's prompts, defaults and validation.
Advanced parameters use a JSON object
with names from the existing namelist database, typed numbers/booleans and quoted
Fortran array lists. As in the terminal wizard, the subsequent generation stage
sets output epochs, IC paths and the selected cooling defaults.

The wizard calls `mkrun.generate_run(ui, write_text)` with an in-memory text sink.
There are no filesystem writes during collection, browsing, validation or
preview. Each generated namelist/CAMB/IC configuration gets a read-only preview
tab. Validation errors disable Save; warnings remain visible. Save lists every
destination and explicitly identifies existing files requiring overwrite
confirmation. A changed destination, symlink, non-file or file larger than 4 MiB
is rejected. Writes are atomic per replaced file, not transactional across the
whole bundle; an I/O failure reports any files already saved.

This generates setup files only: it does not execute CAMB, IC generators,
simulations, schedulers or shell commands. The inherited output schedule and
physics defaults still require a scientific/storage audit before a real run.
The current parameter database does not support every sector offered by the
restored runner (for example FDM when `m_axion` is absent); those selections fail
with a clear message rather than introducing replacement defaults.

Run bounded tests from the root:

```sh
python3 -B -m unittest discover -s patch/cuRamses/aux -p test_ramses_run_gui.py -v
```

Tests exercise all IC branches, CLI/GUI byte equality, selected pre-refactor
comparisons, invalid input, absent Tk/display, canceled saves, overwrite and
concurrent destination changes. The baseline comparison requires the original
`b1d489633822c4ecca2cd9c68cc5b592b4ec25f6:mkrun.py` in Git history.
An additional real-widget test covers Next/Back, advanced editing, read-only
preview and confirmed saving when a display is available; it skips on headless hosts.
