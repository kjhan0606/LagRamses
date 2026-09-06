# Run configuration wizard

From the repository root, run `python3 mkrun.py --mode gui`. The older
`--gui` spelling remains a compatibility alias. With no arguments,
`mkrun.py` retains its terminal wizard. `--help` works without Tkinter or a display.
GUI mode requires Python's Tkinter module and a graphical desktop or forwarded
display; an unavailable dependency/display produces an actionable error and exit
status 2. No packages are downloaded automatically.

Next/Back walks the same prompts and defaults as the terminal wizard, including
DMO/hydro, dark matter and gravity sectors, cosmology, AMR, zoom, IC pipeline,
output epochs and hydro settings. The base and maximum AMR levels are shown in
one side-by-side form, so they can be adjusted together and validated as a pair.
Returning to a question retains its answer; changing it rebuilds subsequent questions.
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
