# P0 Active-Build Output Audit

Status: Draft v0.1  
Date: 2026-08-28

## 1. Inspected active patch files

The audit used only the active patch tree:

- `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/output_amr.kjhan.f90`
- `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/backup_hdf5.f90`
- `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/read_params.jaehyun.f90`

No source file was modified.

## 2. Confirmed facts

1. `outformat='hdf5'` selects the HDF5 checkpoint path.
2. The normal output path writes `hydro`, `part`, `sink`, and, when RT is
   compiled and enabled, a separate `rt` checkpoint component.
3. The HDF5 writer exports AMR topology, one raw `uold` dataset for every
   hydro variable, particle catalogues, and a structured sink catalogue.
4. The sink HDF5 catalogue includes `dMsmbh`, `dMBH_coarse`, `dMEd_coarse`,
   `eps_sink`, `Esave`, spin, position, velocity, and mass.
5. Particle output can include type, birth epoch, and metallicity, allowing a
   stellar-source catalogue after the active particle-type convention is read.
6. `nvar` is a compile-time hydro-variable count. The examined runtime reader
   validates only its minimum size; it does not provide a self-describing map
   from `uold_i` to passive-scalar meanings.

## 3. Evidence not yet available

No complete HDF5 output, `info_*.txt`, `namelist.txt`, or output metadata file
exists under the current Paper-III workspace. Therefore this audit does not
assert either the presence or absence of the following fields in a real HR5
checkpoint:

- RT photon density and flux by group;
- H, He, and H2 chemistry fractions;
- passive-scalar ordering and metallicity slot;
- dust abundance or dust-to-metal prescription;
- exact stellar particle-type and age convention; and
- the physical interval represented by the sink accretion accumulators.

## 4. Required snapshot package for converter development

The first converter test requires one *complete* selected output directory,
not an isolated hydro file. The package must contain:

1. `COMPLETE` marker;
2. HDF5 checkpoint payload or every matching binary component;
3. `info_*.txt` and `namelist.txt`;
4. `compilation.txt` identifying the active source revision; and
5. the sink time-series or a source ledger covering the output interval.

The converter will reject an output lacking a completion marker or one whose
metadata revision cannot be identified.

## 5. Immediate P0 consequence

The standalone benchmark solver and TPU memory work can continue without a
production snapshot. Rasterizer implementation, source SED normalization, and
the first M1--S_N scientific comparison cannot begin until the snapshot package
above is provided or its filesystem location is identified.
