# cuRAMSES Auxiliary Python Scripts

A collection of auxiliary scripts for the cuRAMSES project, organized into five categories: post-processing, generators, analysis plots, visualization, and verification tests.

---

## 1. Post-Processing

### `compute_pk.py`
Post-processing matter power spectrum P(k) and two-point correlation function ξ(r) from RAMSES binary output data. Deposits all matter (DM particles, star particles, gas AMR leaf cells, sinks) onto a uniform grid via NGP, then computes P(k) via FFT and ξ(r) via Hankel transform with Gaussian window.
```bash
python3 compute_pk.py /path/to/output_00005 --ngrid 256 --nkbin 50
python3 compute_pk.py /path/to/output_00005 --no-gas         # Particles only
python3 compute_pk.py /path/to/output_00005 --no-particles   # Gas only
```
- Dependencies: numpy
- Output: `pk_00005.dat` (5 columns: k [h/Mpc], P(k) [(Mpc/h)³], r [Mpc/h], ξ(r), Nmodes)
- Reads RAMSES binary files directly (AMR + hydro + particles)
- Supports hilbert, ksection, bisection orderings

### `plot_density.py`
2D gas density projection from RAMSES output with cell-size-aware AMR rendering.
```bash
python3 plot_density.py /path/to/output_00005
```
- Dependencies: numpy, matplotlib

---

## 2. Generators (Table & Configuration)

### `generate_multi_z_table.py`
Generates a multi-redshift Grackle cooling table used by the Eunha cooling module.
```bash
python3 generate_multi_z_table.py --output grackle_multi_z.bin
```
- Dependencies: numpy, struct, subprocess
- Output: `grackle_multi_z.bin` (98 MB, 10 redshift snapshots)

### `generate_de_table.py`
Generates a dark energy transfer function ratio table R_DE(k,a) using the CAMB Weyl potential method.
```bash
python3 generate_de_table.py --w0=-0.8 --wa=0 --cs2=0.0 -o de_table_wcdm.dat
```
- Dependencies: numpy, camb

### `generate_neutrino_table.py`
Generates a neutrino transfer function ratio table R_nu(k,a) = T_nu/T_cb using CAMB.
```bash
python3 generate_neutrino_table.py neutrino_table.dat
```
- Dependencies: numpy, camb

### `de_linear_test_setup.py`
Batch-generates CAMB transfer functions, MUSIC IC configs, RAMSES namelists, and SLURM scripts for four DE models (LCDM, wCDM, CPL-thaw, CPL-phantom).
```bash
cd test_ksection/run_de_pk_test/
python3 ../../misc/de_linear_test_setup.py
```
- Dependencies: numpy, camb
- Must be run from the target directory

### `ramses_nml_generator.py`
Interactive and CLI tool for generating RAMSES namelist (.nml) files with preset-based configuration.
```bash
python3 ramses_nml_generator.py                                    # Interactive mode
python3 ramses_nml_generator.py --preset hr5_production -o out.nml # Non-interactive
python3 ramses_nml_generator.py --import old.nml -o new.nml        # Import and re-export
python3 ramses_nml_generator.py --list-presets                     # List available presets
```
- Dependencies: none (Python 3.6+ standard library only)
- 5 presets: cosmo_uniform, cosmo_zoomin, restart, gravity_only, hr5_production

### `create_pvar006.py`
Creates `ic_pvar_00006` (zoom geometry scalar) for RAMSES zoom-in simulations. Reads `ic_refmap` at level 7 to build zoom mask, fills zoom sub-levels 8–11 with 1.0.
```bash
python3 create_pvar006.py
```
- Dependencies: numpy

---

## 3. Analysis Plots

### Cooling & Star Formation
| Script | Description | Output |
|--------|-------------|--------|
| `plot_grackle_table.py` | Grackle binary cooling table visualization (net rate + mean molecular weight) | PNG |
| `plot_sfrd.py` | Cosmic SFRD vs redshift with Madau & Dickinson (2014) and observations | `sfrd_comparison.png` |
| `plot_sfrd_comparison.py` | Multi-run SFRD comparison (holdback / no-holdback / different g\*) | PNG |

### Power Spectrum
| Script | Description | Output |
|--------|-------------|--------|
| `compare_pk_camb.py` | Simulation P(k) vs CAMB linear theory comparison | `pk_ratio_*.png`, `growth_factor.png` |

```bash
python3 compare_pk_camb.py --basedir test_ksection/run_de_pk_test/
```

### Performance
| Script | Description | Output |
|--------|-------------|--------|
| `plot_scaling.py` | MPI + OMP scaling test log parser and figure generator | `fig_cosmo1024_mpi_scaling.pdf`, `fig_cosmo1024_omp_scaling.pdf` |
| `plot_strong_scaling.py` | Multi-node MPI strong scaling plot (ncpu = 8, 16, 32, 64) | `fig_strong_scaling.pdf` |

---


## 4. Verification Tests

| Script | Description | Validates |
|--------|-------------|-----------|
| `test_dt_discontinuity.py` | Shows dt-dependent cooling creates spurious T discontinuities at AMR level boundaries | Townsend exact integration superiority |
| `test_equilibrium_convergence.py` | Equilibrium convergence: N-R delayed approach vs Townsend instant | Cooling module accuracy |
| `test_sgs_dissipation.py` | SGS dissipation-only verification against analytical ODE solution | Implicit dissipation formula |

---

## Dependencies

| Package | Required by |
|---------|-------------|
| numpy | Nearly all scripts |
| matplotlib | Most plotting scripts |
| camb | `generate_de_table.py`, `generate_neutrino_table.py`, `de_linear_test_setup.py`, `compare_pk_camb.py` |
| scipy | `ksection_domain_viz.py`, `ksection_progressive_viz.py` |
| Pillow (PIL) | `ksection_domain_viz.py`, `ksection_progressive_viz.py` |
| python-pptx | `make_presentation.py` |
| (none) | `ramses_nml_generator.py`, `compute_pk.py` (numpy only), `create_pvar006.py` |

```bash
pip install numpy matplotlib camb scipy Pillow python-pptx  # Install all
```
