#!/usr/bin/env python3
"""Prepare the five gated L512/N1024 DE production extensions.

The script creates model-specific CAMB transfers, the clustering-DE and
massive-neutrino response tables, LagMUSIC configurations, low-resolution
gate inputs, and production namelists/scripts.  It never submits jobs.
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
from pathlib import Path
import stat
import sys

import numpy as np


OMEGA_M = 0.3111
OMEGA_B = 0.049
OMEGA_L = 1.0 - OMEGA_M
H0 = 67.66
H = H0 / 100.0
N_S = 0.9665
A_S = 2.1e-9
SUM_MNU = 0.06
ZSTART = 49.0
FORCE_PNORM_REFERENCE = 3.209495550129538e-11

MODELS = (
    "cpl_cluster_m09_p02",
    "hs10_m01",
    "nu_lcdm",
    "f6_nu",
    "ede03",
)

DESCRIPTIONS = {
    "cpl_cluster_m09_p02": "clustering CPL, w0=-0.9, wa=0.2, cs2=0",
    "hs10_m01": "lagRamses quasi-static Horndeski, mu0=0.1, m=0.1 h/Mpc",
    "nu_lcdm": "LCDM with one 0.06-eV massive-neutrino species",
    "f6_nu": "Hu-Sawicki F6 with one 0.06-eV massive-neutrino species",
    "ede03": "early dark energy, omega_ede=0.03, zc=3000, w=1/3",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--camb-dir", type=Path, required=True)
    parser.add_argument("--music", type=Path, required=True)
    parser.add_argument("--ramses", type=Path, required=True)
    parser.add_argument("--layout-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def write_text(path: Path, text: str, force: bool, executable: bool = False) -> None:
    if path.exists() and not force:
        if path.read_text() != text:
            raise FileExistsError(f"{path} differs; rerun with --force")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    if executable:
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def load_camb(camb_dir: Path):
    camb_dir = camb_dir.resolve()
    sys.path.insert(0, str(camb_dir))
    camb = importlib.import_module("camb")
    imported = Path(camb.__file__).resolve()
    if camb_dir not in imported.parents:
        raise RuntimeError(f"expected CAMB below {camb_dir}, imported {imported}")
    if not hasattr(camb.model.MuSigmaMGParams, "set_lagramses_horndeski"):
        raise RuntimeError("CAMB build lacks the exact lagRamses Horndeski model")
    return camb


def massive_neutrino_density(camb) -> float:
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=H0,
        ombh2=OMEGA_B * H**2,
        omch2=(OMEGA_M - OMEGA_B) * H**2,
        mnu=SUM_MNU,
        num_massive_neutrinos=1,
        neutrino_hierarchy="degenerate",
    )
    return float(pars.omnuh2 / H**2)


def ede_camb_fraction() -> float:
    """Map lagRamses omega_ede to CAMB's fraction including radiation."""
    ac = 1.0 / 3001.0
    omega_gamma = 2.469e-5 / H**2 * (2.7255 / 2.7255) ** 4
    omega_r = omega_gamma * (1.0 + 0.22710731766 * 3.046)
    e2_no_radiation = OMEGA_M / ac**3 + OMEGA_L
    e2_with_radiation = e2_no_radiation + omega_r / ac**4
    ratio = 0.03 * e2_no_radiation / e2_with_radiation
    return ratio / (1.0 + ratio)


def configure_model(camb, model_name: str):
    omega_nu = massive_neutrino_density(camb) if model_name in {"nu_lcdm", "f6_nu"} else 0.0
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=H0,
        ombh2=OMEGA_B * H**2,
        omch2=(OMEGA_M - OMEGA_B - omega_nu) * H**2,
        mnu=SUM_MNU if omega_nu > 0.0 else 0.0,
        num_massive_neutrinos=1 if omega_nu > 0.0 else 0,
        neutrino_hierarchy="degenerate",
    )
    pars.InitPower.set_params(ns=N_S, As=A_S)
    if model_name == "cpl_cluster_m09_p02":
        pars.set_dark_energy(w=-0.9, wa=0.2, cs2=0.0, dark_energy_model="fluid")
    elif model_name == "hs10_m01":
        pars.MG.set_lagramses_horndeski(mu0=0.1, mass_hmpc=0.1, h=H)
    elif model_name == "f6_nu":
        pars.MG.set_fR(1.0e-6, n=1.0)
    elif model_name == "ede03":
        pars.DarkEnergy = camb.dark_energy.AxionEffectiveFluid()
        pars.DarkEnergy.set_params(
            w_n=1.0 / 3.0,
            fde_zc=ede_camb_fraction(),
            zc=3000.0,
            theta_i=math.pi / 2.0,
        )
    pars.set_matter_power(redshifts=[ZSTART, 0.0], kmax=50.0)
    pars.WantTransfer = True
    pars.NonLinear = camb.model.NonLinear_none
    return pars, omega_nu


def make_transfer(camb, campaign: Path, model_name: str, force: bool) -> dict[str, float | str]:
    transfer_dir = campaign / "transfers"
    transfer_dir.mkdir(parents=True, exist_ok=True)
    output = transfer_dir / f"transfer_{model_name}_z49.dat"
    diagnostics_path = output.with_suffix(".json")
    if output.exists() and diagnostics_path.exists() and not force:
        return json.loads(diagnostics_path.read_text())

    pars, omega_nu = configure_model(camb, model_name)
    results = camb.get_results(pars)
    raw = results.get_matter_transfer_data().transfer_data.copy()
    density_var = "delta_nonu" if omega_nu > 0.0 else "delta_tot"
    density_index = 7 if omega_nu > 0.0 else 6

    kh, redshifts, linear_pk = results.get_linear_matter_power_spectrum(
        var1=density_var,
        var2=density_var,
        hubble_units=True,
        k_hunit=True,
    )
    iz = int(np.argmin(np.abs(redshifts - ZSTART)))
    if abs(float(redshifts[iz]) - ZSTART) > 1.0e-8:
        raise RuntimeError(f"{model_name}: missing z={ZSTART:g} linear spectrum")

    transfer_kh = raw[0, :, 0]
    density_transfer = raw[density_index, :, 0]
    if not np.allclose(transfer_kh, kh, rtol=1.0e-6, atol=0.0):
        raise RuntimeError(f"{model_name}: transfer and P(k) grids differ")
    samples = linear_pk[iz] / (8.0 * math.pi**3 * kh**N_S * density_transfer**2)
    mask = (kh >= 1.0e-3) & (kh <= 10.0)
    pnorm = float(np.mean(samples[mask]))
    pnorm_scatter = float(np.std(samples[mask]) / pnorm)
    if pnorm_scatter > 1.0e-6:
        raise RuntimeError(f"{model_name}: non-constant CAMB-to-MUSIC normalization")
    if abs(pnorm / FORCE_PNORM_REFERENCE - 1.0) > 5.0e-5:
        raise RuntimeError(f"{model_name}: common-As pnorm drift: {pnorm:.16e}")

    if omega_nu > 0.0:
        # LagMUSIC's DMO plugin consumes column 7 (zero-based 6) as the
        # particle-density transfer. Feed it T_cb, not total matter.
        raw[6, :, 0] = raw[7, :, 0]
        # Its vtotal weighting does not know Omega_nu. Rescale CDM and baryon
        # velocity columns so the plugin returns the cb-weighted velocity.
        velocity_rescale = OMEGA_M / (OMEGA_M - omega_nu)
        raw[10, :, 0] *= velocity_rescale
        raw[11, :, 0] *= velocity_rescale

    with output.open("w") as stream:
        for ik in range(raw.shape[1]):
            stream.write(" ".join(f"{raw[iv, ik, 0]:15.8e}" for iv in range(13)) + "\n")

    sigma8 = results.get_sigma8()
    diagnostics: dict[str, float | str] = {
        "A_s": A_S,
        "density_transfer": "cb" if omega_nu > 0.0 else "total",
        "force_pnorm": pnorm,
        "force_pnorm_relative_scatter": pnorm_scatter,
        "omega_nu": omega_nu,
        "sigma8_z0": float(sigma8[-1]),
        "sigma8_zstart": float(sigma8[0]),
        "zstart": ZSTART,
    }
    if model_name == "ede03":
        diagnostics["camb_fde_zc"] = ede_camb_fraction()
        diagnostics["lagramses_omega_ede"] = 0.03
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n")
    return diagnostics


def make_response_tables(campaign: Path, camb_dir: Path, omega_nu: float) -> None:
    aux = Path(__file__).resolve().parents[2] / "cuRamses" / "aux"
    sys.path.insert(0, str(aux))
    from generate_de_table import generate_de_table
    from generate_neutrino_table import generate_neutrino_table

    table_dir = campaign / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    generate_de_table(
        output_file=str(table_dir / "de_cpl_cluster_m09_p02.dat"),
        w0=-0.9,
        wa=0.2,
        cs2_de=0.0,
        omega_m=OMEGA_M,
        omega_b=OMEGA_B,
        h0=H,
        nk=240,
        na=80,
    )
    generate_neutrino_table(
        output_file=str(table_dir / "neutrino_mnu006.dat"),
        omega_b=OMEGA_B,
        omega_m=OMEGA_M,
        h0=H,
        sum_mnu=SUM_MNU,
        n_s=N_S,
        A_s=A_S,
        nk=240,
        na=80,
        a_min=0.01,
        a_max=1.0,
    )
    metadata = {
        "camb_dir": str(camb_dir.resolve()),
        "camb_commit": "abc2fe2",
        "omega_nu": omega_nu,
        "sum_mnu_eV": SUM_MNU,
    }
    (table_dir / "response_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )


def music_config(
    campaign: Path,
    model_name: str,
    diagnostics: dict[str, float | str],
    level: int,
    boxlen: float,
    output_dir: str,
    random_source: str,
) -> str:
    w0, wa = (-0.9, 0.2) if model_name == "cpl_cluster_m09_p02" else (-1.0, 0.0)
    masses = (SUM_MNU, 0.0, 0.0) if model_name in {"nu_lcdm", "f6_nu"} else (0.0, 0.0, 0.0)
    return f"""[setup]
boxlength       = {boxlen}
zstart          = {ZSTART}
levelmin        = {level}
levelmax        = {level}
baryons         = no
use_2LPT        = yes
use_LLA         = no
padding         = 8
force_pnorm     = {float(diagnostics['force_pnorm']):.15e}
vfact_scale     = 1.000000000000000e+00
dmo_velocity_source = transfer
kspace_TF       = yes
slab_solve_unigrid = yes
slab_2lpt_unigrid = yes
lpt2_boost_threads = auto

[cosmology]
Omega_m         = {OMEGA_M}
Omega_L         = {OMEGA_L}
Omega_b         = {OMEGA_B}
H0              = {H0}
sigma_8         = {float(diagnostics['sigma8_z0']):.10f}
nspec           = {N_S}
m_nu1           = {masses[0]}
m_nu2           = {masses[1]}
m_nu3           = {masses[2]}
w0              = {w0}
wa              = {wa}
transfer        = camb_file
transfer_file   = {campaign / 'transfers' / f'transfer_{model_name}_z49.dat'}

[random]
disk_cached     = yes
seed[{level}]        = {random_source}

[output]
format          = grafic2
filename        = {output_dir}

[poisson]
fft_fine        = no
kspace          = yes
accuracy        = 1e-5
pre_smooth      = 3
post_smooth     = 3
smoother        = gs
laplace_order   = 6
grad_order      = 6
"""


def model_controls(campaign: Path, model_name: str, omega_nu: float) -> tuple[str, str]:
    flags: list[str] = []
    blocks: list[str] = []
    if model_name == "cpl_cluster_m09_p02":
        blocks.append(
            f"""&CPL_PARAMS
w0=-0.9
wa=0.2
cs2_de=0.0
de_table='{campaign / 'tables' / 'de_cpl_cluster_m09_p02.dat'}'
/
"""
        )
    elif model_name == "hs10_m01":
        flags.append("use_horndeski=.true.")
        blocks.append("""&HORNDESKI_PARAMS
hs_mu0=0.1
hs_mass=0.1
/
""")
    elif model_name == "nu_lcdm":
        flags.append("use_neutrino=.true.")
    elif model_name == "f6_nu":
        flags.extend(("use_neutrino=.true.", "use_fR=.true."))
        blocks.append("""&FR_PARAMS
fR0=-1.0d-6
fR_n=1
n_iter_fR=6000
fR_eps=1.0d-5
/
""")
    elif model_name == "ede03":
        flags.append("use_ede=.true.")
        blocks.append("""&EDE_PARAMS
omega_ede=0.03
z_ede=3000.0
w_ede=0.3333333333333333
/
""")
    if model_name in {"nu_lcdm", "f6_nu"}:
        blocks.append(
            f"""&NEUTRINO_PARAMS
omega_nu={omega_nu:.15e}
neutrino_table='{campaign / 'tables' / 'neutrino_mnu006.dat'}'
/
"""
        )
    return "\n".join(flags), "\n".join(blocks)


def simulation_namelist(
    campaign: Path,
    model_name: str,
    omega_nu: float,
    level: int,
    ic_dir: str,
    ngridtot: int,
    nparttot: int,
) -> str:
    flags, blocks = model_controls(campaign, model_name, omega_nu)
    de_perturb = ".true." if model_name == "cpl_cluster_m09_p02" else ".false."
    return f"""! {DESCRIPTIONS[model_name]}
! Exact model-specific transfer; common A_s and Gaussian realization.
&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
sink=.false.
nrestart=0
nstepmax=1000000
nsubcycle=1
aexp_step_limit=0.1
ordering='ksection'
memory_balance=.true.
use_fftw=.true.
dump_pk=.false.
de_perturb={de_perturb}
exchange_method='auto'
{flags}
/

&OUTPUT_PARAMS
noutput=4
aout=0.333333333,0.500000000,0.666666667,1.000000000
foutput=100000
match_aout=.true.
/

&COSMO_PARAMS
omega_m={OMEGA_M}
omega_l={OMEGA_L}
omega_b=0.0
h0={H0}
/

&INIT_PARAMS
filetype='grafic'
initfile(1)='../{ic_dir}/level_{level:03d}'
/

&AMR_PARAMS
levelmin={level}
levelmax={level}
nexpand=1
ngridtot={ngridtot}
nparttot={nparttot}
/

&REFINE_PARAMS
m_refine=1*8.0
ivar_refine=0
/

&POISSON_PARAMS
epsilon=1.0d-4
/

{blocks}"""


def production_slurm(campaign: Path, model_name: str, ramses: Path) -> str:
    model_dir = campaign / model_name
    return f"""#!/bin/bash
#SBATCH --job-name=dmo_{model_name}
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --exclude=grammar007,grammar012,grammar023,grammar081,grammar112
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=2
#SBATCH --mem=480G
#SBATCH --time=7-00:00:00
#SBATCH --chdir={model_dir}
#SBATCH --output=run-%j.out
#SBATCH --error=run-%j.err

set -euo pipefail
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M I_MPI_PIN_DOMAIN=omp I_MPI_PIN_ORDER=compact
echo "model={model_name} code=e964a94 extension5 start=$(date --iso-8601=seconds)"
sha256sum {ramses}
mpirun -np "${{SLURM_NTASKS:-32}}" {ramses} run.nml
test -f output_00005/info_00005.txt
echo "model={model_name} end=$(date --iso-8601=seconds)"
"""


def write_layout(args: argparse.Namespace, diagnostics: dict[str, dict[str, float | str]] | None) -> None:
    campaign = args.campaign.resolve()
    control = campaign / "extension5_20260805"
    validation = control / "validation"
    omega_nu = (
        float(diagnostics["nu_lcdm"]["omega_nu"])
        if diagnostics is not None
        else 0.0014
    )
    for index, model_name in enumerate(MODELS):
        prod_dir = campaign / model_name
        prod_dir.mkdir(parents=True, exist_ok=True)
        write_text(
            prod_dir / "run.nml",
            simulation_namelist(
                campaign,
                model_name,
                omega_nu,
                level=10,
                ic_dir=f"ics_{model_name}",
                ngridtot=1_342_177_280,
                nparttot=4_294_967_296,
            ),
            args.force,
        )
        write_text(
            prod_dir / "run-extension5.slurm",
            production_slurm(campaign, model_name, args.ramses.resolve()),
            args.force,
            executable=True,
        )
        val_dir = validation / model_name
        val_dir.mkdir(parents=True, exist_ok=True)
        write_text(
            val_dir / "run.nml",
            simulation_namelist(
                campaign,
                model_name,
                omega_nu,
                level=6,
                ic_dir=f"ics_{model_name}",
                ngridtot=2_000_000,
                nparttot=1_000_000,
            ),
            args.force,
        )
        if diagnostics is not None:
            prod_random = str(campaign / "wnoise_0010.bin")
            write_text(
                campaign / f"music_{model_name}.conf",
                music_config(
                    campaign,
                    model_name,
                    diagnostics[model_name],
                    10,
                    512.0,
                    f"ics_{model_name}",
                    prod_random,
                ),
                args.force,
            )
            val_random = str(20260805) if index == 0 else str(validation / "wnoise_0006.bin")
            write_text(
                validation / f"music_{model_name}.conf",
                music_config(
                    campaign,
                    model_name,
                    diagnostics[model_name],
                    6,
                    64.0,
                    f"ics_{model_name}",
                    val_random,
                ),
                args.force,
            )

    manifest = {
        "models": list(MODELS),
        "common_A_s": A_S,
        "common_n_s": N_S,
        "production_white_noise": str(campaign / "wnoise_0010.bin"),
        "lagcamb_commit": "abc2fe2",
        "lagramses_binary": str(args.ramses.resolve()),
        "diagnostics": diagnostics,
    }
    write_text(
        control / "manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        args.force,
    )


def main() -> int:
    args = parse_args()
    args.campaign = args.campaign.resolve()
    if args.layout_only:
        write_layout(args, None)
        return 0
    camb = load_camb(args.camb_dir)
    diagnostics = {
        model_name: make_transfer(camb, args.campaign, model_name, args.force)
        for model_name in MODELS
    }
    omega_nu = float(diagnostics["nu_lcdm"]["omega_nu"])
    make_response_tables(args.campaign, args.camb_dir, omega_nu)
    write_layout(args, diagnostics)
    print(json.dumps(diagnostics, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
