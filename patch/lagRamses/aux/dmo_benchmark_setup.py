#!/usr/bin/env python3
"""Build a matched-phase DMO benchmark campaign for lagRamses/cuRAMSES.

The default campaign uses model-specific transfer functions from the local
lagCAMB tree when the requested linear model is available. All initial
conditions retain the same random phases. Models without a corresponding
lagCAMB transfer use the LCDM transfer at the high starting redshift and are
identified explicitly in the campaign metadata. The metadata records the
forward species contract separately from the initial transfer selection.

Default validation models
-------------------------
lcdm, f5, f6, n1, n5, sym_a

Example
-------
python3 dmo_benchmark_setup.py \
    --outdir /gpfs/kjhan/Hydro/DE_nonstd/DMO_bench_v1 \
    --make-ics --ic-mode model
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import stat
import subprocess
import sys
from typing import Dict

import numpy as np


OMEGA_M = 0.3111
OMEGA_B = 0.049
H0 = 67.66
N_S = 0.9665
A_S = 2.1e-9

DEFAULT_MUSIC = "/home/kjhan/BACKUP/LagMUSIC/music/build/MUSIC"
DEFAULT_RAMSES = (
    "/home/kjhan/BACKUP/lagRamses-de-nonstd/bin/ramses_dmo_bench3d"
)
DEFAULT_CAMB = "/home/kjhan/BACKUP/CAMB/CAMB"


MODELS: Dict[str, Dict[str, str]] = {
    "lcdm": {
        "description": "LCDM control with every non-standard switch disabled",
        "flags": "",
        "blocks": "",
    },
    "q1": {
        "description": "legacy frozen-field Ratra-Peebles quintessence, alpha=1",
        "flags": "use_quintessence=.true.",
        "blocks": """&QUINT_PARAMS
quint_pot=1
quint_ic_mode=0
quint_alpha=1.0
quint_lambda=1.0
quint_phi_ini=0.01
/
""",
    },
    "phicdm_a01": {
        "description": (
            "Ratra-Peebles phiCDM matter-era tracker, alpha=0.1 "
            "(production benchmark)"
        ),
        "flags": "use_quintessence=.true.",
        "blocks": """&QUINT_PARAMS
quint_pot=1
quint_ic_mode=1
quint_alpha=0.1
quint_lambda=1.0
/
""",
    },
    "cde10": {
        "description": (
            "Effective one-fluid coupled-quintessence DMO model, beta=0.1, "
            "with mass, friction, and fifth-force terms"
        ),
        "flags": """use_quintessence=.true.
use_coupled_de=.true.""",
        "blocks": """&QUINT_PARAMS
quint_pot=1
quint_alpha=1.0
quint_lambda=1.0
quint_phi_ini=0.01
/

&COUPLED_DE_PARAMS
beta_cde=0.1
cde_friction=.true.
cde_vary_mass=.true.
/
""",
    },
    "f5": {
        "description": "Hu-Sawicki f(R), F5, n=1",
        "flags": "use_fR=.true.",
        "blocks": """&FR_PARAMS
fR0=-1.0d-5
fR_n=1
n_iter_fR=20
fR_eps=1.0d-6
/
""",
    },
    "f6": {
        "description": "Hu-Sawicki f(R), F6, n=1",
        "flags": "use_fR=.true.",
        "blocks": """&FR_PARAMS
fR0=-1.0d-6
fR_n=1
n_iter_fR=20
fR_eps=1.0d-6
/
""",
    },
    "n1": {
        "description": "normal-branch DGP, H0*rc=1, Omega_rc=0.25",
        "flags": "use_nDGP=.true.",
        "blocks": """&NDGP_PARAMS
omega_rc=0.25
nDGP_branch=1
n_iter_nDGP=20
nDGP_eps=1.0d-6
/
""",
    },
    "n5": {
        "description": "normal-branch DGP, H0*rc=5, Omega_rc=0.01",
        "flags": "use_nDGP=.true.",
        "blocks": """&NDGP_PARAMS
omega_rc=0.01
nDGP_branch=1
n_iter_nDGP=20
nDGP_eps=1.0d-6
/
""",
    },
    "sym_a": {
        "description": "Symmetron A: a_ssb=0.5, beta=1, L=1 Mpc/h",
        "flags": "use_symmetron=.true.",
        "blocks": """&SYMMETRON_PARAMS
a_ssb=0.5
beta_symmetron=1.0
L_symmetron=1.0
n_iter_symmetron=20
symmetron_eps=1.0d-6
/
""",
    },
    "w09": {
        "description": "smooth constant-w dark energy, w=-0.9",
        "flags": "",
        "music_w0": -0.9,
        "music_wa": 0.0,
        "blocks": """&CPL_PARAMS
w0=-0.9
wa=0.0
cs2_de=1.0
/
""",
    },
    "cpl_m09_p02": {
        "description": "smooth CPL dark energy, w0=-0.9, wa=+0.2",
        "flags": "",
        "music_w0": -0.9,
        "music_wa": 0.2,
        "blocks": """&CPL_PARAMS
w0=-0.9
wa=0.2
cs2_de=1.0
/
""",
    },
    "kess51": {
        "description": "purely kinetic k-essence, X0=0.51",
        "flags": "use_kessence=.true.",
        "blocks": """&KESSENCE_PARAMS
kes_x0=0.51
/
""",
    },
    "kess5001": {
        "description": "purely kinetic k-essence, X0=0.5001",
        "flags": "use_kessence=.true.",
        "blocks": """&KESSENCE_PARAMS
kes_x0=0.5001
/
""",
    },
    "ede03": {
        "description": "early dark energy dispatch test, Omega_EDE=0.03",
        "flags": "use_ede=.true.",
        "blocks": """&EDE_PARAMS
omega_ede=0.03
z_ede=3000.0
w_ede=0.333333333333
/
""",
    },
    "chap075": {
        "description": "generalized Chaplygin gas, As=0.75, alpha=1e-4",
        "flags": "use_chaplygin=.true.",
        "blocks": """&CHAPLYGIN_PARAMS
chaplygin_As=0.75
chaplygin_alpha=1.0d-4
/
""",
    },
    "chap099": {
        "description": "generalized Chaplygin gas, As=0.99, alpha=1e-4",
        "flags": "use_chaplygin=.true.",
        "blocks": """&CHAPLYGIN_PARAMS
chaplygin_As=0.99
chaplygin_alpha=1.0d-4
/
""",
    },
    "rvm001": {
        "description": "running vacuum, nu=1e-3",
        "flags": "use_rvm=.true.",
        "blocks": """&RVM_PARAMS
rvm_nu=1.0d-3
/
""",
    },
    "hs10_m01": {
        "description": "quasi-static Horndeski, mu0=0.1, mass=0.1 h/Mpc",
        "flags": "use_horndeski=.true.",
        "blocks": """&HORNDESKI_PARAMS
hs_mu0=0.1
hs_mass=0.1
/
""",
    },
    "dilaton_a": {
        "description": "environmentally dependent dilaton benchmark",
        "flags": "use_dilaton=.true.",
        "blocks": """&DILATON_PARAMS
beta_dilaton=1.0
L_dilaton=1.0
a0_dilaton=0.5
n_iter_dilaton=20
dilaton_eps=1.0d-6
/
""",
    },
    "gal_tracker": {
        "description": "parameter-free cubic Galileon tracker",
        "flags": "use_galileon=.true.",
        "blocks": """&GALILEON_PARAMS
galileon_tracker=.true.
c2_galileon=-1.0
c3_galileon=1.0
n_iter_galileon=20
galileon_eps=1.0d-6
/
""",
    },
}

DEFAULT_MODELS = ("lcdm", "f5", "f6", "n1", "n5", "sym_a")

# These models have parameter-for-parameter counterparts in the local lagCAMB.
# All other models intentionally fall back to the LCDM high-redshift transfer.
CAMB_MATCHED_MODELS = {
    "lcdm",
    "q1",
    "phicdm_a01",
    "cde10",
    "w09",
    "cpl_m09_p02",
    "kess51",
    "kess5001",
    "chap075",
    "chap099",
    "rvm001",
    "f5",
    "f6",
    "n1",
    "n5",
    "sym_a",
}

SINGLE_FLUID_DYNAMICS_NOTES = {
    "cde10": (
        "lagCAMB couples the scalar field to CDM, while the one-species DMO "
        "forward model applies the effective coupled dynamics to the full "
        "CDM+baryon proxy fluid"
    ),
    "rvm001": (
        "lagCAMB exchanges vacuum energy with CDM, while the one-species DMO "
        "forward model applies the effective running-matter law to the full "
        "CDM+baryon proxy fluid"
    ),
}

SINGLE_FLUID_SPECIES_CONTRACT = {
    "cde10": {
        "linear_exchange_species": "cdm_only",
        "nbody_evolved_species": "single_cdm_baryon_proxy_fluid",
        "nbody_exchange_species": "single_cdm_baryon_proxy_fluid",
        "forward_species_contract_match": False,
        "claim_scope": "effective_one_fluid_response",
        "physical_precision_claim_allowed": False,
        "validation_status": "level_4_internal_effective_one_fluid",
    },
    "rvm001": {
        "linear_exchange_species": "cdm_only",
        "nbody_evolved_species": "single_cdm_baryon_proxy_fluid",
        "nbody_exchange_species": "single_cdm_baryon_proxy_fluid",
        "forward_species_contract_match": False,
        "claim_scope": "effective_one_fluid_response",
        "physical_precision_claim_allowed": False,
        "validation_status": "level_4_internal_effective_one_fluid",
    },
}

SPECIES_BACKGROUND_PREFLIGHT = {
    "file": "BACKGROUND_PREFLIGHT.md",
    "method": "fixed-total-matter comparison of species-resolved and effective one-fluid backgrounds",
    "quantity": "H_species_resolved/H_effective_one_fluid - 1",
    "redshifts": [49.0, 10.0, 2.0, 1.0, 0.5, 0.0],
    "delta_h_over_h": {
        "cde10": [-0.010609, -0.006705, -0.001629, 0.001259, 0.002676, 0.0],
        "rvm001": [0.0008287, 0.0004818, 0.0001681, 0.0000735, 0.00002412, 0.0],
    },
    "interpretation": "approximation_diagnostic_not_correction",
    "derivation_mode": "manually_tabulated_from_direct_background_evaluations",
    "derivation_procedure": (
        "At fixed H0=67.66, Omega_m=0.3111, and Omega_b=0.049, evaluate H(z) "
        "with the lagCAMB CDM-only exchange and with the corresponding "
        "effective law applied to the full proxy matter density, then form "
        "H_species_resolved/H_effective_one_fluid - 1 at the listed redshifts"
    ),
    "derivation_inputs": {
        "H0_km_s_Mpc": H0,
        "Omega_m": OMEGA_M,
        "Omega_b": OMEGA_B,
        "cde10_beta": 0.1,
        "rvm001_nu": 0.001,
    },
}

LCDM_TRANSFER_REASONS = {
    "ede03": "lagRamses EDE dispatch model has no parameter-identical lagCAMB counterpart",
    "hs10_m01": "lagRamses quasi-static mu(k,a) does not match lagCAMB Bellini-alpha Horndeski",
    "dilaton_a": "no parameter-identical lagCAMB dilaton implementation is available",
    "gal_tracker": (
        "no parameter-identical lagCAMB cubic-Galileon transfer implementation "
        "is available"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--models",
        nargs="+",
        choices=tuple(MODELS),
        default=list(DEFAULT_MODELS),
        help="model directories to generate",
    )
    parser.add_argument("--boxlen", type=float, default=500.0, help="Mpc/h")
    parser.add_argument("--levelmin", type=int, default=8)
    parser.add_argument("--levelmax", type=int, default=14)
    parser.add_argument("--zstart", type=float, default=49.0)
    parser.add_argument("--seed", type=int, default=20260725)
    parser.add_argument(
        "--phase-anchor-level",
        type=int,
        help=(
            "generate white noise at this level and restrict it to levelmin; "
            "use one common finest anchor across a resolution ladder"
        ),
    )
    parser.add_argument("--ngridtot", type=int, default=80_000_000)
    parser.add_argument("--nparttot", type=int, default=25_000_000)
    parser.add_argument("--music", default=DEFAULT_MUSIC)
    parser.add_argument("--ramses", default=DEFAULT_RAMSES)
    parser.add_argument(
        "--camb-dir",
        type=Path,
        default=Path(DEFAULT_CAMB),
        help="local lagCAMB source/build directory",
    )
    parser.add_argument(
        "--ic-mode",
        choices=("model", "shared"),
        default="model",
        help="model: lagCAMB-matched transfers; shared: one LCDM transfer for triage",
    )
    parser.add_argument(
        "--dmo-velocity-source",
        choices=("transfer", "density_2lpt"),
        default="transfer",
        help=(
            "transfer: use CAMB vtotal plus MUSIC's 2LPT term (default); "
            "density_2lpt: derive velocities from the density potential"
        ),
    )
    parser.add_argument(
        "--music-tasks",
        type=int,
        default=1,
        help="MPI ranks used while generating ICs",
    )
    parser.add_argument(
        "--music-unigrid-slab",
        action="store_true",
        help=(
            "use LagMUSIC's distributed slab FFT and 2LPT paths; "
            "recommended for uniform 512^3 and larger ICs"
        ),
    )
    parser.add_argument(
        "--reuse-white-noise",
        action="store_true",
        help=(
            "generate the first IC from the seed, then make every later "
            "model read the stored wnoise_LEVEL.bin realization"
        ),
    )
    parser.add_argument(
        "--white-noise-file",
        type=Path,
        help=(
            "read this existing white-noise realization for every generated IC; "
            "the file is never copied or modified"
        ),
    )
    parser.add_argument("--slurm-tasks", type=int, default=32)
    parser.add_argument("--omp-threads", type=int, default=2)
    parser.add_argument("--slurm-memory", default="420G")
    parser.add_argument(
        "--slurm-exclude",
        default="",
        help="comma-separated Slurm nodes to exclude from IC and simulation jobs",
    )
    parser.add_argument(
        "--dump-pk",
        action="store_true",
        help="measure spectra in situ (disabled by default for production memory safety)",
    )
    parser.add_argument(
        "--scalar-iters",
        type=int,
        default=1000,
        help="maximum nonlinear iterations (solvers exit immediately on convergence)",
    )
    parser.add_argument(
        "--scalar-eps",
        type=float,
        default=1.0e-6,
        help="relative residual tolerance for nonlinear scalar solvers",
    )
    parser.add_argument(
        "--aexp-step-limit",
        type=float,
        default=0.1,
        help="maximum fractional expansion-factor change per coarse step",
    )
    parser.add_argument(
        "--output-redshifts",
        nargs="+",
        type=float,
        default=[5.0, 2.0, 1.0, 0.8, 0.5, 0.2, 0.0],
        help="strictly decreasing output redshifts after the initial dump",
    )
    parser.add_argument("--make-ics", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace generated text files; simulation outputs are never removed",
    )
    return parser.parse_args()


def write_text(path: Path, content: str, force: bool, executable: bool = False) -> None:
    if path.exists() and not force:
        old = path.read_text()
        if old != content:
            raise FileExistsError(
                f"{path} exists with different content; use --force to replace it"
            )
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    if executable:
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def file_sha256(path: Path) -> str | None:
    """Return a binary digest when the generated campaign can read the file."""
    path = path.expanduser().resolve()
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_revision(path: Path) -> str | None:
    """Return the source revision that contains path, without changing state."""
    path = path.expanduser().resolve()
    workdir = path if path.is_dir() else path.parent
    result = subprocess.run(
        ["git", "-C", str(workdir), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def git_tracked_worktree_dirty(path: Path) -> bool | None:
    """Report whether tracked source differs from HEAD."""
    path = path.expanduser().resolve()
    workdir = path if path.is_dir() else path.parent
    result = subprocess.run(
        [
            "git",
            "-C",
            str(workdir),
            "status",
            "--porcelain",
            "--untracked-files=no",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    return bool(result.stdout.strip()) if result.returncode == 0 else None


def load_local_camb(camb_dir: Path):
    """Import CAMB from the requested tree, never silently from site-packages."""
    camb_dir = camb_dir.expanduser().resolve()
    if not (camb_dir / "camb" / "__init__.py").is_file():
        raise FileNotFoundError(f"CAMB Python package not found below {camb_dir}")
    sys.path.insert(0, str(camb_dir))
    try:
        camb = importlib.import_module("camb")
    except ImportError as exc:
        raise RuntimeError(f"cannot import local CAMB from {camb_dir}") from exc
    imported = Path(camb.__file__).resolve()
    if camb_dir not in imported.parents:
        raise RuntimeError(f"expected CAMB below {camb_dir}, imported {imported}")
    return camb


def configure_camb_dark_energy(pars, model_name: str, camb) -> None:
    """Apply the lagCAMB parameterization associated with a campaign model."""
    de = camb.dark_energy
    if model_name == "lcdm":
        return
    if model_name == "w09":
        pars.set_dark_energy(w=-0.9, wa=0.0, dark_energy_model="fluid")
    elif model_name == "cpl_m09_p02":
        pars.set_dark_energy(w=-0.9, wa=0.2, dark_energy_model="fluid")
    elif model_name == "q1":
        pars.DarkEnergy = de.TrackerQuintessence()
        pars.DarkEnergy.set_params(
            pot_type=1, ic_mode=0, alpha=1.0, lam=1.0, phi_ini=0.01
        )
    elif model_name == "phicdm_a01":
        pars.DarkEnergy = de.TrackerQuintessence()
        pars.DarkEnergy.set_params(pot_type=1, ic_mode=1, alpha=0.1, lam=1.0)
    elif model_name == "cde10":
        pars.DarkEnergy = de.CoupledQuintessence()
        pars.DarkEnergy.set_params(
            pot_type=1, alpha=1.0, lam=1.0, phi_ini=0.01, beta=0.1
        )
    elif model_name in {"kess51", "kess5001"}:
        pars.DarkEnergy = de.KEssence()
        pars.DarkEnergy.set_params(x0=0.51 if model_name == "kess51" else 0.5001)
    elif model_name in {"chap075", "chap099"}:
        pars.DarkEnergy = de.Chaplygin()
        pars.DarkEnergy.set_params(
            As=0.75 if model_name == "chap075" else 0.99, alpha=1.0e-4
        )
    elif model_name == "rvm001":
        pars.DarkEnergy = de.RunningVacuum()
        pars.DarkEnergy.set_params(nu=1.0e-3)
    elif model_name == "f5":
        pars.MG.set_fR(1.0e-5, n=1.0)
    elif model_name == "f6":
        pars.MG.set_fR(1.0e-6, n=1.0)
    elif model_name == "n1":
        pars.MG.set_nDGP(H0rc=1.0)
    elif model_name == "n5":
        pars.MG.set_nDGP(H0rc=5.0)
    elif model_name == "sym_a":
        # lagRamses stores L in Mpc/h; lagCAMB expects physical Mpc.
        pars.MG.set_symmetron(a_ssb=0.5, beta=1.0, L_Mpc=1.0/(H0/100.0))
    else:
        raise ValueError(f"no exact lagCAMB transfer mapping for {model_name}")


def make_transfer(
    outdir: Path, model_name: str, zstart: float, force: bool, camb
) -> tuple[Path, dict[str, float]]:
    transfer_dir = outdir / "transfers"
    transfer_path = transfer_dir / f"transfer_{model_name}_z{zstart:g}.dat"
    diagnostics_path = transfer_dir / f"transfer_{model_name}_z{zstart:g}.json"
    if transfer_path.exists() and diagnostics_path.exists() and not force:
        return transfer_path, json.loads(diagnostics_path.read_text())

    h = H0 / 100.0
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=H0,
        ombh2=OMEGA_B * h**2,
        omch2=(OMEGA_M - OMEGA_B) * h**2,
        mnu=0.0,
    )
    pars.InitPower.set_params(ns=N_S, As=A_S)
    configure_camb_dark_energy(pars, model_name, camb)
    pars.set_matter_power(redshifts=[zstart, 0.0], kmax=50.0)
    pars.WantTransfer = True
    results = camb.get_results(pars)
    transfer = results.get_matter_transfer_data().transfer_data
    kh, pk_redshifts, linear_pk = results.get_linear_matter_power_spectrum(
        hubble_units=True, k_hunit=True
    )

    # CAMB sorts requested redshifts from high to low; zstart is index zero.
    transfer_path.parent.mkdir(parents=True, exist_ok=True)
    with transfer_path.open("w") as stream:
        for ik in range(transfer.shape[1]):
            stream.write(
                " ".join(f"{transfer[iv, ik, 0]:15.8e}" for iv in range(transfer.shape[0]))
                + "\n"
            )

    sigma8_values = results.get_sigma8()
    zstart_pk_index = int(abs(pk_redshifts - zstart).argmin())
    if not math.isclose(
        float(pk_redshifts[zstart_pk_index]), zstart, rel_tol=0.0, abs_tol=1.0e-8
    ):
        raise RuntimeError(f"CAMB did not return a linear P(k) table at z={zstart}")

    # LagMUSIC defines P(k) = pnorm k^ns T(k)^2 and its diagnostic spectrum is
    # smaller than the conventional CAMB P(k) by 8*pi^3. Derive pnorm directly
    # from CAMB's linear spectrum, so A_s (not sigma8 or MUSIC back-scaling)
    # fixes the absolute IC amplitude.
    transfer_kh = transfer[0, :, 0]
    transfer_total = transfer[6, :, 0]
    if not ((abs(transfer_kh / kh - 1.0)) < 1.0e-6).all():
        raise RuntimeError("CAMB transfer and linear-P(k) grids do not match")
    pnorm_samples = linear_pk[zstart_pk_index] / (
        8.0 * math.pi**3 * kh**N_S * transfer_total**2
    )
    pnorm_mask = (kh >= 1.0e-3) & (kh <= 10.0)
    pnorm = float(pnorm_samples[pnorm_mask].mean())
    pnorm_scatter = float(pnorm_samples[pnorm_mask].std() / pnorm)
    if pnorm_scatter > 1.0e-6:
        raise RuntimeError(
            f"CAMB-to-MUSIC amplitude conversion is not constant: {pnorm_scatter:.3e}"
        )
    diagnostics = {
        "force_pnorm": pnorm,
        "force_pnorm_relative_scatter": pnorm_scatter,
        "sigma8_zstart": float(sigma8_values[0]),
        "sigma8_z0": float(sigma8_values[-1]),
        "zstart": zstart,
        "omega_b": OMEGA_B,
        "omega_cdm": OMEGA_M - OMEGA_B,
        "species_convention": "single_fluid_total_matter",
    }
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n")
    return transfer_path, diagnostics


def velocity_growth_correction(
    model_transfer: Path, lcdm_transfer: Path
) -> tuple[np.ndarray, float, float, float]:
    """Return the mode-by-mode DMO velocity-growth correction to LCDM.

    LagMUSIC's DMO 2LPT path derives velocities from the total-density
    displacement potential and its background (LCDM/w0-wa) growth factor.
    The CAMB Newtonian velocity transfers therefore supply the relative
    correction

        [(T_v/T_delta)_model / (T_v/T_delta)_LCDM](k).

    The ratio is a diagnostic only. LagMUSIC consumes the original CAMB
    ``vtotal`` transfer directly in its opt-in DMO 2LPT velocity solve.
    """
    model = np.loadtxt(model_transfer)
    lcdm = np.loadtxt(lcdm_transfer)
    if model.shape[1] < 12 or lcdm.shape[1] < 12:
        raise RuntimeError("CAMB transfer tables must contain at least 12 columns")

    omega_c = OMEGA_M - OMEGA_B
    model_vtotal = (omega_c * model[:, 10] + OMEGA_B * model[:, 11]) / OMEGA_M
    lcdm_vtotal = (omega_c * lcdm[:, 10] + OMEGA_B * lcdm[:, 11]) / OMEGA_M
    model_growth = model_vtotal / model[:, 6]
    lcdm_growth = lcdm_vtotal / lcdm[:, 6]
    if np.any(model_growth <= 0.0) or np.any(lcdm_growth <= 0.0):
        raise RuntimeError("CAMB velocity-to-density transfer must be positive")
    lcdm_on_model_k_all = np.exp(
        np.interp(
            np.log(model[:, 0]),
            np.log(lcdm[:, 0]),
            np.log(lcdm_growth),
        )
    )
    correction = model_growth / lcdm_on_model_k_all
    mask = (
        (model[:, 0] >= 1.0e-3)
        & (model[:, 0] <= 1.0)
        & (model[:, 0] >= lcdm[:, 0].min())
        & (model[:, 0] <= lcdm[:, 0].max())
    )
    samples = correction[mask]
    if samples.size == 0 or not np.all(np.isfinite(samples)):
        raise RuntimeError("cannot derive a finite velocity-growth correction")

    scale = float(samples.mean())
    relative_scatter = float(samples.std() / scale)
    maximum_deviation = float(np.max(np.abs(samples / scale - 1.0)))
    if scale <= 0.0:
        raise RuntimeError(f"invalid velocity-growth correction: {scale}")
    return correction, scale, relative_scatter, maximum_deviation


def music_config(
    args: argparse.Namespace,
    model_name: str,
    transfer: Path,
    sigma8_z0: float,
    force_pnorm: float,
    vfact_scale: float,
    dmo_velocity_source: str,
    ic_dir: str,
    random_source: str,
) -> str:
    legacy_keys = "LagMUSIC" in str(Path(args.music).expanduser().resolve())
    spectral_key = "nspec" if legacy_keys else "n_s"
    w0_key = "w0" if legacy_keys else "w_0"
    wa_key = "wa" if legacy_keys else "w_a"
    music_w0 = float(MODELS[model_name].get("music_w0", -1.0))
    music_wa = float(MODELS[model_name].get("music_wa", 0.0))
    slab_setup = ""
    poisson_mode = "fft_fine        = yes"
    if args.music_unigrid_slab:
        slab_setup = """kspace_TF       = yes
slab_solve_unigrid = yes
slab_2lpt_unigrid = yes
lpt2_boost_threads = auto
"""
        poisson_mode = """fft_fine        = no
kspace           = yes"""
    return f"""[setup]
boxlength       = {args.boxlen}
zstart          = {args.zstart}
levelmin        = {args.levelmin}
levelmax        = {args.levelmin}
baryons         = no
use_2LPT        = yes
use_LLA         = no
padding         = 8
force_pnorm     = {force_pnorm:.15e}
vfact_scale     = {vfact_scale:.15e}
dmo_velocity_source = {dmo_velocity_source}
{slab_setup}

[cosmology]
Omega_m         = {OMEGA_M}
Omega_L         = {1.0 - OMEGA_M}
Omega_b         = {OMEGA_B}
H0              = {H0}
sigma_8         = {sigma8_z0:.10f}
{spectral_key}  = {N_S}
m_nu1           = 0.0
m_nu2           = 0.0
m_nu3           = 0.0
{w0_key}        = {music_w0:.16g}
{wa_key}        = {music_wa:.16g}
transfer        = camb_file
transfer_file   = {transfer}

[random]
disk_cached     = yes
seed[{args.phase_anchor_level or args.levelmin}] = {random_source}

[output]
format          = grafic2
filename        = {ic_dir}

[poisson]
{poisson_mode}
accuracy        = 1e-5
pre_smooth      = 3
post_smooth     = 3
smoother        = gs
laplace_order   = 6
grad_order      = 6
"""


def effective_model_blocks(
    args: argparse.Namespace, model: Dict[str, str]
) -> str:
    """Return model blocks after applying the requested solver controls."""
    blocks = model["blocks"]
    blocks = blocks.replace("n_iter_fR=20", f"n_iter_fR={args.scalar_iters}")
    blocks = blocks.replace("n_iter_nDGP=20", f"n_iter_nDGP={args.scalar_iters}")
    blocks = blocks.replace(
        "n_iter_symmetron=20", f"n_iter_symmetron={args.scalar_iters}"
    )
    blocks = blocks.replace(
        "n_iter_dilaton=20", f"n_iter_dilaton={args.scalar_iters}"
    )
    blocks = blocks.replace(
        "n_iter_galileon=20", f"n_iter_galileon={args.scalar_iters}"
    )
    blocks = blocks.replace("fR_eps=1.0d-6", f"fR_eps={args.scalar_eps:.8e}")
    blocks = blocks.replace("nDGP_eps=1.0d-6", f"nDGP_eps={args.scalar_eps:.8e}")
    blocks = blocks.replace(
        "symmetron_eps=1.0d-6", f"symmetron_eps={args.scalar_eps:.8e}"
    )
    blocks = blocks.replace(
        "dilaton_eps=1.0d-6", f"dilaton_eps={args.scalar_eps:.8e}"
    )
    blocks = blocks.replace(
        "galileon_eps=1.0d-6", f"galileon_eps={args.scalar_eps:.8e}"
    )
    return blocks


def namelist(
    args: argparse.Namespace, model: Dict[str, str], ic_dir: str, ic_model: str
) -> str:
    outputs_z = args.output_redshifts
    outputs_a = ",".join(f"{1.0 / (1.0 + z):.9f}" for z in outputs_z)
    nlevels = args.levelmax - args.levelmin + 1
    nsubcycle = (
        "1"
        if args.levelmax == args.levelmin
        else f"1,{args.levelmax - args.levelmin}*2"
    )
    flags = model["flags"]
    if flags:
        flags += "\n"
    blocks = effective_model_blocks(args, model)
    return f"""! {model['description']}
! 2LPT IC transfer source = {ic_model}
! All models use the same random phases.
&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
sink=.false.
nrestart=0
nstepmax=1000000
nsubcycle={nsubcycle}
aexp_step_limit={args.aexp_step_limit:.8g}
ordering='ksection'
memory_balance=.true.
use_fftw=.true.
dump_pk={'.true.' if args.dump_pk else '.false.'}
de_perturb=.false.
exchange_method='auto'
{flags}/

&OUTPUT_PARAMS
noutput={len(outputs_z)}
aout={outputs_a}
foutput=100000
match_aout=.true.
/

&COSMO_PARAMS
omega_m={OMEGA_M}
omega_l={1.0 - OMEGA_M}
omega_b=0.0
h0={H0}
/

&INIT_PARAMS
filetype='grafic'
initfile(1)='../{ic_dir}/level_{args.levelmin:03d}'
/

&AMR_PARAMS
levelmin={args.levelmin}
levelmax={args.levelmax}
nexpand=1
ngridtot={args.ngridtot}
nparttot={args.nparttot}
/

&REFINE_PARAMS
m_refine={nlevels}*8.0
ivar_refine=0
/

&POISSON_PARAMS
epsilon=1.0d-4
/

{blocks}"""


def slurm_script(args: argparse.Namespace, model_name: str, model_dir: Path) -> str:
    final_output = len(args.output_redshifts) + 1
    exclude_line = (
        f"#SBATCH --exclude={args.slurm_exclude}\n" if args.slurm_exclude else ""
    )
    return f"""#!/bin/bash
#SBATCH --job-name=dmo_{model_name}
#SBATCH --partition=normal
#SBATCH --nodes=1
{exclude_line}#SBATCH --ntasks={args.slurm_tasks}
#SBATCH --cpus-per-task={args.omp_threads}
#SBATCH --mem={args.slurm_memory}
#SBATCH --time=7-00:00:00
#SBATCH --chdir={model_dir}
#SBATCH --output=run-%j.out
#SBATCH --error=run-%j.err

set -euo pipefail
export OMP_NUM_THREADS={args.omp_threads}
export OMP_STACKSIZE=256M
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact

cd {model_dir}
echo "model={model_name} start=$(date --iso-8601=seconds)"
echo "nodes=${{SLURM_NNODES:-1}} tasks=${{SLURM_NTASKS:-{args.slurm_tasks}}} omp=$OMP_NUM_THREADS"
sha256sum {args.ramses}
mpirun -np "${{SLURM_NTASKS:-{args.slurm_tasks}}}" {args.ramses} run.nml
test -f output_{final_output:05d}/info_{final_output:05d}.txt
echo "model={model_name} end=$(date --iso-8601=seconds)"
"""


def ic_slurm_script(
    args: argparse.Namespace,
    outdir: Path,
    transfer_models: list[str],
    ic_dirs: dict[str, str],
) -> str:
    """Return a memory-qualified production IC job for the generated configs."""
    exclude_line = (
        f"#SBATCH --exclude={args.slurm_exclude}\n" if args.slurm_exclude else ""
    )
    models = " ".join(transfer_models)
    checks = " ".join(
        ("ic_deltab", "ic_poscx", "ic_poscy", "ic_poscz", "ic_velcx", "ic_velcy", "ic_velcz")
    )
    ic_case = "\n".join(
        f'    {name}) ic_dir="{ic_dirs[name]}" ;;' for name in transfer_models
    )
    white_noise_check = ""
    if args.white_noise_file is not None:
        source = args.white_noise_file.expanduser().resolve()
        white_noise_check = f'test -r "{source}"\nsha256sum "{source}"\n'
    return f"""#!/bin/bash
#SBATCH --job-name=dmo_species_ics
#SBATCH --partition=normal
#SBATCH --nodes=1
{exclude_line}#SBATCH --ntasks={args.music_tasks}
#SBATCH --cpus-per-task={args.omp_threads}
#SBATCH --mem={args.slurm_memory}
#SBATCH --time=2-00:00:00
#SBATCH --chdir={outdir}
#SBATCH --output=logs/make-ics-%j.out
#SBATCH --error=logs/make-ics-%j.err

set -euo pipefail
export OMP_NUM_THREADS={args.omp_threads}
export OMP_STACKSIZE=256M
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact

cd {outdir}
mkdir -p logs transfers
echo "IC campaign start=$(date --iso-8601=seconds)"
sha256sum {args.music}
{white_noise_check}for model in {models}; do
  case "$model" in
{ic_case}
    *) echo "unknown IC model: $model" >&2; exit 2 ;;
  esac
  echo "IC model=$model start=$(date --iso-8601=seconds)"
  mpirun -np "${{SLURM_NTASKS:-{args.music_tasks}}}" {args.music} "music_${{model}}.conf"
  level_dir="{outdir}/${{ic_dir}}/level_{args.levelmin:03d}"
  for component in {checks}; do
    test -s "$level_dir/$component"
  done
  for diagnostic in input_powerspec.txt dump_transfer.txt; do
    if [[ -f "$diagnostic" ]]; then
      stem="${{diagnostic%.*}}"
      suffix="${{diagnostic##*.}}"
      mv "$diagnostic" "transfers/${{stem}}_${{model}}.${{suffix}}"
    fi
  done
  echo "IC model=$model bytes=$(du -sb "{outdir}/${{ic_dir}}" | awk '{{print $1}}') end=$(date --iso-8601=seconds)"
done

touch ICS_COMPLETE
echo "IC campaign end=$(date --iso-8601=seconds)"
"""


def submit_chain(models: list[str]) -> str:
    quoted = " ".join(models)
    return f"""#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
previous=""
for model in {quoted}; do
    if [[ -z "$previous" ]]; then
        result=$(sbatch "$model/run.slurm")
    else
        result=$(sbatch --dependency="afterok:$previous" "$model/run.slurm")
    fi
    jobid=${{result##* }}
    echo "$model $jobid"
    previous=$jobid
done
"""


def submit_all(models: list[str]) -> str:
    quoted = " ".join(models)
    return f"""#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
: > submitted_jobs.tsv
for model in {quoted}; do
    jobid=$(cd "$model" && sbatch --parsable run.slurm)
    printf '%s\\t%s\\n' "$model" "$jobid" | tee -a submitted_jobs.tsv
done
"""


def manual_chain(args: argparse.Namespace, models: list[str]) -> str:
    quoted = " ".join(models)
    total_ranks = args.slurm_tasks
    return f"""#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS={args.omp_threads}
export OMP_STACKSIZE=256M
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
for model in {quoted}; do
    echo "model=$model start=$(date --iso-8601=seconds)"
    (
        cd "$model"
        mpirun -np {total_ranks} {args.ramses} run.nml
    ) >"$model/manual.out" 2>"$model/manual.err"
    echo "model=$model end=$(date --iso-8601=seconds)"
done
"""


def manual_all() -> str:
    """Run IC generation and all simulations sequentially without Slurm."""
    return """#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -f ICS_COMPLETE ]]; then
    bash ./make_ics.slurm
fi
exec ./run_manual_chain.sh
"""


def background_preflight() -> str:
    """Describe the quantified species-contract approximation."""
    redshifts = SPECIES_BACKGROUND_PREFLIGHT["redshifts"]
    cde = SPECIES_BACKGROUND_PREFLIGHT["delta_h_over_h"]["cde10"]
    rvm = SPECIES_BACKGROUND_PREFLIGHT["delta_h_over_h"]["rvm001"]
    header = "| model | " + " | ".join(f"z={z:g}" for z in redshifts) + " |"
    rule = "|---|" + "---:|" * len(redshifts)

    def row(name: str, values: list[float]) -> str:
        percentages = " | ".join(f"{100.0 * value:+.5g}%" for value in values)
        return f"| {name} | {percentages} |"

    return f"""# Background species-contract preflight

The particle load follows the standard one-species DMO convention.  lagCAMB
retains `Omega_b={OMEGA_B}` and supplies the mass-weighted CDM+baryon density
and velocity transfers.  RAMSES evolves one collisionless proxy fluid with
`omega_m={OMEGA_M}` and no separate baryon component.

The lagCAMB `cde10` and `rvm001` models assign the dark-sector exchange to CDM
alone.  The DMO forward model applies the corresponding effective dynamics to
the full proxy fluid.  At fixed total matter, the table reports
`H_species_resolved/H_effective_one_fluid - 1`.

{header}
{rule}
{row("cde10", cde)}
{row("rvm001", rvm)}

The values diagnose a deliberate single-fluid dynamics approximation.  They
are not corrections to the transfer, phase, normalization, or numerical
integration.  The LCDM control uses the same total-matter loading convention.
The values were manually tabulated from direct background evaluations at the
parameters recorded in `campaign.json`.
"""


def readme(
    args: argparse.Namespace,
    model_names: list[str],
    ic_models: dict[str, str],
    transfer_diagnostics: dict[str, dict[str, float]],
) -> str:
    dx_base = args.boxlen / (2**args.levelmin)
    dx_fine_kpc = 1000.0 * args.boxlen / (2**args.levelmax)
    kny = 3.141592653589793 * (2**args.levelmin) / args.boxlen
    transfer_lines = []
    dynamics_lines = []
    for name in model_names:
        source = ic_models[name]
        suffix = ""
        if source == "lcdm" and name != "lcdm":
            reason = (
                "shared-IC triage mode"
                if args.ic_mode == "shared"
                else LCDM_TRANSFER_REASONS[name]
            )
            suffix = f" ({reason})"
        transfer_lines.append(f"- `{name}`: `{source}` lagCAMB transfer{suffix}")
        if name in SINGLE_FLUID_DYNAMICS_NOTES:
            dynamics_lines.append(
                f"- `{name}`: {SINGLE_FLUID_DYNAMICS_NOTES[name]}. "
                "The claim is limited to an effective one-fluid response."
            )
    transfer_summary = "\n".join(transfer_lines)
    dynamics_summary = (
        "\n".join(dynamics_lines)
        if dynamics_lines
        else "- No model in this campaign has a recorded one-fluid species-contract qualification."
    )
    lcdm_sigma8 = transfer_diagnostics["lcdm"]["sigma8_z0"]
    return f"""# cuRAMSES DMO benchmark campaign

This is a matched-phase DMO benchmark. The transfer mode is `{args.ic_mode}`.

- Box: {args.boxlen:g} Mpc/h
- Particle load: {2**args.levelmin}^3
- AMR: level {args.levelmin} through {args.levelmax}
- Base spacing: {dx_base:.6f} Mpc/h
- Finest cell: {dx_fine_kpc:.6f} kpc/h comoving
- Particle Nyquist wavenumber: {kny:.6f} h/Mpc
- Start: z={args.zstart:g}, 2LPT
- Seed: {args.seed}
- White-noise phase anchor: level {args.phase_anchor_level or args.levelmin}
- Primordial amplitude: A_s={A_S:.8e}
- CAMB Omega_b: {OMEGA_B:.8f}
- Species convention: single DMO fluid loaded from the mass-weighted CDM+baryon transfer
- RAMSES `omega_b=0` means that no separate baryon particle or fluid is
  evolved. Its `omega_m` still contains the full CDM+baryon density
  represented by one collisionless particle species. It does not mean that
  CAMB was run with zero baryons.
- LCDM sigma8(z=0), diagnostic only: {lcdm_sigma8:.8f}
- Models: {", ".join(model_names)}
- Scalar solver limit/tolerance: {args.scalar_iters} / {args.scalar_eps:.3e}
- Maximum fractional expansion step: {args.aexp_step_limit:.6g}

All models in this campaign use the same random seed and phases. When the
phase anchor is finer than the particle level, LagMUSIC generates the
white-noise realization at the anchor and restricts it to the particle
level. Using one common anchor across a resolution ladder preserves the
shared long-wave realization. Placing the same seed independently at each
particle level does not preserve the realization. LagMUSIC derives
`force_pnorm` directly from lagCAMB's linear P(k,zstart). The common A_s fixes
the absolute amplitude without sigma8 re-normalisation or MUSIC growth
back-scaling.
For DMO 2LPT, `dmo_velocity_source=transfer` constructs the first-order
velocity field from lagCAMB's model-specific `vtotal` transfer and retains
LagMUSIC's high-redshift 2LPT velocity term. LagMUSIC also receives the
model's CPL `w0` and `wa` values when applicable. Its background velocity
factor is consistent with the simulation namelist. In the legacy
`density_2lpt` mode, `vfact_scale` applies the model/LCDM
velocity-to-density growth ratio to the density-derived velocity field.

## Transfer source

{transfer_summary}

## Species-contract qualification

{dynamics_summary}

The file `BACKGROUND_PREFLIGHT.md` quantifies the background difference for
the qualified models.  The values are approximation diagnostics and are not
applied as corrections.

On a dedicated non-Slurm host, run `./run_manual_all.sh` to generate ICs
and execute every simulation sequentially. Under Slurm, submit
`make_ics.slurm` first.  After it creates `ICS_COMPLETE`, run
`./submit_all.sh` on grammar for independent concurrent jobs, or
`./submit_chain.sh` when the models must run sequentially. The
`./run_manual_chain.sh` script is the sequential fallback for a non-Slurm
host. Simulation outputs are written inside each model directory.
"""


def main() -> int:
    args = parse_args()
    if args.levelmax < args.levelmin:
        raise ValueError("levelmax must be >= levelmin")
    if (
        args.phase_anchor_level is not None
        and args.phase_anchor_level < args.levelmin
    ):
        raise ValueError("phase-anchor-level must be >= levelmin")
    if args.music_tasks < 1:
        raise ValueError("music-tasks must be positive")
    if args.aexp_step_limit <= 0.0:
        raise ValueError("aexp-step-limit must be positive")
    if any(z < 0.0 or z >= args.zstart for z in args.output_redshifts):
        raise ValueError("output redshifts must satisfy 0 <= z < zstart")
    if any(
        left <= right
        for left, right in zip(
            args.output_redshifts, args.output_redshifts[1:]
        )
    ):
        raise ValueError("output redshifts must be strictly decreasing")
    if args.slurm_tasks * args.omp_threads > 64:
        raise ValueError("requested Slurm CPU count exceeds one 64-core grammar node")

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    if args.white_noise_file is not None:
        args.white_noise_file = args.white_noise_file.expanduser().resolve()
        if not args.white_noise_file.is_file():
            raise FileNotFoundError(
                f"white-noise realization not found: {args.white_noise_file}"
            )
    model_names = list(dict.fromkeys(args.models))
    species_contract_models = [
        name for name in model_names if name in SINGLE_FLUID_SPECIES_CONTRACT
    ]
    ic_models = {
        name: (
            name
            if args.ic_mode == "model" and name in CAMB_MATCHED_MODELS
            else "lcdm"
        )
        for name in model_names
    }
    transfer_models = list(dict.fromkeys(["lcdm", *ic_models.values()]))
    camb = load_local_camb(args.camb_dir)
    transfer_data = {
        name: make_transfer(outdir, name, args.zstart, args.force, camb)
        for name in transfer_models
    }
    lcdm_transfer = transfer_data["lcdm"][0]
    for name, (transfer, diagnostics) in transfer_data.items():
        _, scale, scatter, maximum_deviation = velocity_growth_correction(
            transfer, lcdm_transfer
        )
        use_velocity_transfer = args.dmo_velocity_source == "transfer"
        diagnostics.update(
            {
                # LagMUSIC now consumes CAMB's vtotal directly in the DMO
                # 2LPT velocity solve.  The model/LCDM growth ratio below is
                # retained as a diagnostic, not collapsed into one scalar.
                "dmo_velocity_mode": (
                    "camb_vtotal"
                    if use_velocity_transfer
                    else "density_2lpt"
                ),
                "dmo_velocity_source": args.dmo_velocity_source,
                "dmo_velocity_transfer_enabled": use_velocity_transfer,
                "vfact_scale": 1.0 if use_velocity_transfer else scale,
                "velocity_growth_ratio_mean": scale,
                "velocity_growth_ratio_relative_scatter": scatter,
                "velocity_growth_ratio_maximum_deviation": maximum_deviation,
                "velocity_growth_ratio_k_range_h_mpc": [1.0e-3, 1.0],
                # Backward-compatible diagnostic names.
                "vfact_scale_relative_scatter": scatter,
                "vfact_scale_maximum_deviation": maximum_deviation,
                "vfact_scale_k_range_h_mpc": [1.0e-3, 1.0],
            }
        )
        diagnostics_path = transfer.with_suffix(".json")
        diagnostics_path.write_text(
            json.dumps(diagnostics, indent=2, sort_keys=True) + "\n"
        )
    transfer_diagnostics = {
        name: diagnostics for name, (_, diagnostics) in transfer_data.items()
    }
    config_paths: dict[str, Path] = {}
    ic_dirs: dict[str, str] = {}
    for name in transfer_models:
        transfer, diagnostics = transfer_data[name]
        ic_dir = "ics_common" if args.ic_mode == "shared" else f"ics_{name}"
        first_transfer_model = transfer_models[0]
        white_noise = outdir / (
            f"wnoise_{args.phase_anchor_level or args.levelmin:04d}.bin"
        )
        if args.white_noise_file is not None:
            random_source = str(args.white_noise_file)
        else:
            random_source = (
                str(white_noise)
                if args.reuse_white_noise and name != first_transfer_model
                else str(args.seed)
            )
        config_path = outdir / f"music_{name}.conf"
        write_text(
            config_path,
            music_config(
                args,
                name,
                transfer,
                diagnostics["sigma8_z0"],
                diagnostics["force_pnorm"],
                diagnostics["vfact_scale"],
                args.dmo_velocity_source,
                ic_dir,
                random_source,
            ),
            args.force,
        )
        config_paths[name] = config_path
        ic_dirs[name] = ic_dir

    for name in model_names:
        model = MODELS[name]
        ic_model = ic_models[name]
        model_dir = outdir / name
        model_dir.mkdir(parents=True, exist_ok=True)
        write_text(
            model_dir / "run.nml",
            namelist(args, model, ic_dirs[ic_model], ic_model),
            args.force,
        )
        write_text(
            model_dir / "run.slurm",
            slurm_script(args, name, model_dir),
            args.force,
            executable=True,
        )

    write_text(
        outdir / "make_ics.slurm",
        ic_slurm_script(args, outdir, transfer_models, ic_dirs),
        args.force,
        executable=True,
    )
    write_text(
        outdir / "submit_chain.sh",
        submit_chain(model_names),
        args.force,
        executable=True,
    )
    write_text(
        outdir / "submit_all.sh",
        submit_all(model_names),
        args.force,
        executable=True,
    )
    write_text(
        outdir / "run_manual_chain.sh",
        manual_chain(args, model_names),
        args.force,
        executable=True,
    )
    write_text(
        outdir / "run_manual_all.sh",
        manual_all(),
        args.force,
        executable=True,
    )
    write_text(
        outdir / "README.md",
        readme(args, model_names, ic_models, transfer_diagnostics),
        args.force,
    )
    if species_contract_models:
        write_text(
            outdir / SPECIES_BACKGROUND_PREFLIGHT["file"],
            background_preflight(),
            args.force,
        )
    model_metadata = {}
    for name in model_names:
        transfer_matches = ic_models[name] == name
        model_metadata[name] = {
            **MODELS[name],
            "blocks": effective_model_blocks(args, MODELS[name]),
            "ic_transfer_model": ic_models[name],
            "ic_transfer_matches_lagcamb_at_zstart": transfer_matches,
            "ic_transfer_match_scope": (
                "mass_weighted_total_matter_density_and_velocity"
                if transfer_matches
                else "shared_or_proxy_transfer"
            ),
        }
        if ic_models[name] != name:
            model_metadata[name]["ic_transfer_note"] = (
                "shared-IC triage mode"
                if args.ic_mode == "shared"
                else LCDM_TRANSFER_REASONS[name]
            )
        if name in SINGLE_FLUID_DYNAMICS_NOTES:
            model_metadata[name]["single_fluid_dynamics_approximation"] = (
                SINGLE_FLUID_DYNAMICS_NOTES[name]
            )
            model_metadata[name].update(SINGLE_FLUID_SPECIES_CONTRACT[name])
            model_metadata[name]["background_species_contract_preflight"] = (
                SPECIES_BACKGROUND_PREFLIGHT["file"]
            )
    metadata = {
        "boxlen_mpc_h": args.boxlen,
        "levelmin": args.levelmin,
        "levelmax": args.levelmax,
        "zstart": args.zstart,
        "seed": args.seed,
        "phase_anchor_level": args.phase_anchor_level or args.levelmin,
        "omega_m": OMEGA_M,
        "omega_b": OMEGA_B,
        "species_convention": "single_fluid_total_matter",
        "H0": H0,
        "n_s": N_S,
        "A_s": A_S,
        "amplitude_normalization": "lagCAMB P(k,zstart) via LagMUSIC force_pnorm",
        "velocity_normalization": (
            "lagCAMB model-specific vtotal transfer with LagMUSIC "
            "model-background vfact and vfact_scale=1"
            if args.dmo_velocity_source == "transfer"
            else
            "density-derived 2LPT velocity with model/LCDM growth ratio "
            "via LagMUSIC vfact_scale"
        ),
        "ic_mode": args.ic_mode,
        "shared_lcdm_ic": args.ic_mode == "shared",
        "transfer_diagnostics": transfer_diagnostics,
        "scalar_iters": args.scalar_iters,
        "scalar_eps": args.scalar_eps,
        "aexp_step_limit": args.aexp_step_limit,
        "output_redshifts": args.output_redshifts,
        "models": model_metadata,
        "camb_dir": str(args.camb_dir.expanduser().resolve()),
        "camb_module": str(Path(camb.__file__).resolve()),
        "camb_version": getattr(camb, "__version__", "unknown"),
        "music": args.music,
        "music_tasks": args.music_tasks,
        "music_unigrid_slab": args.music_unigrid_slab,
        "reuse_white_noise": args.reuse_white_noise or args.white_noise_file is not None,
        "white_noise_file": (
            str(args.white_noise_file) if args.white_noise_file is not None else None
        ),
        "ramses": args.ramses,
        "species_contract_preflight": (
            SPECIES_BACKGROUND_PREFLIGHT if species_contract_models else None
        ),
        "provenance": {
            "lagcamb": {
                "path": str(args.camb_dir.expanduser().resolve()),
                "git_revision": git_revision(args.camb_dir),
                "tracked_worktree_dirty": git_tracked_worktree_dirty(args.camb_dir),
                "omega_b": OMEGA_B,
                "omega_b_meaning": "separate baryon component in the linear calculation",
            },
            "lagmusic": {
                "binary": str(Path(args.music).expanduser().resolve()),
                "binary_sha256": file_sha256(Path(args.music)),
                "source_git_revision": git_revision(Path(args.music)),
                "source_tracked_worktree_dirty": git_tracked_worktree_dirty(
                    Path(args.music)
                ),
            },
            "curamses": {
                "binary": str(Path(args.ramses).expanduser().resolve()),
                "binary_sha256": file_sha256(Path(args.ramses)),
                "omega_b": 0.0,
                "omega_b_meaning": "no separately evolved baryon species. omega_m is the full proxy fluid",
            },
            "campaign_generator": {
                "path": str(Path(__file__).resolve()),
                "source_git_revision": git_revision(Path(__file__)),
                "source_sha256": file_sha256(Path(__file__)),
                "working_tree_dirty_at_generation": git_tracked_worktree_dirty(
                    Path(__file__)
                ),
            },
        },
    }
    write_text(
        outdir / "campaign.json",
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        args.force,
    )

    if args.make_ics:
        music = Path(args.music)
        if not music.is_file():
            raise FileNotFoundError(f"MUSIC binary not found: {music}")
        for name in transfer_models:
            ic_level = outdir / ic_dirs[name] / f"level_{args.levelmin:03d}"
            if ic_level.exists() and not args.force:
                print(f"IC already exists: {ic_level}")
            else:
                subprocess.run(
                    [
                        "mpirun",
                        "-np",
                        str(args.music_tasks),
                        str(music),
                        config_paths[name].name,
                    ],
                    cwd=outdir,
                    check=True,
                )
                for diagnostic_name in ("input_powerspec.txt", "dump_transfer.txt"):
                    diagnostic = outdir / diagnostic_name
                    if diagnostic.exists():
                        suffix = Path(diagnostic_name).suffix
                        archived = (
                            outdir
                            / "transfers"
                            / f"{Path(diagnostic_name).stem}_{name}{suffix}"
                        )
                        diagnostic.replace(archived)

    print(f"campaign ready: {outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
