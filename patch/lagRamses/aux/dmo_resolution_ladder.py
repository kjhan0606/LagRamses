#!/usr/bin/env python3
"""Generate a uniform DMO resolution ladder with consistent physics settings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys


MODELS = ("lcdm", "f5", "f6", "n1", "n5", "sym_a")
RESOURCE_PRESETS = {
    6: {"slurm_tasks": 4, "omp_threads": 4, "memory": "32G", "music_tasks": 1},
    7: {"slurm_tasks": 8, "omp_threads": 4, "memory": "64G", "music_tasks": 2},
    8: {"slurm_tasks": 4, "omp_threads": 4, "memory": "128G", "music_tasks": 8},
    9: {"slurm_tasks": 8, "omp_threads": 2, "memory": "256G", "music_tasks": 16},
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--levels", nargs="+", type=int, default=[6, 7, 8])
    parser.add_argument(
        "--models",
        nargs="+",
        choices=MODELS,
        default=list(MODELS),
        help="model subset to generate; defaults to the full validation set",
    )
    parser.add_argument("--boxlen", type=float, default=64.0)
    parser.add_argument("--zstart", type=float, default=49.0)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument(
        "--phase-anchor-level",
        type=int,
        help=(
            "common white-noise level for all resolutions; defaults to the "
            "finest requested level"
        ),
    )
    parser.add_argument("--aexp-step-limit", type=float, default=0.1)
    parser.add_argument(
        "--output-redshifts",
        nargs="+",
        type=float,
        default=[5.0, 2.0, 1.0, 0.8, 0.5, 0.2, 0.0],
    )
    parser.add_argument("--scalar-iters", type=int, default=6000)
    parser.add_argument("--scalar-eps", type=float, default=1.0e-6)
    parser.add_argument(
        "--ramses",
        default="/home/kjhan/BACKUP/lagRamses-de-nonstd/bin/ramses_cosmo_model_test3d",
    )
    parser.add_argument("--make-ics", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def capacities(level: int) -> tuple[int, int]:
    particles = 2 ** (3 * level)
    # Spatial domain ownership can concentrate particles strongly on a subset
    # of ranks at late times even in a uniform-grid run.  Four slots per
    # particle avoids rank-local exhaustion without changing the dynamics.
    nparttot = max(400_000, 4 * particles)
    # ngridtot is distributed evenly over ranks, while the initial k-section
    # ownership can be substantially imbalanced.  Reserve the same 25% margin
    # as particles against the full cell load; oct-count-based sizing is too
    # small for multi-rank initial-grid construction.
    ngridtot = max(400_000, (5 * particles + 3) // 4)
    return ngridtot, nparttot


def ic_slurm_script(
    root: Path,
    campaign: Path,
    level: int,
    resources: dict[str, object],
    boxlen: float,
    zstart: float,
    seed: int,
    aexp_step_limit: float,
    phase_anchor_level: int,
    ramses: str,
    models: list[str],
    scalar_iters: int,
    scalar_eps: float,
    output_redshifts: list[float],
) -> str:
    model_args = " ".join(models)
    output_args = " ".join(f"{redshift:.15g}" for redshift in output_redshifts)
    return f"""#!/bin/bash
#SBATCH --job-name=ic_L{level}_{2**level}
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks={resources["music_tasks"]}
#SBATCH --cpus-per-task={resources["omp_threads"]}
#SBATCH --mem={resources["memory"]}
#SBATCH --time=1-00:00:00
#SBATCH --chdir={campaign}
#SBATCH --output=make-ics-%j.out
#SBATCH --error=make-ics-%j.err

set -euo pipefail
export OMP_NUM_THREADS={resources["omp_threads"]}
export OMP_STACKSIZE=256M
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
{sys.executable} {Path(__file__).resolve()} --root {root} --levels {level} \
  --models {model_args} \
  --boxlen {boxlen:.15g} --zstart {zstart:.15g} --seed {seed} \
  --aexp-step-limit {aexp_step_limit:.8g} \
  --scalar-iters {scalar_iters} --scalar-eps {scalar_eps:.8e} \
  --output-redshifts {output_args} \
  --phase-anchor-level {phase_anchor_level} --ramses {ramses} --make-ics
"""


def main() -> int:
    args = arguments()
    unsupported = [level for level in args.levels if level not in RESOURCE_PRESETS]
    if unsupported:
        raise ValueError(f"no resource preset for levels: {unsupported}")
    if args.aexp_step_limit <= 0.0:
        raise ValueError("aexp-step-limit must be positive")
    if args.scalar_iters <= 0:
        raise ValueError("scalar-iters must be positive")
    if args.scalar_eps <= 0.0:
        raise ValueError("scalar-eps must be positive")
    phase_anchor_level = (
        args.phase_anchor_level
        if args.phase_anchor_level is not None
        else max(args.levels)
    )
    if phase_anchor_level < max(args.levels):
        raise ValueError(
            "phase-anchor-level must be >= the finest requested level"
        )

    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    setup = Path(__file__).with_name("dmo_benchmark_setup.py")
    records = []
    for level in args.levels:
        resources = dict(RESOURCE_PRESETS[level])
        # LagMUSIC's slab RNG is decomposition-independent at the phase-anchor
        # level, but restricting a finer anchor onto a coarser particle level
        # is not currently reproducible with more than one MUSIC rank. These
        # IC jobs are short, so use the serial RNG path below the anchor until
        # the distributed restriction path is made phase preserving.
        if level < phase_anchor_level:
            resources["music_tasks"] = 1
        ngridtot, nparttot = capacities(level)
        campaign = root / f"L{level}_{2**level:03d}"
        command = [
            sys.executable,
            str(setup),
            "--outdir",
            str(campaign),
            "--models",
            *args.models,
            "--boxlen",
            str(args.boxlen),
            "--levelmin",
            str(level),
            "--levelmax",
            str(level),
            "--zstart",
            str(args.zstart),
            "--seed",
            str(args.seed),
            "--phase-anchor-level",
            str(phase_anchor_level),
            "--ngridtot",
            str(ngridtot),
            "--nparttot",
            str(nparttot),
            "--ramses",
            args.ramses,
            "--music-tasks",
            str(resources["music_tasks"]),
            "--slurm-tasks",
            str(resources["slurm_tasks"]),
            "--omp-threads",
            str(resources["omp_threads"]),
            "--slurm-memory",
            resources["memory"],
            "--scalar-iters",
            str(args.scalar_iters),
            "--scalar-eps",
            str(args.scalar_eps),
            "--aexp-step-limit",
            str(args.aexp_step_limit),
            "--output-redshifts",
            *(str(redshift) for redshift in args.output_redshifts),
            "--ic-mode",
            "model",
        ]
        if args.make_ics:
            command.append("--make-ics")
        if args.force:
            command.append("--force")
        subprocess.run(command, check=True)
        ic_script = campaign / "make_ics.slurm"
        ic_script.write_text(
            ic_slurm_script(
                root,
                campaign,
                level,
                resources,
                args.boxlen,
                args.zstart,
                args.seed,
                args.aexp_step_limit,
                phase_anchor_level,
                args.ramses,
                args.models,
                args.scalar_iters,
                args.scalar_eps,
                args.output_redshifts,
            )
        )
        ic_script.chmod(ic_script.stat().st_mode | 0o110)
        records.append(
            {
                "level": level,
                "particle_load": 2**level,
                "cell_size_mpc_h": args.boxlen / 2**level,
                "particle_nyquist_h_mpc": 3.141592653589793 * 2**level / args.boxlen,
                "ngridtot": ngridtot,
                "nparttot": nparttot,
                "scalar_iters": args.scalar_iters,
                "scalar_eps": args.scalar_eps,
                "aexp_step_limit": args.aexp_step_limit,
                "phase_anchor_level": phase_anchor_level,
                "campaign": str(campaign),
                **resources,
            }
        )

    metadata = {
        "boxlen_mpc_h": args.boxlen,
        "zstart": args.zstart,
        "seed": args.seed,
        "phase_anchor_level": phase_anchor_level,
        "aexp_step_limit": args.aexp_step_limit,
        "output_redshifts": args.output_redshifts,
        "models": args.models,
        "levels": records,
    }
    (root / "resolution_ladder.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
