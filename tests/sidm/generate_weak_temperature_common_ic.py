#!/usr/bin/env python3
"""Generate RAMSES ASCII ICs with x/v RNG draws aligned to the DSMC reference.

The x and velocity RNG stream follows weak_temperature_dsmc.cpp's draw
sequence and call order; cross-language floating-point values are validated
within 8 machine eps rather than bitwise.
Transverse coordinates use a separate SplitMix64 stream so adding y,z for
RAMSES cannot shift the reference x/v draws for later particle IDs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

MASK64 = (1 << 64) - 1
GAMMA64 = 0x9E3779B97F4A7C15
MIX1 = 0xBF58476D1CE4E5B9
MIX2 = 0x94D049BB133111EB
PI = 3.141592653589793238462643383279502884
N_DEFAULT = 131_072
SEED_DEFAULT = 20261027
TRANSVERSE_STREAM_XOR = 0xD1B54A32D192ED03
RHO0 = 1.0
T0 = 1.0
EPSILON = 0.02
K_LAMBDA = 0.25
SIGMA_OVER_M = 1.0
LAMBDA0 = 1.0 / (RHO0 * SIGMA_OVER_M)
BOX = 2.0 * PI / (K_LAMBDA / LAMBDA0)
WAVE_NUMBER = K_LAMBDA / LAMBDA0
REFERENCE_SOURCE_SHA256 = "74b4364c40e8503a75ce7bda68042167aa31077d7d9ff9101b235f01584a55e5"


class SplitMix64:
    def __init__(self, seed: int) -> None:
        self.state = seed & MASK64
        self.has_spare = False
        self.spare = 0.0

    def integer(self) -> int:
        self.state = (self.state + GAMMA64) & MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * MIX1) & MASK64
        z = ((z ^ (z >> 27)) * MIX2) & MASK64
        return z ^ (z >> 31)

    def uniform(self) -> float:
        q = self.integer() >> 13
        return (float(q) + 0.5) * 2.0**-51

    def normal(self) -> float:
        if self.has_spare:
            self.has_spare = False
            return self.spare
        radius = math.sqrt(-2.0 * math.log(self.uniform()))
        angle = 2.0 * PI * self.uniform()
        self.spare = radius * math.sin(angle)
        self.has_spare = True
        return radius * math.cos(angle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_seed(value: str) -> int:
    seed = int(value, 10)
    if not 0 <= seed <= MASK64:
        raise argparse.ArgumentTypeError("seed must be an unsigned 64-bit integer")
    return seed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="new directory; must not exist")
    parser.add_argument("--particles", type=int, default=N_DEFAULT)
    parser.add_argument("--seed", type=parse_seed, default=SEED_DEFAULT)
    parser.add_argument("--energy-matched-control", action="store_true")
    args = parser.parse_args()
    if args.particles < 2:
        parser.error("--particles must be at least 2")
    return args


def generate(output_dir: Path, particle_count: int, seed: int,
             energy_matched_control: bool) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=False)
    part_path = output_dir / "ic_part"
    requested_norm = 1.0 / math.sqrt(1.0 - EPSILON * EPSILON)
    initial_epsilon = 0.0 if energy_matched_control else EPSILON
    initial_temperature = T0 / requested_norm if energy_matched_control else T0
    density_norm = 1.0 / math.sqrt(1.0 - initial_epsilon * initial_epsilon)
    particle_mass = RHO0 * BOX**3 / particle_count
    phase_rng = SplitMix64(seed)
    transverse_rng = SplitMix64(seed ^ TRANSVERSE_STREAM_XOR)

    with part_path.open("x", encoding="ascii", newline="\n") as stream:
        for _particle_id in range(1, particle_count + 1):
            while True:
                x = BOX * phase_rng.uniform()
                cosine = math.cos(WAVE_NUMBER * x)
                density_factor = 1.0 / (
                    density_norm * (1.0 + initial_epsilon * cosine)
                )
                acceptance = (
                    density_factor * density_norm * (1.0 - abs(initial_epsilon))
                )
                if phase_rng.uniform() <= acceptance:
                    break

            local_temperature = initial_temperature * (
                1.0 + initial_epsilon * cosine
            )
            sigma_v = math.sqrt(local_temperature)
            velocity = (sigma_v * phase_rng.normal(),
                        sigma_v * phase_rng.normal(),
                        sigma_v * phase_rng.normal())
            y = BOX * transverse_rng.uniform()
            z = BOX * transverse_rng.uniform()
            centered = (x - BOX / 2.0, y - BOX / 2.0, z - BOX / 2.0)
            row = (*centered, *velocity, particle_mass)
            stream.write(" ".join(format(value, ".17e") for value in row))
            stream.write("\n")

    current = Path(__file__).resolve()
    profile = "energy_matched_uniform_control" if energy_matched_control else "perturbed"
    manifest: dict[str, object] = {
        "status": "GENERATED_REFERENCE_MATCHED_COMMON_ASCII_IC",
        "generator": current.name,
        "generator_sha256": sha256(current),
        "reference_initializer": "weak_temperature_dsmc.cpp::initialize",
        "reference_source_sha256": REFERENCE_SOURCE_SHA256,
        "python_version": platform.python_version(),
        "phase_space_rng": "SplitMix64(seed); reference x, rejection, and vx/vy/vz call order",
        "transverse_rng": "SplitMix64(seed XOR 0xD1B54A32D192ED03); independent y,z draws",
        "ic_seed": seed,
        "particle_count": particle_count,
        "profile": profile,
        "id_rule": "one-based ASCII input line number",
        "ascii_schema": ["x_centered", "y_centered", "z_centered", "vx", "vy", "vz", "mass"],
        "coordinate_rule": "ASCII coordinates are centered; RAMSES and DSMC readers add box_length/2",
        "box_length": BOX,
        "wave_number": WAVE_NUMBER,
        "mean_free_path": LAMBDA0,
        "k_lambda": K_LAMBDA,
        "density_mean": RHO0,
        "temperature_mean_parameter": T0,
        "epsilon": initial_epsilon,
        "density_profile": "rho0/[A*(1+epsilon*cos(k*x))]" if not energy_matched_control else "rho0",
        "temperature_profile": "T0*(1+epsilon*cos(k*x))" if not energy_matched_control else "T0*sqrt(1-epsilon^2)",
        "transverse_sampling": "independent uniform y,z on [0,box_length) from separate RNG stream",
        "velocity_sampling": "reference initializer order; conditional Gaussian vx,vy,vz with variance T(x)",
        "recentering_or_rescaling": False,
        "antithetic_pairing": False,
        "represented_particle_mass": particle_mass,
        "total_represented_mass": particle_count * particle_mass,
        "ic_part_bytes": part_path.stat().st_size,
        "ic_part_sha256": sha256(part_path),
        "reference_match_validation": "NOT_RUN",
        "scientific_admission": "NOT_EVALUATED",
    }
    manifest_path = output_dir / "ic_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="ascii")
    return manifest


def main() -> int:
    args = parse_args()
    try:
        manifest = generate(args.output_dir.resolve(), args.particles, args.seed,
                            args.energy_matched_control)
    except (FileExistsError, OSError, RuntimeError, ValueError) as exc:
        print(f"IC generation failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
