"""Build a provenance-enforced, UVB-free metal-only SNRT thermal atlas."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import sys

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.thermal_atlas import ThermalAtlas, _equilibrium_log_temperature, write_thermal_atlas


GRACKLE_VERSION = "3.4.2-dev"
GRACKLE_REVISION = "f93091ff8456962d7017a5bff7472945a30e3dad"
GRACKLE_DATA_REVISION = "928696482fbe15d9bac4382de6134d95568f099c"
GRACKLE_DATA_REPOSITORY = "https://github.com/grackle-project/grackle_data_files"
SOURCE_COOLING_DATASET = "CoolingRates/Metals/Cooling"


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_scale_factors(path: str | Path) -> np.ndarray:
    values = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            values.append(float(line))
    scale_factor = np.asarray(values, dtype=np.float64)
    if (
        len(scale_factor) < 2
        or not np.isfinite(scale_factor).all()
        or np.any(scale_factor <= 0.0)
        or np.any(np.diff(scale_factor) <= 0.0)
    ):
        raise ValueError("scale-factor file requires at least two strictly increasing positive entries")
    return scale_factor


def read_uvb_free_metal_cooling(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the no-UVB Cloudy metal coefficient table and reject UVB products."""

    with h5py.File(Path(path), "r") as handle:
        if "UVBRates" in handle:
            raise ValueError("metal atlas source contains UV background rates")
        if "CoolingRates/Metals/Heating" in handle:
            raise ValueError("metal atlas source contains a metal heating table")
        if SOURCE_COOLING_DATASET not in handle:
            raise ValueError(f"metal atlas source lacks {SOURCE_COOLING_DATASET}")
        dataset = handle[SOURCE_COOLING_DATASET]
        rank = int(dataset.attrs.get("Rank", -1))
        parameter_name = dataset.attrs.get("Parameter1_Name", b"")
        if isinstance(parameter_name, bytes):
            parameter_name = parameter_name.decode()
        if rank != 2 or parameter_name != "hden":
            raise ValueError("UVB-free metal cooling must be a rank-2 (hden, temperature) table")
        log_density = np.asarray(dataset.attrs.get("Parameter1"), dtype=np.float64)
        temperature = np.asarray(dataset.attrs.get("Temperature"), dtype=np.float64)
        coefficient = np.asarray(dataset, dtype=np.float64)
    if coefficient.shape != (len(log_density), len(temperature)):
        raise ValueError("metal cooling dataset shape does not match its axes")
    if (
        len(log_density) < 2
        or len(temperature) < 2
        or np.any(np.diff(log_density) <= 0.0)
        or np.any(np.diff(temperature) <= 0.0)
        or not np.isfinite(coefficient).all()
        or np.any(coefficient <= 0.0)
    ):
        raise ValueError("metal cooling table has invalid axes or coefficients")
    return log_density, np.log10(temperature), coefficient


def build_atlas(
    source_data: str | Path,
    scale_factor: np.ndarray,
    *,
    generated_utc: str,
) -> ThermalAtlas:
    """Convert Cloudy coefficients to a solar-metallicity volumetric table.

    Grackle revision ``f93091f`` skips the CMB subtraction above 100 T_CMB as
    an optimization. SNRT applies it continuously: this removes the finite
    cutoff step while preserving the same CMB-temperature zero and becoming
    negligible at high temperature. Metallicity is not a table dimension; the
    runtime applies the source model's exact linear Z/Z_solar multiplier.
    """

    log_density, log_temperature, cooling_coefficient = read_uvb_free_metal_cooling(source_data)
    density = 10.0**log_density
    rate = np.empty(
        (len(scale_factor), len(log_density), len(log_temperature)),
        dtype=np.float64,
    )
    log_coefficient = np.log10(cooling_coefficient)
    for time_index, expansion in enumerate(scale_factor):
        cmb_temperature = 2.73 / expansion
        log_cmb_temperature = np.log10(cmb_temperature)
        cmb_coefficient = np.asarray(
            [
                10.0 ** np.interp(log_cmb_temperature, log_temperature, row)
                for row in log_coefficient
            ]
        )
        effective_coefficient = -cooling_coefficient + cmb_coefficient[:, None]
        rate[time_index] = effective_coefficient * density[:, None] ** 2

    shape = rate.shape
    neutral_primordial_mu = 1.0 / (0.76 + 0.24 / 4.0)
    mean_molecular_weight = np.full(shape, neutral_primordial_mu, dtype=np.float64)
    equilibrium_log_temperature = _equilibrium_log_temperature(log_temperature, rate)
    generator = Path(__file__).resolve()
    provenance = {
        "thermal_component": "metal_only",
        "primordial_rates_included": "false",
        "uv_background_included": "false",
        "photoheating_included": "false",
        "metallicity_scaling": "linear_z_solar",
        "metallicity_application": "analytic_runtime_multiplier",
        "rate_sign_convention": "heating_positive_cooling_negative",
        "source_data_name": Path(source_data).name,
        "source_data_sha256": sha256(source_data),
        "source_cooling_dataset": SOURCE_COOLING_DATASET,
        "source_repository": GRACKLE_DATA_REPOSITORY,
        "source_repository_revision": GRACKLE_DATA_REVISION,
        "generator_name": generator.name,
        "generator_version": "2",
        "generator_sha256": sha256(generator),
        "grackle_version": GRACKLE_VERSION,
        "grackle_revision": GRACKLE_REVISION,
        "generated_utc": generated_utc,
        "cmb_metal_floor": "continuous_subtraction_at_tcmb",
        "cmb_floor_reference_deviation": "grackle_2dex_cutoff_removed_for_continuity",
        "mean_molecular_weight_model": "neutral_primordial_initialization_only",
        "equilibrium_temperature_model": "metal_cooling_cmb_floor_only",
    }
    return ThermalAtlas(
        scale_factor=scale_factor,
        log_hydrogen_number_density_cm3=log_density,
        log_temperature_k=log_temperature,
        net_rate_erg_s_cm3=rate,
        mean_molecular_weight=mean_molecular_weight,
        equilibrium_log_temperature_k=equilibrium_log_temperature,
        provenance=provenance,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-data", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--scale-factors", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--generated-utc",
        default=datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        help="explicit UTC provenance timestamp; defaults to the current time",
    )
    args = parser.parse_args()
    source_data = Path(args.source_data).resolve()
    output = Path(args.output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite thermal atlas: {output}")
    actual_sha256 = sha256(source_data)
    if actual_sha256 != args.expected_source_sha256:
        raise ValueError(
            f"source cooling-table checksum mismatch: expected {args.expected_source_sha256}, got {actual_sha256}"
        )
    generated_utc = datetime.fromisoformat(args.generated_utc.replace("Z", "+00:00"))
    if generated_utc.tzinfo is None:
        raise ValueError("generated-utc must include a timezone")
    atlas = build_atlas(
        source_data,
        read_scale_factors(args.scale_factors),
        generated_utc=generated_utc.astimezone(timezone.utc).replace(microsecond=0).isoformat(),
    )
    write_thermal_atlas(output, atlas)
    print(
        "METAL_THERMAL_ATLAS_OK "
        f"output={output} sha256={sha256(output)} source_sha256={actual_sha256} "
        f"shape={atlas.net_rate_erg_s_cm3.shape}"
    )


if __name__ == "__main__":
    main()
