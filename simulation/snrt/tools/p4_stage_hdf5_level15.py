"""Conservatively resample RAMSES HDF5 AMR leaves onto the fixed P4 mesh.

The HDF5 writer exports raw ``uold_N`` arrays, but their physical meaning is
only known from the RAMSES hydro descriptor.  This adapter therefore requires
an explicit JSON field map.  Every mapped field is either volume averaged or
mass weighted, and constants are retained in the metadata as non-production
fallbacks rather than being mistaken for measured snapshot fields.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.snapshot import GridSpec, neutral_primordial_input, write_static_rt_input
from snrt_core.source_ledger import read_photon_source_ledger_csv
from snrt_core.thermal_atlas import read_thermal_atlas


_INFO_VALUE = re.compile(r"^\s*([A-Za-z_]+)\s*=\s*([^!#]+)")
_CHILD_BITS = np.asarray(
    tuple((child & 1, (child >> 1) & 1, (child >> 2) & 1) for child in range(8)), dtype=np.float64
)
_FIELD_UNITS = {
    "code_density",
    "code_momentum_density",
    "code_energy_density",
    "code_velocity",
    "code_pressure",
    "cgs",
    "dimensionless",
    "mass_fraction",
    "metallicity_solar",
}
_AVERAGING_RULES = {"volume", "mass_weighted"}
_SCALAR_FIELDS = (
    "thermal_pressure",
    "total_energy_density",
    "temperature_k",
    "metal_density",
    "metallicity_solar",
    "dust_to_metal",
    "dust_relative_abundance",
    "x_hii",
    "x_heii",
    "x_heiii",
    "x_h2",
)
_PRODUCTION_DATASET_FIELDS = (
    "metallicity_solar",
    "dust_to_metal",
    "x_h2",
    "x_hii",
    "x_heii",
    "x_heiii",
)


def _read_info(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text().splitlines():
        match = _INFO_VALUE.match(line)
        if match is not None:
            try:
                values[match.group(1)] = float(match.group(2).strip().replace("D", "E"))
            except ValueError:
                continue
    required = {"aexp", "unit_l", "unit_d"}
    missing = required - set(values)
    if missing:
        raise ValueError(f"RAMSES info file missing values: {sorted(missing)}")
    if values["unit_l"] <= 0.0 or values["unit_d"] <= 0.0:
        raise ValueError("RAMSES unit_l and unit_d must be positive")
    return values


def _read_field_map(path: Path) -> dict[str, object]:
    """Read and validate the explicit HDF5-to-SNRT variable map."""

    payload = json.loads(path.read_text())
    if payload.get("schema") != "snrt-hdf5-field-map" or payload.get("schema_version") != 1:
        raise ValueError("unsupported SNRT HDF5 field-map schema")
    fields = payload.get("fields")
    if not isinstance(fields, dict):
        raise ValueError("HDF5 field map requires a fields object")

    def validate_source_spec(name: str, spec: object, *, required: bool = False) -> None:
        if spec is None:
            if required:
                raise ValueError(f"HDF5 field map requires {name}")
            return
        if not isinstance(spec, dict):
            raise ValueError(f"field map entry {name} must be an object or null")
        has_dataset = "dataset" in spec
        has_constant = "constant" in spec
        if has_dataset == has_constant:
            raise ValueError(f"field map entry {name} must contain exactly one of dataset or constant")
        if has_dataset and (not isinstance(spec["dataset"], str) or not spec["dataset"]):
            raise ValueError(f"field map entry {name} has an invalid dataset")
        if has_constant:
            try:
                constant = float(spec["constant"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"field map entry {name} has an invalid constant") from exc
            if not np.isfinite(constant):
                raise ValueError(f"field map entry {name} has a non-finite constant")
            if not str(spec.get("reason", "")).strip():
                raise ValueError(f"field map constant {name} requires a reason")
        if "unit" not in spec or spec["unit"] not in _FIELD_UNITS:
            raise ValueError(f"field map entry {name} requires a supported explicit unit")
        if "averaging" not in spec or spec["averaging"] not in _AVERAGING_RULES:
            raise ValueError(f"field map entry {name} requires volume or mass_weighted averaging")
        if spec["unit"] == "mass_fraction":
            try:
                solar_mass_fraction = float(spec["solar_mass_fraction"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"mass_fraction field {name} requires solar_mass_fraction") from exc
            if solar_mass_fraction <= 0.0:
                raise ValueError(f"mass_fraction field {name} requires positive solar_mass_fraction")

    def validate_derived_spec(name: str, spec: object) -> None:
        if not isinstance(spec, dict) or "derive" not in spec:
            raise ValueError(f"derived field map entry {name} requires a derive rule")
        if not isinstance(spec["derive"], str) or not spec["derive"]:
            raise ValueError(f"derived field map entry {name} has an invalid derive rule")
        depends_on = spec.get("depends_on")
        if not isinstance(depends_on, list) or not depends_on or not all(
            isinstance(value, str) and value for value in depends_on
        ):
            raise ValueError(f"derived field map entry {name} requires a non-empty depends_on list")
        if name == "thermal_pressure":
            if spec["derive"] != "ideal_gas_pressure_from_conservative":
                raise ValueError("thermal_pressure has an unsupported derive rule")
            if "gamma" in spec:
                try:
                    gamma = float(spec["gamma"])
                except (TypeError, ValueError) as exc:
                    raise ValueError("derived thermal_pressure gamma must be numeric") from exc
                if gamma <= 1.0:
                    raise ValueError("derived thermal_pressure gamma must be greater than one")
        elif name == "metallicity_solar":
            if spec["derive"] != "metallicity_from_density_and_metal_density":
                raise ValueError("metallicity_solar has an unsupported derive rule")
            try:
                solar_mass_fraction = float(spec["solar_mass_fraction"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("derived metallicity_solar requires solar_mass_fraction") from exc
            if solar_mass_fraction <= 0.0:
                raise ValueError("derived metallicity_solar requires positive solar_mass_fraction")
        else:
            raise ValueError(f"derived field map entry {name} is not supported")

    def validate_spec(name: str, spec: object, *, required: bool = False) -> None:
        if spec is None:
            if required:
                raise ValueError(f"HDF5 field map requires {name}")
            return
        if isinstance(spec, dict) and "derive" in spec:
            if required:
                raise ValueError(f"HDF5 field map requires a source dataset for {name}")
            validate_derived_spec(name, spec)
            return
        validate_source_spec(name, spec, required=required)

    density = fields.get("density")
    validate_spec("density", density, required=True)
    if density["unit"] != "code_density":
        raise ValueError("density must be mapped in code_density units")
    velocity = fields.get("velocity")
    if velocity is not None:
        if not isinstance(velocity, list) or len(velocity) != 3:
            raise ValueError("velocity field map entry must contain three components")
        for index, spec in enumerate(velocity):
            validate_spec(f"velocity[{index}]", spec, required=True)
            if spec["unit"] == "code_momentum_density" and spec.get("quantity") != "momentum_density":
                raise ValueError(
                    f"velocity[{index}] code_momentum_density entries require quantity=momentum_density"
                )
    for name in _SCALAR_FIELDS:
        validate_spec(name, fields.get(name))
    if fields.get("thermal_pressure") is not None and fields.get("temperature_k") is not None:
        raise ValueError("field map must choose either thermal_pressure or temperature_k")
    return payload


def _field_unit_scale(spec: dict[str, object], info: dict[str, float]) -> float:
    unit = spec["unit"]
    if unit in {"cgs", "dimensionless", "metallicity_solar"}:
        return 1.0
    if unit == "code_density":
        return info["unit_d"]
    if unit == "code_momentum_density":
        if "unit_t" not in info or info["unit_t"] <= 0.0:
            raise ValueError("code_momentum_density requires a positive unit_t in the RAMSES info file")
        return info["unit_d"] * info["unit_l"] / info["unit_t"]
    if unit == "code_energy_density":
        if "unit_t" not in info or info["unit_t"] <= 0.0:
            raise ValueError("code_energy_density requires a positive unit_t in the RAMSES info file")
        return info["unit_d"] * (info["unit_l"] / info["unit_t"]) ** 2
    if unit == "code_velocity":
        if "unit_t" not in info or info["unit_t"] <= 0.0:
            raise ValueError("code_velocity field requires a positive unit_t in the RAMSES info file")
        return info["unit_l"] / info["unit_t"]
    if unit == "code_pressure":
        if "unit_t" not in info or info["unit_t"] <= 0.0:
            raise ValueError("code_pressure field requires a positive unit_t in the RAMSES info file")
        return info["unit_d"] * (info["unit_l"] / info["unit_t"]) ** 2
    if unit == "mass_fraction":
        return 1.0 / float(spec["solar_mass_fraction"])
    raise ValueError(f"unsupported field-map unit {unit!r}")


def _read_grid_blocks(dataset: h5py.Dataset, grid_indices: np.ndarray) -> np.ndarray:
    """Read eight RAMSES cell values per selected grid without a full read."""

    result = np.empty((len(grid_indices), 8), dtype=np.float64)
    if dataset.ndim != 1:
        raise ValueError(f"RAMSES HDF5 dataset {dataset.name} must be one-dimensional")
    if len(dataset) < (int(grid_indices[-1]) + 1) * 8:
        raise ValueError(f"RAMSES HDF5 dataset {dataset.name} is shorter than the AMR grid table")
    starts = np.r_[0, np.flatnonzero(np.diff(grid_indices) != 1) + 1]
    stops = np.r_[starts[1:], len(grid_indices)]
    for begin, end in zip(starts, stops, strict=True):
        first_grid = int(grid_indices[begin])
        last_grid = int(grid_indices[end - 1]) + 1
        result[begin:end] = np.asarray(dataset[first_grid * 8 : last_grid * 8]).reshape((-1, 8))
    return result


def _read_hdf5_scalar_attribute(handle: h5py.File, group_name: str, attribute: str) -> float | None:
    """Read a scalar HDF5 attribute without assuming its storage shape."""

    if group_name not in handle or attribute not in handle[group_name].attrs:
        return None
    value = np.asarray(handle[group_name].attrs[attribute])
    if value.size != 1:
        raise ValueError(f"HDF5 attribute {group_name}/{attribute} must be scalar")
    result = float(value.reshape(-1)[0])
    if not np.isfinite(result):
        raise ValueError(f"HDF5 attribute {group_name}/{attribute} is non-finite")
    return result


def _snapshot_backed_field(fields: dict[str, object], name: str, seen: set[str] | None = None) -> bool:
    """Return whether a field is a dataset or a fully dataset-backed derivation."""

    seen = set() if seen is None else seen
    if name in seen:
        return False
    seen.add(name)
    if name == "velocity":
        specs = fields.get("velocity")
        return isinstance(specs, list) and len(specs) == 3 and all(
            isinstance(spec, dict) and "dataset" in spec for spec in specs
        )
    spec = fields.get(name)
    if not isinstance(spec, dict):
        return False
    if "dataset" in spec:
        return True
    if "derive" not in spec:
        return False
    dependencies = spec.get("depends_on")
    if not isinstance(dependencies, list) or not dependencies:
        return False
    return all(_snapshot_backed_field(fields, dependency, seen.copy()) for dependency in dependencies)


def _preflight_hdf5_snapshot(
    path: Path,
    *,
    fields: dict[str, object],
    source_specs: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Validate HDF5 groups and mapped dataset lengths without reading field values."""

    required_datasets = {fields["density"]["dataset"]}
    required_datasets.update(
        spec["dataset"] for spec in source_specs.values() if "dataset" in spec
    )
    levels_summary: list[dict[str, int]] = []
    with h5py.File(path, "r") as handle:
        if "amr" not in handle or "hydro" not in handle:
            raise ValueError("HDF5 snapshot requires /amr and /hydro groups")
        levels = sorted(
            int(name.removeprefix("level_"))
            for name in handle["amr"].keys()
            if name.startswith("level_") and f"hydro/{name}" in handle
        )
        if not levels:
            raise ValueError("HDF5 snapshot has no AMR/hydro levels")
        for level in levels:
            amr = handle[f"amr/level_{level}"]
            hydro = handle[f"hydro/level_{level}"]
            for coordinate in ("xg_1", "xg_2", "xg_3"):
                if coordinate not in amr:
                    raise ValueError(f"AMR level {level} is missing {coordinate}")
            grid_count = len(amr["xg_1"])
            if any(len(amr[coordinate]) != grid_count for coordinate in ("xg_2", "xg_3")):
                raise ValueError(f"AMR level {level} coordinate lengths disagree")
            if "son_flag" not in amr or len(amr["son_flag"]) < grid_count * 8:
                raise ValueError(f"AMR level {level} has an incomplete son_flag dataset")
            missing = sorted(dataset for dataset in required_datasets if dataset not in hydro)
            if missing:
                raise ValueError(f"AMR level {level} is missing mapped datasets: {missing}")
            short = sorted(
                dataset for dataset in required_datasets if len(hydro[dataset]) < grid_count * 8
            )
            if short:
                raise ValueError(f"AMR level {level} has short mapped datasets: {short}")
            levels_summary.append({"level": level, "grid_count": grid_count, "cell_count": grid_count * 8})
        gamma = _read_hdf5_scalar_attribute(handle, "header", "gamma")
        nvar = _read_hdf5_scalar_attribute(handle, "header", "nvar")
        nlevelmax_file = _read_hdf5_scalar_attribute(handle, "amr", "nlevelmax_file")
    return {
        "path": str(path.resolve()),
        "populated_levels": levels,
        "level_summary": levels_summary,
        "header_gamma": gamma,
        "header_nvar": int(nvar) if nvar is not None else None,
        "amr_nlevelmax_file": int(nlevelmax_file) if nlevelmax_file is not None else None,
        "mapped_datasets": sorted(required_datasets),
    }


def _candidate_grid_batches(
    amr: h5py.Group,
    *,
    cell_width_code: float,
    left_edge_code: np.ndarray,
    right_edge_code: np.ndarray,
    chunk_grids: int,
):
    """Yield local AMR grids whose eight-cell boxes intersect the P4 cube."""

    count = len(amr["xg_1"])
    for start in range(0, count, chunk_grids):
        stop = min(start + chunk_grids, count)
        centers = np.column_stack(
            (
                np.asarray(amr["xg_1"][start:stop]),
                np.asarray(amr["xg_2"][start:stop]),
                np.asarray(amr["xg_3"][start:stop]),
            )
        )
        intersects = np.all(centers + cell_width_code > left_edge_code, axis=1)
        intersects &= np.all(centers - cell_width_code < right_edge_code, axis=1)
        local = np.flatnonzero(intersects)
        if len(local):
            yield start + local, centers[local]


def _leaf_overlaps(
    shape: tuple[int, ...] | np.ndarray,
    *,
    center_code: np.ndarray,
    leaf_width_code: float,
    left_edge_code: np.ndarray,
    analysis_cell_width_code: float,
):
    """Yield target-cell indices and overlap fractions for one AMR leaf."""

    shape = np.asarray(shape, dtype=np.int64)
    if shape.shape != (3,) or np.any(shape <= 0):
        raise ValueError("analysis shape must contain three positive dimensions")
    cell_left = center_code - 0.5 * leaf_width_code
    cell_right = center_code + 0.5 * leaf_width_code
    start = np.floor((cell_left - left_edge_code) / analysis_cell_width_code + 1.0e-10).astype(np.int64)
    stop = np.ceil((cell_right - left_edge_code) / analysis_cell_width_code - 1.0e-10).astype(np.int64)
    start = np.maximum(start, 0)
    stop = np.minimum(stop, shape)
    if np.any(stop <= start):
        return
    for index_x in range(start[0], stop[0]):
        for index_y in range(start[1], stop[1]):
            for index_z in range(start[2], stop[2]):
                index = np.asarray((index_x, index_y, index_z))
                analysis_left = left_edge_code + index * analysis_cell_width_code
                analysis_right = analysis_left + analysis_cell_width_code
                overlap = np.minimum(cell_right, analysis_right) - np.maximum(cell_left, analysis_left)
                if np.all(overlap > 0.0):
                    yield (index_x, index_y, index_z), float(np.prod(overlap / analysis_cell_width_code))


def _deposit_leaf(
    density_sum: np.ndarray,
    coverage: np.ndarray,
    *,
    center_code: np.ndarray,
    density_code: float,
    leaf_width_code: float,
    left_edge_code: np.ndarray,
    analysis_cell_width_code: float,
) -> None:
    """Volume-average one AMR leaf cell into every overlapping analysis cell."""

    for index, volume_fraction in _leaf_overlaps(
        density_sum.shape,
        center_code=center_code,
        leaf_width_code=leaf_width_code,
        left_edge_code=left_edge_code,
        analysis_cell_width_code=analysis_cell_width_code,
    ):
        density_sum[index] += density_code * volume_fraction
        coverage[index] += volume_fraction


def _deposit_mapped_field(
    field_sum: np.ndarray,
    field_weight: np.ndarray,
    *,
    index: tuple[int, int, int],
    raw_value: float,
    density_code: float,
    volume_fraction: float,
    averaging: str,
) -> None:
    if not np.isfinite(raw_value):
        raise ValueError("mapped HDF5 field contains a non-finite value")
    weight = density_code if averaging == "mass_weighted" else 1.0
    field_sum[index] += raw_value * weight * volume_fraction
    field_weight[index] += weight * volume_fraction


def _field_description(spec: dict[str, object] | None) -> dict[str, object]:
    if spec is None:
        return {"status": "missing"}
    if "derive" in spec:
        description = {
            "status": "derived",
            "derive": spec["derive"],
            "depends_on": spec["depends_on"],
        }
        for key in ("gamma", "solar_mass_fraction"):
            if key in spec:
                description[key] = float(spec[key])
        return description
    if "dataset" in spec:
        description = {
            "status": "dataset",
            "dataset": spec["dataset"],
            "unit": spec["unit"],
            "averaging": spec["averaging"],
        }
        for key in ("quantity", "source_semantics"):
            if key in spec:
                description[key] = spec[key]
        return description
    return {
        "status": "constant",
        "constant": float(spec["constant"]),
        "reason": spec["reason"],
        "unit": spec["unit"],
        "averaging": spec["averaging"],
    }


def _production_contract_complete(
    field_map: dict[str, object],
    *,
    sources_present: bool,
) -> bool:
    fields = field_map["fields"]
    if not sources_present or not _snapshot_backed_field(fields, "velocity"):
        return False
    if not (
        _snapshot_backed_field(fields, "temperature_k")
        or _snapshot_backed_field(fields, "thermal_pressure")
    ):
        return False
    return all(_snapshot_backed_field(fields, name) for name in _PRODUCTION_DATASET_FIELDS)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True)
    parser.add_argument("--info", required=True)
    parser.add_argument("--zoom-manifest", required=True)
    parser.add_argument("--field-map", required=True)
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument("--source-ledger", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--analysis-level", type=int, default=15)
    parser.add_argument("--scan-chunk-grids", type=int, default=1_000_000)
    parser.add_argument("--metallicity-solar", type=float, default=1.0e-6)
    parser.add_argument(
        "--allow-equilibrium-fallback",
        action="store_true",
        help="allow atlas equilibrium temperature only when no temperature/pressure field is mapped",
    )
    parser.add_argument(
        "--require-production-contract",
        action="store_true",
        help="fail unless all production-required fields are snapshot datasets, not constants",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="validate the HDF5/field-map contract without scanning field values or writing output",
    )
    args = parser.parse_args()
    if args.analysis_level < 0:
        raise ValueError("analysis-level must be non-negative")
    if args.scan_chunk_grids < 1:
        raise ValueError("scan-chunk-grids must be positive")
    if args.metallicity_solar <= 0.0:
        raise ValueError("metallicity-solar must be positive")

    info = _read_info(Path(args.info))
    field_map = _read_field_map(Path(args.field_map))
    fields = field_map["fields"]
    velocity_specs = fields.get("velocity")
    mapped_specs: dict[str, dict[str, object]] = {}
    if velocity_specs is not None:
        for component, spec in enumerate(velocity_specs):
            mapped_specs[f"velocity_{component}"] = spec
    for name in _SCALAR_FIELDS:
        spec = fields.get(name)
        if spec is not None and "derive" not in spec:
            mapped_specs[name] = spec
    for spec in mapped_specs.values():
        _field_unit_scale(spec, info)

    manifest = json.loads(Path(args.zoom_manifest).read_text())
    final = manifest["final"]
    left = np.asarray(final["left_edge_code"], dtype=np.float64)
    shape = np.asarray(final["shape"], dtype=np.int64)
    width = float(final["width_code"])
    if left.shape != (3,) or shape.shape != (3,) or np.any(shape <= 0) or width <= 0.0:
        raise ValueError("invalid P4 final-cube geometry")
    analysis_cell_width = 2.0 ** (-args.analysis_level)
    if not np.allclose(width, analysis_cell_width * shape, rtol=0.0, atol=analysis_cell_width * 1.0e-8):
        raise ValueError("P4 cube does not align with the requested RAMSES analysis level")
    right = left + width
    snapshot_path = Path(args.snapshot)
    if args.preflight_only:
        preflight = _preflight_hdf5_snapshot(snapshot_path, fields=fields, source_specs=mapped_specs)
        atlas = read_thermal_atlas(args.thermal_atlas)
        source_ledger = read_photon_source_ledger_csv(args.source_ledger)
        sources = source_ledger.source_catalog_in_cube(left, width, tuple(int(value) for value in shape))
        if sources is None:
            raise ValueError("photon ledger has no source inside the P4 cube")
        gamma = fields.get("thermal_pressure")
        if isinstance(gamma, dict) and "derive" in gamma and preflight["header_gamma"] is None and "gamma" not in gamma:
            raise ValueError("derived thermal_pressure needs HDF5 header gamma or an explicit map gamma")
        print(
            "P4_HDF5_PREFLIGHT_OK "
            f"levels={preflight['populated_levels']} datasets={len(preflight['mapped_datasets'])} "
            f"sources={len(sources.cell_index)} atlas_shape={atlas.net_rate_erg_s_cm3.shape} "
            f"header_gamma={preflight['header_gamma']}",
            flush=True,
        )
        return
    output_shape = tuple(int(value) for value in shape)
    density_sum = np.zeros(output_shape, dtype=np.float64)
    coverage = np.zeros(output_shape, dtype=np.float64)
    field_sums = {name: np.zeros(output_shape, dtype=np.float64) for name in mapped_specs}
    field_weights = {name: np.zeros(output_shape, dtype=np.float64) for name in mapped_specs}
    level_summary = []
    mass_deposited_code = 0.0

    snapshot_gamma: float | None = None
    with h5py.File(snapshot_path, "r") as handle:
        snapshot_gamma = _read_hdf5_scalar_attribute(handle, "header", "gamma")
        levels = sorted(
            int(name.removeprefix("level_"))
            for name in handle["amr"].keys()
            if name.startswith("level_") and f"hydro/{name}" in handle
        )
        if not levels:
            raise ValueError("HDF5 snapshot has no AMR/hydro levels")
        density_spec = fields["density"]
        for level in levels:
            amr = handle[f"amr/level_{level}"]
            hydro = handle[f"hydro/level_{level}"]
            datasets = {"density": hydro[density_spec["dataset"]]}
            for name, spec in mapped_specs.items():
                if "dataset" in spec:
                    if spec["dataset"] not in hydro:
                        raise ValueError(
                            f"field-map dataset {spec['dataset']!r} is absent from hydro/level_{level}"
                        )
                    datasets[name] = hydro[spec["dataset"]]
            leaf_width = 2.0 ** (-level)
            candidate_grids = 0
            deposited_leaves = 0
            refined_cells = 0
            for grid_indices, grid_centers in _candidate_grid_batches(
                amr,
                cell_width_code=leaf_width,
                left_edge_code=left,
                right_edge_code=right,
                chunk_grids=args.scan_chunk_grids,
            ):
                candidate_grids += len(grid_indices)
                blocks = {name: _read_grid_blocks(dataset, grid_indices) for name, dataset in datasets.items()}
                son_blocks = _read_grid_blocks(amr["son_flag"], grid_indices)
                for grid_position, (grid_center, son_flags) in enumerate(zip(grid_centers, son_blocks, strict=True)):
                    child_centers = grid_center + (_CHILD_BITS - 0.5) * leaf_width
                    for child_position, (child_center, son_flag) in enumerate(zip(child_centers, son_flags, strict=True)):
                        child_left = child_center - 0.5 * leaf_width
                        child_right = child_center + 0.5 * leaf_width
                        if np.any(child_right <= left) or np.any(child_left >= right):
                            continue
                        if son_flag != 0.0:
                            refined_cells += 1
                            continue
                        density_code = float(blocks["density"][grid_position, child_position])
                        if not np.isfinite(density_code) or density_code <= 0.0:
                            raise ValueError(f"invalid density in level {level} grid {int(grid_indices[grid_position])}")
                        leaf_deposited = False
                        for index, volume_fraction in _leaf_overlaps(
                            shape,
                            center_code=child_center,
                            leaf_width_code=leaf_width,
                            left_edge_code=left,
                            analysis_cell_width_code=analysis_cell_width,
                        ):
                            density_sum[index] += density_code * volume_fraction
                            coverage[index] += volume_fraction
                            mass_deposited_code += density_code * volume_fraction
                            for name, spec in mapped_specs.items():
                                raw_value = (
                                    float(spec["constant"])
                                    if "constant" in spec
                                    else float(blocks[name][grid_position, child_position])
                                )
                                _deposit_mapped_field(
                                    field_sums[name],
                                    field_weights[name],
                                    index=index,
                                    raw_value=raw_value,
                                    density_code=density_code,
                                    volume_fraction=volume_fraction,
                                    averaging=spec["averaging"],
                                )
                            leaf_deposited = True
                        deposited_leaves += int(leaf_deposited)
            level_summary.append(
                {
                    "level": level,
                    "candidate_grids": candidate_grids,
                    "deposited_leaf_cells": deposited_leaves,
                    "refined_cells_in_cube": refined_cells,
                }
            )
            print(
                f"P4_HDF5_SCAN level={level} candidates={candidate_grids} "
                f"leaves={deposited_leaves} refined={refined_cells}",
                flush=True,
            )

    if not np.allclose(coverage, 1.0, rtol=0.0, atol=1.0e-8):
        raise ValueError(
            "adaptive AMR leaves do not cover the P4 mesh exactly: "
            f"min={coverage.min():.16g} max={coverage.max():.16g}"
        )
    density_code = density_sum / coverage
    if not np.isfinite(density_code).all() or np.any(density_code <= 0.0):
        raise ValueError("resampled AMR gas density is non-finite or non-positive")
    raw_field_values: dict[str, np.ndarray] = {}
    for name, spec in mapped_specs.items():
        if np.any(field_weights[name] <= 0.0):
            raise ValueError(f"mapped field {name} has cells with no valid averaging weight")
        raw_field_values[name] = field_sums[name] / field_weights[name] * _field_unit_scale(spec, info)

    density_cgs = density_code * info["unit_d"]
    field_values = dict(raw_field_values)
    velocity_components: list[np.ndarray] = []
    if velocity_specs is not None:
        for component, spec in enumerate(velocity_specs):
            name = f"velocity_{component}"
            value = raw_field_values[name]
            if spec.get("quantity") == "momentum_density" or spec["unit"] == "code_momentum_density":
                value = value / np.maximum(density_cgs, np.finfo(np.float64).tiny)
            if np.any(~np.isfinite(value)):
                raise ValueError(f"derived velocity component {component} is non-finite")
            field_values[name] = value
            velocity_components.append(value)

    thermal_spec = fields.get("thermal_pressure")
    if isinstance(thermal_spec, dict) and "derive" in thermal_spec:
        if "total_energy_density" not in raw_field_values or len(velocity_components) != 3:
            raise ValueError(
                "conservative thermal_pressure derivation requires total_energy_density and all three momenta"
            )
        header_gamma = snapshot_gamma
        map_gamma = float(thermal_spec["gamma"]) if "gamma" in thermal_spec else None
        if header_gamma is not None and map_gamma is not None and not np.isclose(
            header_gamma, map_gamma, rtol=0.0, atol=1.0e-10
        ):
            raise ValueError(
                f"field-map gamma {map_gamma:.16g} disagrees with HDF5 header gamma {header_gamma:.16g}"
            )
        gamma = header_gamma if header_gamma is not None else map_gamma
        if gamma is None or gamma <= 1.0:
            raise ValueError("conservative thermal_pressure derivation requires gamma")
        kinetic_energy_density = 0.5 * density_cgs * sum(value**2 for value in velocity_components)
        thermal_pressure = (gamma - 1.0) * (raw_field_values["total_energy_density"] - kinetic_energy_density)
        if np.any(~np.isfinite(thermal_pressure)) or np.any(thermal_pressure <= 0.0):
            raise ValueError(
                "conservative thermal_pressure derivation produced non-positive values: "
                f"min={thermal_pressure.min():.16g}"
            )
        field_values["thermal_pressure"] = thermal_pressure

    metallicity_spec = fields.get("metallicity_solar")
    if isinstance(metallicity_spec, dict) and "derive" in metallicity_spec:
        metal_density_name = str(metallicity_spec.get("metal_density_field", "metal_density"))
        if metal_density_name not in raw_field_values:
            raise ValueError(
                f"derived metallicity_solar requires mapped metal-density field {metal_density_name!r}"
            )
        metal_density = raw_field_values[metal_density_name]
        if np.any(~np.isfinite(metal_density)) or np.any(metal_density < 0.0):
            raise ValueError("mapped metal density is non-finite or negative")
        solar_mass_fraction = float(metallicity_spec["solar_mass_fraction"])
        metallicity = metal_density / np.maximum(density_cgs, np.finfo(np.float64).tiny) / solar_mass_fraction
        field_values["metallicity_solar"] = metallicity

    atlas = read_thermal_atlas(args.thermal_atlas)
    hydrogen_number_density = density_cgs * 0.76 / 1.67262192369e-24
    metallicity = field_values.get("metallicity_solar", np.full(output_shape, args.metallicity_solar))
    if np.any(~np.isfinite(metallicity)) or np.any(metallicity <= 0.0):
        raise ValueError("metallicity_solar must be positive for thermal-atlas interpolation")
    if "temperature_k" in field_values:
        temperature = field_values["temperature_k"]
        thermal_initialization = "mapped temperature_k field"
    elif "thermal_pressure" in field_values:
        temperature = atlas.temperature_from_pressure(
            info["aexp"],
            field_values["thermal_pressure"],
            density_cgs,
            metallicity,
            hydrogen_mass_fraction=0.76,
        )
        thermal_initialization = "mapped thermal_pressure field plus thermal atlas inversion"
    elif args.allow_equilibrium_fallback:
        temperature = atlas.equilibrium_temperature(info["aexp"], hydrogen_number_density, metallicity)
        thermal_initialization = "explicitly allowed thermal-atlas equilibrium fallback"
    else:
        raise ValueError("field map must provide temperature_k or thermal_pressure unless fallback is explicit")
    if np.any(~np.isfinite(temperature)) or np.any(temperature <= 0.0):
        raise ValueError("staged temperature is non-finite or non-positive")

    if "dust_relative_abundance" in field_values:
        dust = field_values["dust_relative_abundance"]
        dust_initialization = "mapped dust_relative_abundance field"
        dust_relative_abundance_origin = "direct"
    elif "dust_to_metal" in field_values and "metallicity_solar" in field_values:
        dust = metallicity * field_values["dust_to_metal"]
        dust_initialization = "derived metallicity_solar times dust_to_metal"
        dust_relative_abundance_origin = "metallicity_solar_times_dust_to_metal"
    else:
        raise ValueError("field map must provide dust_relative_abundance or both metallicity_solar and dust_to_metal")
    if np.any(~np.isfinite(dust)) or np.any(dust < 0.0):
        raise ValueError("dust_relative_abundance must be finite and non-negative")

    ionization_names = ("x_hii", "x_heii", "x_heiii")
    ionization_present = [name in field_values for name in ionization_names]
    if not all(ionization_present):
        raise ValueError("field map must explicitly map the initial H/He ionization fractions")
    x_h2 = field_values.get("x_h2")
    velocity = None
    if velocity_specs is not None:
        velocity = np.stack([field_values[f"velocity_{component}"] for component in range(3)], axis=0)

    source_ledger = read_photon_source_ledger_csv(args.source_ledger)
    sources = source_ledger.source_catalog_in_cube(left, width, output_shape)
    if sources is None:
        raise ValueError("photon ledger has no source inside the P4 cube")
    production_contract_complete = _production_contract_complete(field_map, sources_present=True)
    if args.require_production_contract and not production_contract_complete:
        raise ValueError(
            "production contract is incomplete: all required thermodynamic, kinematic, composition, "
            "ionization, and source fields must be snapshot datasets"
        )
    snapshot = neutral_primordial_input(
        GridSpec(
            cell_width_cm=info["unit_l"] * analysis_cell_width,
            left_edge_cm=left * info["unit_l"],
        ),
        density_cgs,
        temperature,
        dust_relative_abundance=dust,
        sources=sources,
        velocity_cm_s=velocity,
        metallicity_solar=metallicity,
        dust_to_metal=field_values.get("dust_to_metal"),
        dust_relative_abundance_origin=dust_relative_abundance_origin,
        x_h2=x_h2,
        cell_level=np.full(output_shape, args.analysis_level, dtype=np.int16),
        x_hii=field_values["x_hii"],
        x_heii=field_values["x_heii"],
        x_heiii=field_values["x_heiii"],
    )
    if args.require_production_contract:
        snapshot.validate_production_contract(require_sources=True)
    write_static_rt_input(args.output, snapshot)

    field_contract = {"density": _field_description(fields["density"])}
    if velocity_specs is None:
        field_contract["velocity"] = {"status": "missing"}
    else:
        field_contract["velocity"] = [_field_description(spec) for spec in velocity_specs]
    for name in _SCALAR_FIELDS:
        field_contract[name] = _field_description(fields.get(name))
    metadata = {
        "snapshot": str(Path(args.snapshot).resolve()),
        "info": str(Path(args.info).resolve()),
        "field_map": str(Path(args.field_map).resolve()),
        "field_map_schema_version": field_map["schema_version"],
        "field_contract": field_contract,
        "production_contract_complete": production_contract_complete,
        "snapshot_header_gamma": snapshot_gamma,
        "aexp": info["aexp"],
        "analysis_level": args.analysis_level,
        "analysis_rule": "volume-conservative resampling of all intersecting RAMSES AMR leaf cells",
        "resampled_fields": ["density"] + list(mapped_specs),
        "density_mass_conservative": True,
        "populated_amr_levels": levels,
        "level_summary": level_summary,
        "left_edge_code": left.tolist(),
        "shape": shape.tolist(),
        "cell_width_code": analysis_cell_width,
        "cell_width_cm": info["unit_l"] * analysis_cell_width,
        "density_unit_g_cm3": info["unit_d"],
        "mean_n_h_cm3": float(hydrogen_number_density.mean()),
        "max_n_h_cm3": float(hydrogen_number_density.max()),
        "coverage_min": float(coverage.min()),
        "coverage_max": float(coverage.max()),
        "mass_deposited_code": mass_deposited_code,
        "mass_from_resampled_density_code": float(density_code.sum()),
        "mass_relative_error": abs(mass_deposited_code - float(density_code.sum()))
        / max(abs(mass_deposited_code), np.finfo(np.float64).tiny),
        "thermal_initialization": thermal_initialization,
        "dust_initialization": dust_initialization,
        "dust_relative_abundance_origin": dust_relative_abundance_origin,
        "thermal_atlas": str(Path(args.thermal_atlas).resolve()),
        "source_ledger": str(Path(args.source_ledger).resolve()),
        "photon_source_count": int(len(sources.cell_index)),
    }
    metadata_path = Path(args.metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        "P4_HDF5_STAGE_OK "
        f"shape={tuple(shape)} sources={len(sources.cell_index)} mean_n_h={hydrogen_number_density.mean():.6g} "
        f"production_contract={production_contract_complete}",
        flush=True,
    )


if __name__ == "__main__":
    main()
