#!/usr/bin/env python3
"""Validate and compare the first-step nDGP density source on two AMR layouts."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import pathlib
import re
import sys
from dataclasses import dataclass
from typing import Iterable

from compare_cuda_ndgp_outputs import (
    FormatError,
    Particle,
    pinned_particle_ids,
    read_particles,
)
from validate_cuda_ndgp_diagnostic import (
    residuals,
    validate_cpu_only_markers,
    validate_warnings,
)


class SourceDiagnosticError(RuntimeError):
    pass


@dataclass(frozen=True)
class SourceCell:
    rho: float
    rho_tot: float
    refined: bool
    rank: int


Coord = tuple[int, int, int]
SourceMap = dict[Coord, SourceCell]
HEADER = re.compile(
    r"^# schema=ndgp-source-v1 rank=(\d+) level=(\d+) "
    r"nstep=(\d+) icount=(\d+) rho_tot=\s*(\S+)$"
)
FATAL = re.compile(
    r"(?i)MPI_ABORT|forrtl: severe|segmentation fault|error stop|FATAL:|ERROR:"
)
CIC_RTOL = 5.0e-12
CIC_ATOL = 5.0e-12
COMMON_NML_SETTINGS = (
    "cosmo=.true.",
    "pic=.true.",
    "poisson=.true.",
    "hydro=.false.",
    "clumpfind=.false.",
    "sink=.false.",
    "sinkprops=.false.",
    "lightcone=.false.",
    "rt=.false.",
    "aton=.false.",
    "de_perturb=.false.",
    "sidm=.false.",
    "use_nDGP=.true.",
    "use_fR=.false.",
    "use_symmetron=.false.",
    "use_dilaton=.false.",
    "use_galileon=.false.",
    "use_mond=.false.",
    "use_coupled_de=.false.",
    "use_quintessence=.false.",
    "use_kessence=.false.",
    "use_chaplygin=.false.",
    "use_rvm=.false.",
    "use_horndeski=.false.",
    "use_ede=.false.",
    "use_neutrino=.false.",
    "use_sgs=.false.",
    "use_adm=.false.",
    "use_fdm=.false.",
    "use_pbh=.false.",
    "levelmin=5",
    "nstepmax=1",
    "nrestart=0",
    "nexpand=0",
    "nDGP_eps=1.0d-4",
    "scalar_solver_strict=.false.",
    "static=.false.",
    "use_fftw=.false.",
    "mg_merged_rb=.false.",
    "cg_levelmin=999",
    "cic_levelmax=0",
    "gpu_auto_tune=.false.",
    "gpu_hydro=.false.",
    "gpu_poisson=.false.",
    "gpu_fft=.false.",
    "gpu_sink=.false.",
    "gpu_scalar=.false.",
    "gpu_particle=.false.",
)


def read_source_file(path: pathlib.Path) -> tuple[int, int, SourceMap]:
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) < 2:
        raise SourceDiagnosticError(f"{path}: missing header or data")
    match = HEADER.fullmatch(lines[0])
    if not match or lines[1] != "# level x y z rho rho_tot refined":
        raise SourceDiagnosticError(f"{path}: malformed schema header")
    rank, level, nstep, icount = (int(value) for value in match.groups()[:4])
    try:
        header_rho_tot = float(match.group(5).replace("D", "E").replace("d", "e"))
    except ValueError as error:
        raise SourceDiagnosticError(f"{path}: malformed header rho_tot") from error
    if nstep != 0 or icount != 1 or not math.isfinite(header_rho_tot):
        raise SourceDiagnosticError(
            f"{path}: expected nstep=0 icount=1 and finite rho_tot"
        )
    cells: SourceMap = {}
    for line_number, line in enumerate(lines[2:], 3):
        fields = line.split()
        if len(fields) != 7:
            raise SourceDiagnosticError(f"{path}:{line_number}: expected 7 fields")
        try:
            row_level, x, y, z = (int(value) for value in fields[:4])
            rho = float(fields[4].replace("D", "E").replace("d", "e"))
            rho_tot = float(fields[5].replace("D", "E").replace("d", "e"))
            refined_raw = int(fields[6])
        except ValueError as error:
            raise SourceDiagnosticError(f"{path}:{line_number}: malformed field") from error
        if row_level != level or refined_raw not in (0, 1):
            raise SourceDiagnosticError(f"{path}:{line_number}: level/refined mismatch")
        if not math.isfinite(rho) or not math.isfinite(rho_tot):
            raise SourceDiagnosticError(f"{path}:{line_number}: non-finite density")
        if rho < 0.0 or rho_tot <= 0.0:
            raise SourceDiagnosticError(f"{path}:{line_number}: nonpositive density mean")
        if rho_tot != header_rho_tot:
            raise SourceDiagnosticError(f"{path}:{line_number}: rho_tot changed")
        key = (x, y, z)
        if key in cells:
            raise SourceDiagnosticError(f"{path}:{line_number}: duplicate cell {key}")
        cells[key] = SourceCell(rho, rho_tot, bool(refined_raw), rank)
    if not cells:
        raise SourceDiagnosticError(f"{path}: rank owns no level-{level} cells")
    return rank, level, cells


def read_source_level(run_dir: pathlib.Path, level: int) -> SourceMap:
    paths = sorted(run_dir.glob(f"ndgp_source_l{level:03d}_rank*.tsv"))
    if len(paths) != 2:
        raise SourceDiagnosticError(
            f"{run_dir}: expected two level-{level} rank dumps, found {len(paths)}"
        )
    result: SourceMap = {}
    ranks: set[int] = set()
    for path in paths:
        rank, file_level, cells = read_source_file(path)
        suffix = re.search(r"_rank(\d{5})\.tsv$", path.name)
        if file_level != level or suffix is None or int(suffix.group(1)) != rank:
            raise SourceDiagnosticError(f"{path}: filename/header mismatch")
        if rank not in (1, 2) or rank in ranks:
            raise SourceDiagnosticError(f"{path}: rank set is not exactly 1,2")
        ranks.add(rank)
        overlap = result.keys() & cells.keys()
        if overlap:
            raise SourceDiagnosticError(
                f"{path}: cross-rank duplicate cell {next(iter(overlap))}"
            )
        result.update(cells)
    if ranks != {1, 2}:
        raise SourceDiagnosticError(f"{run_dir}: rank set is {sorted(ranks)}")
    return result


def cube(low: int, high: int) -> set[Coord]:
    return set(itertools.product(range(low, high), repeat=3))


def face_neighbors(key: Coord, period: int) -> Iterable[Coord]:
    for axis in range(3):
        for offset in (-1, 1):
            neighbor = list(key)
            neighbor[axis] = (neighbor[axis] + offset) % period
            yield tuple(neighbor)  # type: ignore[return-value]


def validate_geometry(
    uniform_l5: SourceMap, amr_l5: SourceMap, amr_l6: SourceMap
) -> tuple[set[Coord], set[Coord], set[Coord], set[Coord], set[Coord]]:
    expected_l5 = cube(0, 32)
    if set(uniform_l5) != expected_l5 or set(amr_l5) != expected_l5:
        raise SourceDiagnosticError("L5 maps are not exactly the canonical 32^3 box")
    if any(cell.refined for cell in uniform_l5.values()):
        raise SourceDiagnosticError("uniform L5 map unexpectedly contains refined cells")
    refined = {key for key, cell in amr_l5.items() if cell.refined}
    if len(refined) != 16**3:
        raise SourceDiagnosticError(
            f"AMR L5 refined-parent count is {len(refined)}, expected 4096"
        )
    expected_l6 = {
        (2 * x + bx, 2 * y + by, 2 * z + bz)
        for x, y, z in refined
        for bx, by, bz in itertools.product((0, 1), repeat=3)
    }
    if set(amr_l6) != expected_l6 or len(expected_l6) != 32**3:
        raise SourceDiagnosticError("AMR L6 cells are not the eight children of each parent")
    if any(cell.refined for cell in amr_l6.values()):
        raise SourceDiagnosticError("levelmax L6 map unexpectedly contains refined cells")
    interface = {
        key
        for key in expected_l5 - refined
        if any(neighbor in refined for neighbor in face_neighbors(key, 32))
    }
    exterior = expected_l5 - refined - interface
    refined_boundary = {
        key
        for key in refined
        if any(neighbor not in refined for neighbor in face_neighbors(key, 32))
    }
    refined_interior = refined - refined_boundary
    return refined, interface, exterior, refined_boundary, refined_interior


def metric_summary(values: dict[Coord, float]) -> dict[str, object]:
    if not values:
        return {"count": 0, "l1": 0.0, "l2": 0.0, "linf": 0.0, "mean": 0.0}
    ordered = sorted(values.items())
    count = len(ordered)
    signed = [value for _, value in ordered]
    maximum_key, maximum_value = max(ordered, key=lambda item: abs(item[1]))
    return {
        "count": count,
        "l1": math.fsum(abs(value) for value in signed) / count,
        "l2": math.sqrt(math.fsum(value * value for value in signed) / count),
        "linf": abs(maximum_value),
        "argmax": list(maximum_key),
        "argmax_signed": maximum_value,
        "mean": math.fsum(signed) / count,
    }


def rho_tot_value(*maps: SourceMap) -> float:
    values = {cell.rho_tot for mapping in maps for cell in mapping.values()}
    if len(values) != 1:
        raise SourceDiagnosticError(f"rho_tot is not globally consistent: {values}")
    return next(iter(values))


def require_setting(text: str, setting: str, path: pathlib.Path) -> None:
    if text.splitlines().count(setting) != 1:
        raise SourceDiagnosticError(f"{path}: setting is not pinned exactly once: {setting}")


def validate_nml_pair(
    uniform_dir: pathlib.Path, amr_dir: pathlib.Path, cap: int = 1
) -> None:
    uniform_path = uniform_dir / "run.nml"
    amr_path = amr_dir / "run.nml"
    uniform = uniform_path.read_text(encoding="ascii")
    amr = amr_path.read_text(encoding="ascii")
    for path, text in ((uniform_path, uniform), (amr_path, amr)):
        for setting in COMMON_NML_SETTINGS:
            require_setting(text, setting, path)
        require_setting(text, f"n_iter_nDGP={cap}", path)
    expected_ic = uniform_dir.parent / "ic"
    expected_l5 = f"initfile(1)='{expected_ic / 'level_005'}'"
    expected_l6 = f"initfile(2)='{expected_ic / 'level_006'}'"
    require_setting(uniform, expected_l5, uniform_path)
    require_setting(amr, expected_l5, amr_path)
    require_setting(uniform, "levelmax=5", uniform_path)
    require_setting(amr, "levelmax=6", amr_path)
    if "initfile(2)" in uniform or "level_006" in uniform:
        raise SourceDiagnosticError("uniform NML still references level-6 IC")
    require_setting(amr, expected_l6, amr_path)
    ignored = re.compile(r"^(?:levelmax=|initfile\(2\)=)")
    uniform_normalized = [line for line in uniform.splitlines() if not ignored.match(line)]
    amr_normalized = [line for line in amr.splitlines() if not ignored.match(line)]
    if uniform_normalized != amr_normalized:
        raise SourceDiagnosticError(
            "source NMLs differ outside levelmax and the AMR-only initfile(2)"
        )


def validate_log(
    run_dir: pathlib.Path,
    expected_markers: set[tuple[int, int]],
    expected_levels: tuple[int, ...],
    cap: int = 1,
) -> dict[str, list[dict[str, float | int | str]]]:
    text = (run_dir / "run.log").read_text(errors="replace")
    if text.count("Run completed") != 1 or FATAL.search(text):
        raise SourceDiagnosticError(f"{run_dir}: incomplete or fatal runtime")
    if "Some grid are outside initial conditions sub-volume" in text:
        raise SourceDiagnosticError(f"{run_dir}: refinement exceeds the pinned IC")
    validate_cpu_only_markers(text)
    validate_warnings(text, cap)
    legacy_warning = (
        "WARNING: IC header carries no omega_b (legacy grafic); "
        "using namelist omega_b"
    )
    if sum(line.strip() == legacy_warning for line in text.splitlines()) != len(
        expected_levels
    ):
        raise SourceDiagnosticError(f"{run_dir}: legacy IC warning count differs")
    marker_rows = re.findall(
        r"\[NDGP_SOURCE_DIAG\]\s+rank=(\d+)\s+level=(\d+)\s+owned_cells=(\d+)",
        text,
    )
    actual_markers = {(int(rank), int(level)) for rank, level, _ in marker_rows}
    if (
        len(marker_rows) != len(expected_markers)
        or actual_markers != expected_markers
        or any(int(count) <= 0 for _, _, count in marker_rows)
    ):
        raise SourceDiagnosticError(
            f"{run_dir}: source marker rows differ: {marker_rows}"
        )
    if text.count("/ic_particle_ids") != len(expected_levels):
        raise SourceDiagnosticError(f"{run_dir}: particle-ID file read count differs")
    if text.count("init_part: idp set from genetIC ic_particle_ids") != 1:
        raise SourceDiagnosticError(f"{run_dir}: particle-ID positive marker differs")
    for level in expected_levels:
        if not re.search(rf"Level\s+{level}\s+has\s+4096\s+grids", text):
            raise SourceDiagnosticError(f"{run_dir}: initial L{level} topology differs")
    if 6 not in expected_levels and re.search(r"nDGP level\s+6\b|Level\s+6\s+has", text):
        raise SourceDiagnosticError(f"{run_dir}: uniform run entered level 6")
    dmo = [line for line in text.splitlines() if "NaN_CHK_DMO" in line]
    if not dmo:
        raise SourceDiagnosticError(f"{run_dir}: finite diagnostic is absent")
    for line in dmo:
        values = dict(re.findall(r"(rho|phi|f|scalar)=([^\s]+)", line))
        if set(values) != {"rho", "phi", "f", "scalar"} or any(
            int(value) != 0 for value in values.values()
        ):
            raise SourceDiagnosticError(f"{run_dir}: malformed/nonzero finite marker")
    rows = {str(level): residuals(text, level, cap) for level in expected_levels}
    for level, level_rows in rows.items():
        if len(level_rows) != 2:
            raise SourceDiagnosticError(
                f"{run_dir}: expected two level-{level} solver rows, got {len(level_rows)}"
            )
    expected_sequence = (5, 5) if expected_levels == (5,) else (5, 6, 5, 6)
    validate_solver_sequence(text, expected_sequence, run_dir)
    return rows


def validate_solver_sequence(
    text: str, expected_sequence: tuple[int, ...], source: pathlib.Path
) -> None:
    actual_sequence = tuple(
        int(level)
        for level in re.findall(
            r"nDGP level\s+([56])\s+(?:converged in|NOT converged after)", text
        )
    )
    if actual_sequence != expected_sequence:
        raise SourceDiagnosticError(
            f"{source}: solver sequence {actual_sequence} differs from {expected_sequence}"
        )


def validate_amr_report(
    report_path: pathlib.Path, output: pathlib.Path, expected_l6_grids: int
) -> None:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("status") != "PASS":
        raise SourceDiagnosticError(f"{report_path}: canonical report is not PASS")
    if any(
        pathlib.Path(report.get(side, {}).get("path", "")).resolve() != output.resolve()
        for side in ("left", "right")
    ):
        raise SourceDiagnosticError(f"{report_path}: output path mismatch")
    left = report.get("left", {})
    counts = left.get("level_counts", {})
    if int(counts.get("5", 0)) != 4096 or int(counts.get("6", 0)) != expected_l6_grids:
        raise SourceDiagnosticError(f"{report_path}: unexpected topology {counts}")
    owners = [item for item in left.get("inputs", []) if "source_cpu" in item]
    if {int(item["source_cpu"]) for item in owners} != {1, 2} or len(owners) != 2:
        raise SourceDiagnosticError(f"{report_path}: owner CPU set differs")
    for item in owners:
        owner_counts = item.get("owner_level_counts", {})
        if int(owner_counts.get("5", 0)) <= 0:
            raise SourceDiagnosticError(f"{report_path}: a rank owns no L5 grids")
        if expected_l6_grids and int(owner_counts.get("6", 0)) <= 0:
            raise SourceDiagnosticError(f"{report_path}: a rank owns no L6 grids")
        if not expected_l6_grids and int(owner_counts.get("6", 0)) != 0:
            raise SourceDiagnosticError(f"{report_path}: uniform run owns L6 grids")


def validate_particle_values(
    particles: dict[int, Particle], output: pathlib.Path
) -> None:
    for identity, particle in particles.items():
        values = (
            *particle.position,
            *particle.velocity,
            particle.mass,
            particle.potential,
        )
        if not all(math.isfinite(value) for value in values):
            raise SourceDiagnosticError(f"{output}: particle {identity} is non-finite")
        if particle.mass <= 0.0:
            raise SourceDiagnosticError(f"{output}: particle {identity} has nonpositive mass")


def validate_outputs(
    run_dir: pathlib.Path, expected_ids: set[int]
) -> dict[int, Particle]:
    outputs = sorted(path for path in run_dir.glob("output_*") if path.is_dir())
    if [path.name for path in outputs] != ["output_00001", "output_00002"]:
        raise SourceDiagnosticError(f"{run_dir}: output set differs")
    initial_particles: dict[int, Particle] | None = None
    for output in outputs:
        if not (output / "COMPLETE").is_file():
            raise SourceDiagnosticError(f"{output}: COMPLETE is absent")
        particles, per_rank = read_particles(output)
        if set(particles) != expected_ids:
            raise SourceDiagnosticError(
                f"{output}: particle IDs differ; extra={len(set(particles)-expected_ids)} "
                f"missing={len(expected_ids-set(particles))}"
            )
        if len(per_rank) != 2 or any(count <= 0 for count in per_rank.values()):
            raise SourceDiagnosticError(f"{output}: both ranks must own particles")
        validate_particle_values(particles, output)
        if initial_particles is None:
            initial_particles = particles
    if initial_particles is None:
        raise SourceDiagnosticError(f"{run_dir}: initial particles are absent")
    return initial_particles


def cic_project(
    particles: dict[int, Particle], ncell: int, accepted: set[Coord]
) -> SourceMap:
    accum = {key: 0.0 for key in accepted}
    cell_volume_inverse = float(ncell**3)
    for particle in particles.values():
        axes: list[tuple[tuple[int, float], tuple[int, float]]] = []
        for position in particle.position:
            scaled = (position % 1.0) * ncell - 0.5
            left_raw = math.floor(scaled)
            fraction = scaled - left_raw
            axes.append(
                (
                    (left_raw % ncell, 1.0 - fraction),
                    ((left_raw + 1) % ncell, fraction),
                )
            )
        density_scale = particle.mass * cell_volume_inverse
        for x, y, z in itertools.product(*axes):
            key = (x[0], y[0], z[0])
            if key in accum:
                accum[key] += density_scale * x[1] * y[1] * z[1]
    return {
        key: SourceCell(value, 0.0, False, 0) for key, value in accum.items()
    }


def direct_cic_report(
    actual: SourceMap, projected: SourceMap, ncell: int
) -> dict[str, object]:
    if set(actual) != set(projected):
        raise SourceDiagnosticError("direct CIC oracle coordinate set differs")
    differences = {
        key: actual[key].rho - projected[key].rho for key in actual
    }
    metrics = metric_summary(differences)
    metrics["actual_integral"] = math.fsum(
        cell.rho for cell in actual.values()
    ) / ncell**3
    metrics["projected_integral"] = math.fsum(
        cell.rho for cell in projected.values()
    ) / ncell**3
    metrics["within_tolerance"] = all(
        math.isclose(
            actual[key].rho,
            projected[key].rho,
            rel_tol=CIC_RTOL,
            abs_tol=CIC_ATOL,
        )
        for key in actual
    )
    return metrics


def fixture_mass_report(
    uniform_particles: dict[int, Particle], amr_particles: dict[int, Particle]
) -> dict[str, object]:
    uniform_total = math.fsum(particle.mass for particle in uniform_particles.values())
    amr_base = math.fsum(
        particle.mass
        for identity, particle in amr_particles.items()
        if identity <= 32**3
    )
    amr_fine = math.fsum(
        particle.mass
        for identity, particle in amr_particles.items()
        if 32**3 < identity <= 2 * 32**3
    )
    amr_total = math.fsum(particle.mass for particle in amr_particles.values())
    expected = {
        "uniform_total": 1.0,
        "amr_total": 1.0,
        "amr_base": 7.0 / 8.0,
        "amr_fine": 1.0 / 8.0,
    }
    actual = {
        "uniform_total": uniform_total,
        "amr_total": amr_total,
        "amr_base": amr_base,
        "amr_fine": amr_fine,
    }
    checks = {
        name: math.isclose(
            actual[name], value, rel_tol=CIC_RTOL, abs_tol=CIC_ATOL
        )
        for name, value in expected.items()
    }
    return {"actual": actual, "expected": expected, "checks": checks}


def validate_matched_config(path: pathlib.Path) -> dict[str, float]:
    config = json.loads(path.read_text(encoding="utf-8"))
    if (
        config.get("schema") != "matched-coarse-global-k2-a1_16"
        or config.get("profile") != "matched-coarse-v3"
    ):
        raise SourceDiagnosticError(f"{path}: matched profile/schema differs")
    raw = config.get("matched_coarse_cic")
    if not isinstance(raw, dict):
        raise SourceDiagnosticError(f"{path}: matched coarse CIC metrics are absent")
    required = (
        "l1",
        "l2",
        "linf",
        "uniform_integral",
        "matched_integral",
        "hard_ceiling",
    )
    try:
        metrics = {name: float(raw[name]) for name in required}
    except (KeyError, TypeError, ValueError) as error:
        raise SourceDiagnosticError(f"{path}: malformed matched coarse metrics") from error
    if not all(math.isfinite(value) for value in metrics.values()):
        raise SourceDiagnosticError(f"{path}: non-finite matched coarse metric")
    if metrics["hard_ceiling"] != 1.0e-6 or metrics["linf"] > metrics["hard_ceiling"]:
        raise SourceDiagnosticError(f"{path}: matched coarse ceiling is not satisfied")
    for name in ("uniform_integral", "matched_integral"):
        if not math.isclose(metrics[name], 1.0, rel_tol=CIC_RTOL, abs_tol=CIC_ATOL):
            raise SourceDiagnosticError(f"{path}: {name} differs from one")
    return metrics


def rows_converged(
    rows: list[dict[str, float | int | str]], cap: int = 800
) -> bool:
    return len(rows) == 2 and all(
        row["outcome"] == "converged"
        and float(row["residual"]) < 1.0e-4
        and int(row["iterations"]) <= cap
        for row in rows
    )


def assess_matched_solver(
    uniform_solver: dict[str, list[dict[str, float | int | str]]],
    amr_solver: dict[str, list[dict[str, float | int | str]]],
) -> dict[str, object]:
    first_l5 = amr_solver["5"][0]
    solver_pass = (
        first_l5["outcome"] == "converged"
        and float(first_l5["residual"]) < 1.0e-4
        and int(first_l5["iterations"]) <= 800
    )
    return {
        "status": "L5_SOLVER_PASS" if solver_pass else "L5_SOLVER_FAIL",
        "uniform_control_valid": rows_converged(uniform_solver["5"]),
        "first_amr_l5": first_l5,
        "l6_and_second_l5_are_characterization_only": True,
    }


def validate(
    uniform_dir: pathlib.Path,
    amr_dir: pathlib.Path,
    uniform_amr_report: pathlib.Path,
    amr_amr_report: pathlib.Path,
    cap: int = 1,
    profile: str = "projection-v2",
    ic_config: pathlib.Path | None = None,
) -> tuple[dict[str, object], list[str]]:
    matched = profile == "matched-coarse-v3"
    if matched != (cap == 800):
        raise SourceDiagnosticError("matched profile requires cap=800; projection requires cap=1")
    validate_nml_pair(uniform_dir, amr_dir, cap)
    uniform_solver = validate_log(uniform_dir, {(1, 5), (2, 5)}, (5,), cap)
    amr_solver = validate_log(
        amr_dir, {(1, 5), (2, 5), (1, 6), (2, 6)}, (5, 6), cap
    )
    validate_amr_report(uniform_amr_report, uniform_dir / "output_00001", 0)
    validate_amr_report(amr_amr_report, amr_dir / "output_00001", 4096)
    uniform_particles = validate_outputs(uniform_dir, set(range(1, 32**3 + 1)))
    amr_particles = validate_outputs(amr_dir, pinned_particle_ids())

    uniform_l5 = read_source_level(uniform_dir, 5)
    amr_l5 = read_source_level(amr_dir, 5)
    amr_l6 = read_source_level(amr_dir, 6)
    unexpected = sorted(amr_dir.glob("ndgp_source_l???_rank*.tsv"))
    expected_paths = sorted(amr_dir.glob("ndgp_source_l005_rank*.tsv")) + sorted(
        amr_dir.glob("ndgp_source_l006_rank*.tsv")
    )
    if unexpected != sorted(expected_paths):
        raise SourceDiagnosticError(f"{amr_dir}: unexpected source level dump")
    if list(uniform_dir.glob("ndgp_source_l006_rank*.tsv")):
        raise SourceDiagnosticError(f"{uniform_dir}: uniform run dumped L6")

    refined, interface, exterior, refined_boundary, refined_interior = validate_geometry(
        uniform_l5, amr_l5, amr_l6
    )
    all_rho_tot = rho_tot_value(uniform_l5, amr_l5, amr_l6)
    uniform_mean = math.fsum(cell.rho for cell in uniform_l5.values()) / 32**3
    amr_l5_mean = math.fsum(cell.rho for cell in amr_l5.values()) / 32**3
    amr_composite_mean = (
        math.fsum(amr_l5[key].rho for key in set(amr_l5) - refined)
        + math.fsum(cell.rho for cell in amr_l6.values()) / 8.0
    ) / 32**3
    uniform_direct = direct_cic_report(
        uniform_l5, cic_project(uniform_particles, 32, set(uniform_l5)), 32
    )
    amr_l5_direct = direct_cic_report(
        amr_l5, cic_project(amr_particles, 32, set(amr_l5)), 32
    )
    fine_particles = {
        identity: particle
        for identity, particle in amr_particles.items()
        if 32**3 < identity <= 2 * 32**3
    }
    if set(fine_particles) != set(range(32**3 + 1, 2 * 32**3 + 1)):
        raise SourceDiagnosticError("AMR fine-particle ID set differs")
    amr_l6_direct = direct_cic_report(
        amr_l6, cic_project(fine_particles, 64, set(amr_l6)), 64
    )
    direct_reports = {
        "uniform_l5": uniform_direct,
        "amr_l5": amr_l5_direct,
        "amr_l6_active_patch": amr_l6_direct,
    }
    invariant_failures = [
        name
        for name, direct in direct_reports.items()
        if not bool(direct["within_tolerance"])
    ]
    if not math.isclose(uniform_mean, all_rho_tot, rel_tol=CIC_RTOL, abs_tol=CIC_ATOL):
        invariant_failures.append("uniform_l5_mean_vs_rho_tot")
    if not math.isclose(amr_l5_mean, all_rho_tot, rel_tol=CIC_RTOL, abs_tol=CIC_ATOL):
        invariant_failures.append("amr_l5_mean_vs_rho_tot")
    fixture_masses = fixture_mass_report(uniform_particles, amr_particles)
    invariant_failures.extend(
        f"fixture_mass_{name}"
        for name, passed in fixture_masses["checks"].items()
        if not passed
    )

    map_difference = {
        key: uniform_l5[key].rho - amr_l5[key].rho for key in uniform_l5
    }
    restriction_difference: dict[Coord, float] = {}
    for x, y, z in refined:
        child_mean = math.fsum(
            amr_l6[(2 * x + bx, 2 * y + by, 2 * z + bz)].rho
            for bx, by, bz in itertools.product((0, 1), repeat=3)
        ) / 8.0
        restriction_difference[(x, y, z)] = amr_l5[(x, y, z)].rho - child_mean

    map_difference_metrics = metric_summary(map_difference)
    matched_config: dict[str, float] | None = None
    scientific_assessment: dict[str, object] | None = None
    if matched:
        expected_config = uniform_dir.parent / "ic" / "ic_config.json"
        if ic_config is None or ic_config.resolve() != expected_config.resolve():
            raise SourceDiagnosticError("matched IC config path differs from the job IC")
        matched_config = validate_matched_config(ic_config)
        if float(map_difference_metrics["linf"]) > matched_config["hard_ceiling"]:
            invariant_failures.append("runtime_matched_l5_source_linf")
        scientific_assessment = assess_matched_solver(uniform_solver, amr_solver)
        if not bool(scientific_assessment["uniform_control_valid"]):
            invariant_failures.append("uniform_l5_solver_control")

    report: dict[str, object] = {
        "schema": "lagRamses-ndgp-amr-source-diagnostic-v2",
        "status": "DIAGNOSTIC_COMPLETE" if not invariant_failures else "INVARIANT_FAIL",
        "profile": profile,
        "cap": cap,
        "source_map_difference_role": (
            "matched-hard-gate" if matched else "characterization-only"
        ),
        "direct_cic_tolerance": {"relative": CIC_RTOL, "absolute": CIC_ATOL},
        "rho_tot": all_rho_tot,
        "projection_metrics": {
            "uniform_l5_mean": uniform_mean,
            "amr_l5_mean": amr_l5_mean,
            "amr_composite_mean": amr_composite_mean,
            "amr_composite_is_not_a_conservation_oracle": True,
        },
        "fixture_particle_masses": fixture_masses,
        "direct_cic_oracle": direct_reports,
        "matched_generator_metrics": matched_config,
        "scientific_assessment": scientific_assessment,
        "counts": {
            "uniform_l5": len(uniform_l5),
            "amr_l5": len(amr_l5),
            "amr_l6": len(amr_l6),
            "refined_parents": len(refined),
            "interface_exterior_cells": len(interface),
            "exterior_cells": len(exterior),
        },
        "solver_rows": {"uniform": uniform_solver, "amr": amr_solver},
        "uniform_minus_amr_l5": {
            "refined": metric_summary({key: map_difference[key] for key in refined}),
            "refined_boundary": metric_summary(
                {key: map_difference[key] for key in refined_boundary}
            ),
            "refined_interior": metric_summary(
                {key: map_difference[key] for key in refined_interior}
            ),
            "unrefined_interface": metric_summary(
                {key: map_difference[key] for key in interface}
            ),
            "exterior": metric_summary({key: map_difference[key] for key in exterior}),
            "all": map_difference_metrics,
        },
        "amr_l5_minus_mean_children": {
            "refined_all": metric_summary(restriction_difference),
            "refined_boundary": metric_summary(
                {key: restriction_difference[key] for key in refined_boundary}
            ),
            "refined_interior": metric_summary(
                {key: restriction_difference[key] for key in refined_interior}
            ),
        },
        "invariant_failures": invariant_failures,
    }
    return report, invariant_failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("uniform_dir", type=pathlib.Path)
    parser.add_argument("amr_dir", type=pathlib.Path)
    parser.add_argument("--uniform-amr-report", type=pathlib.Path, required=True)
    parser.add_argument("--amr-amr-report", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cap", type=int, choices=(1, 800), default=1)
    parser.add_argument(
        "--profile",
        choices=("projection-v2", "matched-coarse-v3"),
        default="projection-v2",
    )
    parser.add_argument("--ic-config", type=pathlib.Path)
    args = parser.parse_args()
    try:
        report, invariant_failures = validate(
            args.uniform_dir.resolve(),
            args.amr_dir.resolve(),
            args.uniform_amr_report.resolve(),
            args.amr_amr_report.resolve(),
            args.cap,
            args.profile,
            args.ic_config.resolve() if args.ic_config else None,
        )
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        if invariant_failures:
            print(
                "nDGP source diagnostic: invariant failure: "
                + ", ".join(invariant_failures),
                file=sys.stderr,
            )
            return 1
    except (SourceDiagnosticError, FormatError, OSError, ValueError, KeyError, TypeError) as error:
        print(f"nDGP source diagnostic: ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
