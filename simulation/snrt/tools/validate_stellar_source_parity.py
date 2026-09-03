#!/usr/bin/env python3
"""Fail closed unless the production stellar source is the tested source.

This gate is intentionally static.  It parses the production Makefile, the
production source-order sidecar, and the G1/harness scripts, resolves
source/object names under the declared roots, and hashes the production
inputs recorded by a linked-build evidence file.  It does not treat an old
object in ``build/`` as evidence and does not compile or launch RAMSES.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = ROOT / "simulation/snrt/config/stellar_source_identity_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8", errors="replace")


def git_head() -> str:
    return subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, text=True
    ).strip()


def git_is_ancestor(ancestor: str, descendant: str) -> bool:
    """Return whether ``ancestor`` is at or before ``descendant``."""

    if not ancestor or not descendant:
        return False
    result = subprocess.run(
        ("git", "merge-base", "--is-ancestor", ancestor, descendant),
        cwd=ROOT,
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def git_status() -> list[str]:
    return subprocess.check_output(
        ("git", "status", "--short"), cwd=ROOT, text=True
    ).splitlines()


def git_status_without(paths: list[Path]) -> list[str]:
    ignored = {
        str(path.resolve().relative_to(ROOT.resolve())) for path in paths
    }
    retained = []
    for line in git_status():
        relative_path = line[3:].strip() if len(line) >= 4 else line.strip()
        if relative_path not in ignored:
            retained.append(line)
    return retained


def make_assignment(text: str, name: str) -> str:
    """Return one make assignment, including backslash continuations."""

    lines = text.splitlines()
    for index, line in enumerate(lines):
        match = re.match(rf"^\s*{re.escape(name)}\s*(?:\?=|:=|\+=|=)\s*(.*)$", line)
        if not match:
            continue
        parts = [match.group(1)]
        while parts[-1].rstrip().endswith("\\"):
            next_index = index + len(parts)
            if next_index >= len(lines):
                break
            parts[-1] = parts[-1].rstrip()[:-1]
            parts.append(lines[next_index].strip())
        return " ".join(parts)
    return ""


def make_scalar(text: str, name: str) -> str:
    match = re.search(
        rf"(?m)^\s*{re.escape(name)}\s*(?:\?=|:=|\+=|=)\s*([^\s#]+)", text
    )
    return match.group(1) if match else ""


def make_database(makefile_path: Path) -> tuple[bool, str]:
    """Read GNU make's expanded variable database without running recipes."""

    result = subprocess.run(
        ("make", "-C", str(makefile_path.parent), "-pn", "ramses"),
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return result.returncode == 0, result.stdout


def make_database_assignment(text: str, name: str) -> str | None:
    matches = re.findall(
        rf"(?m)^\s*{re.escape(name)}\s*(?:\?=|:=|\+=|=)\s*(.*)$", text
    )
    return matches[-1] if matches else None


def expand_make_variables(value: str, variables: dict[str, str]) -> str:
    expanded = value
    for _ in range(10):
        next_value = re.sub(
            r"\$\(([A-Za-z_][A-Za-z0-9_]*)\)",
            lambda match: variables.get(match.group(1), ""),
            expanded,
        )
        next_value = re.sub(
            r"\$\(filter-out\s+([^,\s]+),([^)]*)\)",
            lambda match: " ".join(
                item
                for item in match.group(2).split()
                if item != match.group(1)
            ),
            next_value,
        )
        if next_value == expanded:
            return next_value
        expanded = next_value
    return expanded


def production_build_input_manifest(
    makefile_path: Path,
    make_db_text: str,
    make_db_variables: dict[str, str],
    assigned_objects: list[str],
    generated_objects: list[str],
) -> dict[str, Any]:
    """Resolve and hash every source selected by the production link inputs."""

    vpath_raw = make_database_assignment(make_db_text, "VPATH") or ""
    vpath = expand_make_variables(vpath_raw, make_db_variables)
    vpath_directories = [
        (makefile_path.parent / directory).resolve()
        for directory in vpath.split(":")
        if directory
    ]
    source_paths: dict[str, str] = {}
    missing_sources: list[str] = []
    for object_name in sorted(set(assigned_objects + ["ramses.o"])):
        if object_name in generated_objects:
            continue
        source_name = fortran_source_for_object(object_name)
        for directory in vpath_directories:
            candidate = directory / source_name
            if candidate.is_file():
                source_paths[object_name] = str(candidate.relative_to(ROOT))
                break
        else:
            missing_sources.append(source_name)

    source_hashes = {
        relative_path: sha256(ROOT / relative_path)
        for relative_path in sorted(set(source_paths.values()))
    }
    tree_digest = hashlib.sha256()
    for relative_path in sorted(source_hashes):
        tree_digest.update(relative_path.encode("utf-8"))
        tree_digest.update(b"\0")
        with (ROOT / relative_path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                tree_digest.update(chunk)
    return {
        "source_paths": source_paths,
        "source_hashes": source_hashes,
        "missing_sources": sorted(set(missing_sources)),
        "tree_sha256": tree_digest.hexdigest(),
    }


def fortran_source_for_object(object_name: str) -> str:
    if not object_name.endswith(".o"):
        raise ValueError(f"not an object name: {object_name}")
    return object_name[:-2] + ".f90"


def source_names_from_runner(
    text: str, array_name: str, extension: str
) -> list[str]:
    match = re.search(
        rf"(?ms)^\s*{re.escape(array_name)}=\(\s*(.*?)^\s*\)", text
    )
    if not match:
        return []
    return re.findall(rf"[A-Za-z0-9_.-]+\.{re.escape(extension)}", match.group(1))


def runner_source_dir(text: str) -> str:
    match = re.search(r"(?m)^\s*SOURCE_DIR=\"([^\"]+)\"", text)
    return match.group(1) if match else ""


def _define_tokens(command: str) -> dict[str, str]:
    defines: dict[str, str] = {}
    for name, value in re.findall(
        r"(?<!\S)-D([A-Za-z_][A-Za-z0-9_]*)(?:=([^\s]+))?", command
    ):
        defines[name] = value or "1"
    return defines


def production_build_log_contract(
    log_text: str,
    required_objects: list[str],
    compile_parameters: dict[str, Any],
    required_compile_flags: list[str],
    embedded_macro: str,
    build_command: str,
    source_prefix: str = "../patch/lagRamses",
    required_optimization_flag: str = "-O3",
    forbidden_compile_flags: tuple[str, ...] = ("-DFDMDEBUG",),
    expected_link_output: str = "ramses_final3d",
    expected_binary_sha256: str | None = None,
) -> dict[str, Any]:
    """Validate the fresh-build claims recorded in the production log."""

    commands: dict[str, list[str]] = {}
    for object_name in required_objects:
        source_name = fortran_source_for_object(object_name)
        pattern = re.compile(
            rf"(?m)^\s*(?P<command>.*\s-c\s+{re.escape(source_prefix)}/{re.escape(source_name)}"
            rf"\s+-o\s+(?:\S+/)?{re.escape(object_name)}\s*)$"
        )
        commands[object_name] = [
            match.group("command").strip() for match in pattern.finditer(log_text)
        ]

    compile_commands = [
        values[0]
        for object_name, values in commands.items()
        if len(values) == 1
    ]
    observations = [
        {
            name: _define_tokens(command).get(name)
            for name in compile_parameters
        }
        for command in compile_commands
    ]
    observed_parameters = observations[0] if observations else {
        name: None for name in compile_parameters
    }
    compile_parameters_consistent = bool(observations) and all(
        observation == observations[0] for observation in observations
    )
    expected_parameters = {
        name: str(value) for name, value in compile_parameters.items()
    }
    compile_parameters_match = (
        compile_parameters_consistent and observed_parameters == expected_parameters
    )
    required_flags_present = bool(compile_commands) and all(
        all(flag in command.split() for flag in required_compile_flags)
        for command in compile_commands
    )
    embedded_yields_disabled = bool(compile_commands) and all(
        embedded_macro not in _define_tokens(command) for command in compile_commands
    )
    optimization_policy_satisfied = bool(compile_commands) and all(
        required_optimization_flag in command.split()
        and all(flag not in command.split() for flag in forbidden_compile_flags)
        for command in compile_commands
    )
    build_command_marker = f"P0_BUILD_COMMAND {build_command}"
    build_command_matches = any(
        line.strip() == build_command_marker for line in log_text.splitlines()
    )
    link_output_matches = any(
        re.search(
            rf"\s-o\s+{re.escape(expected_link_output)}(?:\s|$)", line
        )
        for line in log_text.splitlines()
    )
    binary_sha_matches = re.findall(
        r"(?m)^P0_BINARY_SHA256=([0-9a-f]{64})\s*$", log_text
    )
    logged_binary_sha256 = binary_sha_matches[-1] if binary_sha_matches else None
    binary_sha256_matches = (
        expected_binary_sha256 is None
        or logged_binary_sha256 == expected_binary_sha256
    )
    forced_rebuild = any(
        line.strip() == build_command_marker and " -B " in f" {line.strip()} "
        for line in log_text.splitlines()
    )
    required_objects_compiled = all(
        len(commands[object_name]) == 1 for object_name in required_objects
    )
    return {
        "status": (
            "pass"
            if required_objects_compiled
            and compile_parameters_match
            and required_flags_present
            and embedded_yields_disabled
            and optimization_policy_satisfied
            and build_command_matches
            and link_output_matches
            and binary_sha256_matches
            and forced_rebuild
            else "blocked"
        ),
        "required_objects_compiled": required_objects_compiled,
        "required_object_compile_counts": {
            object_name: len(values) for object_name, values in commands.items()
        },
        "compile_parameters": observed_parameters,
        "compile_parameters_consistent": compile_parameters_consistent,
        "compile_parameters_match": compile_parameters_match,
        "required_compile_flags_present": required_flags_present,
        "embedded_yields_disabled": embedded_yields_disabled,
        "optimization_policy_satisfied": optimization_policy_satisfied,
        "build_command_matches": build_command_matches,
        "link_output_matches": link_output_matches,
        "logged_binary_sha256": logged_binary_sha256,
        "expected_binary_sha256": expected_binary_sha256,
        "binary_sha256_matches": binary_sha256_matches,
        "forced_rebuild": forced_rebuild,
    }


def binary_linkage_from_text(
    nm_text: str, nm_returncode: int, symbol_patterns: list[str]
) -> dict[str, Any]:
    patterns = {
        pattern: any(
            re.search(pattern, line, re.IGNORECASE) is not None
            and re.search(r"\sT\s", line) is not None
            for line in nm_text.splitlines()
        )
        for pattern in symbol_patterns
    }
    return {
        "status": "pass"
        if nm_returncode == 0 and all(patterns.values())
        else "blocked",
        "nm_returncode": nm_returncode,
        "patterns": patterns,
        "text_symbols_required": True,
    }


def binary_linkage_contract(
    binary: Path, symbol_patterns: list[str]
) -> dict[str, Any]:
    """Verify that the production executable exports configured stellar symbols."""

    if not binary.is_file():
        return {
            "status": "blocked",
            "nm_returncode": None,
            "patterns": {pattern: False for pattern in symbol_patterns},
            "text_symbols_required": True,
        }
    try:
        result = subprocess.run(
            ("nm", "--defined-only", str(binary)),
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return {
            "status": "blocked",
            "nm_returncode": None,
            "patterns": {pattern: False for pattern in symbol_patterns},
            "text_symbols_required": True,
        }
    return binary_linkage_from_text(
        result.stdout, result.returncode, symbol_patterns
    )


def production_smoke_contract(
    log_text: str,
    smoke_command: str,
    expected_output: str,
    expected_exit_code: int,
    expected_patch_dir: str,
    expected_repository_head: str,
    expected_binary_sha256: str | None = None,
) -> dict[str, Any]:
    """Validate the no-argument executable startup smoke record."""

    exit_matches = re.findall(r"(?m)^P0_SMOKE_EXIT_CODE=(\d+)\s*$", log_text)
    exit_code = int(exit_matches[-1]) if exit_matches else None
    command_marker = f"P0_SMOKE_COMMAND {smoke_command}"
    command_matches = any(
        line.strip() == command_marker for line in log_text.splitlines()
    )
    expected_output_found = expected_output in log_text
    patch_dir_found = f"patch dir    = {expected_patch_dir}" in log_text
    logged_heads = re.findall(
        r"(?m)^\s*last commit\s+=\s+([0-9a-f]{40})(?:-dirty)?\s*$",
        log_text,
    )
    repository_head_found = any(
        logged_head == expected_repository_head
        or git_is_ancestor(logged_head, expected_repository_head)
        for logged_head in logged_heads
    )
    binary_sha_matches = re.findall(
        r"(?m)^P0_BINARY_SHA256=([0-9a-f]{64})\s*$", log_text
    )
    logged_binary_sha256 = binary_sha_matches[-1] if binary_sha_matches else None
    binary_sha256_matches = (
        expected_binary_sha256 is None
        or logged_binary_sha256 == expected_binary_sha256
    )
    return {
        "status": "pass"
        if command_matches
        and exit_code == expected_exit_code
        and expected_output_found
        and patch_dir_found
        and repository_head_found
        and binary_sha256_matches
        else "blocked",
        "command_matches": command_matches,
        "exit_code": exit_code,
        "expected_exit_code": expected_exit_code,
        "expected_output_found": expected_output_found,
        "patch_dir_found": patch_dir_found,
        "repository_head_found": repository_head_found,
        "logged_repository_head": logged_heads[-1] if logged_heads else None,
        "logged_binary_sha256": logged_binary_sha256,
        "expected_binary_sha256": expected_binary_sha256,
        "binary_sha256_matches": binary_sha256_matches,
    }


def active_fflags(text: str) -> str:
    """Return the non-debug FFLAGS branch, not text from the debug branch."""

    match = re.search(
        r"(?ms)^ifdef\s+FDMDEBUG\s*\n.*?^else\s*\n(?P<active>.*?)^endif\s*$",
        text,
    )
    return match.group("active") if match else ""


def relative_root_from_runner(value: str) -> str:
    if value.startswith("$ROOT/"):
        return value[len("$ROOT/") :]
    return value


def audit(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    production = config["production"]
    g1 = config["g1_runner"]
    harness = config["production_linked_harness"]
    makefile_path = ROOT / production["makefile"]
    runner_path = ROOT / g1["script"]
    makefile = makefile_path.read_text(encoding="utf-8", errors="replace")
    runner = runner_path.read_text(encoding="utf-8", errors="replace")

    patch_root = production["patch_root"]
    declared_patch = make_scalar(makefile, "PATCH")
    expected_patch = str(Path("..") / patch_root)
    # The Makefile is evaluated from bin/, so PATCH=../patch/lagRamses is the
    # expected spelling even though the config stores a repository-relative root.
    declared_patch_repo_relative = str(
        (makefile_path.parent / declared_patch).resolve()
        .relative_to(ROOT.resolve())
    ) if declared_patch else ""

    make_db_ok, make_db_text = make_database(makefile_path)
    make_db_variables: dict[str, str] = {}
    for name, _ in re.findall(
        r"(?m)^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:\?=|:=|\+=|=)\s*(.*)$",
        make_db_text,
    ):
        make_db_variables[name] = make_database_assignment(make_db_text, name) or ""
    assignment_values = {}
    assignment_names = list(
        dict.fromkeys(
            production["object_assignments"]
            + production["stellar_order_assignments"]
        )
    )
    for name in assignment_names:
        raw_value = (
            make_database_assignment(make_db_text, name)
            if make_db_ok
            else None
        )
        if raw_value is None:
            raw_value = make_assignment(makefile, name)
        assignment_values[name] = expand_make_variables(
            raw_value, make_db_variables
        )
    assigned_objects = sorted(
        {
            object_name
            for value in assignment_values.values()
            for object_name in re.findall(r"[A-Za-z0-9_.-]+\.o", value)
        }
    )
    required_objects = production["required_objects"]
    required_sources = [fortran_source_for_object(name) for name in required_objects]
    missing_production_sources = [
        source for source in required_sources if not (ROOT / patch_root / source).is_file()
    ]
    missing_production_objects = [
        name for name in required_objects if name not in assigned_objects
    ]
    build_input_manifest = production_build_input_manifest(
        makefile_path,
        make_db_text,
        make_db_variables,
        assigned_objects,
        production.get("generated_objects", []),
    )
    source_manifest = production["source_manifest"]
    source_manifest_path = ROOT / source_manifest["path"]
    source_manifest_text = (
        source_manifest_path.read_text(encoding="utf-8", errors="replace")
        if source_manifest_path.is_file()
        else ""
    )
    manifest_sources = sorted(
        {
            Path(match).name
            for match in re.findall(
                r"patch/lagRamses/[A-Za-z0-9_.-]+\.f90", source_manifest_text
            )
        }
    )
    manifest_source_order = [
        Path(match).name
        for match in re.findall(
            r"patch/lagRamses/[A-Za-z0-9_.-]+\.f90", source_manifest_text
        )
    ]
    expected_manifest_sources = sorted(
        source for source in required_sources if source.startswith("stellar_")
    )
    source_manifest_matches = (
        source_manifest_path.is_file()
        and manifest_sources == expected_manifest_sources
        and all(source in required_sources for source in manifest_sources)
    )
    manifest_object_order = [source[:-4] + ".o" for source in manifest_source_order]
    production_stellar_object_order = [
        object_name
        for name in production["stellar_order_assignments"]
        for object_name in re.findall(r"[A-Za-z0-9_.-]+\.o", assignment_values[name])
        if object_name.startswith("stellar_")
    ]
    source_manifest_order_matches = (
        manifest_object_order == production_stellar_object_order
    )
    makefile_consumes_source_manifest = bool(
        re.search(
            r"(?m)^\s*include\s+\$\(PATCH\)/stellar_enrichment_sources\.mk\s*$",
            makefile,
        )
    )
    production_link_recipe_complete = bool(
        re.search(
            r"(?m)^ramses:\s*.*\$\(MODOBJ\).*\$\(AMRLIB\).*ramses\.o\s*$",
            makefile,
        )
        and re.search(
            r"(?m)^\s*\$\(F90\).*\$\(MODOBJ\).*\$\(AMRLIB\).*"
            r"-o\s+\$\(EXEC\)\$\(NDIM\)d\s+\$\(LIBS\)\s*$",
            makefile,
        )
    )

    runner_sources = source_names_from_runner(runner, "sources", "f90")
    runner_objects = source_names_from_runner(runner, "objects", "o")
    runner_dir = relative_root_from_runner(runner_source_dir(runner))
    native_root = g1["native_root"]
    missing_native_sources = [
        source for source in runner_sources if not (ROOT / native_root / source).is_file()
    ]
    missing_native_objects = [
        object_name
        for object_name in runner_objects
        if not (ROOT / native_root / fortran_source_for_object(object_name)).is_file()
    ]
    shared_hashes: dict[str, dict[str, str | None]] = {}
    for source in g1["shared_contract_sources"]:
        production_file = ROOT / patch_root / source
        native_file = ROOT / native_root / source
        shared_hashes[source] = {
            "production_sha256": sha256(production_file) if production_file.is_file() else None,
            "native_sha256": sha256(native_file) if native_file.is_file() else None,
        }
    differing_shared_sources = [
        source
        for source, hashes in shared_hashes.items()
        if hashes["production_sha256"] is None
        or hashes["native_sha256"] is None
        or hashes["production_sha256"] != hashes["native_sha256"]
    ]
    absent_shared_sources = [
        source
        for source, hashes in shared_hashes.items()
        if hashes["production_sha256"] is not None
        and hashes["native_sha256"] is None
    ]
    identical_shared_sources = [
        source
        for source, hashes in shared_hashes.items()
        if hashes["production_sha256"] is not None
        and hashes["production_sha256"] == hashes["native_sha256"]
    ]
    content_differing_shared_sources = [
        source
        for source, hashes in shared_hashes.items()
        if hashes["production_sha256"] is not None
        and hashes["native_sha256"] is not None
        and hashes["production_sha256"] != hashes["native_sha256"]
    ]
    shared_contract_source_set = set(g1["shared_contract_sources"])
    shared_contract_partitions = (
        set(absent_shared_sources)
        | set(identical_shared_sources)
        | set(content_differing_shared_sources)
    )
    shared_contract_profile_bounded = (
        all(
            hashes["production_sha256"] is not None
            for hashes in shared_hashes.values()
        )
        and set(absent_shared_sources).issubset(shared_contract_source_set)
        and set(identical_shared_sources).issubset(shared_contract_source_set)
        and set(content_differing_shared_sources).issubset(shared_contract_source_set)
        and not (
            set(absent_shared_sources)
            & set(identical_shared_sources)
        )
        and not (
            set(absent_shared_sources)
            & set(content_differing_shared_sources)
        )
        and not (
            set(identical_shared_sources)
            & set(content_differing_shared_sources)
        )
        and shared_contract_partitions == shared_contract_source_set
    )
    shared_contract_baseline = config["production_linked_harness"].get(
        "shared_contract_baseline", {}
    )
    shared_contract_baseline_matches = (
        sorted(absent_shared_sources)
        == sorted(shared_contract_baseline.get("absent_in_native", []))
        and sorted(identical_shared_sources)
        == sorted(shared_contract_baseline.get("identical_sources", []))
        and sorted(content_differing_shared_sources)
        == sorted(shared_contract_baseline.get("differing_sources", []))
    )

    runner_source_object_match = sorted(runner_sources) == sorted(
        fortran_source_for_object(object_name) for object_name in runner_objects
    )
    native_only_mirror_sources = g1.get("native_only_mirror_sources", [])
    native_only_mirror_sources_declared = (
        all(source in runner_sources for source in native_only_mirror_sources)
        and not set(native_only_mirror_sources).intersection(
            g1["shared_contract_sources"]
        )
        and all(
            source
            not in {
                fortran_source_for_object(object_name)
                for object_name in required_objects
            }
            for source in native_only_mirror_sources
        )
    )
    source_of_truth = config["source_of_truth"]
    source_of_truth_declared = (
        source_of_truth["path"] == patch_root
        and source_of_truth["role"] == "canonical_production_source"
        and source_of_truth["strategy"] == "production_linked_harness"
        and source_of_truth["mirror_role"] == "differential_oracle"
        and source_of_truth["mirror_can_remain_as_differential_oracle"] is True
    )
    compile_parameters = production["compile_parameters"]
    effective_compile_parameters = {
        name: make_scalar(make_db_text, name) if make_db_ok else make_scalar(makefile, name)
        for name in compile_parameters
    }
    compile_parameter_match = all(
        effective_compile_parameters[name] == str(value)
        for name, value in compile_parameters.items()
    )
    required_compile_flags = production["required_compile_flags"]
    active_flags = active_fflags(makefile)
    effective_fflags = expand_make_variables(
        make_database_assignment(make_db_text, "FFLAGS") or active_flags,
        make_db_variables,
    )
    compile_flags_present = all(
        flag in effective_fflags.split() for flag in required_compile_flags
    )
    embedded_macro = production["embedded_yield_policy"]["macro"]
    embedded_yields_default_disabled = (
        not bool(re.search(rf"(?m)^\s*{re.escape(embedded_macro)}\s*(?:\?=|:=|=)", makefile))
        and bool(re.search(rf"(?m)^\s*ifdef\s+{re.escape(embedded_macro)}", makefile))
    )
    compile_policy = production["compile_policy"]
    harness_path = ROOT / harness["script"]
    harness_exists = harness_path.is_file()
    harness_text = harness_path.read_text(encoding="utf-8", errors="replace") if harness_exists else ""
    harness_targets_makefile = bool(
        re.search(
            r"(?m)^\s*make\s+-C\s+\"\$ROOT/bin\"\s+-B\s+ramses(?:\s|$)",
            harness_text,
        )
    )
    harness_object_coverage = all(
        object_name in assigned_objects for object_name in required_objects
    )
    required_stellar_objects = sorted(
        object_name for object_name in required_objects if object_name.startswith("stellar_")
    )
    declared_stellar_objects = sorted(
        object_name for object_name in assigned_objects if object_name.startswith("stellar_")
    )
    production_stellar_objects_match_contract = (
        make_db_ok
        and declared_stellar_objects == required_stellar_objects
    )
    evidence_path = ROOT / harness["evidence_file"]
    evidence: dict[str, Any] | None = None
    evidence_valid = False
    current_source_hashes = {
        source: sha256(ROOT / patch_root / source)
        for source in required_sources
        if (ROOT / patch_root / source).is_file()
    }
    current_harness_hash = sha256(harness_path) if harness_exists else None
    recorder_path = ROOT / harness["recorder"]
    current_validator_hash = sha256(Path(__file__).resolve())
    current_recorder_hash = sha256(recorder_path) if recorder_path.is_file() else None
    expected_binary_path = (ROOT / harness["binary"]).resolve()
    expected_build_log_path = (ROOT / harness["build_log"]).resolve()
    expected_smoke_log_path = (ROOT / harness["smoke_log"]).resolve()
    expected_build_command = harness["build_command"]
    current_linkage_contract = binary_linkage_contract(
        expected_binary_path, harness["linkage_symbol_patterns"]
    )
    current_build_log_contract: dict[str, Any] | None = None
    if expected_build_log_path.is_file():
        try:
            current_build_log_contract = production_build_log_contract(
                expected_build_log_path.read_text(
                    encoding="utf-8", errors="replace"
                ),
                required_objects,
                compile_parameters,
                required_compile_flags,
                embedded_macro,
                expected_build_command,
                expected_patch,
                compile_policy["required_optimization_flag"],
                tuple(compile_policy["forbidden_compile_flags"]),
                harness["link_output"],
                sha256(expected_binary_path)
                if expected_binary_path.is_file()
                else None,
            )
        except OSError:
            current_build_log_contract = None
    current_smoke_contract: dict[str, Any] | None = None
    if expected_smoke_log_path.is_file():
        try:
            current_smoke_contract = production_smoke_contract(
                expected_smoke_log_path.read_text(
                    encoding="utf-8", errors="replace"
                ),
                harness["smoke_command"],
                harness["smoke_expected_output"],
                harness["smoke_expected_exit_code"],
                expected_patch,
                git_head(),
                sha256(expected_binary_path)
                if expected_binary_path.is_file()
                else None,
            )
        except OSError:
            current_smoke_contract = None
    if evidence_path.is_file():
        try:
            evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
            evidence_binary_path = Path(evidence.get("binary_path", ""))
            if not evidence_binary_path.is_absolute():
                evidence_binary_path = ROOT / evidence_binary_path
            evidence_build_log_path = Path(evidence.get("build_log", ""))
            if not evidence_build_log_path.is_absolute():
                evidence_build_log_path = ROOT / evidence_build_log_path
            evidence_binary_path = evidence_binary_path.resolve()
            evidence_build_log_path = evidence_build_log_path.resolve()
            evidence_valid = (
                evidence.get("status") == harness["required_evidence_status"]
                and evidence.get("source_root") == patch_root
                and evidence.get("makefile_sha256") == sha256(makefile_path)
                and sorted(evidence.get("required_objects", []))
                == sorted(required_objects)
                and git_is_ancestor(evidence.get("repository_head", ""), git_head())
                and isinstance(evidence.get("worktree_status"), list)
                and evidence.get("config_sha256") == sha256(config_path)
                and evidence.get("harness_sha256") == current_harness_hash
                and evidence.get("validator_sha256") == current_validator_hash
                and evidence.get("recorder_sha256") == current_recorder_hash
                and evidence.get("runner_sha256") == sha256(runner_path)
                and evidence.get("source_manifest_sha256")
                == (sha256(source_manifest_path) if source_manifest_path.is_file() else None)
                and evidence.get("production_source_hashes") == current_source_hashes
                and evidence.get("build_input_source_hashes")
                == build_input_manifest["source_hashes"]
                and evidence.get("build_input_tree_sha256")
                == build_input_manifest["tree_sha256"]
                and current_build_log_contract is not None
                and current_build_log_contract.get("status") == "pass"
                and evidence.get("build_log_contract") == current_build_log_contract
                and current_linkage_contract.get("status") == "pass"
                and evidence.get("binary_linkage_contract") == current_linkage_contract
                and current_smoke_contract is not None
                and current_smoke_contract.get("status") == "pass"
                and evidence.get("smoke_contract") == current_smoke_contract
                and evidence.get("compile_parameters")
                == current_build_log_contract.get("compile_parameters")
                and evidence.get("required_compile_flags") == required_compile_flags
                and evidence.get("forced_rebuild") is True
                and evidence.get("build_command") == expected_build_command
                and evidence_binary_path == expected_binary_path
                and evidence_binary_path.is_file()
                and evidence.get("binary_sha256") == sha256(evidence_binary_path)
                and evidence_build_log_path == expected_build_log_path
                and evidence_build_log_path.is_file()
                and evidence.get("build_log_sha256") == sha256(evidence_build_log_path)
                and evidence.get("smoke_log") == str(expected_smoke_log_path)
                and expected_smoke_log_path.is_file()
                and evidence.get("smoke_log_sha256") == sha256(expected_smoke_log_path)
            )
        except (OSError, json.JSONDecodeError, TypeError):
            evidence_valid = False
    criteria = {
        "production_makefile_selects_declared_patch": (
            declared_patch_repo_relative == patch_root
            and declared_patch == expected_patch
        ),
        "production_objects_are_listed": not missing_production_objects,
        "production_build_inputs_resolve": not build_input_manifest["missing_sources"],
        "production_stellar_objects_match_contract": production_stellar_objects_match_contract,
        "production_sources_resolve_under_patch": not missing_production_sources,
        "production_source_manifest_exists_and_matches": source_manifest_matches,
        "production_source_manifest_order_matches_make_objects": source_manifest_order_matches,
        "production_makefile_consumes_source_manifest": makefile_consumes_source_manifest,
        "production_link_recipe_is_complete": production_link_recipe_complete,
        "production_shared_contract_profile_bounded": shared_contract_profile_bounded,
        "production_make_database_resolves": make_db_ok,
        "native_runner_source_root_is_declared": runner_dir == native_root,
        "native_runner_sources_resolve": not missing_native_sources,
        "native_runner_objects_resolve": not missing_native_objects,
        "native_runner_sources_match_objects": runner_source_object_match,
        "native_only_mirror_sources_declared": native_only_mirror_sources_declared,
        "production_linked_harness_exists": harness_exists,
        "production_linked_harness_targets_bin_makefile": harness_targets_makefile,
        "production_linked_harness_runs_smoke": bool(
            re.search(
                r"(?m)^\s*\"\$BINARY\"\s+>>\s+\"\$SMOKE_LOG\"\s+2>&1\s*$",
                harness_text,
            )
        ),
        "production_linked_harness_object_coverage": harness_object_coverage,
        "production_linked_build_evidence": evidence_valid,
        "runner_has_required_parity_gate": bool(
            re.search(
                r"(?m)^\s*python3\s+\"\$ROOT/simulation/snrt/tools/validate_stellar_source_parity\.py\"\s+--require-pass\s*$",
                runner,
            )
        ),
        "runner_diagnostic_escape_is_fail_closed": (
            bool(
                re.search(
                    r'(?m)^if \[\[ "\$\{P0_DIAGNOSTIC:-0\}" == 1 \]\];? then$',
                    runner,
                )
            )
            and "G1_NATIVE_DIAGNOSTIC_ONLY" in runner
            and "--require-pass" in runner
        ),
        "source_of_truth_is_canonical_production_tree": source_of_truth_declared,
        "compile_parameters_match_contract": compile_parameter_match,
        "required_compile_flags_present": compile_flags_present,
        "embedded_yields_disabled_by_default": embedded_yields_default_disabled,
    }
    blocking_reasons = [name for name, passed in criteria.items() if not passed]
    return {
        "schema": "snrt-stellar-source-parity-audit-v1",
        "repository_root": str(ROOT),
        "repository_head": git_head(),
        "worktree_status": git_status(),
        "config": str(config_path.relative_to(ROOT)),
        "production": {
            "makefile": production["makefile"],
            "declared_patch": declared_patch,
            "declared_patch_repo_relative": declared_patch_repo_relative,
            "expected_patch_repo_relative": patch_root,
            "assignments": assignment_values,
            "assigned_objects": assigned_objects,
            "required_objects": required_objects,
            "missing_objects": missing_production_objects,
            "missing_sources": missing_production_sources,
            "build_input_manifest": build_input_manifest,
            "source_manifest": source_manifest["path"],
            "source_manifest_order": manifest_source_order,
            "source_manifest_sources": manifest_sources,
            "source_manifest_matches": source_manifest_matches,
            "source_manifest_order_matches": source_manifest_order_matches,
            "production_stellar_object_order": production_stellar_object_order,
            "make_database_resolves": make_db_ok,
            "effective_compile_parameters": effective_compile_parameters,
            "compile_parameters": {
                name: make_scalar(makefile, name)
                for name in production["compile_parameters"]
            },
        },
        "g1_runner": {
            "script": g1["script"],
            "declared_source_dir": runner_source_dir(runner),
            "resolved_source_root": runner_dir,
            "expected_native_root": native_root,
            "sources": runner_sources,
            "objects": runner_objects,
            "missing_sources": missing_native_sources,
            "missing_objects": missing_native_objects,
            "sources_match_objects": runner_source_object_match,
            "native_only_sources": native_only_mirror_sources,
            "native_only_disposition": g1.get("native_only_mirror_disposition", ""),
            "role": source_of_truth["mirror_role"],
        },
        "shared_contract": {
            "hashes": shared_hashes,
            "differing_sources": differing_shared_sources,
            "content_differing_sources": content_differing_shared_sources,
            "absent_in_native": absent_shared_sources,
            "identical_sources": identical_shared_sources,
            "profile_bounded": shared_contract_profile_bounded,
            "baseline_profile_matches": shared_contract_baseline_matches,
            "baseline": shared_contract_baseline,
            "interpretation": "differential diagnostic only; byte identity is not required while the mirror remains an oracle",
        },
        "production_linked_harness": {
            "script": harness["script"],
            "exists": harness_exists,
            "targets_bin_makefile": harness_targets_makefile,
            "required_objects": required_objects,
            "binary": harness["binary"],
            "build_log": harness["build_log"],
            "build_command": harness["build_command"],
            "smoke_log": harness["smoke_log"],
            "smoke_command": harness["smoke_command"],
            "smoke_expected_exit_code": harness["smoke_expected_exit_code"],
            "linkage_symbol_patterns": harness["linkage_symbol_patterns"],
            "current_linkage_contract": current_linkage_contract,
            "current_smoke_contract": current_smoke_contract,
            "evidence_file": harness["evidence_file"],
            "evidence": evidence,
            "evidence_valid": evidence_valid,
            "evidence_validation_inputs": {
                "current_source_hashes": current_source_hashes,
            "current_build_log_contract": current_build_log_contract,
                "current_binary_linkage_contract": current_linkage_contract,
                "current_smoke_contract": current_smoke_contract,
            },
        },
        "compile_contract": {
            "parameters": {
                name: make_scalar(makefile, name)
                for name in compile_parameters
            },
            "required_flags": required_compile_flags,
            "effective_fflags": effective_fflags,
            "required_flags_present": compile_flags_present,
            "embedded_yield_macro": embedded_macro,
            "embedded_yields_default_disabled": embedded_yields_default_disabled,
            "derived_indices_for_metal_delayed_cooling": {
                "inener": int(make_scalar(makefile, "NDIM")) + 3,
                "imetal": int(make_scalar(makefile, "NENER")) + int(make_scalar(makefile, "NDIM")) + 3,
                "idelay": int(make_scalar(makefile, "NENER")) + int(make_scalar(makefile, "NDIM")) + 4,
                "ichem": int(make_scalar(makefile, "NENER")) + int(make_scalar(makefile, "NDIM")) + 5,
            },
        },
        "criteria": criteria,
        "status": "pass" if not blocking_reasons else "blocked",
        "blocking_reasons": blocking_reasons,
        "policy": config["policy"],
        "provenance": {
            "makefile_sha256": sha256(makefile_path),
            "runner_sha256": sha256(runner_path),
            "config_sha256": sha256(config_path),
            "stale_build_objects_are_not_evidence": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()
    payload = audit(args.config)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "STELLAR_SOURCE_PARITY_"
            f"{payload['status'].upper()} "
            f"blocked={','.join(payload['blocking_reasons']) or 'none'}"
        )
    return 0 if payload["status"] == "pass" or not args.require_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
