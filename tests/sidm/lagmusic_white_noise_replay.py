#!/usr/bin/env python3
"""Manual, non-production MPI replay anchored to one external white-noise file."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import struct
import subprocess
from pathlib import Path

import numpy as np

from lagmusic_small_gate import FIELDS, compare_fields, read_field


CASES = {"external_r1": 1, "external_r4": 4}
LEVELS = (4, 5)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def compact_noise_metadata(path: Path) -> dict:
    raw = path.read_bytes()
    if len(raw) < 12:
        raise ValueError("white-noise input is shorter than its compact header")
    dimensions = list(struct.unpack("=3I", raw[:12]))
    cells = int(np.prod(dimensions))
    payload = len(raw) - 12
    if payload == 8 * cells:
        precision = "float64"
    elif payload == 4 * cells:
        precision = "float32"
    else:
        raise ValueError("white-noise input is not a compact MUSIC cube")
    return {"path": str(path.resolve()), "sha256": digest(path), "bytes": len(raw),
            "dimensions": dimensions, "precision": precision}


def replace_assignment(config: str, key: str, value: str) -> str:
    pattern = rf"(?m)^{re.escape(key)}\s*=.*$"
    result, changed = re.subn(pattern, f"{key} = {value}", config)
    if changed != 1:
        raise ValueError(f"expected one assignment for {key}, found {changed}")
    return result


def make_config(template: str, case: Path, white_noise: Path) -> str:
    config = replace_assignment(template, "seed[4]", str(white_noise.resolve()))
    config = replace_assignment(config, "filename", str(case / "ic"))
    setup = "[setup]\nslab_rng_refinement = no\nrng_boost_threads = no\n"
    if "[setup]\n" not in config:
        raise ValueError("MUSIC template has no [setup] section")
    return config.replace("[setup]\n", setup, 1)


def prepare(root: Path, template: Path, binary: Path, white_noise: Path,
            omp_threads: int) -> None:
    if root.exists():
        raise FileExistsError(f"refusing to overwrite existing diagnostic root: {root}")
    if not binary.is_file() or not binary.stat().st_mode & 0o111:
        raise ValueError(f"MUSIC binary is not executable: {binary}")
    if not template.is_file():
        raise ValueError(f"missing template: {template}")
    if not white_noise.is_file():
        raise ValueError(f"missing white-noise input: {white_noise}")
    if omp_threads < 1:
        raise ValueError("OMP thread count must be positive")
    root.mkdir(parents=True)
    source = template.read_text()
    anchor = compact_noise_metadata(white_noise)
    for label in CASES:
        case = root / label
        case.mkdir()
        config = make_config(source, case, white_noise)
        (case / "music.conf").write_text(config)
    manifest = {
        "purpose": "manual external-white-noise MPI replay only",
        "production_release": False,
        "binary": {"path": str(binary.resolve()), "sha256": digest(binary)},
        "template": {"path": str(template.resolve()), "sha256": digest(template)},
        "white_noise_anchor": anchor,
        "cases": {label: {"ranks": ranks,
                            "config_sha256": digest(root / label / "music.conf")}
                  for label, ranks in CASES.items()},
        "required_log_markers": ["Random numbers for level   4 will be read from file.",
                                 "Reading compact MUSIC white noise file"],
        "explicit_switches": {"disk_cached": "yes", "slab_rng_refinement": "no",
                              "rng_boost_threads": "no"},
        "execution": {"omp_num_threads": omp_threads,
                      "omp_dynamic": "FALSE",
                      "reason": "separate MPI reproducibility from OpenMP threading"},
        "scope": "external L4 phase anchor; not a production IC or science release",
    }
    write_json(root / "manifest.json", manifest)


def run(root: Path, binary: Path, timeout_seconds: int) -> None:
    manifest = json.loads((root / "manifest.json").read_text())
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(manifest["execution"]["omp_num_threads"])
    environment["OMP_DYNAMIC"] = manifest["execution"]["omp_dynamic"]
    for label, ranks in CASES.items():
        case = root / label
        config = case / "music.conf"
        if not config.is_file():
            raise FileNotFoundError(f"run prepare first: {config}")
        command = ["nice", "-n", "19", "timeout", str(timeout_seconds),
                   "mpirun", "--bind-to", "none", "-np", str(ranks), str(binary), "music.conf"]
        with (case / "music.log").open("w") as log:
            subprocess.run(command, cwd=case, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True)


def audit(root: Path) -> None:
    manifest = json.loads((root / "manifest.json").read_text())
    anchor = manifest["white_noise_anchor"]
    input_path = Path(anchor["path"])
    if compact_noise_metadata(input_path) != anchor:
        raise ValueError("external white-noise anchor changed after preparation")
    report = {"production_release": False, "phase_anchor_pass": True,
              "mpi_output_consistency_pass": True, "cases": {}, "fields": {}}
    for label in CASES:
        case = root / label
        config = (case / "music.conf").read_text()
        log = (case / "music.log").read_text()
        markers = {marker: marker in log for marker in manifest["required_log_markers"]}
        config_anchor = f"seed[4] = {anchor['path']}" in config
        copied_noise = case / "wnoise_0004.bin"
        local_noise_sha256 = digest(copied_noise) if copied_noise.is_file() else None
        # store_rnd writes the post-load working field, which can include the
        # generator's mean/coordinate convention.  It is not required to be
        # byte-identical to the immutable input.  The independent-rank check
        # below is the relevant cache invariant.
        passed = all(markers.values()) and config_anchor and local_noise_sha256 is not None
        report["cases"][label] = {"log_markers": markers, "config_anchor": config_anchor,
                                  "local_wnoise_sha256": local_noise_sha256,
                                  "anchor_sha256": anchor["sha256"],
                                  "stored_wnoise_matches_input": local_noise_sha256 == anchor["sha256"],
                                  "phase_anchor_pass": passed}
        report["phase_anchor_pass"] &= passed
    local_hashes = {row["local_wnoise_sha256"] for row in report["cases"].values()}
    report["local_wnoise_rank_identical"] = len(local_hashes) == 1 and None not in local_hashes
    report["phase_anchor_pass"] &= report["local_wnoise_rank_identical"]
    for level in LEVELS:
        for field in FIELDS:
            _, one = read_field(root / "external_r1" / "ic" / f"level_{level:03d}" / field)
            _, four = read_field(root / "external_r4" / "ic" / f"level_{level:03d}" / field)
            key = f"L{level}/{field}"
            try:
                metrics = {"pass": True, **compare_fields(one, four)}
            except ValueError as exc:
                delta = one.astype(np.float64) - four.astype(np.float64)
                metrics = {"pass": False, "failure": str(exc),
                           "relative_l2": float(np.linalg.norm(delta)) /
                                          max(float(np.linalg.norm(one.astype(np.float64))), 1e-30),
                           "max_abs": float(np.max(np.abs(delta))), "bitwise_equal": False}
            report["fields"][key] = metrics
            report["mpi_output_consistency_pass"] &= metrics["pass"]
    write_json(root / "white_noise_replay.json", report)
    if not report["phase_anchor_pass"] or not report["mpi_output_consistency_pass"]:
        raise ValueError("external-white-noise replay gate failed; see white_noise_replay.json")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "run", "audit"))
    parser.add_argument("root", type=Path)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--white-noise", type=Path)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--host", default="LagEunha")
    args = parser.parse_args()
    if args.stage == "prepare":
        if not all((args.template, args.binary, args.white_noise)):
            parser.error("prepare requires --template, --binary, and --white-noise")
        prepare(args.root, args.template, args.binary, args.white_noise, args.omp_threads)
    else:
        if socket.gethostname().split(".")[0] != args.host:
            raise RuntimeError(f"manual replay requires {args.host}")
        if args.stage == "run":
            if not args.binary:
                parser.error("run requires --binary")
            run(args.root, args.binary, args.timeout)
        else:
            audit(args.root)


if __name__ == "__main__":
    main()
