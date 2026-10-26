#!/usr/bin/env python3
"""Bounded, manual LagEunha diagnosis; does not alter production sources."""
import argparse
import hashlib
import json
from pathlib import Path
import socket
import struct
import subprocess

import numpy as np
from lagmusic_small_gate import FIELDS, read_field


def metrics(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return {"relative_l2": float(np.linalg.norm(a-b)/max(np.linalg.norm(a), 1e-30)),
            "max_abs": float(np.max(np.abs(a-b))),
            "bitwise_equal": bool(np.array_equal(a, b))}


def noise(path):
    raw = path.read_bytes()
    dims = struct.unpack("=3I", raw[:12])
    if np.prod(dims) > 2**21 or len(raw) != 12+8*np.prod(dims):
        raise ValueError(f"unexpected small noise shape: {path}")
    return np.frombuffer(raw, dtype="=f8", offset=12).reshape(dims)


def periodic_last_write(field, period):
    """Model increasing-order serial writeback followed by periodic readout."""
    indices = [np.arange(n) % period + ((n-1-np.arange(n)%period)//period)*period
               for n in field.shape]
    return field[np.ix_(*indices)]


def followup(root, baseline):
    """Remove averaging for alias proof, then test a non-wrapping geometry."""
    root.mkdir(exist_ok=False)
    exe = baseline / "MUSIC"
    original = (baseline / "z6r1/music.conf").read_text()
    report = {"production_release": False, "pairs": {}}
    for label, levels in (("raw_wrap", (4, 5)), ("nonwrap", (5, 6, 7))):
        paths = []
        for ranks in (1, 4):
            case = root / f"{label}_r{ranks}"
            case.mkdir()
            conf = original.replace(str(baseline / "z6r1/ic"), str(case / "ic"))
            if label == "raw_wrap":
                conf = conf.replace("levelmax = 6", "levelmax = 5")
                conf = conf.replace("kspace_TF = no", "kspace_TF = yes")
                conf = conf.replace("[setup]\n", "[setup]\nrng_boost_threads = no\n")
            else:
                conf = conf.replace("levelmin = 4", "levelmin = 5")
                conf = conf.replace("levelmin_TF = 4", "levelmin_TF = 5")
                conf = conf.replace("levelmax = 6", "levelmax = 7")
                conf = conf.replace("seed[4] = 20260401", "seed[7] = 20260404")
            (case / "music.conf").write_text(conf)
            with (case / "music.log").open("w") as log:
                subprocess.run(["timeout", "180", "mpirun", "-np", str(ranks), str(exe), "music.conf"],
                               cwd=case, stdout=log, stderr=subprocess.STDOUT, check=True)
            paths.append(case)
        pair = {"noise": {}, "fields": {}}
        for level in levels:
            a = noise(paths[0] / f"wnoise_{level:04d}.bin")
            b = noise(paths[1] / f"wnoise_{level:04d}.bin")
            pair["noise"][str(level)] = metrics(a, b)
            pair["noise"][str(level)]["shape"] = list(a.shape)
            if label == "raw_wrap" and level == 5:
                # Serial cube writeback stores the last visit to each periodic
                # cell; slab output instead keeps all extended FFT cells.
                projected = periodic_last_write(b, 32)
                pair["periodic_last_write_projection"] = metrics(a, projected)
            for field in FIELDS:
                _, a = read_field(paths[0] / f"ic/level_{level:03d}" / field)
                _, b = read_field(paths[1] / f"ic/level_{level:03d}" / field)
                pair["fields"][f"{level}/{field}"] = metrics(a, b)
        report["pairs"][label] = pair
        (root / "geometry_report.json").write_text(json.dumps(report, indent=2)+"\n")
        print(label, json.dumps(pair["noise"]), pair.get("periodic_last_write_projection"), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("--followup", action="store_true")
    args = ap.parse_args()
    if socket.gethostname().split(".")[0] != "LagEunha":
        raise RuntimeError("compute probe requires LagEunha")
    if args.followup:
        followup(args.root, args.baseline)
        return
    args.root.mkdir(exist_ok=False)
    exe = args.baseline / "MUSIC"
    conf = (args.baseline / "z6r1/music.conf").read_text()
    cases = {"auto_no_boost": "rng_boost_threads = no\n",
             "legacy_ref_no_boost": "slab_rng_refinement = no\nrng_boost_threads = no\n",
             "legacy_ref_boost": "slab_rng_refinement = no\n"}
    report = {"purpose": "rank-path isolation only", "production_release": False,
              "binary_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(), "cases": {}}
    for label, switches in cases.items():
        case = args.root / label
        case.mkdir()
        config = conf.replace("[setup]\n", "[setup]\n"+switches)
        config = config.replace(str(args.baseline / "z6r1/ic"), str(case / "ic"))
        (case / "music.conf").write_text(config)
        with (case / "music.log").open("w") as log:
            subprocess.run(["timeout", "180", "mpirun", "-np", "4", str(exe), "music.conf"],
                           cwd=case, stdout=log, stderr=subprocess.STDOUT, check=True)
        result = {"config_sha256": hashlib.sha256(config.encode()).hexdigest(), "noise": {}, "fields": {}}
        for level in (4, 5, 6):
            result["noise"][str(level)] = metrics(
                noise(args.baseline / f"z6r1/wnoise_{level:04d}.bin"),
                noise(case / f"wnoise_{level:04d}.bin"))
            for field in FIELDS:
                _, a = read_field(args.baseline / f"z6r1/ic/level_{level:03d}" / field)
                _, b = read_field(case / f"ic/level_{level:03d}" / field)
                result["fields"][f"{level}/{field}"] = metrics(a, b)
        report["cases"][label] = result
        (args.root / "rank_path_report.json").write_text(json.dumps(report, indent=2)+"\n")
        print(label, json.dumps(result["noise"]), flush=True)


if __name__ == "__main__":
    main()
