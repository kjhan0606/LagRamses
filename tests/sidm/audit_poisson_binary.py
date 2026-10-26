#!/usr/bin/env python3
"""Check linked FFTW Poisson code, not just the startup feature banner.

This is a build-evidence gate for an unstripped RAMSES executable. It neither
executes the binary nor proves solver dispatch, numerical accuracy, or a
uniform periodic mesh. Require actual FFT solve evidence in the runtime log
as a separate gate. Static FFTW linkage is supported.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


def parse_symbols(text: str) -> dict[str, str]:
    symbols = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 3:
            _, kind, name = fields
        elif len(fields) == 2:
            kind, name = fields
        else:
            continue
        if len(kind) == 1:
            symbols[name.split("@", 1)[0]] = kind
    return symbols


def classify_symbols(text: str) -> dict:
    symbols = parse_symbols(text)
    routines = sorted(name for name, kind in symbols.items()
                      if kind in "TtWw" and name.rstrip("_").lower() ==
                      "fftw_poisson_solve_uniform")
    plans = sorted(name for name in symbols
                   if re.fullmatch(r"fftw_(?:mpi_)?plan_\w+", name))
    executes = sorted(name for name in symbols
                      if re.fullmatch(r"fftw_(?:mpi_)?execute(?:_\w+)?", name))
    return {
        "fftw_poisson_routines": routines,
        "fftw_plan_symbols": plans,
        "fftw_execute_symbols": executes,
        "regular_code_symbols_present": any(kind in "Tt" for kind in symbols.values()),
        "fftw_binary_evidence_pass": bool(routines and plans and executes),
        "solver_dispatch_verified": False,
        "physics_convergence_verified": False,
        "scope": "linked code evidence only; actual FFT solve and MG runtime gates remain mandatory",
    }


def audit(binary: Path) -> dict:
    binary = binary.resolve(strict=True)
    if not binary.is_file():
        raise ValueError(f"expected an executable file: {binary}")
    # Use the regular table as well as dynamic imports: ldd alone misses
    # statically linked FFTW, and a library alone does not prove RAMSES has
    # its FFTW Poisson entry point. No binary execution is performed.
    regular = subprocess.run(["nm", "-g", "--", str(binary)], text=True,
                             capture_output=True, check=False, timeout=30)
    dynamic = subprocess.run(["nm", "-D", "--", str(binary)], text=True,
                             capture_output=True, check=False, timeout=30)
    if regular.returncode != 0:
        raise ValueError(f"cannot inspect regular symbols: {regular.stderr.strip()}")
    result = classify_symbols(regular.stdout + "\n" + dynamic.stdout)
    digest = hashlib.sha256()
    with binary.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    result.update({
        "binary": str(binary), "binary_sha256": digest.hexdigest(),
        "nm_regular_returncode": regular.returncode,
        "nm_dynamic_returncode": dynamic.returncode,
        "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    })
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--require-fftw", action="store_true")
    args = parser.parse_args()
    if args.json is not None and args.json.exists():
        parser.error("refusing to overwrite a previous binary audit")
    result = audit(args.binary)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    print(text, end="")
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text)
    if args.require_fftw and not result["fftw_binary_evidence_pass"]:
        raise SystemExit("FFTW Poisson binary gate failed (missing or insufficient linked-code evidence)")


if __name__ == "__main__":
    main()
