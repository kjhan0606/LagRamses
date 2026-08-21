#!/usr/bin/env python3
"""Fail-closed census for CUDA cell addressing under block grid-major AMR."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CUDA = ROOT / "patch" / "cuRamses"


def fail(message: str) -> None:
    print(f"CUDA BLOCK INDEX CENSUS: FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def read(name: str) -> str:
    return (CUDA / name).read_text(encoding="utf-8")


def require_count(text: str, pattern: str, expected: int, label: str) -> None:
    count = len(re.findall(pattern, text, flags=re.MULTILINE))
    if count != expected:
        fail(f"{label}: expected {expected}, found {count}")


def main() -> int:
    poisson = read("poisson_cuda_kernels.cu")
    particle = read("particle_cuda_kernels.cu")
    scalar = read("scalar_cuda_kernels.cu")

    require_count(poisson, r"\bamr_cuda_cell_1based\s*\(", 8, "poisson cell sites")
    require_count(particle, r"\bamr_cuda_cell_1based\s*\(", 2, "particle cell sites")
    require_count(particle, r"\bamr_cuda_live_cell_prefix\s*\(", 1, "particle prefix site")
    require_count(scalar, r"\bamr_cuda_cell_1based\s*\(", 2, "scalar cell sites")

    direct_cell = re.compile(
        r"ncoarse\s*\+[^;\n]*(?:icell|ind)[^;\n]*\*\s*ngridmax"
        r"|ncoarse\s*\+[^;\n]*ngridmax[^;\n]*(?:igrid|igr)"
    )
    for name, text in (("poisson", poisson), ("particle", particle), ("scalar", scalar)):
        match = direct_cell.search(text)
        if match:
            fail(f"{name}: legacy cell stride remains: {match.group(0).strip()}")

    require_count(
        poisson,
        r"flag2\s*\[[^\]]+\]\s*/\s*ngridmax",
        2,
        "poisson flag2 packing whitelist",
    )
    if "cudaMemcpy2D" in particle:
        fail("particle upload still uses legacy child-plane cudaMemcpy2D")
    for token in (
        "pm_layout_valid(",
        "amr_cuda_layout_valid(",
        "live_ncell > ncell",
        "(size_t)live_ncell * elemsz",
    ):
        if token not in particle:
            fail(f"particle prefix validation/copy token missing: {token}")

    makefile = (ROOT / "bin" / "Makefile").read_text(encoding="utf-8")
    if not re.search(
        r"VPATH\s*=\s*\$\(PATCH\):\.\./patch/cuda:"
        r"\.\./patch/oct_tree:\.\./patch/cuRamses",
        makefile,
    ):
        fail("unexpected Makefile VPATH; source-winner assumptions need re-audit")
    active_rho = (ROOT / "patch" / "lagRamses" / "rho_fine.kjhan.f90").read_text(
        encoding="utf-8"
    )
    if "call cuda_pm_rho_begin_c" not in active_rho:
        fail("active lagRamses rho source does not dispatch CUDA deposit")
    if "int(amr_block_size,c_int), int(twotondim,c_int)" not in active_rho:
        fail("active lagRamses rho source does not pass B/C")

    for rel in (
        Path("patch/lagRamses/multigrid_fine_commons.f90"),
        Path("patch/cuRamses/multigrid_fine_commons.f90"),
    ):
        source = (ROOT / rel).read_text(encoding="utf-8")
        calls = len(re.findall(
            r"call\s+cuda_mg_(?:gauss_seidel|residual|restrict_execute|interp_execute)_c",
            source,
            flags=re.IGNORECASE,
        ))
        bc = source.count("int(amr_block_size,c_int), int(twotondim,c_int)")
        if calls != 9 or bc != 9:
            fail(f"{rel}: MG ABI calls={calls}, B/C pairs={bc}, expected 9/9")

    print("CUDA BLOCK INDEX CENSUS: PASS (poisson=8 particle=3 scalar=2)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
