from __future__ import annotations

from pathlib import Path
import re


SOURCE = Path(__file__).resolve().parents[1] / "output_fdm.f90"


def _subroutine(source: str, name: str) -> str:
    token = re.escape(name)
    match = re.search(rf"(?im)^\s*subroutine\s+{token}\b", source)
    if match is None:
        raise AssertionError(f"subroutine {name} not found")
    end = re.search(
        rf"(?im)^\s*end\s+subroutine\s+{token}\b", source[match.start() :]
    )
    if end is None:
        raise AssertionError(f"end of subroutine {name} not found")
    return source[match.start() : match.start() + end.end()]


def test_outer_wave_provenance_does_not_start_output_time_ghost_exchange() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    writer = _subroutine(source, "output_fdm_outer_wave_provenance").lower()

    assert not re.search(r"\bcall\s+make_virtual_fine_dp2?\s*\(", writer)
    assert "call fdm_neighbor_cell" in writer
    assert "stencil_loc=stencil_loc+1.0d0" in writer
    assert "complete_current_stencil_fraction" in writer
