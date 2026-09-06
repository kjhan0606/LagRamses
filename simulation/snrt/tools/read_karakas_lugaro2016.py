"""Read pinned KL16 mirror rows without selecting a population or converting ejecta.

This source-specific reader retains commented nodes, initial-composition model
coordinates and source mass discrepancies. It does not run COLIBRE processing.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re

from adapt_g2_candidate_sources import SourceAdapterError

DEFAULT_MATRIX = Path(__file__).resolve().parents[1] / "config/g2_source_selection_matrix_v1.json"
_NUMBER = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
_HEADER = re.compile(rf"Initial mass\s*=\s*({_NUMBER}),\s*Z\s*=\s*({_NUMBER}),\s*Y\s*=\s*({_NUMBER}),\s*M_mix\s*=\s*({_NUMBER})")


def _blocks(text: str, *, initial: bool) -> list[dict]:
    blocks = []
    current = None
    for line_number, raw in enumerate(text.splitlines(), 1):
        match = _HEADER.search(raw)
        if match:
            values = list(map(float, match.groups()))
            if not all(math.isfinite(x) for x in values):
                raise SourceAdapterError("nonfinite KL16 header")
            overshoot = re.search(rf"N_ov\s*=\s*({_NUMBER})", raw)
            current = {
                "source_line": line_number, "source_header": raw,
                "coordinate": dict(zip(("initial_mass_msun", "metallicity_mass_fraction",
                                        "initial_helium_label", "mixing_mass_msun"), values)),
                "overshoot_label": float(overshoot[1]) if overshoot else None,
                "elements_by_atomic_number": {}, "row_comment_flags": [],
            }
            blocks.append(current)
            continue
        if current is None:
            continue
        if "Final mass =" in raw:
            values = list(map(float, re.findall(_NUMBER, raw)))
            if len(values) != 2 or not all(math.isfinite(x) and x >= 0 for x in values):
                raise SourceAdapterError("invalid KL16 returned/final mass")
            current["final_mass_msun"], current["mass_expelled_msun"] = values
            continue
        fields = raw.lstrip("# ").split()
        if len(fields) < 2 or not fields[1].isdigit():
            continue
        if len(fields) != 7:
            raise SourceAdapterError(f"invalid KL16 element row at {line_number}")
        atomic_number = int(fields[1])
        if atomic_number in current["elements_by_atomic_number"]:
            raise SourceAdapterError("duplicate KL16 atomic number")
        values = list(map(float, fields[2:]))
        if not all(math.isfinite(x) for x in values) or any(x < 0 for x in values[-(1 if initial else 2):]):
            raise SourceAdapterError("invalid KL16 source abundance")
        current["elements_by_atomic_number"][atomic_number] = {
            "source_symbol": fields[0], "source_line": line_number,
            "source_numeric_columns": values,
            "mass_fraction": values[-1] if initial else values[-2],
            "gross_mass_msun": None if initial else values[-1],
        }
        current["row_comment_flags"].append(raw.lstrip().startswith("#"))
    for block in blocks:
        flags = block.pop("row_comment_flags")
        if not flags or (any(flags) and not all(flags)):
            raise SourceAdapterError("empty or partly commented KL16 block")
        block["commented_out"] = all(flags)
        if not initial and "mass_expelled_msun" not in block:
            raise SourceAdapterError("KL16 yield block lacks returned mass")
    if not blocks:
        raise SourceAdapterError("no KL16 source blocks")
    return blocks


def read_karakas_lugaro2016(matrix_path: Path = DEFAULT_MATRIX) -> dict:
    matrix = json.loads(Path(matrix_path).read_text())
    candidate = next(c for c in matrix["candidates"] if c["candidate_id"] == "karakas_lugaro2016_agb")
    if candidate["approval_id"] is not None:
        raise SourceAdapterError("KL16 review reader does not implement physical approval")
    base = Path(candidate["source_asset_path"])
    texts = {}
    for suffix in ("007", "014", "030"):
        for prefix in ("yield", "initial", "data"):
            name = f"{prefix}_z{suffix}.txt"
            data = (base / name).read_bytes()
            if hashlib.sha256(data).hexdigest() != candidate["source_file_sha256"].get(name):
                raise SourceAdapterError(f"KL16 fingerprint mismatch: {name}")
            texts[name] = data.decode()
    records, excluded, initial_records = [], [], []
    for suffix in ("007", "014", "030"):
        initial = _blocks(texts[f"initial_z{suffix}.txt"], initial=True)
        for block in initial:
            block["source_file"] = f"initial_z{suffix}.txt"
        initial_records.extend(initial)
        axes = [[float(x) for x in line.split("#", 1)[0].split()]
                for line in texts[f"data_z{suffix}.txt"].splitlines()
                if line.split("#", 1)[0].strip()]
        if len(axes) != 2 or len(axes[0]) != len(axes[1]) or len(set(axes[0])) != len(axes[0]):
            raise SourceAdapterError("KL16 auxiliary mass axes disagree")
        masses = dict(zip(*axes))
        active_masses = []
        for row in _blocks(texts[f"yield_z{suffix}.txt"], initial=False):
            row["source_file"] = f"yield_z{suffix}.txt"
            coord = row["coordinate"]
            if coord["metallicity_mass_fraction"] != float("0." + suffix):
                raise SourceAdapterError("KL16 filename/header metallicity mismatch")
            mass = coord["initial_mass_msun"]
            if row["commented_out"]:
                excluded.append(row)
                continue
            active_masses.append(mass)
            if mass not in masses:
                raise SourceAdapterError("KL16 active mass missing from auxiliary array")
            row["auxiliary_final_mass_msun"] = masses[mass]
            row["auxiliary_minus_header_final_mass_msun"] = masses[mass] - row["final_mass_msun"]
            row["listed_gross_sum_msun"] = math.fsum(e["gross_mass_msun"] for e in row["elements_by_atomic_number"].values())
            row["gross_sum_minus_labelled_expelled_msun"] = row["listed_gross_sum_msun"] - row["mass_expelled_msun"]
            # Full header coordinates are retained: no positional X0 matching or
            # assumption that mixing/overshoot variants are interchangeable.
            matches = [b for b in initial if b["coordinate"] == coord
                       and b["overshoot_label"] == row["overshoot_label"] and not b["commented_out"]]
            row["initial_composition_matching_lines"] = [b["source_line"] for b in matches]
            row["net_yield_diagnostic_msun"] = None
            if len(matches) == 1:
                composition = matches[0]["elements_by_atomic_number"]
                if not set(row["elements_by_atomic_number"]) <= set(composition):
                    raise SourceAdapterError("KL16 initial composition lacks yielded elements")
                row["net_yield_diagnostic_msun"] = {
                    z: e["gross_mass_msun"] - composition[z]["mass_fraction"] * row["mass_expelled_msun"]
                    for z, e in row["elements_by_atomic_number"].items()
                }
            records.append(row)
        if len(active_masses) != len(set(active_masses)) or set(active_masses) != set(masses):
            raise SourceAdapterError("KL16 active yield and auxiliary mass coordinates disagree")
    return {
        "candidate_id": candidate["candidate_id"], "status": "source_rows_review_only",
        "production_ready": False, "canonical_rows_emitted": 0,
        "runtime_activation_allowed": False, "renormalization_applied": False,
        "records": records, "excluded_commented_records": excluded,
        "initial_composition_records": initial_records,
        "net_diagnostic_definition": "gross - source X0 * labelled mass expelled; only unique exact full-header matches, no physical approval",
        "lifetime_yr": None, "injected_energy_erg": None,
    }
