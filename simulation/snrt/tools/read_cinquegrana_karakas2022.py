"""Read the pinned 77-species CK22 mirror, without a decay or ejecta model.

The column convention follows Karakas (2010), section 4, as suggested by
the source headers; the CK22 README calls <X(i)> a final fraction instead.
Preserve that ambiguity and the two inconsistent mass denominators. Never
use this review reader as a production-source approval or normalize its rows.
"""
from __future__ import annotations

import hashlib
import csv
import io
import json
import math
from pathlib import Path
import re
import subprocess

from adapt_g2_candidate_sources import SourceAdapterError

DEFAULT_MATRIX = Path(__file__).resolve().parents[1] / "config/g2_source_selection_matrix_v1.json"
_ATOMIC_NUMBER = {symbol: z for symbol, z in zip(
    ("he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al", "si", "p", "s", "fe", "co", "ni"),
    (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 26, 27, 28))}
_TRACKED = {1: "H", 2: "He", 6: "C", 7: "N", 8: "O", 10: "Ne", 12: "Mg", 14: "Si", 16: "S", 20: "Ca", 26: "Fe"}


def _evolution_rows(reference: dict) -> list[dict]:
    """Published model durations, not a terminal-lumped release prescription."""
    if reference["approval_id"] is not None:
        raise SourceAdapterError("CK22 reader cannot activate a release-time model")
    data = (DEFAULT_MATRIX.parents[1] / reference["path"]).read_bytes()
    if hashlib.sha256(data).hexdigest() != reference["sha256"]:
        raise SourceAdapterError("CK22 evolution table fingerprint mismatch")
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if len(rows) != reference["source_rows"]:
        raise SourceAdapterError("CK22 evolution table row count mismatch")
    for row in rows:
        for key in ("initial_mass_msun", "metallicity_mass_fraction", "initial_helium_label",
                    "stellar_duration_myr", "rgb_duration_myr", "agb_duration_myr"):
            if not math.isfinite(float(row[key])) or float(row[key]) <= 0:
                raise SourceAdapterError("invalid CK22 evolution coordinate or duration")
    return rows


def _pinned_files(base: Path, tree_sha1: str) -> dict[str, bytes]:
    # The immutable tree binds paths AND contents for the sparse raw-yield
    # selection. Its processed HDF5/surface files are not read or adopted.
    listing = subprocess.check_output(["git", "-C", str(base), "ls-tree", "-rz", "--full-tree", tree_sha1])
    files = {}
    for entry in listing.split(b"\0"):
        if not entry:
            continue
        meta, name_bytes = entry.split(b"\t", 1)
        mode, kind, blob = meta.decode().split()
        name = name_bytes.decode()
        if name != "readme_yields.dat" and not re.fullmatch(r"z\d{2}/yields_m[^/]+\.dat", name):
            continue
        path = base / name
        if mode != "100644" or kind != "blob" or path.is_symlink():
            raise SourceAdapterError("CK22 source must contain ordinary data files")
        data = path.read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if actual != blob:
            raise SourceAdapterError(f"CK22 fingerprint mismatch: {name}")
        files[name] = data
    actual_names = {p.relative_to(base).as_posix() for p in base.glob("z*/yields_*.dat")}
    if actual_names | {"readme_yields.dat"} != set(files):
        raise SourceAdapterError("CK22 source file inventory differs from pinned tree")
    return files


def _read_table(text: str) -> dict:
    lines = text.splitlines()
    if " ".join(lines[:2]).split() != ["species", "A", "yield", "mass(i)_lost", "mass(i)_0", "<X(i)>", "X0(i)", "log10(<X(i)>/X0(i))"]:
        raise SourceAdapterError("unrecognized CK22 column header")
    rows = {}
    for line_number, raw in enumerate(lines[2:], 3):
        fields = raw.split()
        if not fields:
            continue
        if len(fields) != 8 or not fields[1].isdigit() or fields[0] in rows:
            raise SourceAdapterError("invalid or duplicate CK22 species row")
        label, mass_number = fields[0], int(fields[1])
        values = list(map(float, fields[2:]))
        if not all(math.isfinite(x) for x in values):
            raise SourceAdapterError("nonfinite CK22 numeric column")
        if label in ("g", "n", "p", "d"):
            if mass_number != (2 if label == "d" else 1):
                raise SourceAdapterError("CK22 special-species mass label mismatch")
            atomic_number = 1 if label in ("p", "d") else None
        elif label in ("al-6", "al*6"):
            if mass_number != 26:
                raise SourceAdapterError("CK22 Al26 state mass label mismatch")
            atomic_number = 13
        else:
            match = re.fullmatch(r"([a-z]+)(\d+)", label)
            if not match or match[1] not in _ATOMIC_NUMBER or int(match[2]) != mass_number:
                raise SourceAdapterError("unrecognized CK22 isotope")
            atomic_number = _ATOMIC_NUMBER[match[1]]
        rows[label] = {
            "source_line": line_number, "source_numeric_tokens": fields[2:],
            "mass_number_label": mass_number, "atomic_number": atomic_number,
            "net_mass_msun": values[0], "gross_mass_msun": values[1],
            "initial_composition_mass_in_wind_msun": values[2],
            "wind_mass_fraction_column": values[3], "initial_mass_fraction": values[4],
            "log_production_factor_column": values[5],
            "net_identity_residual_msun": values[0] - (values[1] - values[2]),
        }
    if len(rows) != 77 or not {"g", "n", "p", "d", "al-6", "al*6"} <= set(rows):
        raise SourceAdapterError("incomplete CK22 77-species table")
    p = rows["p"]
    if p["wind_mass_fraction_column"] <= 0 or p["initial_mass_fraction"] <= 0:
        raise SourceAdapterError("invalid CK22 H denominator")
    return {
        "species_by_source_label": rows,
        "missing_tracked_elements": [name for z, name in _TRACKED.items()
                                     if z not in {r["atomic_number"] for r in rows.values()}],
        "H_gross_over_wind_fraction_msun": p["gross_mass_msun"] / p["wind_mass_fraction_column"],
        "H_initial_wind_mass_over_X0_msun": p["initial_composition_mass_in_wind_msun"] / p["initial_mass_fraction"],
        # Diagnostic only: g is a network bookkeeping species, not measured Ca,
        # Sr or Ba. No assignment of its mass to a physical element is made.
        "listed_gross_sum_including_g_msun": math.fsum(r["gross_mass_msun"] for r in rows.values()),
        "returned_mass_msun": None, "remnant_mass_msun": None,
    }


def read_cinquegrana_karakas2022(matrix_path: Path = DEFAULT_MATRIX) -> dict:
    matrix = json.loads(Path(matrix_path).read_text())
    candidate = next(c for c in matrix["candidates"] if c["candidate_id"] == "cinquegrana_karakas2022_agb")
    if candidate["approval_id"] is not None:
        raise SourceAdapterError("CK22 review reader does not implement physical approval")
    evolution = _evolution_rows(candidate["evolution_reference"])
    files = _pinned_files(Path(candidate["source_asset_path"]), candidate["source_mirror_tree_sha1"])
    if hashlib.sha256(files["readme_yields.dat"]).hexdigest() != candidate["source_file_sha256"]["readme_yields.dat"]:
        raise SourceAdapterError("CK22 README fingerprint mismatch")
    records = []
    for name, data in sorted(files.items()):
        if name == "readme_yields.dat":
            continue
        match = re.fullmatch(r"z(\d{2})/yields_m(\d+(?:\.\d+)?)z\1.dat", name)
        if not match:
            raise SourceAdapterError("unexpected CK22 data file")
        row = _read_table(data.decode())
        row.update(source_file=name, initial_mass_msun=float(match[2]),
                   metallicity_mass_fraction=int(match[1]) / 100)
        # Same evolution family according to CK22 section 3. The yield filename
        # supplies M and Z only: do not claim a full population/Y approval.
        matches = [e for e in evolution if float(e["initial_mass_msun"]) == row["initial_mass_msun"]
                   and float(e["metallicity_mass_fraction"]) == row["metallicity_mass_fraction"]]
        row["evolution_mass_Z_matches"] = matches
        row["evolution_duration_candidate_yr"] = (
            float(matches[0]["stellar_duration_myr"]) * 1e6 if len(matches) == 1 else None)
        row["release_time_yr"] = None
        coord = row["initial_mass_msun"], row["metallicity_mass_fraction"]
        row["evolution_endpoint_caveat"] = None
        if coord in {(8.0, 0.09), (7.0, 0.1), (8.0, 0.1)}:
            row["evolution_endpoint_caveat"] = "KCJ22_section_4.2.1_explicit_early_termination"
        elif coord == (7.5, 0.1):
            row["evolution_endpoint_caveat"] = "KCJ22_carbon_ignition_zero_TP_endpoint_requires_clarification"
        records.append(row)
    excluded = {(7.0, 0.09), (5.5, 0.1)}
    expected = {(m / 2, z / 100) for m in range(2, 17) for z in range(4, 11)} - excluded
    if len(records) != len(expected) or {(r["initial_mass_msun"], r["metallicity_mass_fraction"]) for r in records} != expected:
        raise SourceAdapterError("CK22 model grid differs from published exclusions")
    return {
        "candidate_id": candidate["candidate_id"], "status": "source_rows_review_only",
        "records": records, "excluded_mass_metallicity_nodes": sorted(excluded),
        "evolution_source_rows": evolution,
        "production_ready": False, "canonical_rows_emitted": 0,
        "runtime_activation_allowed": False, "renormalization_applied": False,
        "negative_residue_clipping_applied": False, "decay_applied": False,
        "lifetime_yr": None, "injected_energy_erg": None,
    }
