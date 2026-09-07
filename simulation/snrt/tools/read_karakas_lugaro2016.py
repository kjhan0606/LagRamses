"""Read pinned KL16 rows and normalize selected gross ejecta, preserving raw data.

This source-specific reader retains commented nodes, initial-composition model
coordinates and raw source mass discrepancies. The operator-selected detailed
yield-table mass takes precedence over the conflicting auxiliary array; this
does not approve a complete yield/population model or run COLIBRE processing.
The operator-selected gross normalization uses ALL listed elements, never just
the tracked subset. It does not infer timing, energy or normalized net yields.
"""
from __future__ import annotations

import hashlib
import csv
import io
import json
import math
from pathlib import Path
import re

from adapt_g2_candidate_sources import SourceAdapterError

DEFAULT_MATRIX = Path(__file__).resolve().parents[1] / "config/g2_source_selection_matrix_v1.json"
_NUMBER = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
_HEADER = re.compile(rf"Initial mass\s*=\s*({_NUMBER}),\s*Z\s*=\s*({_NUMBER}),\s*Y\s*=\s*({_NUMBER}),\s*M_mix\s*=\s*({_NUMBER})")
GROSS_NORMALIZATION_POLICY = "all_listed_elements_to_selected_expelled_mass_v1"
TRACKED_ATOMIC_NUMBERS = (1, 2, 6, 7, 8, 10, 12, 14, 16, 20, 26)
LIFETIME_PATH = DEFAULT_MATRIX.parents[1] / 'data/kl16_stellar_lifetimes.csv'
LIFETIME_SHA256 = '62e7fd40ce1d63d9024c1755f8226bb501ce92db95831c945ef18eaa2f80ee13'
FISHLOCK_SHA256 = 'b7379ba6eda1018bfeab527ff9bd17fa57fc4da45b2aa5ff845014c31f7eee00'


def read_fishlock2014(*, lifetime_model: str) -> list[dict]:
    """Explicit low-Z comparison: Fishlock yields plus a Padova lifetime fit.

    Fishlock 2014 Table 1 does NOT give total stellar lifetimes. Raiteri96
    coefficients follow Valiante et al. 2009, equations 3--6 (0905.1691).
    This is a declared cross-model timing approximation, not Monash ages.
    Fishlock section 3 identifies the 7-Msun model as ONe; exclude it.
    Gross mass(i)_lost already includes the residual-envelope assumption
    described in that paper. No KL16-style normalization is applied here.
    """
    if lifetime_model != 'raiteri96_padova':
        raise SourceAdapterError('Fishlock requires explicit raiteri96_padova lifetime approximation')
    matrix = json.loads(DEFAULT_MATRIX.read_text())
    candidate = next(c for c in matrix['candidates'] if c['candidate_id'] == 'karakas_lugaro2016_agb')
    path = Path(candidate['source_asset_path']) / 'yield_z001.txt'
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != FISHLOCK_SHA256:
        raise SourceAdapterError('Fishlock source fingerprint mismatch')
    blocks = []; current = None
    for line in data.decode().splitlines():
        match = re.fullmatch(r'\s*#\s*([\d.]+) Msun, Z = ([\d.]+)\s*', line)
        if match:
            current = dict(mass=float(match[1]), z=float(match[2]), elements={})
            blocks.append(current)
        elif line.strip() and not line.lstrip().startswith('#'):
            fields = line.split()
            if current is None or len(fields) != 8:
                raise SourceAdapterError('malformed Fishlock row')
            atomic_number = int(fields[1]); values = list(map(float, fields[2:]))
            if atomic_number in current['elements'] or not all(map(math.isfinite, values)) or values[1] < 0:
                raise SourceAdapterError('invalid Fishlock element row')
            current['elements'][atomic_number] = values[1]
    expected = [1.,1.25,1.5,2.,2.25,2.5,2.75,3.,3.25,3.5,4.,4.5,5.,5.5,6.,7.]
    if [b['mass'] for b in blocks] != expected or any(b['z'] != .001 for b in blocks):
        raise SourceAdapterError('unexpected Fishlock source coordinates')
    # Independent mass-column check: Fishlock Table 1 core masses are printed
    # to 0.001 Msun. Check against that rounding interval, not against the
    # final total stellar mass (which still includes an unejected envelope).
    source_core_masses = (.667,.649,.646,.661,.673,.709,.746,.792,.843,.857,
                          .883,.908,.938,.972,1.015)
    records = []
    for b, source_core in zip(blocks[:-1], source_core_masses, strict=True):
        m,z,e = b['mass'],b['z'],b['elements']
        if not set(TRACKED_ATOMIC_NUMBERS) <= e.keys() or len(e) < 70:
            raise SourceAdapterError('incomplete Fishlock element payload')
        returned = math.fsum(e.values()); remnant = m-returned
        if not 0 < returned < m or not 0 < remnant < 1.4:
            raise SourceAdapterError('Fishlock CO-AGB mass budget invalid')
        if abs(remnant-source_core) > .0005+1e-12:
            raise SourceAdapterError('Fishlock gross sum inconsistent with published core mass')
        x,y = math.log10(m),math.log10(z)
        a0 = 10.13+.07547*y-.008084*y*y
        a1 = -4.424-.7939*y-.1187*y*y
        a2 = 1.262+.3385*y+.05417*y*y
        lifetime = 10**(a0+a1*x+a2*x*x)
        records.append(dict(coordinate=dict(initial_mass_msun=m,metallicity_mass_fraction=z),
                            overshoot_label='Fishlock2014_original',
                            evolution=dict(stellar_lifetime_yr=lifetime,core_kind='CO',
                                           lifetime_source='Raiteri96_Padova_fit_not_Fishlock_evolution'),
                            selected_ejecta=dict(returned_mass_msun=returned,remnant_mass_msun=remnant,
                                                tracked_ejected_mass_msun=[e[k] for k in TRACKED_ATOMIC_NUMBERS])))
    return records


def attach_kl16_lifetimes(records: list[dict]) -> None:
    """Exact M,Z,Y,overshoot match to K14/KL16 Table 1; no lifetime fit.

    M_mix is a post-processing nucleosynthesis choice, not a second evolution
    calculation. KL16 section 3.1 explicitly defines absent N_ov as no
    overshoot. Do not transfer an overshoot lifetime to a no-overshoot model.
    """
    data = LIFETIME_PATH.read_bytes()
    if hashlib.sha256(data).hexdigest() != LIFETIME_SHA256:
        raise SourceAdapterError('KL16 lifetime data fingerprint mismatch')
    grid = {}
    for line in csv.DictReader(io.StringIO(data.decode())):
        key = tuple(float(line[k]) for k in ('initial_mass_msun', 'metallicity', 'initial_helium', 'overshoot'))
        age = float(line['stellar_duration_myr'])*1e6
        if key in grid or not all(math.isfinite(v) for v in (*key, age)) or age <= 0:
            raise SourceAdapterError('invalid or duplicate KL16 lifetime coordinate')
        grid[key] = (age, line['source'], line['core_kind'])
    matches = []
    for row in records:
        c = row['coordinate']
        key = (c['initial_mass_msun'], c['metallicity_mass_fraction'], c['initial_helium_label'],
               0. if row['overshoot_label'] is None else row['overshoot_label'])
        if key not in grid:
            raise SourceAdapterError(f'no exact KL16 evolution match: {key}')
        age, source, core = grid[key]
        matches.append(dict(stellar_lifetime_yr=age, lifetime_source=source+':Table1',
                            lifetime_data_sha256=LIFETIME_SHA256, core_kind=core,
                            lifetime_convention='source_total_stellar_duration_to_AGB_endpoint'))
    for row, match in zip(records, matches):
        row['evolution'] = match


def normalize_selected_ejecta(row: dict) -> dict:
    """Return a derived payload; never overwrite or renormalize raw fields.

    Apply M_i' = M_expelled * M_i / sum_all(M_i). Initial-composition data
    and old raw net-yield diagnostics must not be mistaken for this payload.
    """
    if row.get("commented_out", True):
        raise SourceAdapterError("cannot normalize an excluded KL16 node")
    gross = {z: e["gross_mass_msun"] for z, e in row["elements_by_atomic_number"].items()}
    if len(gross) != 78 or not set(TRACKED_ATOMIC_NUMBERS) <= gross.keys():
        raise SourceAdapterError("normalization requires the full 78-element KL16 payload")
    initial = row["coordinate"]["initial_mass_msun"]
    expelled = row["selected_mass_expelled_msun"]
    remnant = row["selected_final_mass_msun"]
    if (not all(math.isfinite(x) and x >= 0 for x in (*gross.values(), initial, expelled, remnant))
            or initial <= 0 or expelled <= 0):
        raise SourceAdapterError("invalid KL16 normalization mass")
    if not math.isclose(expelled + remnant, initial, rel_tol=1e-12, abs_tol=0):
        raise SourceAdapterError("selected KL16 mass budget does not close")
    try:
        total = math.fsum(gross.values())
    except OverflowError as error:
        raise SourceAdapterError("KL16 gross sum overflow") from error
    if not math.isfinite(total) or total <= 0:
        raise SourceAdapterError("KL16 gross sum must be finite and positive")
    scale = expelled / total
    if not math.isfinite(scale) or scale <= 0:
        raise SourceAdapterError("invalid KL16 normalization factor")
    selected = {z: mass * scale for z, mass in gross.items()}
    # No element absorbs a residual: preserve every raw ratio, to roundoff.
    selected_sum = math.fsum(selected.values())
    if not math.isclose(selected_sum, expelled, rel_tol=1e-14, abs_tol=0):
        raise SourceAdapterError("normalized KL16 mass budget does not close")
    return {
        "policy": GROSS_NORMALIZATION_POLICY,
        "raw_gross_sum_msun": total,
        "normalization_factor": scale,
        "gross_sum_correction_msun": selected_sum - total,
        "returned_mass_msun": expelled, "remnant_mass_msun": remnant,
        "mass_source": row["selected_mass_source"],
        "gross_mass_msun_by_atomic_number": selected,
        "mass_fraction_by_atomic_number": {z: mass / expelled for z, mass in selected.items()},
        "tracked_ejected_mass_msun": [selected[z] for z in TRACKED_ATOMIC_NUMBERS],
        "untracked_ejecta_msun": math.fsum(v for z, v in selected.items() if z not in TRACKED_ATOMIC_NUMBERS),
        "total_metal_ejecta_msun": math.fsum(v for z, v in selected.items() if z not in (1, 2)),
        "net_yield_msun": None,
        "net_yield_status": "unavailable_initial_composition_normalization_and_model_matching_not_selected",
    }


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


def read_karakas_lugaro2016(matrix_path: Path = DEFAULT_MATRIX, *, include_lifetimes: bool = False) -> dict:
    matrix = json.loads(Path(matrix_path).read_text())
    candidate = next(c for c in matrix["candidates"] if c["candidate_id"] == "karakas_lugaro2016_agb")
    if candidate["approval_id"] is not None:
        raise SourceAdapterError("KL16 review reader does not implement physical approval")
    if candidate.get("selected_gross_normalization_policy") != GROSS_NORMALIZATION_POLICY:
        raise SourceAdapterError("KL16 normalization lacks the matching operator selection")
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
            # Preserve raw fingerprints and both source values. Do not repair
            # the archive. Normalized ejecta live in a separate payload. A printed
            # article value takes precedence if available; Table 7 only prints
            # the 3.5 Msun example, not the disputed 4 Msun node.
            row["selected_final_mass_msun"] = row["final_mass_msun"]
            row["selected_mass_expelled_msun"] = row["mass_expelled_msun"]
            row["selected_mass_source"] = "yield_table_header"
            if mass == 3.5 and coord["metallicity_mass_fraction"] == 0.03:
                row["selected_final_mass_msun"] = 0.727
                row["selected_mass_expelled_msun"] = 2.773
                row["selected_mass_source"] = "KL16_article_table_7"
            row["auxiliary_mass_disagreement_resolved"] = False
            if row["auxiliary_minus_header_final_mass_msun"] != 0.0:
                choice = candidate["source_row_review"]["mass_disagreement"]
                if (choice.get("authoritative_value_selected") is not True or
                        mass != choice["initial_mass_msun"] or
                        coord["metallicity_mass_fraction"] != choice["metallicity_mass_fraction"] or
                        row["final_mass_msun"] != choice["selected_final_mass_msun"] or
                        row["mass_expelled_msun"] != choice["selected_mass_expelled_msun"] or
                        masses[mass] != choice["auxiliary_final_mass_msun"] or
                        choice["selected_source"] != "yield_table_header"):
                    raise SourceAdapterError("KL16 mass disagreement lacks the matching operator selection")
                row["auxiliary_mass_disagreement_resolved"] = True
            if not math.isclose(mass, row["selected_final_mass_msun"] +
                                row["selected_mass_expelled_msun"], rel_tol=0, abs_tol=5e-4):
                raise SourceAdapterError("KL16 selected remnant and ejecta do not close initial mass")
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
            row["selected_ejecta"] = normalize_selected_ejecta(row)
            records.append(row)
        if len(active_masses) != len(set(active_masses)) or set(active_masses) != set(masses):
            raise SourceAdapterError("KL16 active yield and auxiliary mass coordinates disagree")
    if include_lifetimes:
        attach_kl16_lifetimes(records)
    return {
        "candidate_id": candidate["candidate_id"], "status": "source_rows_review_only",
        "production_ready": False, "canonical_rows_emitted": 0,
        "runtime_activation_allowed": False, "renormalization_applied": True,
        "normalization_policy": GROSS_NORMALIZATION_POLICY,
        "source_values_modified": False,
        "records": records, "excluded_commented_records": excluded,
        "initial_composition_records": initial_records,
        "net_diagnostic_definition": "gross - source X0 * labelled mass expelled; only unique exact full-header matches, no physical approval",
        "lifetime_yr": None, "injected_energy_erg": None,
    }
