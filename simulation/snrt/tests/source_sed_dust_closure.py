#!/usr/bin/env python3
"""Test source-SED identity, AGN grouping, and source-bound dust closure."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.dust import read_dust_opacity_metadata  # noqa: E402
from snrt_core.sed import (  # noqa: E402
    EV_ERG,
    SED_QUADRATURE_RELATIVE_TOLERANCE,
    integrate_photon_sed_groups,
    read_lbol_photon_sed,
)
from snrt_core.snapshot import GridSpec, SourceCatalog, neutral_primordial_input, write_static_rt_input  # noqa: E402
from tools import build_draine_dust_opacity  # noqa: E402


def _write_sed(path: Path, edges: np.ndarray) -> None:
    fraction = 1.0 / (edges[-1] - edges[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("energy_ev", "energy_fraction_per_ev"))
        for energy in edges:
            writer.writerow((f"{energy:.17g}", f"{fraction:.17g}"))


def _write_tilted_sed(path: Path, edges: np.ndarray) -> None:
    """Write a distinct source-bound fixture for the stellar-side contract."""

    energy = np.asarray(edges, dtype=np.float64)
    fraction = 1.0 + 0.2 * (energy - energy[0]) / (energy[-1] - energy[0])
    fraction /= np.trapezoid(fraction, energy)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("energy_ev", "energy_fraction_per_ev"))
        writer.writerows((f"{e:.17g}", f"{f:.17g}") for e, f in zip(energy, fraction, strict=True))


def _write_draine(
    path: Path,
    energies: tuple[float, ...] = (1.0, 10.0, 100.0, 1000.0),
    *,
    albedo: float = 0.0,
    cosine: float = 0.0,
    cosine_squared: float = 0.0,
) -> None:
    # A compact table with the same header semantics as the official Draine
    # format.  The opacity is deliberately energy-dependent so the test
    # exercises source weighting rather than only schema loading.
    rows = tuple(
        (
            1.2398419843320026 / energy,
            albedo,
            cosine,
            energy,
            energy * (1.0 - albedo) * 1.0e26,
            cosine_squared,
        )
        for energy in energies
    )
    lines = [
        "1.0e-26 = M_dust per H",
        "100.0 = M_gas/M_dust",
        "# wavelength albedo unused C_ext/H K_abs unused",
        *(" ".join(f"{value:.17g}" for value in row) for row in rows),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_powerlaw_sed(path: Path) -> np.ndarray:
    """Write f_E proportional to E on a grid that omits the 10 eV edge."""

    energies = np.linspace(1.0, 100.0, 4097, dtype=np.float64)
    fraction = 2.0 * energies / (100.0**2 - 1.0)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("energy_ev", "energy_fraction_per_ev"))
        writer.writerows((f"{energy:.17g}", f"{value:.17g}") for energy, value in zip(energies, fraction, strict=True))
    return energies


def _write_powerlaw_draine(path: Path) -> np.ndarray:
    """Write kappa(E)=k0 E^-1/2 on a log grid that omits the 10 eV edge."""

    energies = np.geomspace(1.0, 100.0, 4096, dtype=np.float64)
    opacity = 4.0e-21 * energies**-0.5
    rows = tuple(
        (1.2398419843320026 / energy, 0.0, 0.0, opacity_value, opacity_value / 1.0e-26, 0.0)
        for energy, opacity_value in zip(energies, opacity, strict=True)
    )
    lines = [
        "1.0e-26 = M_dust per H",
        "100.0 = M_gas/M_dust",
        "# wavelength albedo unused C_ext/H K_abs unused",
        *(" ".join(f"{value:.17g}" for value in row) for row in rows),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return energies


def _independent_draine_recompute(
    sed_path: Path,
    draine_path: Path,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Recompute source-weighted opacity without using the production builder."""

    with sed_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        sed_rows = [(float(row["energy_ev"]), float(row["energy_fraction_per_ev"])) for row in reader]
    sed_data = np.asarray(sed_rows, dtype=np.float64)
    draine_rows: list[tuple[float, ...]] = []
    for raw_line in draine_path.read_text(encoding="utf-8").splitlines():
        fields = raw_line.split()
        if len(fields) != 6:
            continue
        try:
            draine_rows.append(tuple(float(field) for field in fields))
        except ValueError:
            continue
    draine_data = np.asarray(draine_rows, dtype=np.float64)
    draine_energy = 1.2398419843320026 / draine_data[:, 0]
    draine_opacity = draine_data[:, 4] * 1.0e-26
    order = np.argsort(draine_energy)
    draine_energy = draine_energy[order]
    draine_opacity = draine_opacity[order]

    weighted_opacity = np.zeros(edges.size - 1, dtype=np.float64)
    weighted_energy = np.zeros_like(weighted_opacity)
    for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        grid = np.linspace(lower, upper, 100001, dtype=np.float64)
        fraction = np.interp(grid, sed_data[:, 0], sed_data[:, 1])
        photons = fraction / (grid * EV_ERG)
        opacity = np.exp(np.interp(np.log(grid), np.log(draine_energy), np.log(draine_opacity)))
        photon_norm = np.trapezoid(photons, grid)
        absorbed_norm = np.trapezoid(photons * opacity, grid)
        weighted_opacity[index] = absorbed_norm / photon_norm
        weighted_energy[index] = np.trapezoid(grid * photons * opacity, grid) / absorbed_norm
    return weighted_opacity, weighted_energy


def _powerlaw_closed_form(edges: np.ndarray, k0: float, exponent: float) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form photon-weighted averages for q(E)=constant, kappa=k0 E^-p."""

    opacity = np.zeros(edges.size - 1, dtype=np.float64)
    energy = np.zeros_like(opacity)
    for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        integral_kappa = k0 * (upper ** (1.0 - exponent) - lower ** (1.0 - exponent)) / (1.0 - exponent)
        integral_energy_kappa = k0 * (upper ** (2.0 - exponent) - lower ** (2.0 - exponent)) / (2.0 - exponent)
        opacity[index] = integral_kappa / (upper - lower)
        energy[index] = integral_energy_kappa / integral_kappa
    return opacity, energy


def _write_candidate(path: Path) -> None:
    path.write_text(
        "source_id,source_kind,x_code,y_code,z_code,bolometric_luminosity_erg_s\n"
        "1,agn,0.5,0.5,0.5,1.0e40\n",
        encoding="utf-8",
    )


def main() -> int:
    with TemporaryDirectory(prefix="source-sed-dust-test-") as directory:
        work = Path(directory)
        edges = np.asarray((1.0, 10.0, 100.0), dtype=np.float64)
        edges_path = work / "edges.txt"
        edges_path.write_text("1\n10\n100\n", encoding="utf-8")
        sed_path = work / "sed.csv"
        _write_sed(sed_path, edges)
        sed = read_lbol_photon_sed(sed_path, expected_bolometric_fraction=1.0)
        moments = integrate_photon_sed_groups(sed, edges, allow_empty_groups=False)
        assert np.isclose(np.sum(moments.group_energy_fraction_per_norm), 1.0)
        assert np.all(moments.group_photon_rate_per_norm > 0.0)
        assert np.allclose(
            moments.photon_weighted_mean_energy_ev,
            [(10.0 - 1.0) / np.log(10.0), (100.0 - 10.0) / np.log(10.0)],
            rtol=SED_QUADRATURE_RELATIVE_TOLERANCE,
        )
        assert moments.quadrature_diagnostics is not None
        assert moments.quadrature_diagnostics.converged
        assert moments.quadrature_diagnostics.maximum_relative_error <= SED_QUADRATURE_RELATIVE_TOLERANCE

        draine_path = work / "draine.all"
        _write_draine(draine_path)
        dust_metadata = build_draine_dust_opacity.build_source_weighted_opacity_metadata(
            draine_path,
            sed_path,
            edges_path,
            expected_bolometric_fraction=1.0,
        )
        assert dust_metadata["schema"] == "snrt_dust_opacity_v2"
        assert dust_metadata["source_sed_identity"] == sed.identity
        assert dust_metadata["group_edges_path"] == str(edges_path.resolve())
        assert dust_metadata["source_table"]["path"] == str(draine_path.resolve())
        assert dust_metadata["builder"]["path"] == str(build_draine_dust_opacity.TOOL_PATH.resolve())
        dust_path = work / "dust.json"
        dust_path.write_text(json.dumps(dust_metadata) + "\n", encoding="utf-8")
        closure = read_dust_opacity_metadata(
            dust_path,
            expected_group_edges_ev=edges,
            expected_source_sed_identity=sed.identity,
            expected_source_sed_sha256=sed.input_sha256,
            require_source_match=True,
        )
        assert closure.schema == "snrt_dust_opacity_v2"
        assert closure.source_sed_identity == sed.identity
        assert closure.group_edges_sha256 == dust_metadata["group_edges_sha256"]
        assert closure.source_table_sha256 == dust_metadata["source_table"]["sha256"]
        assert closure.builder_sha256 == dust_metadata["builder"]["sha256"]
        assert closure.payload_sha256 == dust_metadata["payload_sha256"]
        assert {entry["role"] for entry in dust_metadata["closure_code_manifest"]} == {
            "dust_builder",
            "source_sed",
            "dust_loader",
            "integrity_helper",
        }
        assert np.all(closure.absorption_cross_section_per_h_cm2 > 0.0)

        # The scattering schema is source-bound independently for a second
        # (stellar-side) synthetic SED; it is not an aggregate STAR+AGN
        # closure and must carry a distinct source identity.
        stellar_sed_path = work / "stellar_sed.csv"
        _write_tilted_sed(stellar_sed_path, edges)
        stellar_sed = read_lbol_photon_sed(stellar_sed_path, expected_bolometric_fraction=1.0)
        scattering_draine_path = work / "scattering-draine.all"
        _write_draine(
            scattering_draine_path,
            albedo=0.5,
            cosine=0.4,
            cosine_squared=0.2,
        )
        stellar_scattering_metadata = build_draine_dust_opacity.build_source_weighted_opacity_metadata(
            scattering_draine_path,
            stellar_sed_path,
            edges_path,
            expected_bolometric_fraction=1.0,
            include_scattering=True,
        )
        assert stellar_scattering_metadata["schema"] == "snrt_dust_opacity_v3"
        assert stellar_scattering_metadata["status"] == "candidate_scattering_isotropic"
        assert stellar_scattering_metadata["source_sed_identity"] == stellar_sed.identity
        assert stellar_scattering_metadata["source_sed_identity"] != sed.identity
        assert max(stellar_scattering_metadata["scattering_cross_section_per_h_cm2"]) > 0.0
        stellar_scattering_path = work / "stellar-scattering-dust.json"
        stellar_scattering_path.write_text(
            json.dumps(stellar_scattering_metadata) + "\n", encoding="utf-8"
        )
        stellar_scattering_closure = read_dust_opacity_metadata(
            stellar_scattering_path,
            expected_group_edges_ev=edges,
            expected_source_sed_identity=stellar_sed.identity,
            expected_source_sed_sha256=stellar_sed.input_sha256,
            require_source_match=True,
        )
        assert stellar_scattering_closure.schema == "snrt_dust_opacity_v3"
        assert stellar_scattering_closure.scattering_phase_function == "phase_isotropic_candidate"

        # Independent numerical and closed-form verification on intentionally
        # offset SED/Draine/group samples. Neither calculation uses the
        # production source-weighted integration routine.
        offset_sed_path = work / "offset-powerlaw-sed.csv"
        offset_draine_path = work / "offset-powerlaw-draine.all"
        offset_sed_energy = _write_powerlaw_sed(offset_sed_path)
        offset_draine_energy = _write_powerlaw_draine(offset_draine_path)
        assert not np.any(np.isclose(offset_sed_energy, 10.0, rtol=0.0, atol=1.0e-12))
        assert not np.any(np.isclose(offset_draine_energy, 10.0, rtol=0.0, atol=1.0e-12))
        offset_metadata = build_draine_dust_opacity.build_source_weighted_opacity_metadata(
            offset_draine_path,
            offset_sed_path,
            edges_path,
            expected_bolometric_fraction=1.0,
        )
        independent_opacity, independent_energy = _independent_draine_recompute(
            offset_sed_path, offset_draine_path, edges
        )
        assert np.allclose(
            offset_metadata["absorption_cross_section_per_h_cm2"],
            independent_opacity,
            rtol=2.0e-5,
            atol=1.0e-30,
        )
        assert np.allclose(
            offset_metadata["absorption_weighted_energy_ev"],
            independent_energy,
            rtol=2.0e-5,
            atol=1.0e-12,
        )
        closed_opacity, closed_energy = _powerlaw_closed_form(edges, 4.0e-21, 0.5)
        assert np.allclose(
            offset_metadata["absorption_cross_section_per_h_cm2"],
            closed_opacity,
            rtol=2.0e-4,
            atol=1.0e-30,
        )
        assert np.allclose(
            offset_metadata["absorption_weighted_energy_ev"],
            closed_energy,
            rtol=2.0e-4,
            atol=1.0e-12,
        )

        wrong_edge_hash = dict(dust_metadata)
        wrong_edge_hash["group_edges_sha256"] = "0" * 64
        wrong_edge_path = work / "wrong-edge-hash.json"
        wrong_edge_path.write_text(json.dumps(wrong_edge_hash) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(wrong_edge_path, expected_group_edges_ev=edges)
        except ValueError as error:
            assert "group-edge hash" in str(error)
        else:
            raise AssertionError("tampered group-edge hash was accepted")

        wrong_table_hash = dict(dust_metadata)
        wrong_table_hash["source_table"] = dict(dust_metadata["source_table"])
        wrong_table_hash["source_table"]["sha256"] = "0" * 64
        wrong_table_path = work / "wrong-table-hash.json"
        wrong_table_path.write_text(json.dumps(wrong_table_hash) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(wrong_table_path, expected_group_edges_ev=edges)
        except ValueError as error:
            assert "source-table hash" in str(error)
        else:
            raise AssertionError("tampered Draine source hash was accepted")

        wrong_builder_hash = dict(dust_metadata)
        wrong_builder_hash["builder"] = dict(dust_metadata["builder"])
        wrong_builder_hash["builder"]["sha256"] = "0" * 64
        wrong_builder_path = work / "wrong-builder-hash.json"
        wrong_builder_path.write_text(json.dumps(wrong_builder_hash) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(wrong_builder_path, expected_group_edges_ev=edges)
        except ValueError as error:
            assert "builder hash" in str(error)
        else:
            raise AssertionError("tampered builder hash was accepted")

        wrong_manifest = dict(dust_metadata)
        wrong_manifest["closure_code_manifest"] = [
            entry
            for entry in dust_metadata["closure_code_manifest"]
            if entry["role"] != "source_sed"
        ]
        wrong_manifest_path = work / "wrong-manifest.json"
        wrong_manifest_path.write_text(json.dumps(wrong_manifest) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(wrong_manifest_path, expected_group_edges_ev=edges)
        except ValueError as error:
            assert "closure code manifest roles" in str(error)
        else:
            raise AssertionError("incomplete closure code manifest was accepted")

        wrong_payload = dict(dust_metadata)
        wrong_payload["absorption_cross_section_per_h_cm2"] = [
            float(value) * 1.01 for value in dust_metadata["absorption_cross_section_per_h_cm2"]
        ]
        wrong_payload_path = work / "wrong-payload.json"
        wrong_payload_path.write_text(json.dumps(wrong_payload) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(wrong_payload_path, expected_group_edges_ev=edges)
        except ValueError as error:
            assert "payload hash" in str(error)
        else:
            raise AssertionError("tampered dust payload was accepted")

        wrong_status = dict(dust_metadata)
        wrong_status["status"] = "production_approved"
        wrong_status_path = work / "wrong-status.json"
        wrong_status_path.write_text(json.dumps(wrong_status) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(wrong_status_path, expected_group_edges_ev=edges)
        except ValueError as error:
            assert "binding status" in str(error)
        else:
            raise AssertionError("unrecognized dust status was accepted")

        try:
            read_dust_opacity_metadata(
                dust_path,
                expected_group_edges_ev=edges,
                expected_group_edges_sha256="f" * 64,
            )
        except ValueError as error:
            assert "photon metadata" in str(error)
        else:
            raise AssertionError("mismatched expected group-edge hash was accepted")

        try:
            read_dust_opacity_metadata(
                dust_path,
                expected_group_edges_ev=edges,
                expected_source_sed_identity=sed.identity,
                expected_source_sed_sha256="0" * 64,
                require_source_match=True,
            )
        except ValueError as error:
            assert "input hash" in str(error)
        else:
            raise AssertionError("mismatched expected source hash was accepted")
        try:
            read_dust_opacity_metadata(
                dust_path,
                expected_group_edges_ev=edges,
                expected_source_sed_identity="0" * 64,
                require_source_match=True,
            )
        except ValueError as error:
            assert "identity" in str(error)
        else:
            raise AssertionError("mismatched source SED identity was accepted")

        reference_path = work / "reference.json"
        reference_path.write_text(
            json.dumps(
                {
                    "schema": "snrt_dust_opacity_v1",
                    "group_edges_ev": edges.tolist(),
                    "absorption_cross_section_per_h_cm2": [1.0, 1.0],
                    "absorption_weighted_energy_ev": [3.0, 30.0],
                    "reference_mixture": "test",
                    "opacity_source": "test",
                    "spectral_weighting": "reference",
                }
            )
            + "\n",
            encoding="utf-8",
        )
        try:
            read_dust_opacity_metadata(
                reference_path,
                expected_group_edges_ev=edges,
                expected_source_sed_identity=sed.identity,
                require_source_match=True,
            )
        except ValueError as error:
            assert "reference" in str(error)
        else:
            raise AssertionError("v1 reference dust closure was accepted for a bound SED")

        full_edges = np.asarray((0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500.0, 2000.0, 10000.0))
        full_edges_path = work / "p0_edges.txt"
        full_edges_path.write_text("\n".join(str(value) for value in full_edges) + "\n", encoding="utf-8")
        explicit_sed_path = work / "agn_sed.csv"
        _write_sed(explicit_sed_path, full_edges)
        candidate_path = work / "candidates.csv"
        _write_candidate(candidate_path)
        output_path = work / "agn.csv"
        metadata_path = work / "agn.json"
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(ROOT)
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "p4_build_agn_photon_ledger.py"),
                "--candidates",
                str(candidate_path),
                "--output",
                str(output_path),
                "--metadata-output",
                str(metadata_path),
                "--group-edges",
                str(full_edges_path),
                "--sed-table",
                str(explicit_sed_path),
                "--sed-bolometric-fraction",
                "1.0",
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        agn_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        explicit_sed = read_lbol_photon_sed(explicit_sed_path, expected_bolometric_fraction=1.0)
        assert agn_metadata["source_sed_identity"]
        assert agn_metadata["source_sed_identity"] == explicit_sed.identity
        assert agn_metadata["source_sed_sha256"] == explicit_sed.input_sha256
        assert agn_metadata["payload_sha256"]
        assert {entry["role"] for entry in agn_metadata["closure_code_manifest"]} == {
            "agn_ledger_builder",
            "source_sed",
            "primordial_closure",
            "integrity_helper",
        }
        assert "Sazonov" not in agn_metadata["reference"]
        assert "nu_lnu_at_13p6_ev_over_lbol" not in agn_metadata["normalization"]
        assert all(
            group["sed_supported_interval_ev"] == [float(low), float(high)]
            for group, low, high in zip(
                agn_metadata["groups"], full_edges[:-1], full_edges[1:], strict=True
            )
        )

        # The AGN source-bound side is bound to the explicit AGN SED used by
        # the photon ledger, with the same nine-group edges; it is not the
        # unbound reference control or the stellar fixture above.
        agn_scattering_draine_path = work / "agn-scattering-draine.all"
        _write_draine(
            agn_scattering_draine_path,
            tuple(float(value) for value in full_edges),
            albedo=0.5,
            cosine=0.4,
            cosine_squared=0.2,
        )
        agn_scattering_metadata = build_draine_dust_opacity.build_source_weighted_opacity_metadata(
            agn_scattering_draine_path,
            explicit_sed_path,
            full_edges_path,
            expected_bolometric_fraction=1.0,
            include_scattering=True,
        )
        agn_scattering_path = work / "agn-scattering-dust.json"
        agn_scattering_path.write_text(json.dumps(agn_scattering_metadata) + "\n", encoding="utf-8")
        agn_scattering_closure = read_dust_opacity_metadata(
            agn_scattering_path,
            expected_group_edges_ev=full_edges,
            expected_source_sed_identity=explicit_sed.identity,
            expected_source_sed_sha256=explicit_sed.input_sha256,
            require_source_match=True,
        )
        assert agn_scattering_closure.schema == "snrt_dust_opacity_v3"
        assert agn_scattering_closure.source_sed_identity == explicit_sed.identity
        assert np.any(agn_scattering_closure.scattering_cross_section_per_h_cm2 > 0.0)

        # An unbound v3 fixture exercises the declared missing-array contract
        # without source-bound payload validation masking the error.
        malformed_v3 = dict(agn_scattering_metadata)
        for key in (
            "group_edges_path",
            "group_edges_sha256",
            "source_table",
            "source_sed_identity",
            "source_sed_sha256",
            "source_sed_contract",
            "source_sed_group_energy_fraction_of_lbol",
            "source_sed_quadrature",
            "groups",
            "builder",
            "closure_code_manifest",
            "payload_hash_scheme",
            "payload_sha256",
        ):
            malformed_v3.pop(key, None)
        malformed_v3["status"] = "reference_scattering_control"
        malformed_v3.pop("scattering_angle_cosine_squared")
        malformed_path = work / "malformed-v3-dust.json"
        malformed_path.write_text(json.dumps(malformed_v3) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(malformed_path, expected_group_edges_ev=full_edges)
        except ValueError as error:
            assert "missing dust opacity fields" in str(error)
        else:
            raise AssertionError("malformed v3 sidecar was not rejected as ValueError")
        null_v3 = dict(agn_scattering_metadata)
        for key in (
            "group_edges_path",
            "group_edges_sha256",
            "source_table",
            "source_sed_identity",
            "source_sed_sha256",
            "source_sed_contract",
            "source_sed_group_energy_fraction_of_lbol",
            "source_sed_quadrature",
            "groups",
            "builder",
            "closure_code_manifest",
            "payload_hash_scheme",
            "payload_sha256",
        ):
            null_v3.pop(key, None)
        null_v3["status"] = "reference_scattering_control"
        null_v3["scattering_weighted_energy_ev"] = None
        null_path = work / "null-v3-dust.json"
        null_path.write_text(json.dumps(null_v3) + "\n", encoding="utf-8")
        try:
            read_dust_opacity_metadata(null_path, expected_group_edges_ev=full_edges)
        except ValueError as error:
            assert "v3 scattering fields cannot be null" in str(error)
        else:
            raise AssertionError("null v3 scattering field was accepted")
        assert np.isclose(
            sum(group["energy_fraction_of_lbol"] for group in agn_metadata["groups"]),
            1.0,
        )

        # Check that intrinsic per-Lbol fields and escaped emitted totals stay
        # distinguishable when the escape fraction is not unity.
        quarter_output = work / "agn-quarter.csv"
        quarter_metadata_path = work / "agn-quarter.json"
        quarter_command = list(
            [
                sys.executable,
                str(ROOT / "tools" / "p4_build_agn_photon_ledger.py"),
                "--candidates",
                str(candidate_path),
                "--output",
                str(quarter_output),
                "--metadata-output",
                str(quarter_metadata_path),
                "--group-edges",
                str(full_edges_path),
                "--sed-table",
                str(explicit_sed_path),
                "--sed-bolometric-fraction",
                "1.0",
                "--escape-fraction",
                "0.25",
            ]
        )
        quarter = subprocess.run(
            quarter_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if quarter.returncode != 0:
            raise RuntimeError(quarter.stdout + quarter.stderr)
        quarter_metadata = json.loads(quarter_metadata_path.read_text(encoding="utf-8"))
        assert np.isclose(
            quarter_metadata["groups"][3]["escaped_energy_fraction_of_lbol"],
            0.25 * quarter_metadata["groups"][3]["energy_fraction_of_lbol"],
        )
        assert np.isclose(
            quarter_metadata["groups"][3]["escaped_photon_rate_per_lbol_s_per_erg_s"],
            0.25 * quarter_metadata["groups"][3]["photon_rate_per_lbol_s_per_erg_s"],
        )

        # Exercise the actual P4 runner's v2 binding and its pre-output
        # mismatch rejection, not only the metadata loader in isolation.
        full_edges = np.asarray((0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500.0, 2000.0, 10000.0))
        full_draine_path = work / "draine-full.all"
        _write_draine(full_draine_path, tuple(float(value) for value in full_edges))
        bound_dust_metadata = build_draine_dust_opacity.build_source_weighted_opacity_metadata(
            full_draine_path,
            explicit_sed_path,
            full_edges_path,
            expected_bolometric_fraction=1.0,
        )
        bound_dust_path = work / "bound-dust.json"
        bound_dust_path.write_text(json.dumps(bound_dust_metadata) + "\n", encoding="utf-8")
        photon_rates = np.asarray(
            [group["photon_rate_per_lbol_s_per_erg_s"] for group in agn_metadata["groups"]],
            dtype=np.float64,
        ) * 1.0e35
        static_path = work / "static.h5"
        static = neutral_primordial_input(
            GridSpec(cell_width_cm=1.0e18, left_edge_cm=np.zeros(3)),
            np.full((2, 2, 2), 1.0e-24),
            np.full((2, 2, 2), 1.0e4),
            dust_relative_abundance=1.0,
            sources=SourceCatalog(
                cell_index=np.asarray([[0, 0, 0]]),
                photon_luminosity_s=photon_rates[None, :],
            ),
        )
        write_static_rt_input(static_path, static)
        runner_output = work / "runner.h5"
        runner_command = [
            sys.executable,
            str(ROOT / "tools" / "p4_run_transport_pilot.py"),
            "--input",
            str(static_path),
            "--photon-metadata",
            str(metadata_path),
            "--dust-opacity-metadata",
            str(bound_dust_path),
            "--steps",
            "1",
            "--sn-order",
            "4",
            "--precision",
            "float64",
            "--output",
            str(runner_output),
        ]
        result = subprocess.run(
            runner_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        with h5py.File(runner_output, "r") as handle:
            assert handle.attrs["dust_opacity_schema"] == "snrt_dust_opacity_v2"
            assert handle.attrs["dust_binding_status"] == "candidate_source_sed_matched"
            assert handle.attrs["source_sed_identity"] == explicit_sed.identity
            assert handle.attrs["group_edges_sha256"] == agn_metadata["group_edges_sha256"]
            assert handle.attrs["dust_opacity_metadata_sha256"]
            assert handle.attrs["dust_payload_sha256"] == bound_dust_metadata["payload_sha256"]
            assert handle.attrs["dust_source_table_sha256"] == bound_dust_metadata["source_table"]["sha256"]
            assert handle.attrs["dust_builder_sha256"] == bound_dust_metadata["builder"]["sha256"]

        mismatched_dust = dict(bound_dust_metadata)
        mismatched_dust["source_sed_identity"] = "0" * 64
        mismatched_path = work / "mismatched-dust.json"
        mismatched_path.write_text(json.dumps(mismatched_dust) + "\n", encoding="utf-8")
        rejected_output = work / "rejected.h5"
        bad_command = list(runner_command)
        bad_command[bad_command.index(str(bound_dust_path))] = str(mismatched_path)
        bad_command[bad_command.index(str(runner_output))] = str(rejected_output)
        rejected = subprocess.run(
            bad_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        assert rejected.returncode != 0
        assert not rejected_output.exists()

        bad_photon_metadata = dict(agn_metadata)
        bad_photon_metadata["group_edges_sha256"] = "0" * 64
        bad_photon_path = work / "bad-photon-metadata.json"
        bad_photon_path.write_text(json.dumps(bad_photon_metadata) + "\n", encoding="utf-8")
        bad_photon_output = work / "bad-photon-output.h5"
        bad_photon_command = list(runner_command)
        bad_photon_command[bad_photon_command.index(str(metadata_path))] = str(bad_photon_path)
        bad_photon_command[bad_photon_command.index(str(runner_output))] = str(bad_photon_output)
        rejected_photon = subprocess.run(
            bad_photon_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        assert rejected_photon.returncode != 0
        assert not bad_photon_output.exists()

    print("SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
