#!/usr/bin/env python3
"""Analyze matched FDM wave, hybrid, and fluid RAMSES outputs.

The script reads the native AMR and FDM files.  It deposits FDM leaf-cell
mass onto one common mesh, then measures power spectra, the wave--hybrid
cross-correlation, a density-PDF KS distance, and mass drift from run logs.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import fft as scipy_fft
from scipy.stats import ks_2samp

from compute_pk import FortranReader, get_nchar, read_amr_header, read_info


MASS_RE = re.compile(r"FDM:\s+M_tot=\s*([+\-0-9.EDed]+)")
STEP_RE = re.compile(
    r"Fine step=\s*(\d+).*?\ba=\s*([+\-0-9.EDed]+)"
)


def fortran_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def read_fdm_mode(output_dir: Path) -> tuple[bool, int]:
    candidates = [
        output_dir / "namelist.txt",
        output_dir.parent / "cosmo.nml",
    ]
    texts = [
        path.read_text(errors="replace")
        for path in candidates
        if path.is_file()
    ]

    def scalar(name: str, default: str) -> str:
        pattern = re.compile(
            rf"(?im)^\s*{re.escape(name)}\s*=\s*([^!\n/]+)"
        )
        for text in texts:
            match = pattern.search(text)
            if match:
                return match.group(1).strip()
        return default

    use_hjm = scalar("fdm_use_hjm", ".false.").lower() == ".true."
    first_wave = int(scalar("fdm_first_wave_level", "0"))
    return use_hjm, first_wave


def deposit_worker(task: tuple) -> dict[str, object]:
    (
        label,
        output_name,
        cpu_ids,
        ngrid,
        use_hjm,
        first_wave,
        scratch_name,
        worker_id,
    ) = task
    output_dir = Path(output_name)
    info = read_info(str(output_dir))
    boxlen = float(info["boxlen"])
    nchar = get_nchar(str(output_dir))
    grid = np.zeros(ngrid**3, dtype=np.float64)

    leaf_count = 0
    negative_leaf_count = 0
    leaf_mass = 0.0

    for icpu in cpu_ids:
        amr_path = output_dir / f"amr_{nchar}.out{icpu:05d}"
        fdm_path = output_dir / f"fdm_{nchar}.out{icpu:05d}"
        if not amr_path.is_file() or not fdm_path.is_file():
            raise FileNotFoundError(f"missing CPU {icpu} data in {output_dir}")

        with FortranReader(str(amr_path)) as amr, FortranReader(str(fdm_path)) as fdm:
            header = read_amr_header(amr)
            fdm_ncpu = int(fdm.read_ints()[0])
            fdm_ndim = int(fdm.read_ints()[0])
            fdm_nlevelmax = int(fdm.read_ints()[0])
            fdm_nboundary = int(fdm.read_ints()[0])

            ncpu = int(header["ncpu"])
            ndim = int(header["ndim"])
            nlevelmax = int(header["nlevelmax"])
            nboundary = int(header["nboundary"])
            nx = int(header["nx"])
            numbl = header["numbl"]
            numbb = header["numbb"]

            if (fdm_ncpu, fdm_ndim, fdm_nlevelmax, fdm_nboundary) != (
                ncpu,
                ndim,
                nlevelmax,
                nboundary,
            ):
                raise ValueError(f"AMR/FDM header mismatch in CPU {icpu}")

            twotondim = 2**ndim
            twondim = 2 * ndim

            for ilevel in range(1, nlevelmax + 1):
                dx = boxlen / (nx * 2**ilevel)
                cell_volume = dx**ndim
                fluid_level = use_hjm and ilevel < first_wave

                for ibound in range(1, ncpu + nboundary + 1):
                    if ibound <= ncpu:
                        ncache = int(numbl[ibound - 1, ilevel - 1])
                    else:
                        ncache = int(numbb[ibound - ncpu - 1, ilevel - 1])

                    fdm_level = int(fdm.read_ints()[0])
                    fdm_ncache = int(fdm.read_ints()[0])
                    if fdm_level != ilevel or fdm_ncache != ncache:
                        raise ValueError(
                            f"AMR/FDM block mismatch CPU={icpu} "
                            f"level={ilevel} boundary={ibound}"
                        )

                    if ncache <= 0:
                        continue

                    amr.skip(3)
                    xg = np.empty((ncache, ndim), dtype=np.float64)
                    for axis in range(ndim):
                        xg[:, axis] = amr.read_reals()
                    amr.skip(1)
                    amr.skip(twondim)

                    son = np.empty((ncache, twotondim), dtype=np.int32)
                    for child in range(twotondim):
                        son[:, child] = amr.read_ints()
                    amr.skip(twotondim)
                    amr.skip(twotondim)

                    owner = ibound == icpu
                    for child in range(twotondim):
                        field1 = fdm.read_reals()
                        field2 = fdm.read_reals()
                        if not owner:
                            continue

                        leaf = son[:, child] == 0
                        if not np.any(leaf):
                            continue

                        if fluid_level:
                            density = field1[leaf]
                        else:
                            density = field1[leaf] ** 2 + field2[leaf] ** 2

                        iz = child // 4
                        iy = (child - 4 * iz) // 2
                        ix = child - 2 * iy - 4 * iz
                        offsets = (ix - 0.5, iy - 0.5, iz - 0.5)
                        coords = [
                            xg[leaf, axis] + offsets[axis] * dx
                            for axis in range(ndim)
                        ]
                        indices = [
                            np.floor(coord / boxlen * ngrid).astype(np.int64) % ngrid
                            for coord in coords
                        ]
                        flat = (
                            indices[0] * ngrid * ngrid
                            + indices[1] * ngrid
                            + indices[2]
                        )
                        cell_mass = density * cell_volume
                        np.add.at(grid, flat, cell_mass)

                        leaf_count += int(density.size)
                        negative_leaf_count += int(np.count_nonzero(density <= 0.0))
                        leaf_mass += float(np.sum(cell_mass, dtype=np.float64))

    partial_path = Path(scratch_name) / f"{label}.part{worker_id:03d}.npy"
    np.save(partial_path, grid)
    return {
        "path": str(partial_path),
        "leaf_count": leaf_count,
        "negative_leaf_count": negative_leaf_count,
        "leaf_mass": leaf_mass,
        "cpu_count": len(cpu_ids),
    }


def deposit_arm(
    label: str,
    output_dir: Path,
    ngrid: int,
    workers: int,
    scratch: Path,
) -> tuple[np.ndarray, dict[str, object]]:
    info = read_info(str(output_dir))
    ncpu = int(info["ncpu"])
    use_hjm, first_wave = read_fdm_mode(output_dir)
    cpu_chunks = [
        list(range(1 + worker, ncpu + 1, workers))
        for worker in range(workers)
    ]
    tasks = [
        (
            label,
            str(output_dir),
            cpu_ids,
            ngrid,
            use_hjm,
            first_wave,
            str(scratch),
            worker,
        )
        for worker, cpu_ids in enumerate(cpu_chunks)
        if cpu_ids
    ]

    print(
        f"[{label}] deposit {output_dir} with {len(tasks)} workers "
        f"onto {ngrid}^3",
        flush=True,
    )
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(tasks)) as pool:
        partials = list(pool.map(deposit_worker, tasks))

    grid = np.zeros(ngrid**3, dtype=np.float64)
    for partial in partials:
        path = Path(str(partial["path"]))
        data = np.load(path, mmap_mode="r")
        grid += data
        del data
        path.unlink()

    leaf_count = sum(int(item["leaf_count"]) for item in partials)
    negative_count = sum(int(item["negative_leaf_count"]) for item in partials)
    leaf_mass = sum(float(item["leaf_mass"]) for item in partials)
    grid_mass = float(np.sum(grid, dtype=np.float64))
    summary = {
        "output": str(output_dir),
        "aexp": float(info["aexp"]),
        "redshift": 1.0 / float(info["aexp"]) - 1.0,
        "ncpu": ncpu,
        "use_hjm": use_hjm,
        "first_wave_level": first_wave,
        "leaf_count": leaf_count,
        "negative_leaf_count": negative_count,
        "negative_leaf_fraction": negative_count / leaf_count if leaf_count else math.nan,
        "leaf_mass": leaf_mass,
        "deposited_mass": grid_mass,
        "deposit_roundoff": grid_mass - leaf_mass,
    }
    return grid, summary


def mass_history(log_path: Path, target_a: float, tolerance: float) -> dict[str, object]:
    series: list[tuple[float, float, int]] = []
    pending_step: tuple[float, int] | None = None
    for line in log_path.read_text(errors="replace").splitlines():
        step_match = STEP_RE.search(line)
        if step_match:
            pending_step = (
                fortran_float(step_match.group(2)),
                int(step_match.group(1)),
            )

        # RAMSES prints M_tot after the Fine step that it describes.  Pair it
        # with that pending step instead of carrying it into the next one.
        mass_match = MASS_RE.search(line)
        if mass_match and pending_step is not None:
            series.append(
                (
                    pending_step[0],
                    fortran_float(mass_match.group(1)),
                    pending_step[1],
                )
            )
            pending_step = None

    if not series:
        return {"pass": False, "reason": "no paired mass/step diagnostics"}

    first_a, first_mass, first_step = series[0]
    endpoint = min(series, key=lambda item: abs(item[0] - target_a))
    upto = [item for item in series if item[0] <= endpoint[0] + 1.0e-12]
    drifts = [(item[1] - first_mass) / first_mass for item in upto]
    endpoint_drift = (endpoint[1] - first_mass) / first_mass
    max_abs_drift = max(abs(value) for value in drifts)
    return {
        "log": str(log_path),
        "initial_a": first_a,
        "initial_step": first_step,
        "initial_mass": first_mass,
        "target_a": target_a,
        "matched_a": endpoint[0],
        "matched_step": endpoint[2],
        "matched_mass": endpoint[1],
        "endpoint_fractional_drift": endpoint_drift,
        "max_abs_fractional_drift": max_abs_drift,
        "sample_count": len(upto),
        "tolerance": tolerance,
        "pass": max_abs_drift <= tolerance,
    }


def fft_modes(ngrid: int) -> tuple[np.ndarray, np.ndarray]:
    kx = np.fft.fftfreq(ngrid) * ngrid
    ky = np.fft.fftfreq(ngrid) * ngrid
    kz = np.arange(ngrid // 2 + 1, dtype=np.float64)
    mode = np.sqrt(
        kx[:, None, None] ** 2
        + ky[None, :, None] ** 2
        + kz[None, None, :] ** 2
    )
    wz = np.ones(ngrid // 2 + 1, dtype=np.float64)
    if ngrid > 2:
        wz[1 : ngrid // 2] = 2.0
    weight = np.broadcast_to(wz[None, None, :], mode.shape)
    return mode, weight


def binned_spectra(
    rho: dict[str, np.ndarray],
    box_mpc_h: float,
    ngrid: int,
    nkbin: int,
    workers: int,
) -> dict[str, np.ndarray]:
    transforms = {}
    for label, field in rho.items():
        delta = field.reshape((ngrid, ngrid, ngrid)) - 1.0
        transforms[label] = scipy_fft.rfftn(
            delta,
            workers=workers,
            # The real-space fields are reused by density_pdf below.
            overwrite_x=False,
        )

    mode, symmetry = fft_modes(ngrid)
    edges_mode = np.geomspace(1.0, ngrid / 2.0, nkbin + 1)
    flat_mode = mode.ravel()
    flat_symmetry = symmetry.ravel()
    bins = np.searchsorted(edges_mode, flat_mode, side="right") - 1
    valid = (
        (flat_mode >= edges_mode[0])
        & (flat_mode < edges_mode[-1])
        & (bins >= 0)
        & (bins < nkbin)
    )
    selected_bins = bins[valid]
    selected_weights = flat_symmetry[valid]
    counts = np.bincount(
        selected_bins,
        weights=selected_weights,
        minlength=nkbin,
    )
    mode_sum = np.bincount(
        selected_bins,
        weights=selected_weights * flat_mode[valid],
        minlength=nkbin,
    )
    mean_mode = np.divide(
        mode_sum,
        counts,
        out=np.full(nkbin, np.nan),
        where=counts > 0,
    )

    volume = box_mpc_h**3
    normalization = volume / ngrid**6
    output: dict[str, np.ndarray] = {
        "k": mean_mode * (2.0 * np.pi / box_mpc_h),
        "counts": counts,
        "k_edges": edges_mode * (2.0 * np.pi / box_mpc_h),
    }

    for label, transform in transforms.items():
        power = np.abs(transform.ravel()[valid]) ** 2
        power_sum = np.bincount(
            selected_bins,
            weights=selected_weights * power,
            minlength=nkbin,
        )
        output[f"pk_{label}"] = np.divide(
            power_sum,
            counts,
            out=np.full(nkbin, np.nan),
            where=counts > 0,
        ) * normalization

    cross = np.real(
        transforms["hybrid"].ravel()[valid]
        * np.conjugate(transforms["wave"].ravel()[valid])
    )
    cross_sum = np.bincount(
        selected_bins,
        weights=selected_weights * cross,
        minlength=nkbin,
    )
    output["pk_cross_hybrid_wave"] = np.divide(
        cross_sum,
        counts,
        out=np.full(nkbin, np.nan),
        where=counts > 0,
    ) * normalization
    output["r_hybrid_wave"] = output["pk_cross_hybrid_wave"] / np.sqrt(
        output["pk_hybrid"] * output["pk_wave"]
    )
    output["ratio_hybrid_wave"] = output["pk_hybrid"] / output["pk_wave"]
    output["ratio_fluid_wave"] = output["pk_fluid"] / output["pk_wave"]
    return output


def density_pdf(
    wave: np.ndarray,
    hybrid: np.ndarray,
    lower: float = -4.0,
    upper: float = 4.0,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    samples = {}
    retained = {}
    nonpositive = {}
    for label, field in (("wave", wave), ("hybrid", hybrid)):
        positive = field > 0.0
        nonpositive[label] = 1.0 - float(np.mean(positive))
        log_density = np.log10(field[positive])
        keep = (log_density >= lower) & (log_density <= upper)
        samples[label] = log_density[keep]
        retained[label] = samples[label].size / field.size

    ks = ks_2samp(samples["hybrid"], samples["wave"], method="asymp")
    edges = np.linspace(lower, upper, 257)
    histograms = {
        "logrho_edges": edges,
        "pdf_wave": np.histogram(samples["wave"], bins=edges, density=True)[0],
        "pdf_hybrid": np.histogram(samples["hybrid"], bins=edges, density=True)[0],
    }
    report = {
        "log10_density_range": [lower, upper],
        "ks_distance": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
        "retained_fraction": retained,
        "nonpositive_fraction": nonpositive,
        "sample_count": {key: int(value.size) for key, value in samples.items()},
        "tolerance": 0.05,
        "pass": float(ks.statistic) <= 0.05,
    }
    return report, histograms


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: object) -> object:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--arm",
        action="append",
        nargs=3,
        metavar=("LABEL", "OUTPUT", "LOG"),
        required=True,
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ngrid", type=int, default=256)
    parser.add_argument("--nkbin", type=int, default=48)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--target-a", type=float, default=0.4)
    parser.add_argument("--mass-tol", type=float, default=0.002)
    args = parser.parse_args()

    arms = {
        label: {"output": Path(output), "log": Path(log)}
        for label, output, log in args.arm
    }
    required = {"wave", "hybrid", "fluid"}
    if set(arms) != required:
        raise ValueError(f"arms must be exactly {sorted(required)}")
    for label, arm in arms.items():
        if not (arm["output"] / "COMPLETE").is_file():
            raise FileNotFoundError(f"{label} output lacks COMPLETE")
        if not arm["log"].is_file():
            raise FileNotFoundError(f"{label} log is missing")

    args.output.mkdir(parents=True, exist_ok=True)
    scratch_parent = str(args.scratch) if args.scratch else None
    scratch = Path(tempfile.mkdtemp(prefix="fdm-matched-", dir=scratch_parent))
    fields = {}
    arm_reports = {}
    try:
        for label in ("wave", "hybrid", "fluid"):
            grid, deposit = deposit_arm(
                label,
                arms[label]["output"],
                args.ngrid,
                args.workers,
                scratch,
            )
            mean_mass = float(np.mean(grid))
            if mean_mass <= 0.0:
                raise ValueError(f"{label} deposited mass is not positive")
            fields[label] = grid / mean_mass
            del grid
            np.save(args.output / f"rho_{label}_{args.ngrid}.npy", fields[label])

            history = mass_history(
                arms[label]["log"],
                args.target_a,
                args.mass_tol,
            )
            deposit["mass_history"] = history
            deposit["complete_sha256"] = sha256(arms[label]["output"] / "COMPLETE")
            arm_reports[label] = deposit

        a_values = [float(arm_reports[label]["aexp"]) for label in required]
        if max(a_values) - min(a_values) > 1.0e-8:
            raise ValueError(f"output scale factors do not match: {a_values}")
        info = read_info(str(arms["wave"]["output"]))
        nchar = get_nchar(str(arms["wave"]["output"]))
        with FortranReader(
            str(arms["wave"]["output"] / f"amr_{nchar}.out00001")
        ) as amr:
            box_mpc_h = float(read_amr_header(amr)["boxlen_ini"])

        spectra = binned_spectra(
            fields,
            box_mpc_h,
            args.ngrid,
            args.nkbin,
            args.workers,
        )
        k_nyquist = np.pi * args.ngrid / box_mpc_h
        accepted = (
            np.isfinite(spectra["ratio_hybrid_wave"])
            & np.isfinite(spectra["r_hybrid_wave"])
            & (spectra["k"] <= k_nyquist / 4.0)
            & (spectra["counts"] >= 100.0)
        )
        if not np.any(accepted):
            raise ValueError("no power-spectrum bins satisfy the preregistered range")
        max_pk_error = float(
            np.max(np.abs(spectra["ratio_hybrid_wave"][accepted] - 1.0))
        )
        min_cross = float(np.min(spectra["r_hybrid_wave"][accepted]))
        pk_report = {
            "analysis_grid": args.ngrid,
            "box_mpc_h": box_mpc_h,
            "k_nyquist_h_mpc": k_nyquist,
            "k_accept_max_h_mpc": k_nyquist / 4.0,
            "minimum_modes_per_bin": 100,
            "accepted_bin_count": int(np.count_nonzero(accepted)),
            "accepted_k_min_h_mpc": float(np.min(spectra["k"][accepted])),
            "accepted_k_max_h_mpc": float(np.max(spectra["k"][accepted])),
            "max_abs_hybrid_wave_power_error": max_pk_error,
            "minimum_hybrid_wave_cross_correlation": min_cross,
            "power_tolerance": 0.05,
            "cross_correlation_tolerance": 0.99,
            "power_pass": max_pk_error <= 0.05,
            "cross_correlation_pass": min_cross >= 0.99,
        }

        pdf_report, pdf_arrays = density_pdf(
            fields["wave"],
            fields["hybrid"],
        )
        mass_pass = all(
            bool(arm_reports[label]["mass_history"]["pass"])
            for label in required
        )
        core_pass = (
            bool(pk_report["power_pass"])
            and bool(pk_report["cross_correlation_pass"])
            and bool(pdf_report["pass"])
            and mass_pass
        )
        report = {
            "status": "CORE_PASS_HALO_PENDING" if core_pass else "FAIL",
            "core_pass": core_pass,
            "target_a": args.target_a,
            "output_a": a_values[0],
            "arms": arm_reports,
            "power_and_cross": pk_report,
            "density_pdf": pdf_report,
            "mass_pass": mass_pass,
            "halo_matching": {
                "status": "not_evaluated",
                "required_only_if_resolved_halo_count_at_least": 20,
            },
            "provenance": {
                "script": str(Path(__file__).resolve()),
                "script_sha256": sha256(Path(__file__).resolve()),
                "python": sys.version,
                "numpy": np.__version__,
            },
        }

        np.savez(
            args.output / "spectra_and_pdf.npz",
            **spectra,
            accepted=accepted,
            **pdf_arrays,
        )
        (args.output / "report.json").write_text(
            json.dumps(json_safe(report), indent=2, sort_keys=True) + "\n"
        )
        with (args.output / "power_bins.csv").open("w") as stream:
            stream.write(
                "k_h_mpc,modes,pk_wave,pk_hybrid,pk_fluid,"
                "ratio_hybrid_wave,ratio_fluid_wave,r_hybrid_wave,accepted\n"
            )
            for index in range(args.nkbin):
                stream.write(
                    f"{spectra['k'][index]:.10e},"
                    f"{spectra['counts'][index]:.0f},"
                    f"{spectra['pk_wave'][index]:.10e},"
                    f"{spectra['pk_hybrid'][index]:.10e},"
                    f"{spectra['pk_fluid'][index]:.10e},"
                    f"{spectra['ratio_hybrid_wave'][index]:.10e},"
                    f"{spectra['ratio_fluid_wave'][index]:.10e},"
                    f"{spectra['r_hybrid_wave'][index]:.10e},"
                    f"{int(accepted[index])}\n"
                )

        print(json.dumps(json_safe(report), indent=2, sort_keys=True))
        return 0 if core_pass else 2
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
