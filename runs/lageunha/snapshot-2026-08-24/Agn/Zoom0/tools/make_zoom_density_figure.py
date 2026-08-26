#!/usr/bin/env python3
"""Render a full-box projected DMO density map with a zoom-region panel."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle
from scipy.io import FortranFile
from scipy.ndimage import gaussian_filter


def deposit_cic(
    x: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    bins: int,
    bounds: tuple[float, float, float, float],
    periodic: bool,
) -> np.ndarray:
    """Deposit particles on a cell-centred 2-D mesh with CIC weights."""
    xmin, xmax, ymin, ymax = bounds
    qx = (x - xmin) * bins / (xmax - xmin) - 0.5
    qy = (y - ymin) * bins / (ymax - ymin) - 0.5
    ix = np.floor(qx).astype(np.int64)
    iy = np.floor(qy).astype(np.int64)
    fx = qx - ix
    fy = qy - iy
    grid = np.zeros((bins, bins), dtype=np.float64)

    for dx, wx in ((0, 1.0 - fx), (1, fx)):
        for dy, wy in ((0, 1.0 - fy), (1, fy)):
            gx = ix + dx
            gy = iy + dy
            contribution = weights * wx * wy
            if periodic:
                np.add.at(grid, (gx % bins, gy % bins), contribution)
            else:
                inside = (gx >= 0) & (gx < bins) & (gy >= 0) & (gy < bins)
                np.add.at(
                    grid,
                    (gx[inside], gy[inside]),
                    contribution[inside],
                )
    return grid


def project_part(
    part_path: str,
    full_bins: int,
    zoom_bins: int,
    zoom_bounds: tuple[float, float, float, float],
    shift: tuple[float, float, float],
    assignment: str,
    max_particle_mass: float | None,
) -> tuple[np.ndarray, np.ndarray, float, int, int]:
    path = pathlib.Path(part_path)
    with FortranFile(path, "r") as handle:
        handle.read_record(np.int32)
        ndim = int(handle.read_record(np.int32)[0])
        npart = int(handle.read_record(np.int32)[0])
        for _ in range(5):
            handle.read_record(np.uint8)
        positions = [handle.read_record(np.float64) for _ in range(ndim)]
        for _ in range(ndim):
            handle.read_record(np.float64)
        masses = handle.read_record(np.float64)
    if ndim != 3 or masses.size != npart:
        raise RuntimeError(f"Invalid RAMSES particle records in {path}")

    x_full = (positions[0] + shift[0]) % 1.0
    y_full = (positions[1] + shift[1]) % 1.0
    full_masses = masses

    # The full-volume panel must retain every zoom mass tier.  Restrict only
    # the close-up panel to the finest tier so that coarse particles do not
    # dominate the zoom-region shot noise.
    if max_particle_mass is not None:
        selected = masses <= max_particle_mass
        x_zoom = x_full[selected]
        y_zoom = y_full[selected]
        zoom_masses = masses[selected]
    else:
        x_zoom = x_full
        y_zoom = y_full
        zoom_masses = masses

    xmin, xmax, ymin, ymax = zoom_bounds
    if assignment == "cic":
        full = deposit_cic(
            x_full,
            y_full,
            full_masses,
            full_bins,
            (0.0, 1.0, 0.0, 1.0),
            periodic=True,
        )
        zoom = deposit_cic(
            x_zoom,
            y_zoom,
            zoom_masses,
            zoom_bins,
            zoom_bounds,
            periodic=False,
        )
    else:
        full, _, _ = np.histogram2d(
            x_full,
            y_full,
            bins=full_bins,
            range=((0.0, 1.0), (0.0, 1.0)),
            weights=full_masses,
        )
        zoom, _, _ = np.histogram2d(
            x_zoom,
            y_zoom,
            bins=zoom_bins,
            range=((xmin, xmax), (ymin, ymax)),
            weights=zoom_masses,
        )
    return (
        full,
        zoom,
        float(full_masses.sum()),
        int(zoom_masses.size),
        npart,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--boxlength", type=float, default=128.0)
    parser.add_argument("--full-bins", type=int, default=512)
    parser.add_argument("--zoom-bins", type=int, default=512)
    parser.add_argument("--full-sigma", type=float, default=1.8)
    parser.add_argument("--zoom-sigma", type=float, default=1.2)
    parser.add_argument(
        "--assignment", choices=("ngp", "cic"), default="ngp"
    )
    parser.add_argument(
        "--max-particle-mass",
        type=float,
        default=None,
        help="Keep only particles at or below this RAMSES code-unit mass.",
    )
    parser.add_argument(
        "--zoom-bounds",
        nargs=4,
        type=float,
        default=(0.4259755, 0.5738515, 0.4102990, 0.5966060),
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
    )
    parser.add_argument("--shift", nargs=3, type=float, default=(0.0, 0.0, 0.0))
    parser.add_argument(
        "--marker", nargs=2, type=float, default=(0.4992375, 0.4714)
    )
    parser.add_argument("--title", required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()

    snapshot = args.snapshot.resolve()
    number = snapshot.name.removeprefix("output_")
    files = sorted(snapshot.glob(f"part_{number}.out*"))
    if not files:
        raise SystemExit(f"No particle files found in {snapshot}")
    info = snapshot / f"info_{number}.txt"
    aexp = None
    for line in info.read_text().splitlines():
        if line.split()[:1] == ["aexp"]:
            aexp = float(line.split("=")[1])
            break
    if aexp is None:
        raise SystemExit(f"No aexp in {info}")
    redshift = 1.0 / aexp - 1.0
    display_redshift = max(redshift, 0.0)

    inputs = [
        (
            str(path),
            args.full_bins,
            args.zoom_bins,
            tuple(args.zoom_bounds),
            tuple(args.shift),
            args.assignment,
            args.max_particle_mass,
        )
        for path in files
    ]
    if args.workers == 1:
        parts = [project_part(*item) for item in inputs]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers
        ) as executor:
            futures = [executor.submit(project_part, *item) for item in inputs]
            parts = [future.result() for future in futures]

    full = sum(part[0] for part in parts)
    zoom = sum(part[1] for part in parts)
    total_mass = sum(part[2] for part in parts)
    zoom_particles = sum(part[3] for part in parts)
    total_particles = sum(part[4] for part in parts)
    if zoom_particles == 0 or total_mass <= 0.0:
        raise SystemExit("The particle-mass selection is empty")
    xmin, xmax, ymin, ymax = args.zoom_bounds

    full_ratio = full / (1.0 / args.full_bins**2) / total_mass
    zoom_pixel_area = (xmax - xmin) * (ymax - ymin) / args.zoom_bins**2
    zoom_ratio = zoom / zoom_pixel_area / total_mass
    full_log = np.log10(
        np.clip(
            gaussian_filter(full_ratio, args.full_sigma, mode="wrap"),
            1.0e-12,
            None,
        )
    )
    zoom_log = np.log10(
        np.clip(
            gaussian_filter(zoom_ratio, args.zoom_sigma, mode="nearest"),
            1.0e-12,
            None,
        )
    )

    # Empty pixels are meaningful outside a selected zoom particle tier, but
    # they must not determine the contrast of the populated Lagrangian region.
    full_samples = full_log[full > 0.0][::8]
    zoom_samples = zoom_log[zoom > 0.0][::8]

    def contrast_limits(samples: np.ndarray) -> tuple[float, float]:
        vmin, vmax = np.nanpercentile(samples, (0.5, 99.7))
        if vmax - vmin < 0.08:
            midpoint = 0.5 * (vmin + vmax)
            vmin, vmax = midpoint - 0.04, midpoint + 0.04
        return float(vmin), float(vmax)

    full_vmin, full_vmax = contrast_limits(full_samples)
    zoom_vmin, zoom_vmax = contrast_limits(zoom_samples)

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "font.family": "DejaVu Sans",
        }
    )
    fig = plt.figure(figsize=(14.5, 5.7), constrained_layout=True)
    axes = fig.subplots(1, 2, gridspec_kw={"width_ratios": (1.0, 1.08)})
    box = args.boxlength
    image0 = axes[0].imshow(
        full_log.T,
        origin="lower",
        extent=(0.0, box, 0.0, box),
        cmap="magma",
        vmin=full_vmin,
        vmax=full_vmax,
        interpolation="bilinear",
        rasterized=True,
    )
    # Embed the exact high-resolution array used by the right panel into the
    # full-volume map.  This preserves a smooth all-tier base-level context
    # while making the cyan footprint pixel-for-pixel consistent with the
    # zoom panel.
    if args.max_particle_mass is not None:
        axes[0].imshow(
            zoom_log.T,
            origin="lower",
            extent=(xmin * box, xmax * box, ymin * box, ymax * box),
            cmap="magma",
            vmin=zoom_vmin,
            vmax=zoom_vmax,
            interpolation="bilinear",
            rasterized=True,
            zorder=2,
        )
        axes[0].set_xlim(0.0, box)
        axes[0].set_ylim(0.0, box)
    rectangle = Rectangle(
        (xmin * box, ymin * box),
        (xmax - xmin) * box,
        (ymax - ymin) * box,
        fill=False,
        edgecolor="#3de5ff",
        linewidth=1.8,
        zorder=3,
    )
    axes[0].add_patch(rectangle)
    axes[0].set(
        xlabel=r"$x\ [h^{-1}\,\mathrm{Mpc}]$",
        ylabel=r"$y\ [h^{-1}\,\mathrm{Mpc}]$",
        title="Full volume with embedded HR footprint",
    )
    axes[0].text(
        0.025,
        0.025,
        "DMO: projected dark-matter density",
        transform=axes[0].transAxes,
        color="white",
        fontsize=8.5,
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.55, "pad": 3},
    )

    image1 = axes[1].imshow(
        zoom_log.T,
        origin="lower",
        extent=(xmin * box, xmax * box, ymin * box, ymax * box),
        cmap="magma",
        vmin=zoom_vmin,
        vmax=zoom_vmax,
        interpolation="bilinear",
        rasterized=True,
    )
    axes[1].scatter(
        [args.marker[0] * box],
        [args.marker[1] * box],
        marker="+",
        s=90,
        linewidths=1.5,
        color="#3de5ff",
        label="expected $z=0$ halo centre",
    )
    axes[1].legend(loc="upper right", frameon=True, fontsize=8)
    axes[1].set(
        xlabel=r"$x\ [h^{-1}\,\mathrm{Mpc}]$",
        ylabel=r"$y\ [h^{-1}\,\mathrm{Mpc}]$",
        title=(
            f"Zoom footprint: {(xmax-xmin)*box:.1f}"
            rf"$\times${(ymax-ymin)*box:.1f} $h^{{-1}}$ Mpc"
        ),
    )

    for xy_full, xy_zoom in (
        ((xmax * box, ymin * box), (xmin * box, ymin * box)),
        ((xmax * box, ymax * box), (xmin * box, ymax * box)),
    ):
        fig.add_artist(
            ConnectionPatch(
                xyA=xy_full,
                coordsA=axes[0].transData,
                xyB=xy_zoom,
                coordsB=axes[1].transData,
                color="#3de5ff",
                linewidth=0.7,
                alpha=0.8,
            )
        )

    colorbar0 = fig.colorbar(image0, ax=axes[0], location="right", shrink=0.78)
    colorbar0.set_label(
        r"full: $\log_{10}(\Sigma_{\rm DM}/\overline{\Sigma}_{\rm m})$"
    )
    colorbar1 = fig.colorbar(image1, ax=axes[1], location="right", shrink=0.78)
    colorbar1.set_label(
        r"zoom: $\log_{10}(\Sigma_{\rm DM}/\overline{\Sigma}_{\rm m})$"
    )
    fig.suptitle(
        f"{args.title}  |  a={aexp:.5f}, z={display_redshift:.2f}",
        fontsize=14,
        weight="bold",
    )

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=260)
    fig.savefig(output.with_suffix(".pdf"), dpi=260)
    plt.close(fig)
    metadata = {
        "snapshot": str(snapshot),
        "aexp": aexp,
        "redshift": redshift,
        "particles": total_particles,
        "zoom_selected_particles": zoom_particles,
        "part_files": len(files),
        "assignment": args.assignment,
        "max_particle_mass": args.max_particle_mass,
        "coordinate_shift_box": args.shift,
        "zoom_bounds_box": args.zoom_bounds,
        "expected_z0_halo_center_box": args.marker,
        "full_vmin_log10": full_vmin,
        "full_vmax_log10": full_vmax,
        "zoom_vmin_log10": zoom_vmin,
        "zoom_vmax_log10": zoom_vmax,
        "full_zoom_overlay": args.max_particle_mass is not None,
        "gas_present": False,
    }
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
