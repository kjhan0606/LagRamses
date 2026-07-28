#!/usr/bin/env python3
"""Plot matched base-grid and leaf-AMR density projections."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle
from PIL import Image


def read_meta(path: Path) -> dict[str, str]:
    meta: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            meta[key.strip()] = value.strip()
    return meta


def numbers(meta: dict[str, str], key: str) -> list[float]:
    return [float(value) for value in meta[key].split()]


def add_distance_bar(
    axis: plt.Axes,
    length: float,
    label: str,
    side: str = "right",
) -> None:
    """Place a physical-distance bar inside an image panel."""
    xmin, xmax = axis.get_xlim()
    ymin, ymax = axis.get_ylim()
    width = xmax - xmin
    height = ymax - ymin
    if side == "left":
        x0 = xmin + 0.055 * width
        x1 = x0 + length
    else:
        x1 = xmax - 0.055 * width
        x0 = x1 - length
    y = ymin + 0.070 * height
    halo = [
        path_effects.Stroke(linewidth=5.2, foreground="black"),
        path_effects.Normal(),
    ]
    bar = axis.plot(
        [x0, x1],
        [y, y],
        color="white",
        linewidth=3.2,
        solid_capstyle="butt",
        zorder=10,
    )[0]
    bar.set_path_effects(halo)
    for x in (x0, x1):
        endcap = axis.plot(
            [x, x],
            [y - 0.012 * height, y + 0.012 * height],
            color="white",
            linewidth=2.2,
            solid_capstyle="butt",
            zorder=10,
        )[0]
        endcap.set_path_effects(halo)
    text = axis.text(
        0.5 * (x0 + x1),
        y + 0.026 * height,
        label,
        color="white",
        fontsize=12,
        ha="center",
        va="bottom",
        zorder=10,
    )
    text.set_path_effects(
        [
            path_effects.Stroke(linewidth=3.0, foreground="black"),
            path_effects.Normal(),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument(
        "--contrast",
        choices=("fixed", "auto"),
        default="fixed",
        help="Use the historical fixed display range or percentile limits.",
    )
    args = parser.parse_args()

    prefix = args.prefix
    output = args.output or prefix.with_name(prefix.name + "_matched_amr_projection.png")
    meta = read_meta(prefix.with_name(prefix.name + "_meta.txt"))
    base_shape = [int(value) for value in meta["base_shape"].split()]
    leaf_shape = [int(value) for value in meta["leaf_shape"].split()]
    lbox = float(meta["lbox_mpc_h"])
    zmin, zmax = numbers(meta, "slab_z_mpc_h")
    xmin, xmax, ymin, ymax = numbers(meta, "zoom_xy_mpc_h")
    mean_column = float(meta["mean_column_code"])
    redshift = float(meta["redshift"])
    aexp = float(meta["aexp"])
    levelmin = int(meta["levelmin"])
    levelmax = int(meta["levelmax"])

    base = np.fromfile(prefix.with_name(prefix.name + "_base.bin"), dtype="<f8")
    leaf = np.fromfile(prefix.with_name(prefix.name + "_leaf.bin"), dtype="<f8")
    base = base.reshape(tuple(base_shape), order="F").T
    leaf = leaf.reshape(tuple(leaf_shape), order="F").T
    base_log = np.log10(np.maximum(base / mean_column, 1.0e-12))
    leaf_log = np.log10(np.maximum(leaf / mean_column, 1.0e-12))

    zoom_aspect = (xmax - xmin) / (ymax - ymin)
    fig = plt.figure(figsize=(10.35, 6.0))
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.0, zoom_aspect),
        wspace=0.012,
    )
    ax_full = fig.add_subplot(grid[0, 0])
    ax_zoom = fig.add_subplot(grid[0, 1])
    fig.subplots_adjust(
        left=0.003,
        right=0.997,
        bottom=0.003,
        top=0.997,
    )

    if args.contrast == "auto":
        samples = np.concatenate(
            (
                base_log[np.isfinite(base_log)].ravel()[::16],
                leaf_log[np.isfinite(leaf_log)].ravel()[::16],
            )
        )
        vmin, vmax = np.nanpercentile(samples, (0.5, 99.6))
        if vmax - vmin < 0.25:
            midpoint = 0.5 * (vmin + vmax)
            vmin, vmax = midpoint - 0.125, midpoint + 0.125
    else:
        vmin, vmax = -0.42, 1.55
    image_full = ax_full.imshow(
        base_log,
        origin="lower",
        extent=(0.0, lbox, 0.0, lbox),
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax_zoom.imshow(
        leaf_log,
        origin="lower",
        extent=(xmin, xmax, ymin, ymax),
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )

    cyan = "#43d8f2"
    rect = Rectangle(
        (xmin, ymin),
        xmax - xmin,
        ymax - ymin,
        fill=False,
        edgecolor=cyan,
        linewidth=1.6,
    )
    ax_full.add_patch(rect)
    for y0, y1 in ((ymax, ymax), (ymin, ymin)):
        fig.add_artist(
            ConnectionPatch(
                xyA=(xmax, y0),
                coordsA=ax_full.transData,
                xyB=(xmin, y1),
                coordsB=ax_zoom.transData,
                color=cyan,
                linewidth=0.8,
                alpha=0.9,
            )
        )

    for axis in (ax_full, ax_zoom):
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_xlabel("")
        axis.set_ylabel("")
        axis.tick_params(
            axis="both",
            which="both",
            bottom=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )
        for spine in axis.spines.values():
            spine.set_visible(False)
    ax_full.set_xlim(0.0, lbox)
    ax_full.set_ylim(0.0, lbox)
    ax_zoom.set_xlim(xmin, xmax)
    ax_zoom.set_ylim(ymin, ymax)
    add_distance_bar(
        ax_full,
        20.0,
        r"$20\ h^{-1}\,\mathrm{Mpc}$",
        side="left",
    )
    add_distance_bar(
        ax_zoom,
        5.0,
        r"$5\ h^{-1}\,\mathrm{Mpc}$",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output,
        dpi=args.dpi,
        facecolor="white",
        transparent=False,
    )
    with Image.open(output) as image:
        image.convert("RGB").save(output, dpi=(args.dpi, args.dpi))
    print(output.resolve())


if __name__ == "__main__":
    main()
