"""Deterministic high-density subvolume selection for P4 staging."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HighDensityRegion:
    """A periodic scout-grid cube whose mean density is globally maximal."""

    start_index: tuple[int, int, int]
    center_index: tuple[int, int, int]
    window_shape: tuple[int, int, int]
    scout_shape: tuple[int, int, int]
    mean_density: float
    peak_density: float

    @property
    def left_edge_code(self) -> tuple[float, float, float]:
        return tuple(start / size for start, size in zip(self.start_index, self.scout_shape))

    @property
    def width_code(self) -> tuple[float, float, float]:
        return tuple(width / size for width, size in zip(self.window_shape, self.scout_shape))


def select_high_density_region(
    mass_density_g_cm3: np.ndarray,
    window_shape: tuple[int, int, int],
) -> HighDensityRegion:
    """Select the periodic cube with maximum mean gas density.

    The fixed window size prevents a single unresolved density spike from
    defining the target.  On a uniform scout grid, maximizing mean density is
    identical to maximizing contained gas mass.
    """

    density = np.asarray(mass_density_g_cm3, dtype=np.float64)
    if density.ndim != 3 or any(size == 0 for size in density.shape):
        raise ValueError("mass_density_g_cm3 must be a non-empty three-dimensional array")
    if not np.isfinite(density).all() or np.any(density < 0.0):
        raise ValueError("mass_density_g_cm3 must be finite and non-negative")
    if len(window_shape) != 3 or any(width < 1 for width in window_shape):
        raise ValueError("window_shape must contain three positive lengths")
    if any(width > size for width, size in zip(window_shape, density.shape)):
        raise ValueError("window_shape must fit inside the scout grid")
    kernel = np.zeros(density.shape, dtype=np.float64)
    indices = tuple((-np.arange(width)) % size for width, size in zip(window_shape, density.shape))
    kernel[np.ix_(*indices)] = 1.0
    mean_density = np.fft.ifftn(np.fft.fftn(density) * np.fft.fftn(kernel)).real
    mean_density /= np.prod(window_shape)
    start = tuple(int(index) for index in np.unravel_index(np.argmax(mean_density), density.shape))
    center = tuple((index + width // 2) % size for index, width, size in zip(start, window_shape, density.shape))
    return HighDensityRegion(
        start_index=start,
        center_index=center,
        window_shape=window_shape,
        scout_shape=tuple(int(size) for size in density.shape),
        mean_density=float(mean_density[start]),
        peak_density=float(np.max(density)),
    )
