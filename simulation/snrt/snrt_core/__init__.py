"""Static-grid multigroup S_N transport primitives for TPU execution."""

from .quadrature import s4_quadrature
from .chemistry import HydrogenUpdate, advance_hydrogen, hydrogen_absorption
from .coupling import build_hydrogen_radiation_step
from .benchmarks import StromgrenProblem, StromgrenState, build_stromgren_runner, make_stromgren_problem
from .ledger import PhotonLedger, energy_from_photons, photon_ledger
from .sources import PointSources, deposit_point_sources
from .transport import (
    TransportConfig,
    advance_explicit,
    build_explicit_step,
    cfl_number,
    initial_intensity,
    radiation_moments,
)

__all__ = [
    "PointSources",
    "TransportConfig",
    "HydrogenUpdate",
    "PhotonLedger",
    "StromgrenProblem",
    "StromgrenState",
    "advance_explicit",
    "advance_hydrogen",
    "build_explicit_step",
    "build_hydrogen_radiation_step",
    "build_stromgren_runner",
    "cfl_number",
    "deposit_point_sources",
    "energy_from_photons",
    "hydrogen_absorption",
    "initial_intensity",
    "make_stromgren_problem",
    "photon_ledger",
    "radiation_moments",
    "s4_quadrature",
]
