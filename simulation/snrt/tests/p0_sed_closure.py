"""SED-weighted group-opacity and photoelectron-closure checks."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.primordial import (
    H_I_FIT,
    _verner_cross_section_numpy,
    group_spectral_closure_from_metadata,
    sed_weighted_group_closure,
)


def main() -> None:
    edges = np.asarray((13.6, 24.59, 54.42, 500.0, 2000.0), dtype=np.float64)
    energy = np.unique(
        np.concatenate([np.geomspace(low, high, 1025) for low, high in zip(edges[:-1], edges[1:], strict=True)])
    )
    photon_spectrum = energy**-2.0
    closure = sed_weighted_group_closure(edges, energy, photon_spectrum)

    group = 3
    selected = (energy >= edges[group]) & (energy <= edges[group + 1])
    photons = np.trapezoid(photon_spectrum[selected], energy[selected])
    expected_sigma = np.trapezoid(
        photon_spectrum[selected] * _verner_cross_section_numpy(energy[selected], H_I_FIT),
        energy[selected],
    ) / photons
    assert np.isclose(closure.cross_sections.hydrogen_i[group], expected_sigma, rtol=1.0e-12, atol=0.0)
    assert closure.photoelectron_excess_energy_ev[0, group] > 0.0
    assert not np.isclose(
        closure.cross_sections.hydrogen_i[group],
        _verner_cross_section_numpy(np.asarray([closure.photon_weighted_energy_ev[group]]), H_I_FIT)[0],
        rtol=1.0e-3,
        atol=0.0,
    )

    metadata = json.loads((PROJECT_ROOT / "data/p4_pilot_agn_photon_ledger.json").read_text())
    restored = group_spectral_closure_from_metadata(metadata)
    assert restored.cross_sections.hydrogen_i.shape == (5,)
    assert restored.photoelectron_excess_energy_ev.shape == (3, 5)
    assert np.all(np.isfinite(restored.cross_sections.hydrogen_i))
    print("P0_SED_CLOSURE_OK synthetic_groups=4 metadata_groups=5 species=3")


if __name__ == "__main__":
    main()
