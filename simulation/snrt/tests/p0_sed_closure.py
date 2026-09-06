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
    reference_energy = np.geomspace(edges[group], edges[group + 1], 131073)
    reference_spectrum = np.interp(reference_energy, energy, photon_spectrum)
    photons = np.trapezoid(reference_spectrum, reference_energy)
    expected_sigma = np.trapezoid(
        reference_spectrum * _verner_cross_section_numpy(reference_energy, H_I_FIT),
        reference_energy,
    ) / photons
    assert np.isclose(closure.cross_sections.hydrogen_i[group], expected_sigma, rtol=1.0e-6, atol=0.0)
    assert closure.photoelectron_excess_energy_ev[0, group] > 0.0
    assert not np.isclose(
        closure.cross_sections.hydrogen_i[group],
        _verner_cross_section_numpy(np.asarray([closure.photon_weighted_energy_ev[group]]), H_I_FIT)[0],
        rtol=1.0e-3,
        atol=0.0,
    )

    threshold_edges = np.asarray((11.2, 13.6, 24.59, 54.42, 100.0))
    threshold_energy = np.unique(
        np.concatenate(
            [
                np.geomspace(low, high, 257)
                for low, high in zip(
                    threshold_edges[:-1], threshold_edges[1:], strict=True
                )
            ]
        )
    )
    threshold_closure = sed_weighted_group_closure(
        threshold_edges, threshold_energy, threshold_energy**-2.0
    )
    assert threshold_closure.cross_sections.hydrogen_i[0] == 0.0
    assert np.all(threshold_closure.cross_sections.helium_i[:2] == 0.0)
    assert np.all(threshold_closure.cross_sections.helium_ii[:3] == 0.0)

    metadata = json.loads(
        (PROJECT_ROOT / "data/p4_pilot_agn_photon_ledger.json").read_text()
    )
    restored = group_spectral_closure_from_metadata(metadata)
    assert restored.cross_sections.hydrogen_i.shape == (9,)
    assert restored.photoelectron_excess_energy_ev.shape == (3, 9)
    assert np.all(np.isfinite(restored.cross_sections.hydrogen_i))
    print(
        "P0_SED_CLOSURE_OK synthetic_groups=4 metadata_groups=9 species=3 "
        "subthreshold_opacity=zero"
    )


if __name__ == "__main__":
    main()
