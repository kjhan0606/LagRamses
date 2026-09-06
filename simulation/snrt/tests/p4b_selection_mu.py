"""P4-B high-density selection and primordial mu-table checks.

Run with: JAX_PLATFORMS=cpu .venv/bin/python tests/p4b_selection_mu.py
"""

from __future__ import annotations

import numpy as np

from snrt_core.mu_table import neutral_primordial_mu_table, primordial_mean_molecular_weight
from snrt_core.selection import select_high_density_region


def main() -> None:
    density = np.ones((16, 16, 16))
    density[5:9, 6:10, 7:11] = 100.0
    region = select_high_density_region(density, (4, 4, 4))
    assert region.start_index == (5, 6, 7)
    assert region.center_index == (7, 8, 9)
    assert region.mean_density == 100.0

    neutral_mu = float(primordial_mean_molecular_weight(0.0, 0.0, 0.0))
    ionized_mu = float(primordial_mean_molecular_weight(1.0, 0.0, 1.0))
    assert np.isclose(neutral_mu, 1.2195121951219512)
    assert np.isclose(ionized_mu, 0.5882352941176471)
    table = neutral_primordial_mu_table(np.array([1.0, 8.0]), np.array([-8.0, 4.0]))
    temperature = table.temperature_from_pressure(
        pressure_dyn_cm2=np.array([1.380649e-12]),
        mass_density_g_cm3=np.array([1.67262192369e-24]),
        n_h_cm3=np.array([1.0]),
    )
    assert np.isclose(temperature[0], neutral_mu * 1.0e4)
    print(f"P4B_SELECTION_MU_OK center={region.center_index} neutral_mu={neutral_mu:.6f}")


if __name__ == "__main__":
    main()
