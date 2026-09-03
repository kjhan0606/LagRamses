"""Small provenance helpers shared by review-only F-P2 report builders."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]

PROMOTION_REQUIRED_FIELDS = (
    "selected_model_or_population_mixture",
    "minimum_delay_gyr",
    "maximum_delay_gyr",
    "events_per_initial_msun_and_imf_conversion",
    "decay_convention_and_horizon_yr",
    "isotope_to_project_element_policy",
    "returned_mass_msun_per_event",
    "terminal_remnant_msun_per_event",
    "wd_remnant_debit_policy",
    "energy_erg_per_event",
    "momentum_vector_convention",
    "momentum_deposition_policy",
    "event_realization_policy",
    "snia_thermal_coupling",
    "metallicity_dependence",
    "source_warning_quarantine_policy",
    "portable_artifact_paths",
    "source_commit_binding",
    "named_approval_id",
)


def project_relative(path: Path) -> str:
    """Return a stable repository-relative display path when possible."""

    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(resolved)
