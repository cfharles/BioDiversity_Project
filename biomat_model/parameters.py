"""Ecological parameters and scenario settings for the biomat restoration model."""

from __future__ import annotations

from dataclasses import dataclass


SEED = 42
GRID_SIZE = 100
CELL_SIZE_METERS = 100
CELL_AREA_HA = 1.0
SIM_YEARS = 30
DT = 1.0
PILOT_LATITUDE_RANGE = (37.5, 40.0)
PILOT_REGION_LABEL = "Araucania / Biobio, Chile"


@dataclass(frozen=True)
class ScenarioConfig:
    """Deployment schedule for one restoration scenario."""

    name: str
    slug: str
    annual_deployment: tuple[int, ...]


ECOLOGY = {
    "native": {
        # Blended native growth rate for Araucaria araucana + slower understory recovery.
        "r": 0.05,
        # Carrying capacity in normalized vegetation cover units.
        "K": 1.0,
        # Competitive pressure from invasives without protection.
        "alpha_from_invasive": 1.8,
        # Competitive pressure from invasives while the biomat blocks light.
        "alpha_from_invasive_biomat": 0.2,
        # Heavy-seeded native species disperse mainly to adjacent cells.
        "spread_radius": 1,
        "spread_rate": 0.075,
        # Mature native cells become stronger suppressors once canopy forms.
        "alpha_on_invasive_early": 0.4,
        "alpha_on_invasive_late": 0.9,
        "maturity_threshold": 0.5,
        "restoration_source_threshold": 0.6,
    },
    "invasive": {
        # Fast invasive growth representing Pinus radiata / Eucalyptus / Acacia blend.
        "r": 0.20,
        "K": 1.0,
        # Seed rain is scaled down to the 100 m grid resolution.
        "spread_radius": 2,
        "spread_rate": 0.12,
        # Seed banks and root fragments can recolonize cleared cells.
        "regrowth_probability": 0.15,
        "regrowth_increment": 0.08,
    },
    "biomat": {
        # 36-48 month design life, with degradation beginning after year 3.
        "active_years": 3,
        "fade_years": 2,
        "native_growth_boost": 2.5,
        "invasive_removal_fraction": 0.80,
        # Residual suppression while the mat is active.
        "suppression": 0.05,
    },
    "initial_conditions": {
        "invasive_dominant_fraction": 0.85,
        "native_patch_fraction": 0.05,
        "bare_fraction": 0.10,
        "invasive_cover_range": (0.70, 0.90),
        "native_patch_range": (0.10, 0.30),
        "bare_native_range": (0.00, 0.12),
        "bare_invasive_range": (0.00, 0.12),
    },
}


SCENARIOS = (
    ScenarioConfig(
        name="No Intervention",
        slug="baseline",
        annual_deployment=tuple(0 for _ in range(SIM_YEARS)),
    ),
    ScenarioConfig(
        name="Moderate Deployment",
        slug="moderate",
        annual_deployment=tuple([50] * 5 + [100] * 10 + [150] * 15),
    ),
    ScenarioConfig(
        name="Aggressive Deployment",
        slug="aggressive",
        annual_deployment=tuple([100] * 5 + [200] * 10 + [300] * 15),
    ),
)
