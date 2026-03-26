"""Ecological parameters and scenario settings for the biomat restoration model."""

from __future__ import annotations

from dataclasses import dataclass


SEED = 42
GRID_SIZE = 100
CELL_SIZE_METERS = 100
CELL_AREA_HA = 1.0
SIM_YEARS = 75
DT = 1.0
PILOT_LATITUDE_RANGE = (37.5, 40.0)
PILOT_REGION_LABEL = "Araucania / Biobio, Chile"


@dataclass(frozen=True)
class ScenarioConfig:
    """Deployment schedule for one restoration scenario."""

    name: str
    slug: str
    annual_deployment: tuple[int, ...]
    deployment_centers: tuple[tuple[int, int], ...] = ((50, 50),)
    # "radial" = expand outward from center (default)
    # "encirclement" = form a hollow ring at ring_radius to encircle invasive core
    deployment_mode: str = "radial"
    ring_radius: int = 28  # cells from center for encirclement mode


ECOLOGY = {
    "native": {
        # Blended native growth rate for Araucaria araucana + slower understory recovery.
        "r": 0.05,
        # Carrying capacity in normalized vegetation cover units.
        "K": 1.0,
        # Competitive pressure from invasives without protection.
        "alpha_from_invasive": 1.1,
        # Competitive pressure from invasives while the biomat blocks light.
        "alpha_from_invasive_biomat": 0.2,
        # Heavy-seeded native species disperse mainly to adjacent cells.
        "spread_radius": 1,
        "spread_rate": 0.075,
        # Mature native cells become stronger suppressors once canopy forms.
        "alpha_on_invasive_early": 0.3,
        "alpha_on_invasive_late": 0.85,
        "maturity_threshold": 0.5,
        "restoration_source_threshold": 0.6,
    },
    "invasive": {
        # Fast invasive growth representing Pinus radiata / Eucalyptus / Acacia blend.
        "r": 0.18,
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
        "invasive_dominant_fraction": 0.05,
        "native_patch_fraction": 0.01,
        "bare_fraction": 0.94,
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
        name="Light Deployment",
        slug="light",
        annual_deployment=tuple([25] * SIM_YEARS),
    ),
    ScenarioConfig(
        name="Moderate Deployment",
        slug="moderate",
        annual_deployment=tuple([75] * SIM_YEARS),
    ),
    ScenarioConfig(
        name="Aggressive Deployment",
        slug="aggressive",
        annual_deployment=tuple([125] * SIM_YEARS),
    ),
    ScenarioConfig(
        name="Multi-Point Deployment",
        slug="multipoint",
        annual_deployment=tuple([25] * SIM_YEARS),
        deployment_centers=(
            (50, 50),   # center
            (25, 75),   # top-left
            (75, 75),   # top-right
            (25, 25),   # bottom-left
            (75, 25),   # bottom-right
        ),
    ),
    ScenarioConfig(
        name="Encirclement Ring",
        slug="encirclement",
        # 75 ha/year deployed in a hollow circle to surround an invasive core
        annual_deployment=tuple([75] * SIM_YEARS),
        deployment_centers=((50, 50),),  # ring centered on the grid
        deployment_mode="encirclement",
        ring_radius=28,  # ~28-cell radius circle ≈ 175-ha encircled area
    ),
)

