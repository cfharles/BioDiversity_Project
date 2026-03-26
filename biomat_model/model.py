"""Core simulation engine for the hybrid Lotka-Volterra + cellular automata restoration model."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import csv
import numpy as np

from .parameters import CELL_AREA_HA, DT, ECOLOGY, GRID_SIZE, SEED, SIM_YEARS, ScenarioConfig
from .visualization import create_gif, plot_scenario_comparison, plot_state, plot_trajectory_summary


def _neighbor_signal(field: np.ndarray, offsets: list[tuple[int, int]]) -> np.ndarray:
    """Aggregate neighboring cover using wrapped-free shifted copies and edge masking."""

    signal = np.zeros_like(field)
    rows, cols = field.shape
    for dr, dc in offsets:
        shifted = np.roll(field, shift=(dr, dc), axis=(0, 1))

        if dr > 0:
            shifted[:dr, :] = 0.0
        elif dr < 0:
            shifted[dr:, :] = 0.0

        if dc > 0:
            shifted[:, :dc] = 0.0
        elif dc < 0:
            shifted[:, dc:] = 0.0

        signal += shifted

    return signal / max(len(offsets), 1)


def _radius_offsets(radius: int) -> list[tuple[int, int]]:
    offsets: list[tuple[int, int]] = []
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            if dr == 0 and dc == 0:
                continue
            if max(abs(dr), abs(dc)) <= radius:
                offsets.append((dr, dc))
    return offsets


NATIVE_OFFSETS = _radius_offsets(ECOLOGY["native"]["spread_radius"])
INVASIVE_OFFSETS = _radius_offsets(ECOLOGY["invasive"]["spread_radius"])


def initialize_grid(seed: int = SEED, grid_size: int = GRID_SIZE) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create the initial native, invasive, and biomat-age grids."""

    rng = np.random.default_rng(seed)
    native = np.zeros((grid_size, grid_size), dtype=float)
    invasive = np.zeros((grid_size, grid_size), dtype=float)
    biomat_age = np.full((grid_size, grid_size), -1, dtype=int)

    initial = ECOLOGY["initial_conditions"]
    random_field = rng.random((grid_size, grid_size))
    invasive_mask = random_field < initial["invasive_dominant_fraction"]
    native_patch_mask = (random_field >= initial["invasive_dominant_fraction"]) & (
        random_field < initial["invasive_dominant_fraction"] + initial["native_patch_fraction"]
    )
    bare_mask = ~(invasive_mask | native_patch_mask)

    invasive[invasive_mask] = rng.uniform(*initial["invasive_cover_range"], size=invasive_mask.sum())
    native[invasive_mask] = rng.uniform(0.0, 0.05, size=invasive_mask.sum())

    native[native_patch_mask] = rng.uniform(*initial["native_patch_range"], size=native_patch_mask.sum())
    invasive[native_patch_mask] = rng.uniform(0.05, 0.20, size=native_patch_mask.sum())

    native[bare_mask] = rng.uniform(*initial["bare_native_range"], size=bare_mask.sum())
    invasive[bare_mask] = rng.uniform(*initial["bare_invasive_range"], size=bare_mask.sum())

    # Remnant forests along edges/corners act as realistic seed reservoirs.
    edge_width = 6
    edge_mask = np.zeros_like(native, dtype=bool)
    edge_mask[:edge_width, :edge_width] = True
    edge_mask[:edge_width, -edge_width:] = True
    edge_mask[-edge_width:, :edge_width] = True
    edge_mask[-edge_width:, -edge_width:] = True
    edge_mask[:4, 30:55] = True
    edge_mask[-4:, 45:70] = True

    native[edge_mask] = np.maximum(native[edge_mask], rng.uniform(0.45, 0.75, size=edge_mask.sum()))
    invasive[edge_mask] = np.minimum(invasive[edge_mask], rng.uniform(0.02, 0.20, size=edge_mask.sum()))

    return np.clip(native, 0.0, 1.0), np.clip(invasive, 0.0, 1.0), biomat_age


def deploy_biomats(
    native: np.ndarray,
    invasive: np.ndarray,
    biomat_age: np.ndarray,
    year: int,
    scenario: ScenarioConfig,
) -> int:
    """Deploy biomats according to the scenario's deployment mode.

    radial        — expand outward from each center (default behaviour).
    encirclement  — form a hollow ring at ring_radius around the center,
                    targeting cells closest to the ring circumference so the
                    biomat encircles an invasive-dominated core.
    """

    quota_per_center = scenario.annual_deployment[year - 1]
    if quota_per_center <= 0:
        return 0

    row_idx, col_idx = np.indices(native.shape)
    active = biomat_age >= 0
    priority = (~active) & ((invasive > native) | (native < ECOLOGY["native"]["maturity_threshold"]))

    total_deployed = 0
    already_chosen = np.zeros(native.shape, dtype=bool)

    for center_r, center_c in scenario.deployment_centers:
        distances = np.sqrt((row_idx - center_r) ** 2 + (col_idx - center_c) ** 2)

        if scenario.deployment_mode == "encirclement":
            # Sort by |distance - ring_radius| so cells on the ring circumference
            # are filled first, creating a hollow-circle encirclement pattern.
            ring_dist = np.abs(distances - scenario.ring_radius)
            sort_key = ring_dist
        else:
            sort_key = distances

        candidate_indices = np.flatnonzero(priority & ~already_chosen)
        if candidate_indices.size == 0:
            continue
        ordered = candidate_indices[np.argsort(sort_key.flat[candidate_indices])]
        chosen = ordered[:quota_per_center]
        already_chosen.flat[chosen] = True
        invasive.flat[chosen] *= 1.0 - ECOLOGY["biomat"]["invasive_removal_fraction"]
        native.flat[chosen] = np.maximum(native.flat[chosen], 0.12)
        biomat_age.flat[chosen] = 0
        total_deployed += int(chosen.size)

    return total_deployed


def update_dynamics(native: np.ndarray, invasive: np.ndarray, biomat_age: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Apply one Euler step of the local Lotka-Volterra competition model."""

    params_native = ECOLOGY["native"]
    params_invasive = ECOLOGY["invasive"]
    params_biomat = ECOLOGY["biomat"]

    active_mask = (biomat_age >= 0) & (biomat_age < params_biomat["active_years"])
    fading_mask = (biomat_age >= params_biomat["active_years"]) & (
        biomat_age < params_biomat["active_years"] + params_biomat["fade_years"]
    )

    fade_progress = np.zeros_like(native)
    fade_progress[fading_mask] = (
        biomat_age[fading_mask] - params_biomat["active_years"] + 1
    ) / params_biomat["fade_years"]

    alpha_ni = np.full_like(native, params_native["alpha_from_invasive"])
    alpha_ni[active_mask] = params_native["alpha_from_invasive_biomat"]
    alpha_ni[fading_mask] = (
        params_native["alpha_from_invasive_biomat"]
        + fade_progress[fading_mask]
        * (params_native["alpha_from_invasive"] - params_native["alpha_from_invasive_biomat"])
    )

    native_growth = np.full_like(native, params_native["r"])
    native_growth[active_mask] *= params_biomat["native_growth_boost"]
    native_growth[fading_mask] *= 1.0 + (params_biomat["native_growth_boost"] - 1.0) * (1.0 - fade_progress[fading_mask])

    # Natives become stronger suppressors as stands mature.
    maturity = np.clip((native - 0.2) / 0.5, 0.0, 1.0)
    alpha_in = params_native["alpha_on_invasive_early"] + maturity * (
        params_native["alpha_on_invasive_late"] - params_native["alpha_on_invasive_early"]
    )

    native_delta = native_growth * native * (1.0 - (native + alpha_ni * invasive) / params_native["K"]) * DT
    invasive_delta = (
        params_invasive["r"] * invasive * (1.0 - (invasive + alpha_in * native) / params_invasive["K"]) * DT
    )

    suppression = np.zeros_like(invasive)
    suppression[active_mask] = params_biomat["suppression"]
    suppression[fading_mask] = params_biomat["suppression"] * (1.0 - fade_progress[fading_mask])

    native = np.clip(native + native_delta, 0.0, 1.0)
    invasive = np.clip(invasive + invasive_delta - suppression, 0.0, 1.0)

    # Enforce physical constraint: N + I cannot exceed 1.0 in any cell.
    # Where they overflow, scale both down proportionally — preserving the
    # competitive ratio between species while respecting total coverage limits.
    total = native + invasive
    overflow = total > 1.0
    native[overflow] /= total[overflow]
    invasive[overflow] /= total[overflow]

    return native, invasive


def spread_seeds(native: np.ndarray, invasive: np.ndarray, biomat_age: np.ndarray, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """Spread native and invasive cover through cellular-automata style neighborhood coupling."""

    params_native = ECOLOGY["native"]
    params_invasive = ECOLOGY["invasive"]
    active_biomat = (biomat_age >= 0) & (biomat_age < ECOLOGY["biomat"]["active_years"])

    native_sources = np.where(native > 0.4, native, 0.0)
    boosted_sources = np.where(native > params_native["restoration_source_threshold"], native * 1.25, native_sources)
    invasive_sources = np.where(invasive > 0.3, invasive, 0.0)

    native_signal = _neighbor_signal(boosted_sources, NATIVE_OFFSETS)
    invasive_signal = _neighbor_signal(invasive_sources, INVASIVE_OFFSETS)

    native_random = rng.uniform(0.85, 1.15, size=native.shape)
    invasive_random = rng.uniform(0.85, 1.15, size=invasive.shape)

    native_increment = params_native["spread_rate"] * native_signal * (1.0 - native) * native_random
    invasive_increment = params_invasive["spread_rate"] * invasive_signal * (1.0 - invasive) * invasive_random
    invasive_increment[active_biomat] *= 0.25

    bare_cleared = (invasive < 0.08) & (~active_biomat)
    regrowth_draw = rng.random(invasive.shape)
    regrowth = np.where(
        bare_cleared & (regrowth_draw < params_invasive["regrowth_probability"]),
        params_invasive["regrowth_increment"],
        0.0,
    )

    native = np.clip(native + native_increment, 0.0, 1.0)
    invasive = np.clip(invasive + invasive_increment + regrowth, 0.0, 1.0)
    return native, invasive


def _year_stats(
    year: int,
    native: np.ndarray,
    invasive: np.ndarray,
    biomat_age: np.ndarray,
    cumulative_biomats: int,
) -> Dict[str, float]:
    restored_cells = int(np.count_nonzero(native > 0.5))
    biodiversity_proxy = float(np.count_nonzero(native > 0.3) / native.size)
    return {
        "year": year,
        "mean_native_cover": float(native.mean()),
        "mean_invasive_cover": float(invasive.mean()),
        "cells_restored": restored_cells,
        "biomats_deployed": cumulative_biomats,
        "restored_area_ha": restored_cells * CELL_AREA_HA,
        "biodiversity_proxy": biodiversity_proxy,
        "active_biomats": int(np.count_nonzero((biomat_age >= 0) & (biomat_age < ECOLOGY["biomat"]["active_years"]))),
    }


def _write_stats_csv(stats: List[Dict[str, float]], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stats[0].keys()))
        writer.writeheader()
        writer.writerows(stats)


def run_simulation(
    scenario: ScenarioConfig,
    base_output_dir: Path,
    make_gif: bool = True,
) -> Dict[str, object]:
    """Run one scenario and save yearly frames, CSV summaries, and a scenario-level plot."""

    rng = np.random.default_rng(SEED)
    native, invasive, biomat_age = initialize_grid(seed=SEED)

    frame_dir = base_output_dir / "output_frames" / scenario.slug
    results_dir = base_output_dir / "results"
    frame_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    stats: List[Dict[str, float]] = []
    frame_paths: List[Path] = []
    cumulative_biomats = 0

    sim_years = len(scenario.annual_deployment)
    for year in range(0, sim_years + 1):
        stats.append(_year_stats(year, native, invasive, biomat_age, cumulative_biomats))
        frame_path = frame_dir / f"{scenario.slug}_year_{year:02d}.png"
        plot_state(native, invasive, biomat_age, stats, scenario, year, frame_path)
        frame_paths.append(frame_path)

        current = stats[-1]
        print(
            f"[{scenario.name}] Year {year}/{sim_years} - "
            f"Native: {current['mean_native_cover'] * 100:5.1f}% | "
            f"Invasive: {current['mean_invasive_cover'] * 100:5.1f}% | "
            f"Restored cells: {current['cells_restored']:4d}"
        )

        if year == sim_years:
            break

        deployed = deploy_biomats(native, invasive, biomat_age, year + 1, scenario)
        cumulative_biomats += deployed
        native, invasive = update_dynamics(native, invasive, biomat_age)
        native, invasive = spread_seeds(native, invasive, biomat_age, rng)

        biomat_age = np.where(biomat_age >= 0, biomat_age + 1, biomat_age)
        biomat_age = np.where(
            biomat_age >= ECOLOGY["biomat"]["active_years"] + ECOLOGY["biomat"]["fade_years"],
            -1,
            biomat_age,
        )

    csv_path = results_dir / f"{scenario.slug}_yearly_stats.csv"
    _write_stats_csv(stats, csv_path)

    summary_path = results_dir / f"{scenario.slug}_trajectory.png"
    plot_trajectory_summary(stats, scenario, summary_path)

    gif_path = None
    if make_gif:
        gif_path = results_dir / f"{scenario.slug}_animation.gif"
        create_gif(frame_paths, gif_path)

    return {
        "scenario": scenario,
        "stats": stats,
        "frame_paths": frame_paths,
        "csv_path": csv_path,
        "summary_path": summary_path,
        "gif_path": gif_path,
        "final_native": native.copy(),
        "final_invasive": invasive.copy(),
        "final_biomat_age": biomat_age.copy(),
    }


def run_all_scenarios(base_output_dir: Path, scenarios: Tuple[ScenarioConfig, ...]) -> List[Dict[str, object]]:
    """Run the full scenario suite and create a final comparison figure."""

    results = [run_simulation(scenario, base_output_dir) for scenario in scenarios]
    comparison_path = base_output_dir / "results" / "scenario_comparison.png"
    plot_scenario_comparison(results, comparison_path)
    return results
