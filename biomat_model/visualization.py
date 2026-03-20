"""Plotting helpers for yearly biomat simulation maps and summary figures."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

try:
    import imageio.v2 as imageio
except ModuleNotFoundError:  # pragma: no cover - optional export dependency
    imageio = None
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from .parameters import GRID_SIZE, PILOT_LATITUDE_RANGE, PILOT_REGION_LABEL, SIM_YEARS, ScenarioConfig


def _state_rgb(native: np.ndarray, invasive: np.ndarray) -> np.ndarray:
    """Build a custom RGB map that emphasizes ecological state instead of raw values."""

    rgb = np.zeros(native.shape + (3,), dtype=float)
    bare = (native < 0.2) & (invasive < 0.2)
    healthy_native = native > 0.7
    recovering_native = (native > invasive) & ~healthy_native
    invasive_dom = invasive > native

    rgb[bare] = np.array([0.55, 0.42, 0.24])
    rgb[healthy_native] = np.array([0.07, 0.42, 0.16])
    rgb[recovering_native] = np.array([0.67, 0.80, 0.30])
    rgb[invasive_dom] = np.stack(
        [
            np.full_like(invasive[invasive_dom], 0.75),
            0.22 + 0.35 * (1.0 - invasive[invasive_dom]),
            0.12 + 0.12 * (1.0 - invasive[invasive_dom]),
        ],
        axis=-1,
    )
    return np.clip(rgb, 0.0, 1.0)


def _draw_chile_inset(ax: plt.Axes) -> None:
    chile_outline = np.array(
        [
            [0.56, 0.96],
            [0.51, 0.90],
            [0.54, 0.84],
            [0.47, 0.77],
            [0.51, 0.69],
            [0.44, 0.59],
            [0.49, 0.51],
            [0.40, 0.42],
            [0.46, 0.31],
            [0.37, 0.20],
            [0.42, 0.08],
            [0.60, 0.10],
            [0.57, 0.22],
            [0.64, 0.35],
            [0.59, 0.49],
            [0.66, 0.62],
            [0.60, 0.77],
            [0.64, 0.90],
        ]
    )
    polygon = Polygon(chile_outline, closed=True, facecolor="#d9ead3", edgecolor="#284b63", linewidth=1.2)
    ax.add_patch(polygon)

    pilot_y = 0.60 - (np.mean(PILOT_LATITUDE_RANGE) - 37.5) * 0.08
    ax.scatter([0.56], [pilot_y], s=28, color="#c44536", zorder=5)
    ax.text(0.08, 0.02, PILOT_REGION_LABEL, fontsize=8, color="#203040")
    ax.set_xlim(0.25, 0.75)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Chile Pilot Site", fontsize=9)
    ax.axis("off")


def plot_state(
    native: np.ndarray,
    invasive: np.ndarray,
    biomat_age: np.ndarray,
    stats: List[Dict[str, float]],
    scenario: ScenarioConfig,
    year: int,
    output_path: Path,
) -> None:
    """Create one yearly figure with a spatial state map and three time-series panels."""

    years = [row["year"] for row in stats]
    native_series = [row["mean_native_cover"] * 100.0 for row in stats]
    invasive_series = [row["mean_invasive_cover"] * 100.0 for row in stats]
    restored_series = [row["cells_restored"] for row in stats]
    biomat_series = [row["biomats_deployed"] for row in stats]
    biodiversity_series = [row["biodiversity_proxy"] * 100.0 for row in stats]

    fig = plt.figure(figsize=(15, 9), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, width_ratios=[1.35, 1.0])

    ax_map = fig.add_subplot(gs[:, 0])
    rgb = _state_rgb(native, invasive)
    ax_map.imshow(rgb, origin="lower")
    active_biomat = (biomat_age >= 0) & (biomat_age < 3)
    if np.any(active_biomat):
        ax_map.contour(active_biomat.astype(float), levels=[0.5], colors=["#3f88c5"], linewidths=0.7, origin="lower")

    ax_map.set_title(f"{scenario.name} Spatial State")
    ax_map.set_xlabel("Grid Easting (100 m cells)")
    ax_map.set_ylabel("Grid Northing (100 m cells)")
    ax_map.set_xlim(-0.5, GRID_SIZE - 0.5)
    ax_map.set_ylim(-0.5, GRID_SIZE - 0.5)

    inset = ax_map.inset_axes([0.76, 0.04, 0.20, 0.32])
    _draw_chile_inset(inset)

    ax_cover = fig.add_subplot(gs[0, 1])
    ax_cover.plot(years, native_series, color="#1b7f3b", linewidth=2.5, label="Native cover")
    ax_cover.plot(years, invasive_series, color="#c44536", linewidth=2.5, label="Invasive cover")
    ax_cover.set_ylabel("Cover (%)")
    ax_cover.set_xlim(0, SIM_YEARS)
    ax_cover.set_ylim(0, 100)
    ax_cover.grid(alpha=0.25)
    ax_cover.legend(loc="upper right", frameon=False)

    ax_restoration = fig.add_subplot(gs[1, 1])
    ax_restoration.plot(years, biomat_series, color="#3f88c5", linewidth=2.2, label="Cumulative biomats")
    ax_restoration.plot(years, restored_series, color="#8fba52", linewidth=2.2, label="Restored cells")
    ax_restoration.set_ylabel("Cells")
    ax_restoration.set_xlim(0, SIM_YEARS)
    ax_restoration.grid(alpha=0.25)
    ax_restoration.legend(loc="upper left", frameon=False)

    ax_biodiversity = fig.add_subplot(gs[2, 1])
    ax_biodiversity.plot(years, biodiversity_series, color="#7a5cfa", linewidth=2.2)
    ax_biodiversity.fill_between(years, biodiversity_series, color="#7a5cfa", alpha=0.15)
    ax_biodiversity.set_ylabel("Cells >30% native (%)")
    ax_biodiversity.set_xlabel("Simulation year")
    ax_biodiversity.set_xlim(0, SIM_YEARS)
    ax_biodiversity.set_ylim(0, 100)
    ax_biodiversity.grid(alpha=0.25)

    fig.suptitle(f"Biomat Restoration Model - Year {year} / {SIM_YEARS}", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def plot_trajectory_summary(stats: List[Dict[str, float]], scenario: ScenarioConfig, output_path: Path) -> None:
    """Create a compact scenario-level summary chart for the full 30-year trajectory."""

    years = [row["year"] for row in stats]
    native_series = [row["mean_native_cover"] * 100.0 for row in stats]
    invasive_series = [row["mean_invasive_cover"] * 100.0 for row in stats]
    restored_series = [row["restored_area_ha"] for row in stats]
    biodiversity_series = [row["biodiversity_proxy"] * 100.0 for row in stats]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    axes[0].plot(years, native_series, color="#1b7f3b", linewidth=2.5)
    axes[0].plot(years, invasive_series, color="#c44536", linewidth=2.5)
    axes[0].set_title("Vegetation Cover")
    axes[0].set_ylabel("Cover (%)")
    axes[0].set_xlabel("Year")
    axes[0].grid(alpha=0.25)

    axes[1].plot(years, restored_series, color="#8fba52", linewidth=2.5)
    axes[1].set_title("Restored Area")
    axes[1].set_ylabel("Hectares")
    axes[1].set_xlabel("Year")
    axes[1].grid(alpha=0.25)

    axes[2].plot(years, biodiversity_series, color="#7a5cfa", linewidth=2.5)
    axes[2].fill_between(years, biodiversity_series, color="#7a5cfa", alpha=0.15)
    axes[2].set_title("Biodiversity Proxy")
    axes[2].set_ylabel("% Cells >30% Native")
    axes[2].set_xlabel("Year")
    axes[2].grid(alpha=0.25)

    fig.suptitle(f"{scenario.name} - 30 Year Restoration Trajectory", fontsize=15, fontweight="bold")
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def plot_scenario_comparison(results: List[Dict[str, object]], output_path: Path) -> None:
    """Compare the three scenarios with year-30 maps and native-cover trajectories."""

    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 0.85])

    for idx, result in enumerate(results):
        scenario = result["scenario"]
        ax = fig.add_subplot(gs[0, idx])
        ax.imshow(_state_rgb(result["final_native"], result["final_invasive"]), origin="lower")
        active_biomat = (result["final_biomat_age"] >= 0) & (result["final_biomat_age"] < 3)
        if np.any(active_biomat):
            ax.contour(active_biomat.astype(float), levels=[0.5], colors=["#3f88c5"], linewidths=0.7, origin="lower")
        ax.set_title(f"{scenario.name}\nYear 30")
        ax.set_xticks([])
        ax.set_yticks([])

    ax_lines = fig.add_subplot(gs[1, :])
    colors = ["#c44536", "#f0a202", "#1b7f3b"]
    for color, result in zip(colors, results):
        stats = result["stats"]
        years = [row["year"] for row in stats]
        native_series = [row["mean_native_cover"] * 100.0 for row in stats]
        ax_lines.plot(years, native_series, linewidth=2.8, color=color, label=result["scenario"].name)

    ax_lines.set_title("Native Cover Trajectory Comparison")
    ax_lines.set_xlabel("Year")
    ax_lines.set_ylabel("Mean native cover (%)")
    ax_lines.set_xlim(0, SIM_YEARS)
    ax_lines.set_ylim(0, 100)
    ax_lines.grid(alpha=0.25)
    ax_lines.legend(frameon=False, ncol=3, loc="upper center")

    fig.suptitle("Biomat Restoration Scenarios", fontsize=17, fontweight="bold")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def create_gif(frame_paths: List[Path], output_path: Path, duration: float = 0.7) -> None:
    """Build an animated GIF from yearly PNG frames."""

    if imageio is None:
        return

    images = [imageio.imread(path) for path in frame_paths]
    imageio.mimsave(output_path, images, duration=duration, loop=0)
