#!/usr/bin/env python3
"""
Interactive simulation of post-fire / post-invasive restoration in the
Chilean Winter Rainfall-Valdivian Forest region using biomat technology.

Run locally:
    streamlit run forest_restoration_app.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

MISSING_PACKAGES: list[str] = []

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None  # type: ignore[assignment]
    MISSING_PACKAGES.append("pandas")

try:
    import plotly.graph_objects as go
except ModuleNotFoundError:
    go = None  # type: ignore[assignment]
    MISSING_PACKAGES.append("plotly")

try:
    import streamlit as st
except ModuleNotFoundError:
    st = None  # type: ignore[assignment]
    MISSING_PACKAGES.append("streamlit")


# -----------------------------------------------------------------------------
# Optional GIS imports (real-data-ready mode)
# -----------------------------------------------------------------------------
try:
    import geopandas as gpd  # type: ignore
except Exception:
    gpd = None


# -----------------------------------------------------------------------------
# Placeholder paths for future real Chile data integration
# -----------------------------------------------------------------------------
STUDY_AREA_PATH = Path("study_area.geojson")
BURN_SCARS_PATH = Path("burn_scars.geojson")
INVASIVE_OCCURRENCES_PATH = Path("invasive_occurrences.csv")
INVASIVE_SUITABILITY_PATH = Path("invasive_suitability.tif")


def ensure_runtime_dependencies() -> None:
    """Fail early with a clear setup message if UI dependencies are missing."""
    if not MISSING_PACKAGES:
        return

    pkgs = " ".join(sorted(MISSING_PACKAGES))
    raise SystemExit(
        "Missing required package(s): "
        f"{', '.join(sorted(MISSING_PACKAGES))}\n"
        "Install dependencies with:\n"
        f"  python3 -m pip install {pkgs}\n"
        "Then run the app with:\n"
        "  streamlit run forest_restoration_app.py"
    )


@dataclass
class ModelParams:
    """Ecological model parameters."""

    r: float = 0.08  # natural regrowth rate
    alpha: float = 0.20  # biomat-assisted establishment
    beta: float = 0.12  # invasive competition pressure on natives
    m: float = 0.03  # mortality / post-fire stress
    g: float = 0.10  # invasive regrowth
    c: float = 0.55  # clearing effectiveness at treatment year
    d: float = 0.10  # native suppression of invasives
    native_spread: float = 0.05  # neighborhood native spread coefficient
    invasive_spread: float = 0.06  # neighborhood invasive recolonization
    drought_prob: float = 0.0  # optional stochastic drought
    fire_recurrence_prob: float = 0.0  # optional fire recurrence


# -----------------------------------------------------------------------------
# Core spatial functions
# -----------------------------------------------------------------------------
def hex_cell_geometry(side_length_m: float = 0.10) -> tuple[float, float]:
    """Return (area_per_hex_cell_m2, cells_per_m2)."""
    area = (3.0 * math.sqrt(3.0) / 2.0) * side_length_m**2
    return area, 1.0 / area


def neighborhood_mean(arr: np.ndarray) -> np.ndarray:
    """Moore-neighborhood mean using 8 neighbors with wrap-free approximation.

    We use edge-padded rolls and average over 8 surrounding cells.
    """
    padded = np.pad(arr, pad_width=1, mode="edge")
    neighbors_sum = (
        padded[:-2, :-2]
        + padded[:-2, 1:-1]
        + padded[:-2, 2:]
        + padded[1:-1, :-2]
        + padded[1:-1, 2:]
        + padded[2:, :-2]
        + padded[2:, 1:-1]
        + padded[2:, 2:]
    )
    return neighbors_sum / 8.0


def update_cell(
    n: float,
    i: float,
    b: float,
    n_neigh: float,
    i_neigh: float,
    treatment_pulse: float,
    params: ModelParams,
    drought_multiplier: float,
) -> tuple[float, float]:
    """Single-cell update for native and invasive cover.

    Kept explicit for readability and future rule customization.
    """
    n_next = (
        n
        + params.r * drought_multiplier * n * (1.0 - n)
        + params.alpha * b * (1.0 - n)
        - params.beta * i * n
        - params.m * n
        + params.native_spread * n_neigh * (1.0 - n)
    )

    i_next = (
        i
        + params.g * i * (1.0 - i)
        - params.c * treatment_pulse
        - params.d * n * i
        + params.invasive_spread * i_neigh * (1.0 - i) * (1.0 - b)
    )

    return float(np.clip(n_next, 0.0, 1.0)), float(np.clip(i_next, 0.0, 1.0))


def derive_state(
    native_cover: np.ndarray,
    invasive_cover: np.ndarray,
    ready_mask: np.ndarray,
    biomat_mask: np.ndarray,
) -> np.ndarray:
    """Convert continuous variables into categorical state codes.

    States:
      0 degraded/empty
      1 invasive-dominated
      2 burned-cleared and ready
      3 biomat installed
      4 recovering native
      5 mature native
    """
    state = np.zeros(native_cover.shape, dtype=np.int8)

    state[invasive_cover >= 0.6] = 1
    state[ready_mask] = 2
    state[biomat_mask] = 3
    state[native_cover >= 0.30] = 4
    state[native_cover >= 0.70] = 5

    return state


def load_or_generate_landscape(
    mode: str = "demo",
    grid_shape: tuple[int, int] = (60, 60),
    random_seed: int = 8,
) -> dict[str, np.ndarray]:
    """Create synthetic landscape (demo) or fallback in real-data-ready mode."""
    rng = np.random.default_rng(random_seed)

    if mode == "real":
        # Real-data-ready path: if files and geopandas are available, this is where
        # loading/rasterization should occur. For now we gracefully fallback.
        if gpd is not None and STUDY_AREA_PATH.exists():
            # Placeholder for future rasterization workflow.
            # Keeping deterministic synthetic baseline to ensure the app always runs.
            pass

    rows, cols = grid_shape
    y, x = np.mgrid[0:rows, 0:cols]

    # Build smooth synthetic invasion hotspots and burned/cleared patches.
    hotspot1 = np.exp(-(((x - 18) ** 2 + (y - 22) ** 2) / (2 * 8.5**2)))
    hotspot2 = np.exp(-(((x - 42) ** 2 + (y - 38) ** 2) / (2 * 7.0**2)))
    invasive = np.clip(0.18 + 0.72 * (0.6 * hotspot1 + 0.7 * hotspot2), 0.0, 1.0)

    burn_center = np.exp(-(((x - 33) ** 2 + (y - 18) ** 2) / (2 * 6.0**2)))
    burned = burn_center > 0.35

    native = np.clip(0.25 - 0.20 * invasive + 0.05 * rng.random((rows, cols)), 0.0, 1.0)
    native[burned] *= 0.12
    invasive[burned] *= 0.20

    biomat = np.zeros((rows, cols), dtype=bool)
    ready = burned.copy()  # burned/cleared and ready for restoration

    return {
        "native_cover": native,
        "invasive_cover": invasive,
        "biomat_present": biomat,
        "ready_mask": ready,
        "years_since_treatment": np.zeros((rows, cols), dtype=np.int16),
    }


def select_treatment_cells(
    invasive_cover: np.ndarray,
    ready_mask: np.ndarray,
    fraction: float,
    strategy: str,
    seed: int,
) -> np.ndarray:
    """Select cells to treat based on strategy and requested fraction."""
    rng = np.random.default_rng(seed)
    n_cells = invasive_cover.size
    k = max(1, int(fraction * n_cells))

    if strategy == "Burned/ready first":
        score = ready_mask.astype(float) * 10.0 + invasive_cover
    elif strategy == "Highest invasion first":
        score = invasive_cover.copy()
    else:  # Mixed priority
        score = 0.65 * invasive_cover + 0.35 * ready_mask.astype(float)

    flat_idx = np.argpartition(score.ravel(), -k)[-k:]
    mask = np.zeros_like(invasive_cover, dtype=bool)
    mask.ravel()[flat_idx] = True

    # Add slight randomness so neighboring ties do not look overly rigid.
    jitter = rng.random(invasive_cover.shape) < 0.02
    return np.logical_or(mask, np.logical_and(jitter, score > np.percentile(score, 75)))


def apply_treatment(
    landscape: dict[str, np.ndarray],
    treatment_mask: np.ndarray,
    remove_invasives: bool,
    install_biomat: bool,
) -> dict[str, np.ndarray]:
    """Apply user restoration actions to selected cells."""
    out = {k: v.copy() for k, v in landscape.items()}

    if remove_invasives:
        out["invasive_cover"][treatment_mask] *= 0.25
        out["ready_mask"][treatment_mask] = True

    if install_biomat:
        out["biomat_present"][treatment_mask] = True
        out["ready_mask"][treatment_mask] = False
        out["years_since_treatment"][treatment_mask] = 0

    return out


def update_landscape(
    native_cover: np.ndarray,
    invasive_cover: np.ndarray,
    biomat_present: np.ndarray,
    treatment_mask: np.ndarray,
    year: int,
    params: ModelParams,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Update all cells for one year using neighborhood-coupled dynamics."""
    n_neigh = neighborhood_mean(native_cover)
    i_neigh = neighborhood_mean(invasive_cover)

    # Drought years reduce natural regrowth.
    drought_multiplier = 0.65 if rng.random() < params.drought_prob else 1.0

    # Optional recurrent fire: pushes some cells back to degraded condition.
    fire_mask = rng.random(native_cover.shape) < params.fire_recurrence_prob

    b = biomat_present.astype(float)
    treatment_pulse = treatment_mask.astype(float) if year == 0 else np.zeros_like(b)

    # Vectorized update from the same equation implemented in update_cell.
    n_next = (
        native_cover
        + params.r * drought_multiplier * native_cover * (1.0 - native_cover)
        + params.alpha * b * (1.0 - native_cover)
        - params.beta * invasive_cover * native_cover
        - params.m * native_cover
        + params.native_spread * n_neigh * (1.0 - native_cover)
    )

    i_next = (
        invasive_cover
        + params.g * invasive_cover * (1.0 - invasive_cover)
        - params.c * treatment_pulse
        - params.d * native_cover * invasive_cover
        + params.invasive_spread * i_neigh * (1.0 - invasive_cover) * (1.0 - b)
    )

    n_next = np.clip(n_next, 0.0, 1.0)
    i_next = np.clip(i_next, 0.0, 1.0)

    if np.any(fire_mask):
        n_next[fire_mask] *= 0.10
        i_next[fire_mask] *= 0.50

    return n_next, i_next


def run_simulation(
    landscape: dict[str, np.ndarray],
    treatment_mask: np.ndarray,
    years: int,
    params: ModelParams,
    seed: int,
) -> dict[str, Any]:
    """Run yearly simulation and return full trajectories."""
    rng = np.random.default_rng(seed)

    native = landscape["native_cover"].copy()
    invasive = landscape["invasive_cover"].copy()
    biomat = landscape["biomat_present"].copy()
    ready = landscape["ready_mask"].copy()

    native_ts = [native.copy()]
    invasive_ts = [invasive.copy()]
    state_ts = [derive_state(native, invasive, ready, biomat)]

    years_since_treatment = landscape["years_since_treatment"].copy()

    for yr in range(years):
        native, invasive = update_landscape(
            native_cover=native,
            invasive_cover=invasive,
            biomat_present=biomat,
            treatment_mask=treatment_mask,
            year=yr,
            params=params,
            rng=rng,
        )

        treated_cells = np.logical_or(treatment_mask, biomat)
        years_since_treatment[treated_cells] += 1

        # Cells with biomat are no longer in "ready" state.
        ready = np.logical_and(ready, ~biomat)

        native_ts.append(native.copy())
        invasive_ts.append(invasive.copy())
        state_ts.append(derive_state(native, invasive, ready, biomat))

    return {
        "native_ts": np.stack(native_ts),
        "invasive_ts": np.stack(invasive_ts),
        "state_ts": np.stack(state_ts),
        "biomat_mask": biomat,
        "years_since_treatment": years_since_treatment,
    }


def compute_metrics(
    native_ts: np.ndarray,
    invasive_ts: np.ndarray,
    biomat_mask: np.ndarray,
    grid_cell_area_m2: float,
) -> pd.DataFrame:
    """Compute yearly summary metrics."""
    years = native_ts.shape[0] - 1
    idx = np.arange(years + 1)

    pct_native_gt_07 = (native_ts > 0.70).mean(axis=(1, 2)) * 100.0
    pct_invasive_dom = (invasive_ts > 0.60).mean(axis=(1, 2)) * 100.0
    mean_native = native_ts.mean(axis=(1, 2))
    mean_invasive = invasive_ts.mean(axis=(1, 2))
    total_treated_area = biomat_mask.sum() * grid_cell_area_m2

    df = pd.DataFrame(
        {
            "year": idx,
            "pct_native_gt_07": pct_native_gt_07,
            "pct_invasive_dom": pct_invasive_dom,
            "mean_native": mean_native,
            "mean_invasive": mean_invasive,
            "total_treated_area_m2": total_treated_area,
        }
    )

    return df


def first_year_metric_at_least(metric: pd.Series, threshold: float) -> int | None:
    """Return first year index where metric >= threshold, else None."""
    hits = np.where(metric.values >= threshold)[0]
    if hits.size == 0:
        return None
    return int(hits[0])


def render_map(state_grid: np.ndarray, title: str) -> go.Figure:
    """Render one map frame with intuitive ecological colors."""
    colorscale = [
        [0.00, "#5a5a5a"],  # 0 degraded/empty (dark gray)
        [0.20, "#c0392b"],  # 1 invasive-dominated (red)
        [0.40, "#2f2f2f"],  # 2 burned/cleared-ready (black-ish)
        [0.60, "#3f51b5"],  # 3 biomat installed (blue)
        [0.80, "#7ccf6b"],  # 4 recovering native (light green)
        [1.00, "#1b5e20"],  # 5 mature native (dark green)
    ]

    fig = go.Figure(
        data=[
            go.Heatmap(
                z=state_grid,
                zmin=0,
                zmax=5,
                colorscale=colorscale,
                colorbar=dict(
                    title="State",
                    tickvals=[0, 1, 2, 3, 4, 5],
                    ticktext=[
                        "Degraded",
                        "Invasive",
                        "Burned-ready",
                        "Biomat",
                        "Recovering",
                        "Mature",
                    ],
                ),
                showscale=True,
            )
        ]
    )

    fig.update_layout(
        title=title,
        xaxis=dict(showgrid=False, visible=False),
        yaxis=dict(showgrid=False, visible=False, scaleanchor="x", scaleratio=1),
        margin=dict(l=10, r=10, t=45, b=10),
        height=560,
    )
    return fig


def animate_simulation(state_ts: np.ndarray) -> go.Figure:
    """Create animated map with built-in play/pause and year slider."""
    years = state_ts.shape[0] - 1

    colorscale = [
        [0.00, "#5a5a5a"],
        [0.20, "#c0392b"],
        [0.40, "#2f2f2f"],
        [0.60, "#3f51b5"],
        [0.80, "#7ccf6b"],
        [1.00, "#1b5e20"],
    ]

    frames = []
    for yr in range(years + 1):
        frames.append(
            go.Frame(
                data=[
                    go.Heatmap(
                        z=state_ts[yr],
                        zmin=0,
                        zmax=5,
                        colorscale=colorscale,
                        showscale=False,
                    )
                ],
                name=str(yr),
            )
        )

    fig = go.Figure(
        data=[
            go.Heatmap(
                z=state_ts[0],
                zmin=0,
                zmax=5,
                colorscale=colorscale,
                colorbar=dict(
                    title="State",
                    tickvals=[0, 1, 2, 3, 4, 5],
                    ticktext=[
                        "Degraded",
                        "Invasive",
                        "Burned-ready",
                        "Biomat",
                        "Recovering",
                        "Mature",
                    ],
                ),
            )
        ],
        frames=frames,
    )

    steps = [
        {
            "label": str(yr),
            "method": "animate",
            "args": [[str(yr)], {"mode": "immediate", "frame": {"duration": 0, "redraw": True}}],
        }
        for yr in range(years + 1)
    ]

    fig.update_layout(
        title="Landscape recovery animation",
        xaxis=dict(showgrid=False, visible=False),
        yaxis=dict(showgrid=False, visible=False, scaleanchor="x", scaleratio=1),
        height=560,
        margin=dict(l=10, r=10, t=45, b=10),
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "x": 0.0,
                "y": 1.08,
                "buttons": [
                    {
                        "label": "Play",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "frame": {"duration": 240, "redraw": True},
                                "fromcurrent": True,
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                    {
                        "label": "Pause",
                        "method": "animate",
                        "args": [[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False}}],
                    },
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "currentvalue": {"prefix": "Year: "},
                "pad": {"t": 35},
                "steps": steps,
            }
        ],
    )

    return fig


# -----------------------------------------------------------------------------
# Streamlit app
# -----------------------------------------------------------------------------
def app() -> None:
    ensure_runtime_dependencies()
    st.set_page_config(page_title="Chilean Forest Restoration Simulator", layout="wide")
    st.title("Chilean Forest Restoration Simulator")
    st.caption(
        "Post-fire / post-invasive recovery in the Chilean Winter Rainfall-Valdivian Forest region"
    )

    # Biomat geometry metrics
    cell_area_m2, cells_per_m2 = hex_cell_geometry(side_length_m=0.10)
    st.info(
        f"Hex biomat cell area: **{cell_area_m2:.6f} m²** | "
        f"Cells (and seeds) per m²: **{cells_per_m2:.2f}**"
    )

    with st.sidebar:
        st.header("Simulation Controls")
        mode = st.radio("Data mode", ["demo", "real-data-ready"], index=0)
        years = st.slider("Simulation years", min_value=30, max_value=50, value=40)
        grid_size = st.slider("Grid size", min_value=30, max_value=90, value=60, step=10)
        grid_cell_area_m2 = st.number_input(
            "Landscape cell area (m²)", min_value=1.0, max_value=5000.0, value=100.0, step=10.0
        )

        st.subheader("Treatment")
        treatment_fraction = st.slider("Fraction of landscape targeted", 0.05, 0.90, 0.30, 0.05)
        strategy = st.selectbox(
            "Targeting strategy",
            ["Burned/ready first", "Highest invasion first", "Mixed priority"],
            index=0,
        )

        st.subheader("Ecological Parameters")
        r = st.slider("Natural regrowth r", 0.01, 0.20, 0.08, 0.01)
        alpha = st.slider("Biomat effectiveness alpha", 0.00, 0.40, 0.20, 0.01)
        beta = st.slider("Invasive competition beta", 0.00, 0.30, 0.12, 0.01)
        m = st.slider("Mortality m", 0.00, 0.20, 0.03, 0.01)
        g = st.slider("Invasive regrowth g", 0.00, 0.30, 0.10, 0.01)
        c = st.slider("Treatment effectiveness c", 0.00, 1.00, 0.55, 0.01)

        st.subheader("Optional Disturbance")
        drought_prob = st.slider("Drought probability", 0.00, 0.50, 0.00, 0.01)
        fire_prob = st.slider("Fire recurrence probability", 0.00, 0.20, 0.00, 0.01)

    params = ModelParams(
        r=r,
        alpha=alpha,
        beta=beta,
        m=m,
        g=g,
        c=c,
        drought_prob=drought_prob,
        fire_recurrence_prob=fire_prob,
    )

    if "base_landscape" not in st.session_state:
        st.session_state.base_landscape = load_or_generate_landscape(
            mode="demo",
            grid_shape=(grid_size, grid_size),
            random_seed=8,
        )

    if st.button("Reset"):
        selected_mode = "real" if mode == "real-data-ready" else "demo"
        st.session_state.base_landscape = load_or_generate_landscape(
            mode=selected_mode,
            grid_shape=(grid_size, grid_size),
            random_seed=8,
        )
        st.session_state.sim_result = None
        st.session_state.treatment_mask = None

    base = st.session_state.base_landscape

    if mode == "real-data-ready":
        st.warning(
            "Real-data-ready mode is scaffolded: if GIS layers are missing, demo landscape is used. "
            "Use placeholder files: study_area.geojson, burn_scars.geojson, "
            "invasive_occurrences.csv, invasive_suitability.tif"
        )

    treatment_mask = select_treatment_cells(
        invasive_cover=base["invasive_cover"],
        ready_mask=base["ready_mask"],
        fraction=treatment_fraction,
        strategy=strategy,
        seed=11,
    )

    c1, c2, c3 = st.columns(3)
    remove_clicked = c1.button("Remove invasives")
    biomat_clicked = c2.button("Apply biomat")
    run_clicked = c3.button("Run simulation")

    if remove_clicked or biomat_clicked:
        st.session_state.base_landscape = apply_treatment(
            landscape=st.session_state.base_landscape,
            treatment_mask=treatment_mask,
            remove_invasives=remove_clicked,
            install_biomat=biomat_clicked,
        )
        base = st.session_state.base_landscape

    if run_clicked:
        st.session_state.treatment_mask = treatment_mask
        st.session_state.sim_result = run_simulation(
            landscape=base,
            treatment_mask=treatment_mask,
            years=years,
            params=params,
            seed=17,
        )

    left, right = st.columns([1.6, 1.0])

    with left:
        if st.session_state.get("sim_result") is not None:
            sim_result = st.session_state.sim_result
            st.plotly_chart(animate_simulation(sim_result["state_ts"]), use_container_width=True)
        else:
            initial_state = derive_state(
                base["native_cover"],
                base["invasive_cover"],
                base["ready_mask"],
                base["biomat_present"],
            )
            st.plotly_chart(render_map(initial_state, "Initial landscape"), use_container_width=True)

    with right:
        st.subheader("Live Metrics")
        if st.session_state.get("sim_result") is None:
            st.write("Run simulation to view time-series metrics.")
        else:
            sim_result = st.session_state.sim_result
            metrics = compute_metrics(
                native_ts=sim_result["native_ts"],
                invasive_ts=sim_result["invasive_ts"],
                biomat_mask=sim_result["biomat_mask"],
                grid_cell_area_m2=grid_cell_area_m2,
            )

            selected_year = st.slider("Year (metrics view)", 0, years, 0)
            row = metrics.iloc[selected_year]

            years_to_50 = first_year_metric_at_least(metrics["pct_native_gt_07"], 50.0)
            years_to_80 = first_year_metric_at_least(metrics["pct_native_gt_07"], 80.0)

            st.metric("Current year", f"{int(row['year'])}")
            st.metric("% area native > 0.7", f"{row['pct_native_gt_07']:.1f}%")
            st.metric("% invasive-dominated", f"{row['pct_invasive_dom']:.1f}%")
            st.metric("Total treated area", f"{row['total_treated_area_m2']:.0f} m²")
            st.metric("Mean native cover", f"{row['mean_native']:.3f}")

            if years_to_50 is None:
                st.write("Estimated years to 50% restored area: not reached")
            else:
                st.write(f"Estimated years to 50% restored area: {years_to_50}")

            if years_to_80 is None:
                st.write("Estimated years to 80% restored area: not reached")
            else:
                st.write(f"Estimated years to 80% restored area: {years_to_80}")

            fig_native = go.Figure()
            fig_native.add_trace(
                go.Scatter(x=metrics["year"], y=metrics["mean_native"], mode="lines", name="Mean native")
            )
            fig_native.add_vline(x=selected_year, line_width=1, line_dash="dash")
            fig_native.update_layout(
                title="Mean Native Cover vs Time",
                xaxis_title="Year",
                yaxis_title="Mean native cover",
                height=250,
                margin=dict(l=10, r=10, t=35, b=10),
            )
            st.plotly_chart(fig_native, use_container_width=True)

            fig_invasive = go.Figure()
            fig_invasive.add_trace(
                go.Scatter(
                    x=metrics["year"], y=metrics["mean_invasive"], mode="lines", name="Mean invasive"
                )
            )
            fig_invasive.add_vline(x=selected_year, line_width=1, line_dash="dash")
            fig_invasive.update_layout(
                title="Mean Invasive Cover vs Time",
                xaxis_title="Year",
                yaxis_title="Mean invasive cover",
                height=250,
                margin=dict(l=10, r=10, t=35, b=10),
            )
            st.plotly_chart(fig_invasive, use_container_width=True)

    st.subheader("Ecological Model")
    st.markdown(
        """
        **Native cover update**

        `N[t+1] = N[t] + r*N[t]*(1-N[t]) + alpha*B*(1-N[t]) - beta*I[t]*N[t] - m*N[t] + neighborhood_native_spread`

        **Invasive cover update**

        `I[t+1] = I[t] + g*I[t]*(1-I[t]) - c*treatment - d*N[t]*I[t] + neighborhood_invasive_recolonization`

        - Logistic regrowth models natural recovery in available space.
        - Biomat term models seed-assisted establishment in treated cells.
        - Competition and mortality reduce native persistence under stress.
        - Neighborhood terms emulate spatial spread/recolonization.
        """
    )

    st.caption(
        "Color scheme: red=invasive, dark gray=degraded, near-black=burned/ready, "
        "blue=biomat-treated, light green=recovering native, dark green=mature native."
    )

    st.subheader("Real GIS Integration Notes")
    st.code(
        """# Replace synthetic demo data with real Chile layers:
# 1) Put files in project root:
#    - study_area.geojson
#    - burn_scars.geojson
#    - invasive_occurrences.csv
#    - invasive_suitability.tif
# 2) In load_or_generate_landscape(mode='real'), rasterize polygons/points to grid.
# 3) Map burn scars -> ready_mask, invasion layer -> invasive_cover,
#    and baseline forest layer -> native_cover.
# 4) Keep array dimensions consistent with simulation update functions.
# 5) If some layers are missing, fallback logic can still run with partial inputs.
""",
        language="python",
    )


if __name__ == "__main__":
    ensure_runtime_dependencies()
    app()
