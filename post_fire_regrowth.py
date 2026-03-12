#!/usr/bin/env python3
"""Simulate post-fire native vegetation recovery with and without biomat support."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def simulate_regrowth(
    years: int,
    n0: float,
    r: float,
    alpha: float,
    m: float,
    biomat_intensity: float,
) -> np.ndarray:
    """Run a discrete-time vegetation cover simulation.

    Parameters
    ----------
    years : int
        Number of years to simulate.
    n0 : float
        Initial native vegetation cover fraction (0 to 1).
    r : float
        Natural annual regrowth rate (logistic regeneration strength).
    alpha : float
        Biomat-assisted establishment rate (seed + biomaterial effect).
    m : float
        Mortality / environmental stress rate.
    biomat_intensity : float
        Fractional biomat deployment intensity M (0 to 1).

    Returns
    -------
    np.ndarray
        Native vegetation cover for each year from 0..years.
    """
    n = np.zeros(years + 1, dtype=float)
    n[0] = np.clip(n0, 0.0, 1.0)

    for t in range(years):
        growth_natural = r * n[t] * (1.0 - n[t])
        growth_biomat = alpha * biomat_intensity * (1.0 - n[t])
        losses = m * n[t]

        n[t + 1] = n[t] + growth_natural + growth_biomat - losses
        n[t + 1] = np.clip(n[t + 1], 0.0, 1.0)

    return n


def first_year_reaching(series: np.ndarray, threshold: float) -> int | None:
    """Return the first year index where series reaches threshold, or None."""
    reached = np.where(series >= threshold)[0]
    if reached.size == 0:
        return None
    return int(reached[0])


def main() -> None:
    # Model parameters (ecological interpretation):
    # r: natural regeneration speed as cover increases.
    # alpha: added establishment from seed-containing biomat cells.
    # m: annual stress/mortality reducing existing cover.
    r = 0.08
    alpha = 0.20
    m = 0.03
    years = 50
    n0 = 0.05

    # Regular hexagonal cell geometry (one seed per cell).
    s = 0.10  # side length in meters
    a_cell = (3.0 * np.sqrt(3.0) / 2.0) * (s**2)
    cells_per_m2 = 1.0 / a_cell

    print("Hexagonal biomat cell geometry")
    print(f"  Side length s: {s:.2f} m")
    print(f"  Cell area A_cell: {a_cell:.6f} m^2")
    print(f"  Cells per m^2 (and seeds per m^2): {cells_per_m2:.2f}")
    print()

    scenarios = {
        "No biomat (M=0.0)": 0.0,
        "Partial biomat (M=0.4)": 0.4,
        "Full biomat (M=0.8)": 0.8,
    }

    time = np.arange(years + 1)
    results: dict[str, np.ndarray] = {}

    for name, m_intensity in scenarios.items():
        results[name] = simulate_regrowth(
            years=years,
            n0=n0,
            r=r,
            alpha=alpha,
            m=m,
            biomat_intensity=m_intensity,
        )

    print("Final native vegetation cover after 50 years")
    for name, series in results.items():
        print(f"  {name}: {series[-1]:.4f}")
    print()

    thresholds = [0.50, 0.80]
    for name, series in results.items():
        print(name)
        for threshold in thresholds:
            year_hit = first_year_reaching(series, threshold)
            pct = int(threshold * 100)
            if year_hit is None:
                print(f"  Never reaches {pct}% cover within {years} years.")
            else:
                print(f"  Reaches {pct}% cover in year {year_hit}.")
        print()

    plt.figure(figsize=(9, 5.5))
    for name, series in results.items():
        plt.plot(time, series, linewidth=2.2, label=name)

    plt.xlabel("Time (years)")
    plt.ylabel("Native vegetation cover fraction")
    plt.title("Post-fire native vegetation recovery with biomat-assisted establishment")
    plt.ylim(0.0, 1.02)
    plt.xlim(0, years)
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Ecological meaning summary:
    # - Logistic regrowth term r*N*(1-N) captures natural regeneration limits.
    # - Biomat term alpha*M*(1-N) adds seed-assisted establishment in open/degraded space.
    # - Mortality term m*N represents ongoing post-fire stress and seedling loss.


if __name__ == "__main__":
    main()
