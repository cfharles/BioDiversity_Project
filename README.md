# Biodegradable Biomat Restoration Model

A hybrid **Lotka-Volterra + Cellular Automata** simulation of native vegetation recovery against invasive species pressure in the Chilean Araucanía / Biobío region (~38°S), with biodegradable biomat deployment as the restoration intervention.

---

## Background

The Chilean Winter Rainfall–Valdivian Forests hotspot is one of the world's most threatened biodiversity zones. Invasive species — primarily *Pinus radiata* (Monterey pine), *Eucalyptus globulus* (blue gum), and *Acacia dealbata* (silver wattle) — have displaced native communities including the keystone *Araucaria araucana* (monkey-puzzle tree) and understory species such as *Quillaja saponaria*, *Cryptocarya alba*, and *Peumus boldus*.

The proposed intervention is a **biodegradable biomat**: a hexagonal cellulose mat made from timber industry waste, placed around native seedlings to suppress invasive competition and protect seedlings during their vulnerable establishment phase.

---

## Model Architecture

The simulation runs on a **100 × 100 spatial grid**, where each cell represents a 100 m × 100 m patch (1 hectare). The pilot area covers approximately 10 km × 10 km in the Andean foothills of Araucanía.

Each cell stores three state variables:

- `N` — native vegetation cover fraction [0, 1]
- `I` — invasive vegetation cover fraction [0, 1]
- `biomat_age` — years since biomat was deployed (−1 if never deployed)

The constraint `N + I ≤ 1` is enforced at every time step. The remainder `1 − N − I` represents bare ground, litter, or non-vegetated surface.

Each yearly time step applies three operations in sequence:

1. **Biomat deployment** — scheduled cells receive a biomat, invasives are cleared
2. **Local Lotka-Volterra dynamics** — competition equations update N and I per cell
3. **Seed dispersal** — cellular automata spread seeds to neighbouring cells

---

## The Mathematics

### Lotka-Volterra Competition (per cell, per year)

**Without biomat:**
$$
N_{t+1} = N_t + r_N N_t \left(1 - \frac{N_t + \alpha_{NI} I_t}{K_N}\right) dt
$$

$$
I_{t+1} = I_t + r_I I_t \left(1 - \frac{I_t + \alpha_{IN} N_t}{K_I}\right) dt
$$

**With active biomat (years 0–3 post-deployment):**
$$
N_{t+1} = N_t + (r_N \cdot \mathrm{boost}) N_t \left(1 - \frac{N_t + \alpha_{NI,\mathrm{reduced}} I_t}{K_N}\right) dt
$$

$$
I_{t+1} = I_t + r_I I_t \left(1 - \frac{I_t + \alpha_{IN} N_t}{K_I}\right) dt - \mathrm{suppression}
$$

Where:
- `rN`, `rI` = intrinsic growth rates
- `KN`, `KI` = carrying capacities (both = 1.0, normalised cover fractions)
- `α_NI` = competitive effect of invasives on natives (reduced by biomat)
- `α_IN` = competitive effect of natives on invasives (increases as natives mature)

### Physical Validity Constraint
After each update, cells where `N + I > 1` are rescaled proportionally:
$$
N \leftarrow \frac{N}{N + I}
$$

$$
I \leftarrow \frac{I}{N + I}
$$
This preserves the competitive ratio between species while enforcing total coverage limits.

### Seed Dispersal (per year)
$$
\mathrm{invasive\_increment} = \mathrm{spread\_rate}_I \cdot \mathrm{neighbourhood\_avg}(I) \cdot (1 - I) \cdot \mathrm{noise}
$$

$$
\mathrm{native\_increment} = \mathrm{spread\_rate}_N \cdot \mathrm{neighbourhood\_avg}(N) \cdot (1 - N) \cdot \mathrm{noise}
$$
Where `neighbourhood_avg` averages source cells within a circular radius, and `noise` is a uniform random factor (±15–40%) representing variable weather and animal disperser activity.

---

## Parameters

All parameters are set in `parameters.py`. Values were informed by reading the literature on Araucaria and Pinus radiata ecology in Chile, then tested and tuned through repeated model runs to produce dynamics consistent with observed invasion timelines and restoration outcomes.

---

### Grid and Simulation

| Parameter | Value |
|---|---|
| Grid size | 100 × 100 cells |
| Cell size | 100 m × 100 m (1 ha per cell) |
| Total pilot area | 10,000 ha |
| Time step | 1 year |
| Simulation duration | 75 years |

---

### Native Species
*(Combined guild: Araucaria araucana + Quillaja saponaria + Cryptocarya alba + Peumus boldus)*

| Parameter | Value | What it represents |
|---|---|---|
| Growth rate `rN` | 0.05 /year | Natives gain ~5% cover per year in ideal conditions — reflects how slowly *Araucaria* establishes |
| Carrying capacity `KN` | 1.0 | Maximum possible cover (normalised fraction) |
| `alpha_NI` without biomat | 1.8 | Each unit of invasive cover removes 1.8 units of native growing space — invasives punch above their weight |
| `alpha_NI` with biomat | 0.2 | Biomat physically blocks invasive competition, dropping pressure to near zero |
| `alpha_IN` young | 0.4 | Young native seedlings provide little suppression of invasives |
| `alpha_IN` mature | 0.9 | Closed native canopy strongly suppresses invasive regeneration |
| Spread radius | 1 cell (100 m) | *Araucaria* seeds are heavy (~10 g) and only travel short distances |
| Spread rate | 0.075 /year | Probability of successful colonisation of a neighbouring cell per year |

---

### Invasive Species
*(Combined guild: Pinus radiata + Eucalyptus globulus + Acacia dealbata)*

| Parameter | Value | What it represents |
|---|---|---|
| Growth rate `rI` | 0.20 /year | Invasives gain ~20% cover per year — 4× faster than natives, consistent with *Pinus radiata* growing up to 2 m/year |
| Carrying capacity `KI` | 1.0 | Maximum possible cover (normalised fraction) |
| Spread radius | 2 cells (200 m) | Wind-dispersed pine seeds reach further than native seeds |
| Spread rate | 0.12 /year | Invasive colonisation rate per neighbouring source cell per year |
| Regrowth probability | 0.15 /year | 15% annual chance of re-establishment from soil seed banks after clearing |
| Regrowth increment | 0.08 | Cover fraction re-added per regrowth event |

---

### Biomat

| Parameter | Value | What it represents |
|---|---|---|
| Active years | 3 | Full protection for 36 months before the mat begins degrading |
| Fade years | 2 | Mat fully absorbed into soil by year 5 (months 48–60) |
| Native growth boost | 2.5× | Biomat nutrients and mycorrhizal inoculants boost native growth rate during active phase |
| Invasive removal fraction | 0.80 | 80% of invasive cover cleared manually at installation |
| Suppression rate | 0.05 /year | Annual fractional reduction in invasive cover while the mat is physically present |

---

### Starting Conditions

| Parameter | Value | What it represents |
|---|---|---|
| Invasive-dominated cells | 85% | Starting state — heavily invaded Araucanía landscape |
| Native patch cells | 5% | Remnant native forest fragments that survived the invasion front |
| Bare ground cells | 10% | Degraded or cleared ground |
| Invasive cover in invaded cells | 70–90% | Range of invasive density across the dominated area |
| Native cover in native patches | 10–30% | Fragmented, suppressed remnant patches |

---

## Scenarios

All intervention scenarios use a single constant deployment rate across the full 75-year simulation, prioritising cells where invasive cover exceeds native cover. Three deployment modes are supported:

- **Radial** — biomats expand outward from one or more center points
- **Multi-point** — biomats deploy simultaneously from five distributed centers across the grid
- **Encirclement** — biomats form a hollow ring at a fixed radius around a center point, surrounding an invasive-dominated core and squeezing inward over time

| Scenario | Annual Rate | Total over 75 years | Mode |
|---|---|---|---|
| No Intervention | 0 ha/yr | — | — |
| Light Deployment | 25 ha/yr | 1,875 ha | Radial (1 center) |
| Moderate Deployment | 75 ha/yr | 5,625 ha | Radial (1 center) |
| Aggressive Deployment | 125 ha/yr | 9,375 ha | Radial (1 center) |
| Multi-Point Deployment | 25 ha/yr per center (125 ha/yr total) | 9,375 ha | Radial (5 centers) |
| Encirclement Ring | 75 ha/yr | 5,625 ha | Hollow circle, radius 28 cells (~175 ha enclosed) |

---

## Simulation Results

Results below are from the previous scenario configuration and will be updated after re-running the model with the new constant-rate scenarios. Re-run with `python3 -m biomat_model.run_simulation` to regenerate.

---

## How To Run

From the project root directory:

```bash
python3 -m biomat_model.run_simulation
```

This runs the 75-year simulation across all six scenarios.

Dependencies: `numpy`, `matplotlib`, `imageio`. Install with:

```bash
pip install numpy matplotlib imageio pillow
```

---

## Outputs

| Path | Description |
|---|---|
| `output_frames/<scenario>/` | 76 PNG frames (year 0–75) showing spatial vegetation map + time series |
| `results/<scenario>_animation.gif` | Animated GIF of the 75-year spread |
| `results/<scenario>_yearly_stats.csv` | Year-by-year statistics |
| `results/<scenario>_trajectory.png` | 3-panel summary chart |
| `results/scenario_comparison.png` | Side-by-side year-75 maps for all scenarios |

---

## Limitations and Caveats

- This is a **proof-of-concept model** for presentation purposes. It is designed to be ecologically plausible, not site-calibrated or validated against field data.
- Several parameters (spread rates, alpha values, growth boost) are **model calibrations** with no direct published citation. They are tuned to produce dynamics consistent with qualitative descriptions in the literature, not fitted to measured data.
- The model treats each species guild (native, invasive) as a **single functional type**, ignoring differences between *Araucaria* and understory species, or between *Pinus*, *Eucalyptus*, and *Acacia*.
- No **climate, topography, or soil** variables are included. The pilot area is treated as spatially homogeneous except for initial conditions.
- Biomat parameters (growth boost, suppression, degradation timeline) are based on **design specifications**, not published field trials of this specific product.
- The model does not account for **fire**, which is a major driver of both invasive spread and native regeneration in this region.

---

## References

- Bustamante, R.O. & Simonetti, J.A. (2005). Is *Pinus radiata* invading the native vegetation in central Chile? *Biological Invasions*, 7(2), 243–249.
- Pauchard, A. & Alaback, P.B. (2004). Influence of elevation, land use, and landscape context on patterns of alien plant invasions in protected areas of south-central Chile. *Conservation Biology*, 18(1), 238–248.
- Richardson, D.M. & Higgins, S.I. (1998). Pines as invaders in the southern hemisphere. In: *Ecology and Biogeography of Pinus*. Cambridge University Press, pp. 450–473.
- Veblen, T.T. et al. (1995). Snowstorm damage and the dynamics of a southern Andean *Nothofagus* forest. *Journal of Ecology*, 83(5), 911–922.
- Vilà, M. et al. (2011). Ecological impacts of invasive alien plants: a meta-analysis. *Ecology Letters*, 14(7), 702–708.
