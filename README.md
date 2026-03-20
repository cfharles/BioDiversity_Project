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
```
N(t+1) = N(t) + rN · N(t) · (1 − (N(t) + α_NI · I(t)) / KN) · dt
I(t+1) = I(t) + rI · I(t) · (1 − (I(t) + α_IN · N(t)) / KI) · dt
```

**With active biomat (years 0–3 post-deployment):**
```
N(t+1) = N(t) + (rN × boost) · N(t) · (1 − (N(t) + α_NI_reduced · I(t)) / KN) · dt
I(t+1) = I(t) + rI · I(t) · (1 − (I(t) + α_IN · N(t)) / KI) · dt − suppression
```

Where:
- `rN`, `rI` = intrinsic growth rates
- `KN`, `KI` = carrying capacities (both = 1.0, normalised cover fractions)
- `α_NI` = competitive effect of invasives on natives (reduced by biomat)
- `α_IN` = competitive effect of natives on invasives (increases as natives mature)

### Physical Validity Constraint
After each update, cells where `N + I > 1` are rescaled proportionally:
```
N ← N / (N + I)
I ← I / (N + I)
```
This preserves the competitive ratio between species while enforcing total coverage limits.

### Seed Dispersal (per year)
```
invasive_increment = spread_rate_I · neighbourhood_avg(I) · (1 − I) · noise
native_increment   = spread_rate_N · neighbourhood_avg(N) · (1 − N) · noise
```
Where `neighbourhood_avg` averages source cells within a circular radius, and `noise` is a uniform random factor (±15–40%) representing variable weather and animal disperser activity.

---

## Parameters and Justifications

All parameters are set in `parameters.py`. Below is a full account of every parameter, its value, and the basis for that value. Where no published citation exists, this is stated explicitly.

---

### Grid and Simulation

| Parameter | Value | Justification |
|---|---|---|
| Grid size | 100 × 100 cells | Chosen to balance spatial resolution and computational cost. No citation — modelling decision. |
| Cell size | 100 m × 100 m | Standard resolution for landscape-scale vegetation models. No citation — modelling decision. |
| Time step (dt) | 1 year | Appropriate for slow-growing native species; standard for decadal vegetation models. No citation — modelling decision. |
| Simulation duration | 30 years | Covers one full Araucaria establishment cycle and multiple biomat deployment waves. No citation — project decision. |
| Random seed | 42 | Fixed for reproducibility. No citation. |

---

### Native Species Parameters
*(Combined guild: Araucaria araucana + Quillaja saponaria + Cryptocarya alba + Peumus boldus)*

| Parameter | Value | Justification and Citation |
|---|---|---|
| Intrinsic growth rate `rN` | 0.05 /year | *Araucaria araucana* is one of the world's slowest-growing conifers, taking 5–10 years to exceed grass height. Understory shrubs are faster but similarly suppressed under invasive canopy. The blended value of 0.05 is consistent with observed Araucaria regeneration dynamics. **Veblen, T.T. et al. (1995). "Snowstorm damage and the dynamics of a southern Andean *Nothofagus* forest." *Journal of Ecology*, 83(5), 911-922** documents Andean forest recovery rates. Growth rate also consistent with **Lara, A. et al. (2009)** work on slow-growing Chilean conifers, though the specific value 0.05 is a model calibration, not a directly measured parameter. |
| Carrying capacity `KN` | 1.0 | Normalised cover fraction. Mature Valdivian forests can achieve 80–90% canopy closure. **No citation for this specific normalisation — modelling convention.** |
| `alpha_NI` (no biomat) | 1.8 | Competitive effect of invasives on natives. Value > 1 indicates invasives displace natives from shared resources disproportionately — a well-documented outcome in invaded Araucaria zones. **Pauchard, A. & Alaback, P.B. (2004). "Influence of elevation, land use, and landscape context on patterns of alien plant invasions along roadsides in protected areas of south-central Chile." *Conservation Biology*, 18(1), 238-248** documents strong suppression of native recruitment under invasive pine. The specific value 1.8 is a model calibration, not directly measured. |
| `alpha_NI` (with biomat) | 0.2 | Represents near-elimination of invasive competitive pressure when the biomat's blackout layer is active. The design specification states the mat blocks >99% of photosynthetically active radiation (PAR) to invasive seedlings. **No published citation for the specific value — derived from biomat design specification.** |
| Native spread radius | 1 cell (100 m) | *Araucaria* piñones are large, heavy seeds (3–4 cm, ~10 g) dispersed by gravity and corvids (*Enicognathus ferrugineus*, slender-billed parakeet) over short distances only. **Veblen, T.T. (1982). "Growth patterns of *Chusquea* bamboo in the understory of Chilean *Nothofagus* forests and their influences on tree regeneration." *Biotropica*, 14(2), 116-124** and **Burns, B.R. (1993). "Fire-induced dynamics of *Araucaria araucana* – *Nothofagus antarctica* forest in the Andes of southern Chile." *Journal of Biogeography*, 20, 669-685** both document limited Araucaria seed dispersal range. |
| Native spread rate | 0.075 /year | Rate of successful colonisation per neighbouring cell per year. **No published citation — calibrated by model tuning to produce plausible rates of natural regeneration from refugia.** |
| `alpha_IN` early | 0.4 | Competitive suppression of invasives by young, sparse native vegetation. Low because seedlings and small trees provide little canopy shading. **No published citation for specific value — consistent with general Lotka-Volterra competition theory and the weak-competitor status of native seedlings documented in Vilà, M. et al. (2011). "Ecological impacts of invasive alien plants: a meta-analysis of their effects on species, communities and ecosystems." *Ecology Letters*, 14(7), 702-708.** |
| `alpha_IN` mature | 0.9 | Competitive suppression of invasives by established native canopy. High because mature *Araucaria* + dense understory produces strong shading (leaf area index > 4). **No published citation for specific value — calibrated to reflect documented suppression of invasive regeneration under closed native canopy, consistent with Bustamante, R.O. & Simonetti, J.A. (2005). "Is *Pinus radiata* invading the native vegetation in central Chile? Demographic responses in a fragmented forest." *Biological Invasions*, 7(2), 243-249.** |
| Maturity threshold | 0.5 | Native cover level above which `alpha_IN` begins approaching its maximum. **No published citation — modelling decision.** |
| Restoration source threshold | 0.6 | Native cover level above which a cell acts as a strong seed source. **No published citation — modelling decision.** |

---

### Invasive Species Parameters
*(Combined guild: Pinus radiata + Eucalyptus globulus + Acacia dealbata)*

| Parameter | Value | Justification and Citation |
|---|---|---|
| Intrinsic growth rate `rI` | 0.20 /year | *Pinus radiata* grows up to 2 m/year in Chilean conditions and can reach 25 m in 20 years. *Acacia dealbata* colonises disturbed areas extremely rapidly. **Richardson, D.M. & Higgins, S.I. (1998). "Pines as invaders in the southern hemisphere." In: Richardson, D.M. (ed.) *Ecology and Biogeography of Pinus*. Cambridge University Press, pp. 450-473** documents this growth behaviour. The specific value 0.20 is a model calibration consistent with observed cover expansion rates. |
| Carrying capacity `KI` | 1.0 | Invasive monocultures can reach near-complete ground cover. **No citation for this specific normalisation — modelling convention.** |
| Invasive spread radius | 2 cells (200 m) | *Pinus radiata* seed dispersal has been documented up to 3 km in some studies, but effective colonisation of new patches typically occurs over shorter distances. Scaled to the 100 m grid resolution. **Nathan, R. et al. (2002). "Mechanisms of long-distance dispersal of seeds by wind." *Nature*, 418, 409-413** documents pine seed dispersal mechanisms. The value of 2 cells is a conservative grid-scaled calibration, not a directly measured parameter. |
| Invasive spread rate | 0.12 /year | Rate of invasive colonisation per source cell per year. **No published citation — calibrated to produce an invasive takeover rate consistent with observed Pinus invasion timelines in Araucanía (~20 years to landscape dominance from plantation edges).** |
| Regrowth probability | 0.15 /year | Annual probability of invasive regrowth from persistent soil seed banks after clearing. *Pinus* and *Acacia* seed banks can remain viable for 5–20+ years. **No specific citation for the 0.15 value — consistent with general seed bank persistence literature (e.g. Thompson, K. et al. (1997). *The soil seed banks of north west Europe*. Cambridge University Press) but not calibrated from Chilean-specific data.** |
| Regrowth increment | 0.08 | Cover fraction added per regrowth event. **No published citation — modelling decision.** |

---

### Biomat Parameters

| Parameter | Value | Justification and Citation |
|---|---|---|
| Active years | 3 years | The biomat provides full protection for 36 months before degradation begins. **Source: biomat design specification (cellulose degradation rate in humid temperate conditions). No independent published citation.** |
| Fade years | 2 years | Biomat fully absorbed into soil by year 5 (48–60 months). **Source: biomat design specification. No independent published citation.** |
| Native growth boost | 2.5× | Multiplier on native growth rate during active biomat protection, reflecting embedded slow-release nutrients and mycorrhizal inoculants. The Global Trees Campaign has reported ~90% 10-year seedling survival with structured protection vs. ~10–20% without. **No published citation for the specific 2.5× multiplier — calibrated to produce a survival improvement consistent with reported field outcomes. Mycorrhizal effects on native plant growth consistent with Hoeksema, J.D. et al. (2010). "A meta-analysis of context-dependency in plant response to inoculation with mycorrhizal fungi." *Ecology Letters*, 13(3), 394-407**, though that paper does not study this specific system. |
| Invasive removal fraction | 0.80 | 80% of invasive cover is removed by manual/mechanical clearing at installation. **No published citation — based on project protocol assumption. No field trial data from this specific intervention exists.** |
| Suppression rate | 0.05 /year | Annual fractional reduction in invasive cover from the physical barrier while the biomat is active. **No published citation — modelling decision.** |

---

### Initial Conditions

| Parameter | Value | Justification and Citation |
|---|---|---|
| Invasive-dominated fraction | 0.85 | ~85% of the pilot area assumed to be invasive-dominated at the start. Consistent with reported land cover in heavily invaded Araucanía hillsides. **No specific citation for this exact figure — calibrated to represent a severely degraded landscape.** |
| Native patch fraction | 0.05 | ~5% of cells initialised as remnant native patches. Represents isolated forest fragments that survived the invasion front. **No specific citation — modelling decision.** |
| Bare fraction | 0.10 | ~10% of cells initialised as bare or heavily degraded ground. **No specific citation — modelling decision.** |
| Invasive cover range | 0.70–0.90 | Cover range in invasive-dominated cells. **No specific citation — calibrated estimate.** |
| Native patch cover range | 0.10–0.30 | Cover range in remnant native patch cells. **No specific citation — calibrated estimate.** |

---

## Scenarios

| Scenario | Deployment Rate |
|---|---|
| No Intervention | 0 cells/year — invasive spread continues unchecked |
| Moderate Deployment | 50 cells/yr (years 1–5) → 100 cells/yr (years 6–15) → 150 cells/yr (years 16–30) |
| Aggressive Deployment | 100 cells/yr (years 1–5) → 200 cells/yr (years 6–15) → 300 cells/yr (years 16–30) |

Biomats deploy in expanding rings from the central pilot point (cell [50, 50]), prioritising cells where invasive cover exceeds native cover.

---

## Simulation Results (30-year outcomes after physical validity fix)

| Scenario | Native Cover Y0 | Native Cover Y30 | Invasive Cover Y30 | Restored Cells Y30 |
|---|---|---|---|---|
| No Intervention | 5.6% | 2.5% | 97.7% | 87 ha |
| Moderate Deployment | 5.6% | 7.9% | 86.2% | 589 ha |
| Aggressive Deployment | 5.6% | 13.3% | 74.6% | 1,098 ha |

---

## How To Run

From the project root directory:

```bash
python3 -m biomat_model.run_simulation
```

Dependencies: `numpy`, `matplotlib`, `imageio`. Install with:

```bash
pip install numpy matplotlib imageio pillow
```

---

## Outputs

| Path | Description |
|---|---|
| `output_frames/<scenario>/` | 31 PNG frames (one per year) showing spatial vegetation map + time series |
| `results/<scenario>_animation.gif` | Animated GIF of the 30-year spread |
| `results/<scenario>_yearly_stats.csv` | Year-by-year statistics (native cover, invasive cover, restored cells, biomats deployed) |
| `results/<scenario>_trajectory.png` | 3-panel summary chart for each scenario |
| `results/scenario_comparison.png` | Side-by-side year-30 maps and trajectory overlay for all 3 scenarios |

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

The following are the published sources referenced in the parameter justifications above. Parameters marked "no published citation" in the table above are not supported by any of these references.

- Bustamante, R.O. & Simonetti, J.A. (2005). Is *Pinus radiata* invading the native vegetation in central Chile? Demographic responses in a fragmented forest. *Biological Invasions*, 7(2), 243–249.
- Burns, B.R. (1993). Fire-induced dynamics of *Araucaria araucana* – *Nothofagus antarctica* forest in the Andes of southern Chile. *Journal of Biogeography*, 20, 669–685.
- Hoeksema, J.D. et al. (2010). A meta-analysis of context-dependency in plant response to inoculation with mycorrhizal fungi. *Ecology Letters*, 13(3), 394–407.
- Nathan, R. et al. (2002). Mechanisms of long-distance dispersal of seeds by wind. *Nature*, 418, 409–413.
- Pauchard, A. & Alaback, P.B. (2004). Influence of elevation, land use, and landscape context on patterns of alien plant invasions along roadsides in protected areas of south-central Chile. *Conservation Biology*, 18(1), 238–248.
- Richardson, D.M. & Higgins, S.I. (1998). Pines as invaders in the southern hemisphere. In: Richardson, D.M. (ed.) *Ecology and Biogeography of Pinus*. Cambridge University Press, pp. 450–473.
- Thompson, K. et al. (1997). *The Soil Seed Banks of North West Europe*. Cambridge University Press.
- Veblen, T.T. (1982). Growth patterns of *Chusquea* bamboo in the understory of Chilean *Nothofagus* forests and their influences on tree regeneration. *Biotropica*, 14(2), 116–124.
- Veblen, T.T. et al. (1995). Snowstorm damage and the dynamics of a southern Andean *Nothofagus* forest. *Journal of Ecology*, 83(5), 911–922.
- Vilà, M. et al. (2011). Ecological impacts of invasive alien plants: a meta-analysis of their effects on species, communities and ecosystems. *Ecology Letters*, 14(7), 702–708.
