# README: Supporting Data and R Code

**Manuscript title:** *Did the Latitudinal Diversity Gradient hold during the Late Palaeozoic Ice Age? Using trace fossils from continental environments to assess diversity and arthropod size patterns*

**Authors:** Gabriel E. B. de Barros, Nicholas J. Minter, Daniel Sedorko, Mírian L. A. F. Pacheco

---

## Overview

This repository contains all R scripts and data files necessary to reproduce the statistical analyses, figures, and paleogeographic reconstructions presented in the manuscript. The code is fully documented and organised to facilitate replication, verification, and extension by other researchers.

**Data files:**

- `Database.xlsx` – Complete ichnological database (occurrences, measurements, metadata)
- `temp_o2.xlsx` – Temperature and oxygen model reconstructions (Royer et al., 2014; Song et al., 2019; Scotese et al., 2021; Mills et al., 2023)
- `280Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv`, `300Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv`, `320Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv`, `340Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv`, `360Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv` – Pohl et al. (2022) Köppen–Geiger climate classification data for each time slice

**R scripts:**

- `diversity_analysis.R` – Diversity, disparity, and latitudinal gradient analyses
- `measurements_analysis.R` – Body size analyses and environmental correlations
- `paleomaps.R` – Paleogeographic map reconstruction and visualisation

**Output:** The scripts generate all figures, statistical summaries, and model outputs reported in the manuscript.

---

## 1. Requirements and Installation

All analyses were performed in **R version 4.6.1** using **RStudio 2026.7.1.147**. The following packages are required:

```r
library(readxl)      # Read Excel files
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(purrr)       # Functional programming
library(ggplot2)     # Data visualisation
library(vegan)       # Community ecology (diversity, PERMANOVA)
library(iNEXT)       # Coverage-based rarefaction
library(mgcv)        # Generalised additive models
library(xgboost)     # Gradient boosting
library(caret)       # Machine learning utilities
library(lme4)        # Linear mixed-effects models
library(lmerTest)    # p-values for lme4
library(MuMIn)       # R² for mixed models
library(performance) # Model performance metrics (ICC, R²)
library(ggeffects)   # Predictions from mixed models
library(reshape2)    # Data reshaping (correlation heatmaps)
library(FSA)         # Dunn's test for multiple comparisons
library(sf)          # Spatial data handling
library(rgplates)    # Plate tectonic reconstructions
library(geojsonsf)   # GeoJSON to Simple Feature conversion
```

---

## 2. File Descriptions

### 2.1 Data Files

#### `Database.xlsx`

This is the primary database, containing six sheets:

| Sheet | Description |
| :--- | :--- |
| **Diversity** | Summary diversity data per deposit |
| **Measurements** | Arthropod trace fossil external widths | 
| **Ichnotaxa** | Ichnotaxonomic reference list |
| **Deposits** | Presence–absence matrix of ichnogenera |
| **Summary 1** | Summary of the data divided by unit |
| **Summary 2** | Summary of the data divided by different variables |
| **Lat** | Ichnodisparity data by latitude |
| **References** | Reference list of all the used works to fill the diversity sheet |
| **Legend** | Legend of the architectural designs and the institutions |

#### `temp_o2.xlsx`

Contains palaeotemperature and palaeo‑oxygen reconstructions used in the correlation analyses.

| Sheet | Columns | Source |
| :--- | :--- | :--- |
| **Temp** | `Scotese` (global temperature); `Song` (seawater); `Geocarbsulf` (continental) | Scotese et al. (2021); Song et al. (2019); Royer et al. (2014) |
| **O2** | `Royer` (atmospheric); `Mils` (atmospheric); `Song` (oceanic dissolved) | Royer et al. (2014); Mills et al. (2023); Song et al. (2019) |

#### `*Ma.csv` (Pohl et al., 2022)

Köppen–Geiger climate classification data for each time slice (280, 300, 320, 340, 360 Ma).

---

### 2.2 R Scripts

#### `diversity_analysis.R`

**Purpose:** Performs all diversity, disparity, and latitudinal gradient analyses, including coverage‑standardised richness, alpha/beta diversity, XGBoost feature importance, and temporal trends.

**Workflow:**

1. Loads data from `Database.xlsx`
2. **Coverage‑based rarefaction** (Chao & Jost, 2012) for climate zones, settings, categories, and time slices using `iNEXT`
3. **XGBoost** feature importance with coverage‑standardised richness as response
4. **Latitudinal diversity patterns** using absolute paleolatitude (LM, GAM, Spearman)
5. **Setting‑specific latitudinal patterns** (Alluvial, Coastal, Lacustrine)
6. **Alpha diversity** (Shannon, Simpson) and **Beta diversity** (Jaccard, PERMANOVA)
7. **Temporal trends** by climate zone and depositional setting

**Outputs:** Figure 2 (diversity plots), Figure 3 (XGBoost and latitudinal models), Appendix D (Supplementary figure S1).

#### `measurements_analysis.R`

**Purpose:** Analyses arthropod trace fossil external widths in relation to palaeolatitude, time, temperature, and oxygen.

**Workflow:**

1. Loads and log‑transforms external width data (`logEW = log(EW)`)
2. **Summary statistics** by zone, setting, category
3. **Kruskal‑Wallis** and **Dunn's post‑hoc** tests
4. **Latitudinal patterns** (global and by time slice: 320, 300, 280 Ma) with LM and GAM
5. **Temporal trends** through geological time
6. **Correlations** between size and temperature/oxygen models (Spearman)
7. **Environmental partitioning**: subaqueous (coastal marine) vs. non‑subaqueous (terrestrial/transitional)
8. **Linear mixed‑effects models** (`lmer`) with random intercepts for Setting, Zone, and Map

**Outputs:** Figure 4 (size by zone/setting/category, latitudinal patterns), Figure 5 (temporal trends, temperature/oxygen correlations, partitioned correlations), Appendix D (summary tables S2–S4).

#### `paleomaps.R`

**Purpose:** Generates paleogeographic maps with Köppen–Geiger climate zones and reconstructed sample localities.

**Workflow:**

1. Reconstructs GPS coordinates using `rgplates` (Scotese, 2016 model)
2. Loads Pohl et al. (2022) climate data for each time slice
3. Plots coastlines, climate zones, and sample localities
4. Saves maps as ggplot2 objects

**Outputs:** Figure 1 (paleogeographic maps with climate zones and localities).

**Note:** Requires an active internet connection for GPlates Web Service access.

---

## 3. Data Availability and Citation

- **Primary database:** `Database.xlsx` (included)
- **Climate data:** Pohl et al. (2022) – see `*Ma.csv` files
- **Temperature/oxygen data:** See `temp_o2.xlsx` and references in manuscript

**Original data sources to cite:**
- Scotese, C.R., 2016. PALEOMAP PaleoAtlas for GPlates and the PaleoData Plotter Program. PALEOMAP Proj. https://doi.org/10.13140/RG2.2.34367.00166
- Pohl, A. et al., 2022. Dataset of Phanerozoic continental climate and Köppen–Geiger climate classes. *Data in Brief* 43, 108424. https://doi.org/10.1016/j.dib.2022.108424
- Royer, D.L. et al., 2014. Error analysis of CO2 and O2 estimates from the long-term geochemical model GEOCARBSULF. *American Journal of Science* 314, 1259–1283. https://doi.org/10.2475/09.2014.01
- Song, H. et al., 2019. Seawater Temperature and Dissolved Oxygen over the Past 500 Million Years. *Journal of Earth Science* 30, 236–243. https://doi.org/10.1007/s12583-018-1002-2
- Scotese, C.R. et al., 2021. Phanerozoic paleotemperatures: The earth’s changing climate during the last 540 million years. *Earth-Science Reviews* 215, 103503. https://doi.org/10.1016/j.earscirev.2021.103503
- Mills, B.J.W. et al., 2023. Evolution of Atmospheric O2 Through the Phanerozoic, Revisited. *Annual Review of Earth and Planetary Sciences* 51, 253–276. https://doi.org/10.1146/annurev-earth-032320-095425

---

## 4. Contact

For questions, comments, or bug reports, please contact:

**Gabriel E. B. de Barros**
Email: gbareabarros@gmail.com
Laboratory of Paleobiology and Astrobiology, Federal University of São Carlos, Brazil
