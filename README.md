# County-Level Grain Production Correlation and Cooperative Mergers

R code and data for:

> Hassanizadeh, B., & Kenkel, P. (2026). County Level Grain Production Correlation and
> Implications for Cooperative Mergers. *Journal of Cooperatives*, 42, 133–165.

Grain cooperatives that merge across counties pool their production risk. The paper
measures how correlated grain production is between Oklahoma and Kansas counties
(2013–2023), and estimates how that correlation depends on east–west and north–south
distance, crop mix, adjacency and county size. Merging with a partner whose production
is less correlated diversifies a cooperative's volume risk.

## Method in brief

1. Total grain production (corn, soybean, sorghum, wheat; million bushels) is built per
   county-year from USDA NASS data. Counties with fewer than two years of positive
   production are dropped.
2. Log production is detrended with a county-specific linear trend; the residuals are
   correlated pairwise over the years both counties report.
3. The correlation `rho` is regressed on `x_dir` and `y_dir` (east–west and north–south
   distance, 100-mile units), `ratio_r` (crop-diversity ratio), `adj` (shared border)
   and `small_prod` (bottom-quartile producer), by OLS and by a Bayesian model (`brms`)
   whose priors are centred on the OLS estimates.

## Repository layout

```
R/
  PairWise_final.R   all county pairs: OLS + Bayesian regression
  Baseline_final.R   one baseline county vs. all others: OLS + Bayesian regression
  Map_Full.R         map of correlation with a baseline county
  Map_D_Full.R       maps of crop-diversity ratio, total production, small producers
data/
  Data.xlsx          county production, one sheet per crop (USDA NASS Quick Stats)
  Points.xls         county centroids for Kansas and Oklahoma
  Adj.csv            county adjacency (polygon neighbours)
figures/             maps used in the paper
```

## Requirements

R (tested with 4.5) and the packages:

```r
install.packages(c("readxl", "dplyr", "readr", "writexl",
                   "ggplot2", "sf", "tigris", "purrr", "ggnewscale", "viridis", "ggspatial",
                   "brms"))
# cmdstanr (backend for brms):
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## How to run

All paths are relative to the repository root. Open `grain-production-correlation.Rproj`
in RStudio (or `setwd()` to this folder), then run any script, e.g.

```r
source("R/PairWise_final.R")
```

| Script | Output |
|---|---|
| `PairWise_final.R` | OLS and Bayesian fits printed to the console |
| `Baseline_final.R` | same, for the county set in `baseline_county`; also `output/Regression_Data_with_County_Names.xlsx` |
| `Map_Full.R` | correlation map in `figures/` (baseline set at the top of the script) |
| `Map_D_Full.R` | `diversity_ratio_map.jpg`, `total_production_map.jpg`, `small_production_map.jpg` in `figures/` |

The map scripts download county boundaries with `tigris`, so they need an internet
connection on the first run. The Bayesian step (20,000 iterations × 4 chains) takes
several minutes.

## Data sources

- Production: USDA National Agricultural Statistics Service, Quick Stats — county-level
  production in bushels for corn, soybean, sorghum and wheat, Kansas and Oklahoma, 2013–2023.
- County centroids and boundaries: U.S. Census Bureau TIGER/Line (boundaries fetched at run time).

## License

Code is released under the MIT License (see `LICENSE`). The data are U.S. government
public-domain sources. If you use this code, please cite the paper above.
