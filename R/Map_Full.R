# =============================================================================
# Correlation Map — FINAL VERSION USING EXACT Map_D.R STYLE
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(tigris)
library(purrr)
library(ggnewscale)
library(viridis)
library(ggspatial)

# =============================================================================
# 1. SETTINGS
# =============================================================================
main_dir           <- "data"      # input folder, relative to the repository root
output_dir         <- "figures"
baseline_county    <- "BARTON"          # File name / Legend county name / correlation with ... county 
baseline_statefp   <- 20                # 20 = Kansas, 40 = Oklahoma
min_years_required <- 2

options(tigris_use_cache = TRUE)
set.seed(456)

# =============================================================================
# 2. LOAD CROP DATA
# =============================================================================
Corn    <- read_excel(file.path(main_dir, "Data.xlsx"), sheet = "Corn") %>%
  mutate(Value = as.numeric(as.character(Value))) %>%
  group_by(County, `State ANSI`, Year) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

Soybean <- read_excel(file.path(main_dir, "Data.xlsx"), sheet = "Soybean") %>%
  mutate(Value = as.numeric(as.character(Value))) %>%
  group_by(County, `State ANSI`, Year) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

Sorghum <- read_excel(file.path(main_dir, "Data.xlsx"), sheet = "Sorghum") %>%
  mutate(Value = as.numeric(as.character(Value))) %>%
  group_by(County, `State ANSI`, Year) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

Wheat   <- read_excel(file.path(main_dir, "Data.xlsx"), sheet = "Wheat") %>%
  mutate(Value = as.numeric(as.character(Value))) %>%
  group_by(County, `State ANSI`, Year) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

# =============================================================================
# 3. LOAD POINTS & BUILD PRODUCTION MATRIX
# =============================================================================
Shapefile <- read_excel(file.path(main_dir, "Points.xls"))
Shapefile$NAME <- toupper(Shapefile$NAME)
Shapefile <- Shapefile %>% filter(STATEFP %in% c(40, 20))

years <- 2013:2023
N <- nrow(Shapefile)
# -------------------------------------------------------------------
# 3. Initialize matrices with NA
# -------------------------------------------------------------------
prod_matrix   <- matrix(NA, nrow = N, ncol = 11, dimnames = list(NULL, years))
Ratio_matrix  <- matrix(NA, nrow = N, ncol = 11, dimnames = list(NULL, years))
valid_years   <- integer(N)  # number of years with positive production

# -------------------------------------------------------------------
# 4. Fill matrices — only when total production > 0
# -------------------------------------------------------------------
for (i in seq_len(N)) {
  cnty <- Shapefile$NAME[i]
  st   <- Shapefile$STATEFP[i]
  
  for (t in seq_along(years)) {
    yr <- years[t]
    
    C  <- Corn$Value[Corn$County == cnty & Corn$`State ANSI` == st & Corn$Year == yr]
    C  <- ifelse(length(C) == 0 || is.na(C), 0, C / 1e6)
    
    Sr <- Sorghum$Value[Sorghum$County == cnty & Sorghum$`State ANSI` == st & Sorghum$Year == yr]
    Sr <- ifelse(length(Sr) == 0 || is.na(Sr), 0, Sr / 1e6)
    
    Sy <- Soybean$Value[Soybean$County == cnty & Soybean$`State ANSI` == st & Soybean$Year == yr]
    Sy <- ifelse(length(Sy) == 0 || is.na(Sy), 0, Sy / 1e6)
    
    W  <- Wheat$Value[Wheat$County == cnty & Wheat$`State ANSI` == st & Wheat$Year == yr]
    W  <- ifelse(length(W) == 0 || is.na(W), 0, W / 1e6)
    
    total <- C + Sr + Sy + W
    
    if (total > 0) {
      valid_years[i]     <- valid_years[i] + 1
      prod_matrix[i, t]  <- total
      Ratio_matrix[i, t] <- (C + Sr + Sy) / total
    }
  }
}

# =============================================================================
# 4. FILTER + PAIRWISE OVERLAP CORRELATION
# =============================================================================
observed_years <- rowSums(!is.na(prod_matrix))
sufficient     <- observed_years >= min_years_required
valid_idx      <- which(sufficient)

log_prod <- log(prod_matrix[valid_idx, , drop = FALSE] + 1e-4)
z_prod   <- scale(log_prod)

baseline_global_idx <- which(Shapefile$NAME == toupper(baseline_county) &
                               Shapefile$STATEFP == baseline_statefp)
baseline_valid_idx <- which(valid_idx == baseline_global_idx)

correlations_valid <- numeric(length(valid_idx))
correlations_valid[baseline_valid_idx] <- 1

for (j in seq_along(valid_idx)) {
  if (j == baseline_valid_idx) next
  overlap <- !is.na(prod_matrix[valid_idx[baseline_valid_idx], ]) &
    !is.na(prod_matrix[valid_idx[j], ])
    correlations_valid[j] <- cor(z_prod[baseline_valid_idx, overlap],
                                 z_prod[j, overlap])
}

correlations_all <- rep(NA_real_, N)
correlations_all[valid_idx] <- correlations_valid
length(correlations_all)
# =============================================================================
# 5. SPATIAL DATA — IDENTICAL TO Map_D.R
# =============================================================================
full_data <- data.frame(
  FID          = Shapefile$FID,
  NAME         = Shapefile$NAME,
  X            = Shapefile$X,
  Y            = Shapefile$Y,
  Correlation  = correlations_all,
  Missing_Data = !sufficient,
  stringsAsFactors = FALSE
)

full_sf <- st_as_sf(full_data, coords = c("X", "Y"), crs = 4326)

county_poly <- map_df(c("OK", "KS"), ~ counties(.x, cb = TRUE, resolution = "20m")) %>%
  st_transform(crs = 4326) %>%
  mutate(NAME = toupper(NAME))

map_data <- county_poly %>%
  st_join(full_sf, join = st_nearest_feature) %>%
  left_join(select(full_data, FID, Correlation, Missing_Data), by = "FID") %>%
  mutate(
    Correlation  = coalesce(Correlation.x,  Correlation.y),
    Missing_Data = coalesce(Missing_Data.x, Missing_Data.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Special fill: white for missing data, red for baseline
map_data$Special_Fill <- ifelse(map_data$Missing_Data, "Missing Data", NA)
baseline_fid <- Shapefile$FID[baseline_global_idx]
map_data$Special_Fill[map_data$FID == baseline_fid] <- "Baseline County"

# =============================================================================
# 6. PLOT — EXACTLY LIKE Map_D.R
# =============================================================================
county_label <- tools::toTitleCase(tolower(baseline_county))   # "BARTON" -> "Barton"

p <- ggplot() +
  geom_sf(data = filter(map_data, is.na(Special_Fill)),
          aes(fill = Correlation), color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D",
                       name = paste0("Correlation with\n", county_label, " County (2013–2023)"),
                       breaks = seq(-1, 1, 0.5),
                       limits = c(-1, 1)) +
  new_scale_fill() +
  geom_sf(data = filter(map_data, !is.na(Special_Fill)),
          aes(fill = Special_Fill), color = "black", size = 0.3) +
  scale_fill_manual(values = c("Baseline County" = "red",
                               "Missing Data" = "white"),
                    name = "",
                    labels = c(paste0(county_label, " County (baseline)"), "Missing Data")) +
  geom_sf(data = map_data, fill = NA, color = "black", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.2, unit = "km") +
  coord_sf(xlim = c(-103.0, -94.0), ylim = c(34, 40)) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial"))

print(p)

# Save

ggsave(file.path(output_dir, paste0("Correlation_", county_label, "_Map_Final.jpg")),
       plot = p,
       width = 8.5, height = 6, dpi = 300, units = "in", bg = "white")

cat("\nCorrelation map saved with full Map_D.R style to:\n", output_dir, "\n")

# =============================================================================
# FINAL SUMMARY STATISTICS
# =============================================================================
mean_prod  <- rowMeans(prod_matrix, na.rm = TRUE)
mean_ratio <- rowMeans(Ratio_matrix, na.rm = TRUE)
q1         <- quantile(mean_prod[sufficient], 0.25, na.rm = TRUE)   # bottom 25% among included counties
small_prod <- mean_prod <= q1

cat(sprintf("Total counties in KS + OK (raw Points.xls):               %4d\n", N))
cat(sprintf("Counties with ≥%d years of positive production (included): %4d\n", min_years_required, sum(sufficient)))

cat("BASELINE COUNTY SUMMARY:", baseline_county, "COUNTY,\n")
cat(sprintf("Mean total grain production (2013–2023):      %6.3f million bushels\n",
            mean_prod[baseline_global_idx]))
cat(sprintf("Mean crop diversity ratio (ratio_r):          %6.3f\n",
            mean_ratio[baseline_global_idx]))
cat(sprintf("Small production county (bottom 25%%)?         %s\n",
            ifelse(small_prod[baseline_global_idx], "YES", "NO")))
cat(sprintf("Valid production years (out of 11):           %2d years\n",
            valid_years[baseline_global_idx]))

