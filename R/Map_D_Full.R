library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(tigris)
library(purrr)
library(ggnewscale)
library(viridis)
library(ggspatial)

# Define paths
main_dir <- "data"   # input folder, relative to the repository root
output_dir <- "figures"
missed <-2

set.seed(456)

# Load and aggregate crop data — DO NOT force NA → 0
Corn <- read_excel(file.path(main_dir, "Data.xlsx"), sheet = "Corn") %>%
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

Wheat <- read_excel(file.path(main_dir, "Data.xlsx"), sheet = "Wheat") %>%
  mutate(Value = as.numeric(as.character(Value))) %>%
  group_by(County, `State ANSI`, Year) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

# Load point file (county centroids)
Shapefile <- read_excel(file.path(main_dir, "Points.xls"))
Shapefile$NAME <- toupper(Shapefile$NAME)
Shapefile <- Shapefile %>% filter(STATEFP %in% c(40, 20))  # OK=40, KS=20

years <- 2013:2023
N <- nrow(Shapefile)

# Initialize matrices with NA (critical change)
prod_matrix  <- matrix(NA, nrow = N, ncol = length(years), dimnames = list(1:N, years))
Ratio_matrix <- matrix(NA, nrow = N, ncol = length(years), dimnames = list(1:N, years))
valid_years_count <- integer(N)  # Count of years with positive production

# Fill matrices — only record data when total production > 0
for (i in 1:N) {
  County <- Shapefile$NAME[i]
  State  <- Shapefile$STATEFP[i]
  
  for (t in seq_along(years)) {
    yr <- years[t]
    
    C  <- Corn$Value[Corn$County == County & Corn$`State ANSI` == State & Corn$Year == yr]
    C  <- ifelse(length(C) == 0 || is.na(C), 0, C / 1e6)
    
    Sr <- Sorghum$Value[Sorghum$County == County & Sorghum$`State ANSI` == State & Sorghum$Year == yr]
    Sr <- ifelse(length(Sr) == 0 || is.na(Sr), 0, Sr / 1e6)
    
    Sy <- Soybean$Value[Soybean$County == County & Soybean$`State ANSI` == State & Soybean$Year == yr]
    Sy <- ifelse(length(Sy) == 0 || is.na(Sy), 0, Sy / 1e6)
    
    W  <- Wheat$Value[Wheat$County == County & Wheat$`State ANSI` == State & Wheat$Year == yr]
    W  <- ifelse(length(W) == 0 || is.na(W), 0, W / 1e6)
    
    total <- C + Sr + Sy + W
    
    if (total>0  ) {
      valid_years_count[i] <- valid_years_count[i] + 1
      prod_matrix[i, t]  <- total 
      Ratio_matrix[i, t] <- (C + Sr + Sy)/ total
    }
  }
}

# Define missing-data counties: < 2 years with positive production
missing_data <- valid_years_count < missed

# Compute means only over observed years
mean_prod  <- rowMeans(prod_matrix, na.rm = TRUE)
mean_prod[is.nan(mean_prod)] <- NA
mean_ratio <- rowMeans(Ratio_matrix, na.rm = TRUE)
mean_ratio[is.nan(mean_ratio)] <- NA

# Small production: bottom 25% among counties with ≥6 valid years
included <- !missing_data & !is.na(mean_prod)
prod_mean_included <- mean_prod[included]
q1 <- quantile(prod_mean_included, 0.25, na.rm = TRUE)

small_prod <- rep(FALSE, N)
small_prod[included] <- mean_prod[included] <= q1

# County identifiers
IDs <- 1:N
county_names <- Shapefile$NAME[IDs]
county_fids  <- Shapefile$FID[IDs]
X <- Shapefile$X[IDs]
Y <- Shapefile$Y[IDs]

# Fetch county polygons
options(tigris_use_cache = TRUE)
county_polygons <- map_df(c("OK", "KS"), ~ counties(.x, cb = TRUE, resolution = "20m")) %>%
  st_transform(crs = 4326)

# Prepare point-based county data
county_data <- data.frame(
  NAME = county_names,
  FID = county_fids,
  X = X,
  Y = Y,
  Missing_Data = missing_data,
  stringsAsFactors = FALSE
)

county_sf <- st_as_sf(county_data, coords = c("X", "Y"), crs = 4326)

# Spatial join: match polygons to nearest point (by FID)
county_polygons$NAME <- toupper(county_polygons$NAME)
map_data <- county_polygons %>%
  st_join(county_sf, join = st_nearest_feature) %>%
  left_join(select(county_data, FID, NAME), by = "FID") %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Attach computed variables
map_data$Diversity_Ratio  <- mean_ratio[match(map_data$FID, county_fids)]
map_data$Total_Production <- mean_prod[match(map_data$FID, county_fids)]
map_data$Missing_Data     <- missing_data[match(map_data$FID, county_fids)]
map_data$Small_Prod       <- small_prod[match(map_data$FID, county_fids)]

cat("Counties with <6 years of positive production (shown white):", sum(map_data$Missing_Data, na.rm = TRUE), "\n")

# Unified Special_Fill for all maps
map_data$Special_Fill <- ifelse(map_data$Missing_Data, "Missing Data", NA)

# Map 1: Diversity Ratio
p_diversity <- ggplot() +
  geom_sf(data = filter(map_data, is.na(Special_Fill)),
          aes(fill = Diversity_Ratio), color = "black", size = 0.3) +
  scale_fill_viridis(option = "D", name = "Crop Diversity Ratio\n(2013–2023)", 
                     breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  new_scale_fill() +
  geom_sf(data = filter(map_data, !is.na(Special_Fill)),
          aes(fill = Special_Fill), color = "black", size = 0.3) +
  scale_fill_manual(values = c("Missing Data" = "white"),
                    name = "", labels = "Missing Data") +
  geom_sf(data = map_data, fill = NA, color = "black", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.2, unit = "km") +
  coord_sf(xlim = c(-103.0, -94.0), ylim = c(34, 40)) +
  theme_minimal() +
  theme(legend.position = "right", panel.grid = element_blank(),
        axis.title = element_blank(), text = element_text(family = "Arial"))

ggsave(file.path(output_dir, "diversity_ratio_map.jpg"), p_diversity,
       width = 8.5, height = 6, dpi = 300, units = "in")

# Map 2: Total Production
p_total <- ggplot() +
  geom_sf(data = filter(map_data, is.na(Special_Fill)),
          aes(fill = Total_Production), color = "black", size = 0.3) +
  scale_fill_viridis(option = "D", name = "Mean Total Grain Production\n(Million Bushels, 2013–2023)",
                     breaks = seq(0, 30, 5), limits = c(0, 30)) +
  new_scale_fill() +
  geom_sf(data = filter(map_data, !is.na(Special_Fill)),
          aes(fill = Special_Fill), color = "black", size = 0.3) +
  scale_fill_manual(values = c("Missing Data" = "white"),
                    name = "", labels = "Missing Data") +
  geom_sf(data = map_data, fill = NA, color = "black", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.2, unit = "km") +
  coord_sf(xlim = c(-103.0, -94.0), ylim = c(34, 40)) +
  theme_minimal() +
  theme(legend.position = "right", panel.grid = element_blank(),
        axis.title = element_blank(), text = element_text(family = "Arial"))

ggsave(file.path(output_dir, "total_production_map.jpg"), p_total,
       width = 8.5, height = 6, dpi = 300, units = "in")

# Map 3: Small vs Large Producers
map_data$Fill_Cat <- ifelse(map_data$Missing_Data, "Missing Data",
                            ifelse(map_data$Small_Prod, "Small Production", "Large Production"))

p_small <- ggplot() +
  geom_sf(data = map_data, aes(fill = Fill_Cat), color = "black", size = 0.3) +
  scale_fill_manual(values = c("Small Production" = "yellow",
                               "Large Production" = "darkblue",
                               "Missing Data" = "white"),
                    name = "Production Category",
                    labels = c("Small Production" = "Small Production (Bottom 25%)",
                               "Large Production" = "Large Production",
                               "Missing Data" = "Missing Data")) +
  geom_sf(data = map_data, fill = NA, color = "black", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.2, unit = "km") +
  coord_sf(xlim = c(-103.0, -94.0), ylim = c(34, 40)) +
  theme_minimal() +
  theme(legend.position = "right", panel.grid = element_blank(),
        axis.title = element_blank(), text = element_text(family = "Arial"))

ggsave(file.path(output_dir, "small_production_map.jpg"), p_small,
       width = 8.5, height = 6, dpi = 300, units = "in")

# Done
cat("All three maps saved to:", output_dir, "\n")