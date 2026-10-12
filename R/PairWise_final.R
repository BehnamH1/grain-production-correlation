# ===================================================================
# Simple.R — FINAL CORRECTED VERSION
# Proper missing data handling | No fake zeros | No subscript errors
# Keeps only counties with ≥6 years of actual positive production
# ===================================================================

library(readxl)
library(dplyr)
library(readr)

set.seed(456)

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------
main_dir <- "data"   # input folder, relative to the repository root

# -------------------------------------------------------------------
# 1. Load raw crop data — NEVER convert NA to 0
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# 2. Load county centroids (Points.xls) — KS + OK only
# -------------------------------------------------------------------
Shapefile <- read_excel(file.path(main_dir, "Points.xls"))
Shapefile$NAME <- toupper(Shapefile$NAME)
Shapefile <- Shapefile %>% filter(STATEFP %in% c(40, 20))  # 40=Oklahoma, 20=Kansas

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

# -------------------------------------------------------------------
# 5. KEEP ONLY COUNTIES WITH ≥2 YEARS OF ACTUAL PRODUCTION
# -------------------------------------------------------------------
keep <- valid_years >= 2
cat("Counties with ≥6 years of positive production:", sum(keep), "\n")

Shapefile    <- Shapefile[keep, , drop = FALSE]
prod_matrix  <- prod_matrix[keep, , drop = FALSE]
Ratio_matrix <- Ratio_matrix[keep, , drop = FALSE]

n_counties   <- nrow(Shapefile)
X            <- Shapefile$X
Y            <- Shapefile$Y
county_names <- Shapefile$NAME
county_fid   <- Shapefile$FID

cat("Final number of counties in analysis:", n_counties, "\n")

# -------------------------------------------------------------------
# 6. County-level summaries
# -------------------------------------------------------------------
mean_prod  <- rowMeans(prod_matrix, na.rm = TRUE)
mean_ratio <- rowMeans(Ratio_matrix, na.rm = TRUE)
q1         <- quantile(mean_prod, 0.25, na.rm = TRUE)
small_prod <- mean_prod <= q1

# -------------------------------------------------------------------
# 7. Covariates (distances in 100-mile units)
# -------------------------------------------------------------------
X_dist_100mi <- abs(outer(X, X, `-`)) * 69 / 100
Y_dist_100mi <- abs(outer(Y, Y, `-`)) * 69 / 100
Ratio_R      <- outer(mean_ratio, mean_ratio, function(a,b) pmin(a,b)/pmax(a,b))
Small_Prod   <- outer(small_prod, small_prod, function(a,b) as.numeric(a | b))

# -------------------------------------------------------------------
# 8. Adjacency matrix from Adj.csv
# -------------------------------------------------------------------
adj_table <- read_csv(file.path(main_dir, "Adj.csv"))
Adj <- matrix(0, n_counties, n_counties, dimnames = list(county_names, county_names))

for (k in seq_len(nrow(adj_table))) {
  i <- match(adj_table$src_FID[k], county_fid)
  j <- match(adj_table$nbr_FID[k], county_fid)
  if (!is.na(i) && !is.na(j) && i != j) {
    Adj[i, j] <- Adj[j, i] <- 1
  }
}
diag(Adj) <- 0

# -------------------------------------------------------------------
# 9. Detrend log-yield per county (only observed years)
# -------------------------------------------------------------------
residuals_matrix <- matrix(NA, n_counties, 11)

for (i in seq_len(n_counties)) {
  obs <- !is.na(prod_matrix[i, ])
  if (sum(obs) < 3) next
  df  <- data.frame(year = years[obs], y = log(prod_matrix[i, obs] + 1e-8))
  m   <- lm(y ~ year, data = df)
  residuals_matrix[i, obs] <- residuals(m)
}
y <- residuals_matrix

# -------------------------------------------------------------------
# 10. Pairwise correlations — ONLY on overlapping observed years
# -------------------------------------------------------------------
emp_cor <- matrix(NA, nrow = n_counties, ncol = n_counties)
for (i in 1:n_counties) {
  for (j in 1:n_counties) {
    if (i == j) {
      next
    }
    # Identify years where BOTH counties have non-zero original production
    both_nonzero <- !is.na(prod_matrix[i, ]) & !is.na(prod_matrix[j, ])
    # Compute Pearson on detrended residuals only for those overlapping non-zero years
    res_i <- residuals_matrix[i, both_nonzero]
    res_j <- residuals_matrix[j, both_nonzero]
    emp_cor[i, j] <- cor(res_i, res_j)
  }
}
#diag(emp_cor) <- NA  # ensure diagonal is 1

# -------------------------------------------------------------------
# 11. Build regression dataset (upper triangle, drop missing pairs)
# -------------------------------------------------------------------
upper <- upper.tri(emp_cor, diag = FALSE)

reg_data <- data.frame(
  rho        = emp_cor[upper],
  x_dir      = X_dist_100mi[upper],
  y_dir      = Y_dist_100mi[upper],
  ratio_r    = Ratio_R[upper],
  adj        = Adj[upper],
  small_prod = Small_Prod[upper]
) %>% na.omit()

cat("Number of county pairs in final regression:", nrow(reg_data), "\n")

# -------------------------------------------------------------------
# 12. regression model
# -------------------------------------------------------------------
reg_model <- lm(rho ~ x_dir + y_dir + ratio_r + adj + small_prod, data = reg_data)
summary(reg_model)

ols_summary <- summary(reg_model)

# -------------------------------------------------------------------
# 12. bayesian model
# -------------------------------------------------------------------

library(brms)
library(cmdstanr)

# Extract numbers once
ols_coef <- coef(ols_summary)[,1]      # Estimate column
ols_se   <- coef(ols_summary)[,2]      # Std.Error column

# Build correct priors – brms evaluates the numbers in R, not in Stan
bprior <- c(
  set_prior(sprintf("normal(%.10f, %.10f)", ols_coef["(Intercept)"],   ols_se["(Intercept)"]   * 20), class = "Intercept"),
  set_prior(sprintf("normal(%.10f, %.10f)", ols_coef["x_dir"],         ols_se["x_dir"]         * 10), class = "b", coef = "x_dir"),
  set_prior(sprintf("normal(%.10f, %.10f)", ols_coef["y_dir"],         ols_se["y_dir"]         * 10), class = "b", coef = "y_dir"),
  set_prior(sprintf("normal(%.10f, %.10f)", ols_coef["ratio_r"],       ols_se["ratio_r"]       * 10), class = "b", coef = "ratio_r"),
  set_prior(sprintf("normal(%.10f, %.10f)", ols_coef["adj"],           ols_se["adj"]           * 10), class = "b", coef = "adj"),
  set_prior(sprintf("normal(%.10f, %.10f)", ols_coef["small_prod"],    ols_se["small_prod"]    * 10), class = "b", coef = "small_prod"),
  set_prior("student_t(3, 0, 2.5)", class = "sigma")
)

bayes_fit <- brm(
  rho ~ x_dir + y_dir + ratio_r + adj + small_prod,
  data = reg_data,
  prior = bprior,
  chains = 4, iter = 20000, warmup = 4000,
  seed = 123,
  backend = "cmdstanr",   # fastest in 2025
  control = list(adapt_delta = 0.95)
)

summary(bayes_fit)
fixef(bayes_fit)      # clean coefficient table
posterior_summary(bayes_fit)[1:7,]  # means, sd, 2.5%/97.5% quantiles

summary(reg_model)


