#!/usr/bin/env Rscript

# ============================================================
# Validation of modeled coal PM2.5 against AQS observations
# (1980 and 1990)
# ============================================================
#
# Purpose
#   1. Match AQS monitoring sites to the nearest HyADS grid cell.
#   2. Evaluate the spatial agreement between modeled coal PM2.5
#      and observed particulate measurements.
#   3. Generate the evaluation statistics used in manuscript
#      Table 1 and Figure 5.
#
# Main manuscript outputs
#
#   Table 1
#     Spatial agreement between modeled coal PM2.5 and
#     observed sulfate concentrations.
#       - 1980 Sulfate (TSP)
#       - 1990 Sulfate (TSP)
#
#   Figure 5
#       (a) 1980 Sulfate (TSP) vs modeled coal PM2.5
#       (b) 1990 Sulfate (TSP) vs modeled coal PM2.5
#
# Supplementary outputs
#
#   Figures
#       S1  1990 PM2.5 mass
#       S2  1990 PM10
#       S3  1990 Total Suspended Particulate (TSP)
#       S4  1980 Total Suspended Particulate (TSP)
#
#   Supplementary Table
#       Complete evaluation statistics for PM2.5, PM10, TSP,
#       and Sulfate (TSP).
#
# Evaluation statistics
#   - Number of monitoring sites
#   - Pearson correlation
#   - Spearman correlation
#   - Coefficient of determination (R2)
#   - Mean bias (Modeled - Observed)
#   - Mean absolute error (MAE)
#   - Root mean squared error (RMSE)
#
# Notes
#   Coal PM2.5 is defined as modeled secondary sulfate PM2.5
#   attributable to SO2 emissions from coal-fired power plants.
#
#   AQS monitoring sites are matched to the nearest HyADS grid cell.
#
#   Scatter plots display:
#       x-axis = Observed annual mean concentration (AQS)
#       y-axis = Modeled decade-mean coal PM2.5 (HyADS)
#
#   Distances greater than MAX_DIST_M are excluded from the
#   evaluation to avoid unrealistic monitor-grid matches.
#
# Output files (out_dir)
#
#   Tables
#       evaluation_summary_all_metrics_distfiltered.csv
#       Table1_main_sulfate_only.csv
#       TableS_AQS_all_metrics_distfiltered.csv
#
#   Main Figure
#       Fig5_main_sulfate_1980_1990_polished.pdf
#
#   Individual manuscript panels
#       Fig5a_1980_Sulfate_TSP_polished.pdf
#       Fig5b_1990_Sulfate_TSP_polished.pdf
#
#   Supplementary figures
#       Scatter_1980_TSP_mass_distfiltered.pdf
#       Scatter_1990_PM25_mass_distfiltered.pdf
#       Scatter_1990_PM10_mass_distfiltered.pdf
#       Scatter_1990_TSP_mass_distfiltered.pdf
#
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(ggplot2)
})

# -----------------------------
# User inputs
# -----------------------------

p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

base_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"
eval_dir <- file.path(base_dir, "evaluation")

out_dir <- file.path(eval_dir, "aqs_hyads_eval_revised")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

MAX_DIST_M <- 100000  # 100 km; for sensitivity testing, try 50000

aqs_files <- list(
  `1980` = file.path(eval_dir, "annual_conc_by_monitor_1980.csv"),
  `1990` = file.path(eval_dir, "annual_conc_by_monitor_1990.csv")
)

hyads_files <- list(
  `1980` = file.path(base_dir, "grids_pm25_total_1980.fst"),
  `1990` = file.path(base_dir, "grids_pm25_total_1990.fst")
)

stopifnot(all(file.exists(unlist(aqs_files))))
stopifnot(all(file.exists(unlist(hyads_files))))

cat("Output directory:\n", out_dir, "\n\n")
cat("Distance threshold:", MAX_DIST_M / 1000, "km\n\n")

# -----------------------------
# Helper functions
# -----------------------------

read_aqs_annual <- function(f) {
  dt <- fread(f)
  setnames(dt, gsub(" ", "_", names(dt)))

  req <- c(
    "State_Code", "County_Code", "Site_Num",
    "Latitude", "Longitude",
    "Parameter_Name", "Parameter_Code",
    "Arithmetic_Mean", "Year"
  )
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    stop("AQS file missing columns: ", paste(miss, collapse = ", "))
  }

  dt[, State_Code := as.integer(State_Code)]
  dt[, County_Code := as.integer(County_Code)]
  dt[, Site_Num := as.integer(Site_Num)]
  dt[, Latitude := as.numeric(Latitude)]
  dt[, Longitude := as.numeric(Longitude)]
  dt[, Arithmetic_Mean := as.numeric(Arithmetic_Mean)]
  dt[, Year := as.integer(Year)]

  dt
}

read_hyads_grid <- function(f, p4s) {
  dt <- fst::read_fst(f, as.data.table = TRUE)

  req <- c("x", "y", "vals.out")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    stop("HyADS fst missing columns: ", paste(miss, collapse = ", "))
  }

  dt[, x := as.numeric(x)]
  dt[, y := as.numeric(y)]
  dt[, vals.out := as.numeric(vals.out)]

  st_as_sf(dt, coords = c("x", "y"), crs = st_crs(p4s), remove = FALSE)
}

make_site_id <- function(state_code, county_code, site_num) {
  paste0(
    sprintf("%02d", as.integer(state_code)), "_",
    sprintf("%03d", as.integer(county_code)), "_",
    sprintf("%04d", as.integer(site_num))
  )
}

to_monitor_sf <- function(dt) {
  if (nrow(dt) == 0) {
    stop("No rows passed to to_monitor_sf(). Check AQS parameter filter.")
  }

  pm <- copy(dt)
  pm[, site_id := make_site_id(State_Code, County_Code, Site_Num)]

  pm <- pm[
    is.finite(Longitude) &
      is.finite(Latitude) &
      is.finite(Arithmetic_Mean)
  ]

  # Collapse repeated records to one observation per monitoring site.
  pm_site <- pm[
    ,
    .(
      lon = mean(Longitude, na.rm = TRUE),
      lat = mean(Latitude, na.rm = TRUE),
      obs = mean(Arithmetic_Mean, na.rm = TRUE),
      year = Year[1],
      Parameter_Name = Parameter_Name[1],
      Parameter_Code = Parameter_Code[1],
      n_aqs_records = .N
    ),
    by = site_id
  ]

  pm_site <- pm_site[
    is.finite(lon) &
      is.finite(lat) &
      is.finite(obs)
  ]

  st_as_sf(pm_site, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

match_nearest_hyads <- function(obs_sf_ll, hyads_sf, p4s) {
  obs_sf <- st_transform(obs_sf_ll, st_crs(p4s))
  idx <- st_nearest_feature(obs_sf, hyads_sf)

  obs_sf$mod <- hyads_sf$vals.out[idx]
  obs_sf$hyads_x <- hyads_sf$x[idx]
  obs_sf$hyads_y <- hyads_sf$y[idx]
  obs_sf$match_dist_m <- as.numeric(
    st_distance(obs_sf, hyads_sf[idx, ], by_element = TRUE)
  )

  obs_sf
}

calc_metrics <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  obs <- obs[ok]
  mod <- mod[ok]

  n <- length(obs)
  if (n < 10) {
    return(list(
      n_obs = n,
      mean_observed = NA_real_,
      mean_modeled = NA_real_,
      mean_bias = NA_real_,
      mae = NA_real_,
      rmse = NA_real_,
      normalized_mean_bias = NA_real_,
      pearson_r = NA_real_,
      spearman_rho = NA_real_,
      r2 = NA_real_,
      lm_intercept = NA_real_,
      lm_slope = NA_real_
    ))
  }

  diff <- mod - obs
  pearson_r <- suppressWarnings(cor(obs, mod, method = "pearson"))
  spearman_rho <- suppressWarnings(cor(obs, mod, method = "spearman"))
  lm_fit <- suppressWarnings(lm(mod ~ obs))
  lm_coef <- coef(lm_fit)

  list(
    n_obs = n,
    mean_observed = mean(obs, na.rm = TRUE),
    mean_modeled = mean(mod, na.rm = TRUE),
    mean_bias = mean(diff, na.rm = TRUE),
    mae = mean(abs(diff), na.rm = TRUE),
    rmse = sqrt(mean(diff^2, na.rm = TRUE)),
    normalized_mean_bias = sum(diff, na.rm = TRUE) / sum(obs, na.rm = TRUE),
    pearson_r = pearson_r,
    spearman_rho = spearman_rho,
    r2 = pearson_r^2,
    lm_intercept = unname(lm_coef[1]),
    lm_slope = unname(lm_coef[2])
  )
}

safe_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# -----------------------------
# Evaluation specification
# -----------------------------
# Important: "PM2.5 STP" is the AQS Parameter_Name.
# In the manuscript, report it as PM2.5 mass, not TSP.

eval_specs <- list(
  list(
    pair_id = "1990_PM25_mass",
    year = 1990,
    match_type = "exact",
    parameter_name = "PM2.5 STP",
    parameter_pattern = NA_character_,
    aqs_metric = "PM2.5 mass",
    size_fraction = "PM2.5",
    aqs_parameter_name_for_table = "PM2.5 STP"
  ),
  list(
    pair_id = "1990_PM10_mass",
    year = 1990,
    match_type = "exact",
    parameter_name = "PM10 Total 0-10um STP",
    parameter_pattern = NA_character_,
    aqs_metric = "PM10 mass",
    size_fraction = "PM10",
    aqs_parameter_name_for_table = "PM10 Total 0-10um STP"
  ),
  list(
    pair_id = "1990_TSP_mass",
    year = 1990,
    match_type = "exact",
    parameter_name = "Suspended particulate (TSP)",
    parameter_pattern = NA_character_,
    aqs_metric = "total suspended particulate matter",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Suspended particulate (TSP)"
  ),
  list(
    pair_id = "1990_Sulfate_TSP",
    year = 1990,
    match_type = "regex",
    parameter_name = NA_character_,
    parameter_pattern = "^Sulfate \\(TSP\\)",
    aqs_metric = "sulfate",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Sulfate (TSP)"
  ),
  list(
    pair_id = "1980_TSP_mass",
    year = 1980,
    match_type = "exact",
    parameter_name = "Suspended particulate (TSP)",
    parameter_pattern = NA_character_,
    aqs_metric = "total suspended particulate matter",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Suspended particulate (TSP)"
  ),
  list(
    pair_id = "1980_Sulfate_TSP",
    year = 1980,
    match_type = "regex",
    parameter_name = NA_character_,
    parameter_pattern = "^Sulfate \\(TSP\\)",
    aqs_metric = "sulfate",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Sulfate (TSP)"
  )
)

# -----------------------------
# Plot functions
# -----------------------------

make_scatter_plot <- function(plot_dt, spec) {
  m <- calc_metrics(plot_dt$obs, plot_dt$mod)

  ann <- sprintf(
    "n = %d\nPearson r = %.2f\nSpearman rho = %.2f\nR2 = %.2f",
    m$n_obs, m$pearson_r, m$spearman_rho, m$r2
  )

  xlab_txt <- paste0(
    "Observed annual mean ", spec$aqs_metric,
    " (", spec$size_fraction, ", AQS, ", spec$year, "; 
