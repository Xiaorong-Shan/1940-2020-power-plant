#!/usr/bin/env Rscript

# ============================================================
# Validation of modeled coal PM2.5 against AQS observations
# (1980 and 1990)
#
# Purpose
#   1. Match AQS monitoring sites to the nearest HyADS grid cell.
#   2. Evaluate spatial agreement between modeled coal PM2.5 and
#      observed particulate measurements.
#   3. Generate manuscript Table 1, Figure 5, and supplementary outputs.
#
# Main manuscript outputs
#   Table 1
#     Spatial agreement between modeled coal PM2.5 and observed sulfate
#     concentrations:
#       - 1980 Sulfate (TSP)
#       - 1990 Sulfate (TSP)
#
#   Figure 5
#     Monitor-level comparisons for:
#       - 1980 Sulfate (TSP) vs modeled coal PM2.5
#       - 1990 Sulfate (TSP) vs modeled coal PM2.5
#     Note: panel labels such as (a) and (b) are NOT drawn inside the plots.
#
# Supplementary outputs
#   Complete evaluation statistics for PM2.5, PM10, TSP, and sulfate.
#
# Notes
#   Coal PM2.5 is modeled secondary sulfate PM2.5 attributable to SO2
#   emissions from coal-fired power plants.
#
#   Scatter plots display:
#     x-axis = observed annual mean concentration (AQS)
#     y-axis = modeled decade-mean coal PM2.5 (HyADS)
#
#   Monitoring sites farther than MAX_DIST_M from the nearest HyADS grid cell
#   are excluded to avoid unrealistic monitor-grid matches.
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

out_dir <- file.path(eval_dir, "aqs_hyads_eval_final_no_panel_labels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

MAX_DIST_M <- 100000  # 100 km

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

# -----------------------------
# Helper functions
# -----------------------------
read_aqs_annual <- function(f) {
  dt <- fread(f)
  setnames(dt, gsub(" ", "_", names(dt)))

  req <- c(
    "State_Code", "County_Code", "Site_Num", "Latitude", "Longitude",
    "Parameter_Name", "Parameter_Code", "Arithmetic_Mean", "Year"
  )
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("AQS file missing columns: ", paste(miss, collapse = ", "))

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
  if (length(miss) > 0) stop("HyADS fst missing columns: ", paste(miss, collapse = ", "))

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
  if (nrow(dt) == 0) stop("No rows passed to to_monitor_sf(). Check AQS parameter filter.")

  pm <- copy(dt)
  pm[, site_id := make_site_id(State_Code, County_Code, Site_Num)]
  pm <- pm[is.finite(Longitude) & is.finite(Latitude) & is.finite(Arithmetic_Mean)]

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

  pm_site <- pm_site[is.finite(lon) & is.finite(lat) & is.finite(obs)]
  st_as_sf(pm_site, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

match_nearest_hyads <- function(obs_sf_ll, hyads_sf, p4s) {
  obs_sf <- st_transform(obs_sf_ll, st_crs(p4s))
  idx <- st_nearest_feature(obs_sf, hyads_sf)

  obs_sf$mod <- hyads_sf$vals.out[idx]
  obs_sf$hyads_x <- hyads_sf$x[idx]
  obs_sf$hyads_y <- hyads_sf$y[idx]
  obs_sf$match_dist_m <- as.numeric(st_distance(obs_sf, hyads_sf[idx, ], by_element = TRUE))
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
      r2 = NA_real_
    ))
  }

  diff <- mod - obs
  pearson_r <- suppressWarnings(cor(obs, mod, method = "pearson"))
  spearman_rho <- suppressWarnings(cor(obs, mod, method = "spearman"))

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
    r2 = pearson_r^2
  )
}

safe_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

# -----------------------------
# Plot functions
# -----------------------------
make_scatter_plot <- function(plot_dt, spec) {
  m <- calc_metrics(plot_dt$obs, plot_dt$mod)
  
  ann <- sprintf(
    "n = %d\nPearson r = %.2f\nSpearman r = %.2f\nR² = %.2f",
    m$n_obs, m$pearson_r, m$spearman_rho, m$r2
  )
  
  # Short axis labels to avoid truncation in PDF/RStudio viewer
  xlab_expr <- expression(
    paste("Observed sulfate (", mu, "g ", m^{-3}, ")")
  )
  
  ylab_expr <- expression(
    paste("Modeled coal PM"[2.5], " (", mu, "g ", m^{-3}, ")")
  )
  
  ggplot(plot_dt, aes(x = obs, y = mod)) +
    geom_point(alpha = 0.45, size = 1.35) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, formula = y ~ x) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.15,
      label = ann,
      size = 3.8,
      lineheight = 0.95
    ) +
    labs(
      x = xlab_expr,
      y = ylab_expr
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# -----------------------------
# Evaluation specifications
# -----------------------------
eval_specs <- list(
  list(
    pair_id = "1980_Sulfate_TSP",
    year = 1980,
    match_type = "regex",
    parameter_name = NA_character_,
    parameter_pattern = "^Sulfate \\(TSP\\)",
    aqs_metric = "sulfate",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Sulfate (TSP)",
    main_figure = TRUE
  ),
  list(
    pair_id = "1990_Sulfate_TSP",
    year = 1990,
    match_type = "regex",
    parameter_name = NA_character_,
    parameter_pattern = "^Sulfate \\(TSP\\)",
    aqs_metric = "sulfate",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Sulfate (TSP)",
    main_figure = TRUE
  ),
  list(
    pair_id = "1980_TSP_mass",
    year = 1980,
    match_type = "exact",
    parameter_name = "Suspended particulate (TSP)",
    parameter_pattern = NA_character_,
    aqs_metric = "total suspended particulate matter",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Suspended particulate (TSP)",
    main_figure = FALSE
  ),
  list(
    pair_id = "1990_PM25_mass",
    year = 1990,
    match_type = "exact",
    parameter_name = "PM2.5 STP",
    parameter_pattern = NA_character_,
    aqs_metric = "PM2.5 mass",
    size_fraction = "PM2.5",
    aqs_parameter_name_for_table = "PM2.5 STP",
    main_figure = FALSE
  ),
  list(
    pair_id = "1990_PM10_mass",
    year = 1990,
    match_type = "exact",
    parameter_name = "PM10 Total 0-10um STP",
    parameter_pattern = NA_character_,
    aqs_metric = "PM10 mass",
    size_fraction = "PM10",
    aqs_parameter_name_for_table = "PM10 Total 0-10um STP",
    main_figure = FALSE
  ),
  list(
    pair_id = "1990_TSP_mass",
    year = 1990,
    match_type = "exact",
    parameter_name = "Suspended particulate (TSP)",
    parameter_pattern = NA_character_,
    aqs_metric = "total suspended particulate matter",
    size_fraction = "TSP",
    aqs_parameter_name_for_table = "Suspended particulate (TSP)",
    main_figure = FALSE
  )
)

# -----------------------------
# Load data
# -----------------------------
cat("Loading AQS files...\n")
aqs_dt <- lapply(aqs_files, read_aqs_annual)

cat("Loading HyADS grid files...\n")
hyads_sf <- lapply(hyads_files, read_hyads_grid, p4s = p4s)

# -----------------------------
# Run evaluations
# -----------------------------
summary_list <- list()
plot_data_list <- list()
plot_list <- list()

for (spec in eval_specs) {
  cat("\nEvaluating:", spec$pair_id, "\n")
  yy <- as.character(spec$year)
  dt_year <- aqs_dt[[yy]]
  
  if (spec$match_type == "exact") {
    dt_sub <- dt_year[Parameter_Name == spec$parameter_name]
  } else if (spec$match_type == "regex") {
    dt_sub <- dt_year[grepl(spec$parameter_pattern, Parameter_Name)]
  } else {
    stop("Unknown match_type: ", spec$match_type)
  }
  
  cat("  Raw AQS rows:", nrow(dt_sub), "\n")
  if (nrow(dt_sub) == 0) next
  
  obs_sf <- to_monitor_sf(dt_sub)
  matched_sf <- match_nearest_hyads(obs_sf, hyads_sf[[yy]], p4s)
  
  plot_dt <- as.data.table(st_drop_geometry(matched_sf))
  plot_dt <- plot_dt[is.finite(obs) & is.finite(mod)]
  
  cat("  Before distance filter:", nrow(plot_dt), "\n")
  cat("  Match distance summary:\n")
  print(summary(plot_dt$match_dist_m))
  
  plot_dt <- plot_dt[match_dist_m <= MAX_DIST_M]
  cat("  After distance filter:", nrow(plot_dt), "\n")
  
  plot_dt[, pair_id := spec$pair_id]
  plot_dt[, aqs_metric := spec$aqs_metric]
  plot_dt[, size_fraction := spec$size_fraction]
  
  plot_data_list[[spec$pair_id]] <- plot_dt
  
  fwrite(
    plot_dt,
    file.path(out_dir, paste0("plot_data_", safe_file_name(spec$pair_id), "_distfiltered.csv"))
  )
  
  m <- calc_metrics(plot_dt$obs, plot_dt$mod)
  
  summary_list[[spec$pair_id]] <- data.table(
    pair_id = spec$pair_id,
    year = spec$year,
    aqs_metric = spec$aqs_metric,
    size_fraction = spec$size_fraction,
    aqs_parameter_name = spec$aqs_parameter_name_for_table,
    n_sites = uniqueN(plot_dt$site_id),
    n_obs = m$n_obs,
    mean_observed = m$mean_observed,
    mean_modeled_coal_pm25 = m$mean_modeled,
    mean_bias = m$mean_bias,
    mae = m$mae,
    rmse = m$rmse,
    normalized_mean_bias = m$normalized_mean_bias,
    pearson_r = m$pearson_r,
    spearman_rho = m$spearman_rho,
    r2 = m$r2,
    median_match_dist_m = median(plot_dt$match_dist_m, na.rm = TRUE),
    max_match_dist_m = max(plot_dt$match_dist_m, na.rm = TRUE)
  )
  
  p <- make_scatter_plot(plot_dt, spec)
  plot_list[[spec$pair_id]] <- p
  
  ggsave(
    file.path(out_dir, paste0("Scatter_", safe_file_name(spec$pair_id), "_distfiltered.pdf")),
    p,
    width = 6.8,
    height = 5.2
  )
}

summary_dt <- rbindlist(summary_list, fill = TRUE)
setorder(summary_dt, year, aqs_metric, size_fraction)

fwrite(summary_dt, file.path(out_dir, "evaluation_summary_all_metrics_distfiltered.csv"))
print(summary_dt)

# -----------------------------
# Main Table 1: sulfate only
# -----------------------------
table1_main <- summary_dt[
  aqs_metric == "sulfate" & size_fraction == "TSP",
  .(
    year,
    observed_metric = "Sulfate",
    size_fraction,
    n_sites,
    pearson_r = round(pearson_r, 2),
    spearman_r = round(spearman_rho, 2),
    r2 = round(r2, 2)
  )
][order(year)]

fwrite(table1_main, file.path(out_dir, "Table1_main_sulfate_only.csv"))
print(table1_main)

# -----------------------------
# SI Table: all metrics
# -----------------------------
table_si <- summary_dt[
  ,
  .(
    year,
    aqs_metric,
    size_fraction,
    n_sites,
    mean_observed = round(mean_observed, 3),
    mean_modeled_coal_pm25 = round(mean_modeled_coal_pm25, 3),
    mean_bias = round(mean_bias, 3),
    mae = round(mae, 3),
    rmse = round(rmse, 3),
    pearson_r = round(pearson_r, 3),
    spearman_r = round(spearman_rho, 3),
    r2 = round(r2, 3),
    median_match_dist_km = round(median_match_dist_m / 1000, 1),
    max_match_dist_km = round(max_match_dist_m / 1000, 1)
  )
]

fwrite(table_si, file.path(out_dir, "TableS_AQS_all_metrics_distfiltered.csv"))
print(table_si)

# -----------------------------
# Main Figure 5: sulfate only
# No panel labels or years inside plots.
# Add (a)/(b) and year labels later in Word/Illustrator if needed.
# -----------------------------
if (requireNamespace("gridExtra", quietly = TRUE)) {
  p1980 <- plot_list[["1980_Sulfate_TSP"]]
  p1990 <- plot_list[["1990_Sulfate_TSP"]]
  
  pdf(
    file.path(out_dir, "Figure5_main_no_panel_labels.pdf"),
    width = 14,
    height = 5.2
  )
  
  gridExtra::grid.arrange(
    p1980,
    p1990,
    ncol = 2
  )
  
  dev.off()
  cat("Saved Figure5_main_no_panel_labels.pdf\n")
} else {
  warning("Package 'gridExtra' is not installed. Individual sulfate panels were saved, but the combined Figure 5 was not created.")
}

cat("\nDone. Outputs saved to:\n", out_dir, "\n")
