#!/usr/bin/env Rscript
###############################################################################
# Script: PWE_county_1940_1990_full.R
# Author: Xiaorong Shan
#
# Goal:
#   Compute county-level coal PM2.5 (HyADS, 36-km grid) and population-weighted
#   PM2.5 exposure (PWE) for 1940–1990 (decadal).
#
# Inputs:
#   (1) HyADS gridded PM2.5 .fst files at 36-km resolution, Albers (x,y,vals.out)
#       Expected filenames: grids_pm25_total_YYYY.fst
#   (2) NHGIS county population time series (wide), e.g. A00AA1940 ... A00AA1990
#
# Outputs (written to out_dir):
#   - county_pm25_pop_YYYY.csv                 (county-level results each decade)
#   - county_pm25_pop_all_years.csv            (stacked county results)
#   - national_PWE_1940_1990.csv               (national PWE time series)
#   - state_PWE_1940_1990.csv                  (state PWE time series)
#
# Visualization (PDF):
#   MAIN FIG (recommended for manuscript):
#     - FIG_MAIN_burden_log_maps_1940_1990.pdf  (log-scaled PM2.5×pop maps)
#     - FIG_MAIN_national_PWE_trend_1940_1990.pdf (PWE trend line)
#
# SI FIG (recommended for supplement):
#     - FIG_SI_pm25_maps_1940_1990.pdf          (county mean PM2.5 maps)
#
# Notes:
#   - For mapping "human exposure intensity", we visualize log(1 + PM2.5×pop).
#     This is unitless and intended for comparative visualization only.
#   - For quantitative reporting, use PWE (μg/m^3).
###############################################################################

suppressPackageStartupMessages({
  library(raster)
  library(sf)
  library(USAboundaries)
  library(data.table)
  library(ggplot2)
  library(fst)
  library(methods)
})

# ===================== USER SETTINGS =====================
years <- c(1940, 1950, 1960, 1970, 1980, 1990)

base_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"
pop_file <- "/scratch/xshan2/R_Code/powerplant/data/population/nhgis0012_csv/nhgis0012_ts_nominal_county.csv"

# Expected HyADS filenames:
fst_template <- "grids_pm25_total_%d.fst"

# Output folder (feel free to rename)
out_dir <- file.path(base_dir, "pwe_outputs_decades_FULL_hopperSafe")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Albers Equal-Area projection used by HyADS grids
p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

cat("Output directory:\n", out_dir, "\n\n")

stopifnot(file.exists(pop_file))

# ===================== HELPERS =====================

# Format county GEOID to 5-digit string (STATEFP(2)+COUNTYFP(3))
fmt_geoid5 <- function(x) {
  x <- as.character(x)
  x <- gsub("\\D", "", x)
  x <- sprintf("%05d", as.integer(x))
  x
}

# Convert .fst (x,y,vals.out) to RasterLayer and assign CRS
make_raster_albers <- function(fst_path, p4s) {
  dt <- as.data.table(fst::read_fst(fst_path))
  if (!all(c("x", "y", "vals.out") %in% names(dt))) {
    stop("FST must include x, y, vals.out. Found: ", paste(names(dt), collapse = ", "))
  }
  r <- raster::rasterFromXYZ(dt[, .(x, y, vals.out)])
  raster::crs(r) <- p4s
  r
}

# Robust area-weighted mean using raster::extract(weights=TRUE)
# Handles different extract return formats across raster versions.
area_weighted_mean_by_polygon <- function(r, polys_sf) {
  polys_sp <- methods::as(polys_sf, "Spatial")
  v_list <- raster::extract(r, polys_sp, weights = TRUE, na.rm = TRUE)

  v <- sapply(v_list, function(x) {
    if (is.null(x)) return(NA_real_)

    if (is.data.frame(x) || is.matrix(x)) {
      if (nrow(x) == 0) return(NA_real_)
      vals <- x[, 1]

      # Try common weight column names
      cand <- c("weight", "weights", "w", "wt", "frac", "fraction")
      wcol <- cand[cand %in% colnames(x)][1]

      if (!is.na(wcol)) {
        wts <- x[, wcol]
        return(stats::weighted.mean(vals, wts, na.rm = TRUE))
      }

      # Fallback: if there are >=2 columns, the 2nd column is often weights in [0,1]
      if (ncol(x) >= 2) {
        wts <- x[, 2]
        if (!all(is.na(wts)) && max(wts, na.rm = TRUE) <= 1.0001) {
          return(stats::weighted.mean(vals, wts, na.rm = TRUE))
        }
      }

      # Final fallback: simple mean (keeps pipeline running)
      return(mean(vals, na.rm = TRUE))
    }

    if (is.numeric(x)) return(mean(x, na.rm = TRUE))

    NA_real_
  })

  as.numeric(v)
}

# Cap extreme values for visualization without requiring extra packages
cap_at <- function(x, upper) {
  x2 <- x
  x2[is.na(x2)] <- 0
  x2[x2 > upper] <- upper
  x2
}

# ===================== 1) LOAD COUNTY GEOMETRY (ONCE) =====================
counties_sf <- USAboundaries::us_counties()
counties_sf <- counties_sf[, !duplicated(names(counties_sf))]

keep_cols <- intersect(c("geoid", "state_abbr", "stusps", "state_name", "name", "geometry"), names(counties_sf))
counties_sf <- counties_sf[, keep_cols]

# Ensure state_abbr exists
if (!("state_abbr" %in% names(counties_sf)) && ("stusps" %in% names(counties_sf))) {
  counties_sf$state_abbr <- counties_sf$stusps
}

# CONUS only
counties_sf <- counties_sf[!(counties_sf$state_abbr %in% c("AK", "HI", "PR")), ]

# Fix GEOID early (critical for reliable merges)
counties_sf$geoid <- fmt_geoid5(counties_sf$geoid)

# Reproject to match HyADS grid CRS
counties_sf <- sf::st_transform(counties_sf, crs = p4s)

cat("CONUS counties loaded:", nrow(counties_sf), "\n\n")

# ===================== 2) LOAD POPULATION (WIDE -> LONG FOR TARGET YEARS) =====================
pop_ts <- fread(pop_file)

pop_cols <- paste0("A00AA", years)
cols_needed <- c("STATEFP", "COUNTYFP", pop_cols)
stopifnot(all(cols_needed %in% names(pop_ts)))

pop_long <- melt(
  pop_ts[, ..cols_needed],
  id.vars = c("STATEFP", "COUNTYFP"),
  variable.name = "VAR",
  value.name = "pop_year"
)

pop_long[, year := as.integer(sub("A00AA", "", VAR))]
pop_long[, geoid := paste0(sprintf("%02d", as.integer(STATEFP)),
                           sprintf("%03d", as.integer(COUNTYFP)))]
pop_long[, geoid := fmt_geoid5(geoid)]
pop_dt <- pop_long[year %in% years, .(geoid, year, pop_year)]

cat("Population rows (all decades):", nrow(pop_dt), "\n\n")

# ===================== 3) MAIN LOOP: COUNTY PM25 + PWE =====================
pwe_national <- list()
pwe_state    <- list()
county_all   <- list()

for (yy in years) {

  cat("====== Processing year:", yy, "======\n")

  fst_path <- file.path(base_dir, sprintf(fst_template, yy))
  if (!file.exists(fst_path)) {
    warning("Missing file: ", fst_path, " (skipping year)")
    pwe_national[[as.character(yy)]] <- data.frame(year = yy, PWE = NA_real_)
    next
  }

  # 3.1 Read raster
  r_yy <- make_raster_albers(fst_path, p4s)

  # quick sanity
  rmin <- raster::cellStats(r_yy, stat = "min", na.rm = TRUE)
  rmax <- raster::cellStats(r_yy, stat = "max", na.rm = TRUE)
  cat("  Raster min/max:", signif(rmin, 6), "/", signif(rmax, 6), "\n")

  # 3.2 County area-weighted PM2.5 means
  cat("  Extracting county weighted means (this can be slow) ...\n")
  pm25_county <- area_weighted_mean_by_polygon(r_yy, counties_sf)

  cat("  County pm25 NA%:", round(mean(is.na(pm25_county)) * 100, 4),
      " | min/max:", signif(min(pm25_county, na.rm = TRUE), 6),
      "/", signif(max(pm25_county, na.rm = TRUE), 6), "\n")

  # 3.3 Build county table (no geometry)
  ctab <- sf::st_drop_geometry(counties_sf)

  dt <- data.table(
    geoid = fmt_geoid5(ctab$geoid),
    state_abbr = as.character(ctab$state_abbr),
    state_name = as.character(ctab$state_name),
    county_name = as.character(ctab$name),
    year  = yy,
    pm25  = as.numeric(pm25_county)
  )

  # 3.4 Merge population for the same year
  pop_year <- pop_dt[year == yy]
  dt <- merge(dt, pop_year, by = c("geoid", "year"), all.x = TRUE)

  # 3.5 Compute burden and contributions (for reporting & maps)
  # Keep only valid rows for PWE (pm25 and population both present)
  dt_calc <- dt[!is.na(pm25) & !is.na(pop_year) & pop_year > 0]

  if (nrow(dt_calc) == 0) {
    warning("Year ", yy, ": no valid rows after merge/filter; PWE is NA")
    pwe_val <- NA_real_
  } else {
    pwe_val <- sum(dt_calc$pm25 * dt_calc$pop_year) / sum(dt_calc$pop_year)
  }

  # Add burden/contrib for county outputs
  dt_calc[, burden  := pm25 * pop_year]
  dt_calc[, contrib := (pm25 * pop_year) / sum(pop_year)]

  cat("  Rows after merge+filter:", nrow(dt_calc), "\n")
  cat("  National PWE:", signif(pwe_val, 8), "\n\n")

  # 3.6 Save county CSV for this year
  out_csv <- file.path(out_dir, sprintf("county_pm25_pop_%d.csv", yy))
  fwrite(dt_calc, out_csv)

  # 3.7 Store for combined outputs
  county_all[[as.character(yy)]] <- dt_calc
  pwe_national[[as.character(yy)]] <- data.frame(year = yy, PWE = pwe_val)

  # 3.8 State-level PWE
  # PWE_state = sum(pm25*pop) / sum(pop) within each state
  st <- dt_calc[, .(
    PWE_state = sum(pm25 * pop_year) / sum(pop_year),
    pop_total = sum(pop_year)
  ), by = .(year, state_abbr, state_name)]
  pwe_state[[as.character(yy)]] <- st
}

# ===================== 4) WRITE COMBINED TABLES =====================
county_all_dt <- rbindlist(county_all, fill = TRUE)
fwrite(county_all_dt, file.path(out_dir, "county_pm25_pop_all_years.csv"))

national_df <- rbindlist(pwe_national, fill = TRUE)
fwrite(national_df, file.path(out_dir, "national_PWE_1940_1990.csv"))

state_df <- rbindlist(pwe_state, fill = TRUE)
fwrite(state_df, file.path(out_dir, "state_PWE_1940_1990.csv"))

cat("Saved combined outputs:\n")
cat("  - county_pm25_pop_all_years.csv\n")
cat("  - national_PWE_1940_1990.csv\n")
cat("  - state_PWE_1940_1990.csv\n\n")

# ===================== 5) MAIN FIG: NATIONAL PWE TREND (μg/m^3) =====================
p_trend <- ggplot(national_df, aes(x = year, y = PWE)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 13) +
  labs(
    title = "National population-weighted coal PM2.5 exposure (1940–1990)",
    x = "Year",
    y = expression(PWE~(mu*g/m^3))
  )

ggsave(file.path(out_dir, "FIG_MAIN_national_PWE_trend_1940_1990.pdf"),
       p_trend, width = 7.6, height = 4.4, dpi = 300)

# ===================== 6) MAPS (MAIN + SI) =====================
# Join to county geometry for mapping
# For visualization, fill missing counties as 0 (lowest color), per your request.
map_dt <- copy(county_all_dt)
map_dt[, geoid := fmt_geoid5(geoid)]
map_dt[, year := factor(year, levels = years, labels = years)]

# Compute map variables
map_dt[, pm25_plot := pm25]
map_dt[is.na(pm25_plot), pm25_plot := 0]

map_dt[, burden := pm25_plot * pop_year]
map_dt[is.na(burden), burden := 0]

map_dt[, burden_log := log1p(burden)]

# Global caps for consistent color across decades
p99_pm25 <- as.numeric(quantile(map_dt$pm25_plot, 0.99, na.rm = TRUE))
if (!is.finite(p99_pm25) || p99_pm25 <= 0) p99_pm25 <- max(map_dt$pm25_plot, na.rm = TRUE)

p99_bur <- as.numeric(quantile(map_dt$burden_log, 0.99, na.rm = TRUE))
if (!is.finite(p99_bur) || p99_bur <= 0) p99_bur <- max(map_dt$burden_log, na.rm = TRUE)

map_dt[, pm25_cap := cap_at(pm25_plot, p99_pm25)]
map_dt[, burden_log_cap := cap_at(burden_log, p99_bur)]

# Merge with sf (keeps geometry)
map_sf <- merge(counties_sf, as.data.frame(map_dt[, .(geoid, year, pm25_cap, burden_log_cap)]),
                by = "geoid", all.x = TRUE)

# Replace any remaining NA in sf after merge as 0 (lowest)
map_sf$pm25_cap[is.na(map_sf$pm25_cap)] <- 0
map_sf$burden_log_cap[is.na(map_sf$burden_log_cap)] <- 0

# --- SI FIG: county mean PM2.5 maps (capped at global P99) ---
p_pm25 <- ggplot(map_sf) +
  geom_sf(aes(fill = pm25_cap), color = "grey25", linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "County mean coal PM2.5 (capped at global P99), 1940–1990",
    fill = expression(paste("PM"[2.5], " (capped)"))
  ) +
  scale_fill_gradient(low = "black", high = "yellow", limits = c(0, p99_pm25))

ggsave(file.path(out_dir, "FIG_SI_pm25_maps_1940_1990.pdf"),
       p_pm25, width = 14, height = 4.6, dpi = 300)

# --- MAIN FIG: exposure burden maps (log-scaled, capped at global P99) ---
p_burden <- ggplot(map_sf) +
  geom_sf(aes(fill = burden_log_cap), color = "grey25", linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "Population exposure burden to coal PM2.5 (log-scaled; capped at global P99), 1940–1990",
    fill = "log(1 + PM2.5×pop)\n(capped)"
  ) +
  scale_fill_gradient(low = "black", high = "yellow", limits = c(0, p99_bur))

ggsave(file.path(out_dir, "FIG_MAIN_burden_log_maps_1940_1990.pdf"),
       p_burden, width = 14, height = 4.6, dpi = 300)

cat("✅ Done. All outputs saved to:\n", out_dir, "\n")
