#!/usr/bin/env Rscript
###############################################################################
# Script: PWE_county_1940_1990.R
# Author: Xiaorong Shan
#
# Goal:
#   Compute county-level coal PM2.5 (HyADS, 36-km grid) and population-weighted
#   PM2.5 exposure (PWE) for 1940–1990 (decadal), and generate figures
#   WITHOUT capping and WITHOUT white (NA) counties on maps (force NA->0).
#
# Inputs:
#   (1) HyADS gridded PM2.5 .fst files (x,y,vals.out), Albers:
#         grids_pm25_total_YYYY.fst
#   (2) NHGIS county population time series (wide):
#         A00AA1940 ... A00AA1990
#
# Outputs (written to out_dir):
#   - county_pm25_pop_YYYY.csv
#   - county_pm25_pop_all_years.csv
#   - national_PWE_1940_1990.csv
#   - state_PWE_1940_1990.csv
#
# Figures (PDF):
#   MAIN:
#     - FIG_MAIN_contrib_log1p_scaled_1990median_1940_1990.pdf
#     - FIG_MAIN_national_PWE_trend_1940_1990.pdf
#   SI:
#     - FIG_SI_pm25_maps_1940_1990.pdf
#
# Notes:
#   - No capping/truncation is applied.
#   - Any missing counties on maps are forced to 0 (lowest color), no white.
#   - Contribution map uses monotone transform for readability:
#       log(1 + contrib / median_1990)
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
fst_template <- "grids_pm25_total_%d.fst"

out_dir <- file.path(base_dir, "pwe_outputs_decades_FULL_NOcap_NOwhite")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

cat("Output directory:\n", out_dir, "\n\n")
stopifnot(file.exists(pop_file))

# ===================== HELPERS =====================
fmt_geoid5 <- function(x) {
  x <- as.character(x)
  x <- gsub("\\D", "", x)
  x <- sprintf("%05d", as.integer(x))
  x
}

make_raster_albers <- function(fst_path, p4s) {
  dt <- as.data.table(fst::read_fst(fst_path))
  if (!all(c("x", "y", "vals.out") %in% names(dt))) {
    stop("FST must include x, y, vals.out. Found: ", paste(names(dt), collapse = ", "))
  }
  r <- raster::rasterFromXYZ(dt[, .(x, y, vals.out)])
  raster::crs(r) <- p4s
  r
}

area_weighted_mean_by_polygon <- function(r, polys_sf) {
  polys_sp <- methods::as(polys_sf, "Spatial")
  v_list <- raster::extract(r, polys_sp, weights = TRUE, na.rm = TRUE)

  v <- sapply(v_list, function(x) {
    if (is.null(x)) return(NA_real_)

    if (is.data.frame(x) || is.matrix(x)) {
      if (nrow(x) == 0) return(NA_real_)
      vals <- x[, 1]

      cand <- c("weight", "weights", "w", "wt", "frac", "fraction")
      wcol <- cand[cand %in% colnames(x)][1]

      if (!is.na(wcol)) {
        wts <- x[, wcol]
        return(stats::weighted.mean(vals, wts, na.rm = TRUE))
      }

      if (ncol(x) >= 2) {
        wts <- x[, 2]
        if (!all(is.na(wts)) && max(wts, na.rm = TRUE) <= 1.0001) {
          return(stats::weighted.mean(vals, wts, na.rm = TRUE))
        }
      }

      return(mean(vals, na.rm = TRUE))
    }

    if (is.numeric(x)) return(mean(x, na.rm = TRUE))
    NA_real_
  })

  as.numeric(v)
}

# ===================== 1) LOAD COUNTY GEOMETRY =====================
counties_sf <- USAboundaries::us_counties()
counties_sf <- counties_sf[, !duplicated(names(counties_sf))]

keep_cols <- intersect(c("geoid", "state_abbr", "stusps", "state_name", "name", "geometry"), names(counties_sf))
counties_sf <- counties_sf[, keep_cols]

if (!("state_abbr" %in% names(counties_sf)) && ("stusps" %in% names(counties_sf))) {
  counties_sf$state_abbr <- counties_sf$stusps
}

# CONUS only
counties_sf <- counties_sf[!(counties_sf$state_abbr %in% c("AK", "HI", "PR")), ]
counties_sf$geoid <- fmt_geoid5(counties_sf$geoid)
counties_sf <- sf::st_transform(counties_sf, crs = p4s)

cat("CONUS counties loaded:", nrow(counties_sf), "\n\n")

# ===================== 2) LOAD POPULATION (WIDE -> LONG) =====================
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
pop_long[, pop_year := as.numeric(pop_year)]
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
    pwe_national[[as.character(yy)]] <- data.frame(year = yy, PWE = NA_real_, nat_pop = NA_real_)
    next
  }

  r_yy <- make_raster_albers(fst_path, p4s)

  # Hopper-safe sanity (no raster::range)
  rmin <- raster::cellStats(r_yy, stat = "min", na.rm = TRUE)
  rmax <- raster::cellStats(r_yy, stat = "max", na.rm = TRUE)
  cat("  Raster min/max:", signif(rmin, 6), "/", signif(rmax, 6), "\n")

  cat("  Extracting county weighted means (this can be slow) ...\n")
  pm25_county <- area_weighted_mean_by_polygon(r_yy, counties_sf)

  cat("  County pm25 NA%:", round(mean(is.na(pm25_county)) * 100, 4),
      " | min/max:", signif(min(pm25_county, na.rm = TRUE), 6),
      "/", signif(max(pm25_county, na.rm = TRUE), 6), "\n")

  ctab <- sf::st_drop_geometry(counties_sf)

  dt <- data.table(
    geoid = fmt_geoid5(ctab$geoid),
    state_abbr = as.character(ctab$state_abbr),
    state_name = as.character(ctab$state_name),
    county_name = as.character(ctab$name),
    year  = yy,
    pm25  = as.numeric(pm25_county)
  )

  pop_year <- pop_dt[year == yy]
  dt <- merge(dt, pop_year, by = c("geoid", "year"), all.x = TRUE)

  # Valid rows for PWE math
  dt_calc <- dt[!is.na(pm25) & is.finite(pm25) &
                  !is.na(pop_year) & is.finite(pop_year) & pop_year > 0]

  if (nrow(dt_calc) == 0) {
    warning("Year ", yy, ": no valid rows after merge/filter; PWE is NA")
    pwe_val <- NA_real_
    nat_pop <- NA_real_
  } else {
    nat_pop <- sum(dt_calc$pop_year)
    pwe_val <- sum(dt_calc$pm25 * dt_calc$pop_year) / nat_pop
  }

  # Add burden and contribution to national PWE (μg/m^3)
  if (is.finite(nat_pop) && nat_pop > 0) {
    dt_calc[, burden := pm25 * pop_year]
    dt_calc[, contrib_nat := burden / nat_pop]  # μg/m^3
  } else {
    dt_calc[, burden := NA_real_]
    dt_calc[, contrib_nat := NA_real_]
  }

  cat("  Rows after merge+filter:", nrow(dt_calc), "\n")
  cat("  National PWE:", signif(pwe_val, 8), "\n\n")

  fwrite(dt_calc, file.path(out_dir, sprintf("county_pm25_pop_%d.csv", yy)))

  county_all[[as.character(yy)]] <- dt_calc
  pwe_national[[as.character(yy)]] <- data.frame(year = yy, PWE = pwe_val, nat_pop = nat_pop)

  # State-level PWE
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

cat("Saved combined outputs in:", out_dir, "\n\n")

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

# ===================== 6) MAIN MAP: COUNTY CONTRIBUTION TO NATIONAL PWE =====================
# Force NA->0 so there are no white counties.
cat("Creating MAIN map (no white): contribution to national PWE (log1p scaled to 1990 median)...\n")

map_dt <- copy(county_all_dt)
map_dt[, geoid := fmt_geoid5(geoid)]
map_dt[, year := as.integer(year)]

ref_med_1990 <- as.numeric(median(map_dt[year == 1990]$contrib_nat, na.rm = TRUE))
if (!is.finite(ref_med_1990) || ref_med_1990 <= 0) ref_med_1990 <- 1e-12

map_dt[, fill_log1p := log1p(contrib_nat / ref_med_1990)]
map_dt[!is.finite(fill_log1p) | is.na(fill_log1p), fill_log1p := 0]

map_dt2 <- as.data.frame(map_dt[, .(geoid, year, fill_log1p)])

map_sf <- merge(counties_sf, map_dt2, by = "geoid", all.x = TRUE)

# ---- HARD ZERO-FILL (no white counties) ----
map_sf$year <- factor(as.integer(as.character(map_sf$year)), levels = years)
map_sf$fill_log1p <- as.numeric(map_sf$fill_log1p)
map_sf$fill_log1p[!is.finite(map_sf$fill_log1p) | is.na(map_sf$fill_log1p)] <- 0

p_main_map <- ggplot(map_sf) +
  geom_sf(aes(fill = fill_log1p), color = "grey25", linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(size = 15)
  ) +
  labs(
    title = "County contribution to national population-weighted coal PM2.5 (scaled to 1990 median; log1p), 1940–1990",
    fill  = "log(1 + contrib / median_1990)\n(contrib in μg/m³)"
  ) +
  scale_fill_gradient(
  low = "black",
  high = "yellow",
  limits = c(0, p99_bur),
  na.value = "black",
  guide = guide_colorbar(na.translate = FALSE)
)

ggsave(file.path(out_dir, "FIG_MAIN_contrib_log1p_scaled_1990median_1940_1990.pdf"),
       p_main_map, width = 14, height = 4.6, dpi = 300)

cat("Reference median_1990 contrib_nat (μg/m^3):", signif(ref_med_1990, 8), "\n\n")

# ===================== 7) SI MAP: COUNTY MEAN PM2.5 (LINEAR) =====================
# Force NA->0 so there are no white counties.
cat("Creating SI map (no white): county mean PM2.5 (linear, no cap)...\n")

pm25_dt <- copy(county_all_dt)
pm25_dt[, geoid := fmt_geoid5(geoid)]
pm25_dt[, year := as.integer(year)]
pm25_dt[, pm25_plot := as.numeric(pm25)]
pm25_dt[!is.finite(pm25_plot) | is.na(pm25_plot), pm25_plot := 0]

pm25_sf <- merge(counties_sf, as.data.frame(pm25_dt[, .(geoid, year, pm25_plot)]),
                 by = "geoid", all.x = TRUE)

# ---- HARD ZERO-FILL (no white counties) ----
pm25_sf$year <- factor(as.integer(as.character(pm25_sf$year)), levels = years)
pm25_sf$pm25_plot <- as.numeric(pm25_sf$pm25_plot)
pm25_sf$pm25_plot[!is.finite(pm25_sf$pm25_plot) | is.na(pm25_sf$pm25_plot)] <- 0

p_pm25 <- ggplot(pm25_sf) +
  geom_sf(aes(fill = pm25_plot), color = "grey25", linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(size = 15)
  ) +
  labs(
    title = "County mean coal PM2.5 (HyADS; 36-km grid aggregated to counties), 1940–1990",
    fill = expression(PM[2.5]~(mu*g/m^3))
  ) +
  scale_fill_gradient(
  low = "black",
  high = "yellow",
  na.value = "black",
  guide = guide_colorbar(na.translate = FALSE)
)

ggsave(file.path(out_dir, "FIG_SI_pm25_maps_1940_1990.pdf"),
       p_pm25, width = 14, height = 4.6, dpi = 300)

cat("✅ Done. All outputs saved to:\n", out_dir, "\n")
