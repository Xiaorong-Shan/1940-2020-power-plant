#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

###############################################################################
# Script: PWE_county_1940_1990_with_Fig2B.R
#
# Purpose:
#   1. Compute county-level coal PM2.5 exposure and population-weighted exposure
#      for 1940–1990.
#   2. Generate:
#        - National PWE trend
#        - County contribution map
#        - County PM2.5 map
#        - Fig. 2b: top contributing states and counties bar chart
#
# Important:
#   - No capping is applied to original calculation.
#   - Map color limits are lowered only for visualization.
#   - Values above visual maximum are shown using the top color.
###############################################################################

suppressPackageStartupMessages({
  library(raster)
  library(sf)
  library(USAboundaries)
  library(data.table)
  library(ggplot2)
  library(fst)
  library(methods)
  library(scales)
})

# ===================== USER SETTINGS =====================
years <- c(1940, 1950, 1960, 1970, 1980, 1990)

base_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"

pop_file <- "/scratch/xshan2/R_Code/powerplant/data/population/nhgis0013_ts_nominal_county.csv"

fst_template <- "grids_pm25_total_%d.fst"

out_dir <- file.path(base_dir, "pwe_outputs_decades_FULL_NOcap_NOwhite")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# Visualization-only percentile limits
main_color_q <- 0.95
pm25_color_q <- 0.95

# Fig. 2b settings
fig2b_year <- 1990
top_n <- 10

cat("Output directory:\n", out_dir, "\n\n")
cat("Population file:\n", pop_file, "\n\n")

stopifnot(file.exists(pop_file))

# ===================== HELPER FUNCTIONS =====================

fmt_geoid5 <- function(x) {
  x <- as.character(x)
  x <- gsub("\\D", "", x)
  x <- sprintf("%05d", as.integer(x))
  x
}

make_raster_albers <- function(fst_path, p4s) {
  dt <- as.data.table(fst::read_fst(fst_path))

  if (!all(c("x", "y", "vals.out") %in% names(dt))) {
    stop("FST must include x, y, vals.out.")
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

get_plot_max <- function(x, q = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]

  if (length(x) == 0) return(1)

  x_pos <- x[x > 0]

  if (length(x_pos) > 0) {
    out <- as.numeric(quantile(x_pos, probs = q, na.rm = TRUE, names = FALSE))
  } else {
    out <- max(x, na.rm = TRUE)
  }

  if (!is.finite(out) || out <= 0) out <- 1
  out
}

# ===================== 1) LOAD COUNTY GEOMETRY =====================

counties_sf <- USAboundaries::us_counties()
counties_sf <- counties_sf[, !duplicated(names(counties_sf))]

keep_cols <- intersect(
  c("geoid", "state_abbr", "stusps", "state_name", "name", "geometry"),
  names(counties_sf)
)

counties_sf <- counties_sf[, keep_cols]

if (!("state_abbr" %in% names(counties_sf)) && ("stusps" %in% names(counties_sf))) {
  counties_sf$state_abbr <- counties_sf$stusps
}

if (!("state_name" %in% names(counties_sf))) {
  counties_sf$state_name <- NA_character_
}

if (!("name" %in% names(counties_sf))) {
  counties_sf$name <- NA_character_
}

# CONUS only
counties_sf <- counties_sf[!(counties_sf$state_abbr %in% c("AK", "HI", "PR")), ]
counties_sf$geoid <- fmt_geoid5(counties_sf$geoid)
counties_sf <- sf::st_transform(counties_sf, crs = p4s)

cat("CONUS counties loaded:", nrow(counties_sf), "\n\n")

# ===================== 2) LOAD POPULATION DATA =====================

pop_ts <- fread(pop_file)

pop_cols <- paste0("A00AA", years)
cols_needed <- c("STATEFP", "COUNTYFP", pop_cols)

missing_cols <- setdiff(cols_needed, names(pop_ts))

if (length(missing_cols) > 0) {
  stop(
    "Population file is missing required columns:\n",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are:\n",
    paste(names(pop_ts), collapse = ", ")
  )
}

pop_long <- melt(
  pop_ts[, ..cols_needed],
  id.vars = c("STATEFP", "COUNTYFP"),
  variable.name = "VAR",
  value.name = "pop_year"
)

pop_long[, year := as.integer(sub("A00AA", "", VAR))]
pop_long[, geoid := paste0(
  sprintf("%02d", as.integer(STATEFP)),
  sprintf("%03d", as.integer(COUNTYFP))
)]
pop_long[, geoid := fmt_geoid5(geoid)]
pop_long[, pop_year := as.numeric(pop_year)]

pop_dt <- pop_long[year %in% years, .(geoid, year, pop_year)]

cat("Population rows:", nrow(pop_dt), "\n\n")

# ===================== 3) COMPUTE COUNTY PM2.5 AND PWE =====================

pwe_national <- list()
pwe_state <- list()
county_all <- list()

for (yy in years) {

  cat("====== Processing year:", yy, "======\n")

  fst_path <- file.path(base_dir, sprintf(fst_template, yy))

  if (!file.exists(fst_path)) {
    warning("Missing file: ", fst_path)
    next
  }

  r_yy <- make_raster_albers(fst_path, p4s)

  rmin <- raster::cellStats(r_yy, stat = "min", na.rm = TRUE)
  rmax <- raster::cellStats(r_yy, stat = "max", na.rm = TRUE)

  cat("  Raster min/max:", signif(rmin, 6), "/", signif(rmax, 6), "\n")
  cat("  Extracting county weighted means...\n")

  pm25_county <- area_weighted_mean_by_polygon(r_yy, counties_sf)

  cat(
    "  County pm25 NA%:", round(mean(is.na(pm25_county)) * 100, 4),
    " | min/max:", signif(min(pm25_county, na.rm = TRUE), 6),
    "/", signif(max(pm25_county, na.rm = TRUE), 6), "\n"
  )

  ctab <- sf::st_drop_geometry(counties_sf)

  dt <- data.table(
    geoid = fmt_geoid5(ctab$geoid),
    state_abbr = as.character(ctab$state_abbr),
    state_name = as.character(ctab$state_name),
    county_name = as.character(ctab$name),
    year = yy,
    pm25 = as.numeric(pm25_county)
  )

  pop_year <- pop_dt[year == yy]

  dt <- merge(
    dt,
    pop_year,
    by = c("geoid", "year"),
    all.x = TRUE
  )

  dt_calc <- dt[
    !is.na(pm25) & is.finite(pm25) &
      !is.na(pop_year) & is.finite(pop_year) & pop_year > 0
  ]

  if (nrow(dt_calc) == 0) {
    warning("No valid county-population rows for year ", yy)
    next
  }

  nat_pop <- sum(dt_calc$pop_year)
  pwe_val <- sum(dt_calc$pm25 * dt_calc$pop_year) / nat_pop

  dt_calc[, burden := pm25 * pop_year]
  dt_calc[, contrib_nat := burden / nat_pop]

  cat("  Rows after merge+filter:", nrow(dt_calc), "\n")
  cat("  National PWE:", signif(pwe_val, 8), "\n\n")

  fwrite(dt_calc, file.path(out_dir, sprintf("county_pm25_pop_%d.csv", yy)))

  county_all[[as.character(yy)]] <- dt_calc

  pwe_national[[as.character(yy)]] <- data.frame(
    year = yy,
    PWE = pwe_val,
    nat_pop = nat_pop
  )

  st <- dt_calc[, .(
    PWE_state = sum(pm25 * pop_year) / sum(pop_year),
    pop_total = sum(pop_year)
  ), by = .(year, state_abbr, state_name)]

  pwe_state[[as.character(yy)]] <- st
}

# ===================== 4) SAVE COMPUTED TABLES =====================

county_all_dt <- rbindlist(county_all, fill = TRUE)
national_df <- rbindlist(pwe_national, fill = TRUE)
state_df <- rbindlist(pwe_state, fill = TRUE)

if (nrow(county_all_dt) == 0) {
  stop("No county-level output was generated. Check FST file paths and population merge.")
}

fwrite(county_all_dt, file.path(out_dir, "county_pm25_pop_all_years.csv"))
fwrite(national_df, file.path(out_dir, "national_PWE_1940_1990.csv"))
fwrite(state_df, file.path(out_dir, "state_PWE_1940_1990.csv"))

cat("Saved combined output tables.\n\n")

# ===================== 5) FIGURE: NATIONAL PWE TREND =====================

p_trend <- ggplot(national_df, aes(x = year, y = PWE)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Year",
    y = expression(PWE~(mu*g/m^3))
  )

ggsave(
  file.path(out_dir, "FIG_MAIN_national_PWE_trend_1940_1990.pdf"),
  p_trend,
  width = 7.6,
  height = 4.4,
  dpi = 300
)

# ===================== 6) FIGURE: COUNTY CONTRIBUTION MAP =====================

cat("Creating county contribution map...\n")

county_all_dt <- as.data.table(county_all_dt)
county_all_dt[, geoid := fmt_geoid5(geoid)]
county_all_dt[, year := as.integer(year)]
county_all_dt[, contrib_nat := as.numeric(contrib_nat)]
county_all_dt[!is.finite(contrib_nat) | is.na(contrib_nat), contrib_nat := 0]

panel <- CJ(
  geoid = unique(counties_sf$geoid),
  year = years,
  unique = TRUE
)

panel <- merge(
  panel,
  county_all_dt[, .(geoid, year, contrib_nat)],
  by = c("geoid", "year"),
  all.x = TRUE
)

panel[is.na(contrib_nat) | !is.finite(contrib_nat), contrib_nat := 0]

med_1990 <- panel[year == 1990 & contrib_nat > 0, median(contrib_nat, na.rm = TRUE)]

if (!is.finite(med_1990) || med_1990 <= 0) {
  med_1990 <- 1e-12
}

panel[, fill_log1p := log1p(contrib_nat / med_1990)]
panel[is.na(fill_log1p) | !is.finite(fill_log1p), fill_log1p := 0]
panel[, year := factor(year, levels = years, labels = years)]

map_sf <- merge(
  counties_sf,
  as.data.frame(panel[, .(geoid, year, fill_log1p)]),
  by = "geoid",
  all.x = FALSE,
  all.y = TRUE
)

map_sf$fill_log1p <- as.numeric(map_sf$fill_log1p)
map_sf$fill_log1p[is.na(map_sf$fill_log1p) | !is.finite(map_sf$fill_log1p)] <- 0

main_plot_max <- get_plot_max(map_sf$fill_log1p, q = main_color_q)

p_main_map <- ggplot(map_sf) +
  geom_sf(
    aes(fill = fill_log1p),
    color = NA
  ) +
  facet_wrap(~ year, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 8, lineheight = 0.9),
    legend.text = element_text(size = 7),
    panel.spacing.x = grid::unit(0.15, "lines"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  labs(
    fill = "Contrib.\nto PWE"
  ) +
  scale_fill_gradient(
    low = "black",
    high = "yellow",
    limits = c(0, main_plot_max),
    oob = scales::squish,
    na.value = "black",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

ggsave(
  file.path(out_dir, "FIG_MAIN_contrib_log1p_scaled_1990median_1940_1990.pdf"),
  p_main_map,
  width = 14,
  height = 4.6,
  dpi = 300
)

cat("Reference median_1990 contrib_nat:", signif(med_1990, 8), "\n")
cat("Main map visual max:", signif(main_plot_max, 8), "\n\n")

# ===================== 7) FIGURE: COUNTY MEAN PM2.5 MAP =====================

cat("Creating county PM2.5 map...\n")

pm25_dt <- copy(county_all_dt)
pm25_dt[, geoid := fmt_geoid5(geoid)]
pm25_dt[, year := as.integer(year)]
pm25_dt[, pm25_plot := as.numeric(pm25)]
pm25_dt[!is.finite(pm25_plot) | is.na(pm25_plot), pm25_plot := 0]

panel_pm25 <- CJ(
  geoid = unique(counties_sf$geoid),
  year = years,
  unique = TRUE
)

panel_pm25 <- merge(
  panel_pm25,
  pm25_dt[, .(geoid, year, pm25_plot)],
  by = c("geoid", "year"),
  all.x = TRUE
)

panel_pm25[is.na(pm25_plot) | !is.finite(pm25_plot), pm25_plot := 0]
panel_pm25[, year := factor(year, levels = years, labels = years)]

pm25_sf <- merge(
  counties_sf,
  as.data.frame(panel_pm25[, .(geoid, year, pm25_plot)]),
  by = "geoid",
  all.x = FALSE,
  all.y = TRUE
)

pm25_sf$pm25_plot <- as.numeric(pm25_sf$pm25_plot)
pm25_sf$pm25_plot[is.na(pm25_sf$pm25_plot) | !is.finite(pm25_sf$pm25_plot)] <- 0

pm25_plot_max <- get_plot_max(pm25_sf$pm25_plot, q = pm25_color_q)

p_pm25 <- ggplot(pm25_sf) +
  geom_sf(
    aes(fill = pm25_plot),
    color = NA
  ) +
  facet_wrap(~ year, nrow = 1, drop = TRUE) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 8, lineheight = 0.9),
    legend.text = element_text(size = 7),
    panel.spacing.x = grid::unit(0.15, "lines"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  labs(
    fill = "PM2.5\n(ug/m3)"
  ) +
  scale_fill_gradient(
    low = "black",
    high = "yellow",
    limits = c(0, pm25_plot_max),
    oob = scales::squish,
    na.value = "black",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

ggsave(
  file.path(out_dir, "FIG_SI_pm25_maps_1940_1990.pdf"),
  p_pm25,
  width = 14,
  height = 4.6,
  dpi = 300
)

cat("PM2.5 map visual max:", signif(pm25_plot_max, 8), "\n\n")

# ===================== 8) FIG. 2B: TOP CONTRIBUTING STATES AND COUNTIES =====================

cat("Creating Fig. 2b: top contributing states and counties...\n")

fig2b_dt <- copy(county_all_dt)
fig2b_dt[, year := as.integer(year)]
fig2b_dt[, contrib_nat := as.numeric(contrib_nat)]
fig2b_dt[!is.finite(contrib_nat) | is.na(contrib_nat), contrib_nat := 0]

# ---- Top contributing states ----
top_states <- fig2b_dt[
  year == fig2b_year,
  .(
    contribution = sum(contrib_nat, na.rm = TRUE)
  ),
  by = .(state_abbr)
][
  order(-contribution)
][
  1:top_n
]

top_states[, label := state_abbr]
top_states[, type := "States"]

# ---- Top contributing counties ----
top_counties <- fig2b_dt[
  year == fig2b_year,
  .(
    contribution = sum(contrib_nat, na.rm = TRUE)
  ),
  by = .(county_name, state_abbr)
][
  order(-contribution)
][
  1:top_n
]

top_counties[, label := paste0(county_name, ", ", state_abbr)]
top_counties[, type := "Counties"]

# ---- Combine states and counties into one plotting table ----
fig2b_plot_dt <- rbind(
  top_states[, .(type, label, contribution)],
  top_counties[, .(type, label, contribution)]
)

fig2b_plot_dt[, label := factor(label, levels = rev(unique(label)))]

# ---- Plot Fig. 2b without title ----
p_fig2b <- ggplot(fig2b_plot_dt, aes(x = label, y = contribution)) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(~ type, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    x = NULL,
    y = expression("Contribution to national PWE"~(mu*g/m^3))
  )

ggsave(
  file.path(out_dir, paste0("FIG_2B_top_contributing_states_counties_", fig2b_year, ".pdf")),
  p_fig2b,
  width = 10,
  height = 5,
  dpi = 300
)

fwrite(
  top_states,
  file.path(out_dir, paste0("FIG_2B_top_states_", fig2b_year, ".csv"))
)

fwrite(
  top_counties,
  file.path(out_dir, paste0("FIG_2B_top_counties_", fig2b_year, ".csv"))
)

cat("Fig. 2b saved.\n\n")

# ===================== 9) SAVE VISUAL SCALE LIMITS =====================

plot_scale_limits <- data.table(
  figure = c("MAIN_contribution_map", "SI_pm25_map"),
  visual_percentile = c(main_color_q, pm25_color_q),
  scale_min = c(0, 0),
  scale_max = c(main_plot_max, pm25_plot_max),
  note = c(
    "Visualization only; values above scale_max are shown with top color using scales::squish.",
    "Visualization only; values above scale_max are shown with top color using scales::squish."
  )
)

fwrite(
  plot_scale_limits,
  file.path(out_dir, "plot_visual_scale_limits.csv")
)

cat("All figures and outputs saved in:\n", out_dir, "\n")
